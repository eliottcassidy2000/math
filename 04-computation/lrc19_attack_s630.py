#!/usr/bin/env python3
"""
S630 — attack on LRC(19).  Convention: n=19 runners, 18 speeds, gap 1/19.
LRC(19) <= C'(19): if 19|v for some speed, M(S)>1/19  (THM-398).
2n-1 = 37 is PRIME with 2 a primitive root => the UNRAMIFIED case (THM-407): the favorable side.
A good multiplier at shell 37 gives M>=2/37 > 1/19 (loose).  Crux: the residual configs that FAIL
the shell-37 dodge (residues cover all 18 +-pairs of (Z/37)*) — are they all loose via B' or a
slightly higher shell?  If yes, C'(19) closes.
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, sys
sys.path.insert(0, '04-computation')
from lrc_fast_s628 import gap_shells, good_mult_fast, gap_brute

N = 19; LVL = Fr(1, N)

def shell_witness(S, Mmax):
    """best (value, m, a) over shells m<=Mmax, a in (Z/m)* (the actual witness, not just value)."""
    best = (Fr(0), None, None)
    for m in range(2, Mmax+1):
        for a in range(1, m//2+1):
            if gcd(a, m) != 1: continue
            md = min(min((a*v) % m, m-((a*v) % m)) for v in S)
            val = Fr(md, m)
            if val > best[0]: best = (val, m, a)
    return best

def dodge37(S):
    if any(v % 37 == 0 for v in S): return None
    return good_mult_fast(S, 37)

def dominant(S):
    """Lemma B: some runner > (n-1)*max(others) => loose."""
    S = sorted(S)
    return S[-1] > (N-1) * (S[-2] if len(S) > 1 else 0)

if __name__ == "__main__":
    print("LRC(19): 2n-1=37 prime, 2 primitive root (unramified). good mult @37 => M>=2/37>1/19.\n")

    print("[A] Near-AP multiple-of-19 configs — the crux (do they pass shell-37 or need a higher shell?)")
    AP = list(range(1, 19))   # {1,...,18}, the tight AP (M=1/19)
    probes = [
        ("AP {1..18} (tight, no mult of 19)", AP),
        ("{1..17,19} = AP w/ 18->19", list(range(1, 18))+[19]),
        ("{1..18,19}\\{1}? {2..19}", list(range(2, 20))),
        ("{1..17,38}", list(range(1, 18))+[38]),
        ("{1..16,18,19}", list(range(1, 17))+[18, 19]),
    ]
    for name, S in probes:
        S = sorted(set(S))
        if len(S) != 18:
            print(f"   {name}: (size {len(S)}, skip)"); continue
        g = reduce(gcd, S)
        d37 = dodge37(S)
        val, m, a = shell_witness(S, 74)        # search shells up to 2*37
        mult19 = any(v % 19 == 0 for v in S)
        print(f"   {name}: mult19={mult19} gcd={g} dodge@37={d37} "
              f"| best shell witness M>={val}={float(val):.5f} at {a}/{m} (>1/19? {val>LVL})")

    print("\n[B] Among multiple-of-19 configs, find shell-37-dodge FAILURES that are NOT dominant")
    import random; rng = random.Random(0)
    fails = []; tested = 0
    for _ in range(60000):
        # build 18 distinct speeds, one == 0 mod 19, smallish
        S = {19*rng.randint(1, 2)}
        while len(S) < 18: S.add(rng.randint(1, 40))
        S = tuple(sorted(S))
        if len(S) != 18: continue
        if reduce(gcd, S) != 1: continue
        tested += 1
        if dodge37(S) is not None: continue      # loose via shell 37
        if dominant(S): continue                 # loose via Lemma B
        fails.append(S)
    print(f"   tested {tested} multiple-of-19 configs; shell-37-fail AND non-dominant: {len(fails)}")
    for S in fails[:6]:
        val, m, a = shell_witness(S, 74)
        print(f"     {S}  -> best higher-shell witness M>={val}={float(val):.5f} at {a}/{m} (>1/19? {val>LVL})")

    print("\n[C] the structural failure set: residues mod 37 cover all 18 +-pairs.")
    S = sorted(set(list(range(1, 18))+[19]))
    res = sorted(v % 37 for v in S)
    pairs = set(frozenset((u, 37-u)) for u in range(1, 19))
    hit = set(frozenset((r, 37-r)) for r in res if r != 0)
    print(f"   {tuple(S)} residues mod37={res};  pairs hit {len(hit)}/18 (cover all => dodge fails)")
