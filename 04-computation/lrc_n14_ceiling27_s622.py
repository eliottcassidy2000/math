#!/usr/bin/env python3
"""
S622 — stress-test the ceiling:  does EVERY loose config at n=14 escape via a
twisted-shell clock a/m with m <= 27 = 2n-1 (residue band-distance >= 2 for 15<=m<=27)?
A no-escape config (within m<=27) would be the true residual / candidate counterexample.
Adversarial search: hill-climb to MAXIMISE the minimal escaping m.
"""
from math import gcd
import itertools, random

def escapes(S, m):
    for a in range(1, m):
        if gcd(a, m) != 1: continue
        if all(14*min((a*s) % m, m-((a*s) % m)) > m for s in S):
            return a
    return None

def min_escape_m(S, mmax=80):
    for m in range(2, mmax+1):
        a = escapes(S, m)
        if a is not None: return m, a
    return None, None

def is_loose_via_anyclock(S):
    """quick necessary check that S is loose (some 1-clock or the structure); we instead
       rely on min_escape_m finding any beat-1/14 witness, which proves looseness directly."""
    return min_escape_m(S)[0] is not None

def random_mult14(rng, lo=1, hi=40):
    while True:
        S = {14*rng.randint(1,2)}
        for m in rng.sample([2,3,4,5,6,7,8,9,10,11,12,13], rng.randint(7,12)):
            S.add(m*rng.randint(1,3))
        while len(S) < 13: S.add(rng.randint(lo,hi))
        S = tuple(sorted(S))
        if len(S)==13:
            g=0
            for x in S: g=gcd(g,x)
            if g==1: return S

if __name__ == "__main__":
    rng = random.Random(7)
    worst = (0, None)      # (minimal escaping m, config)
    over27 = []            # configs whose minimal escape needs m>27
    N = 40000
    for _ in range(N):
        S = random_mult14(rng)
        m, a = min_escape_m(S)
        if m is None:
            over27.append((S, ">80")); continue
        if m > 27: over27.append((S, m))
        if m > worst[0]: worst = (m, S)
    print(f"searched {N} adversarial multiple-of-14 configs.")
    print(f"WORST (largest minimal escaping m): m={worst[0]}  S={worst[1]}")
    print(f"configs needing m>27: {len(over27)}  {over27[:5]}")

    # hill-climb from the worst to push m up
    cur_m, cur = worst
    for it in range(3000):
        S = list(cur); i = rng.randrange(13)
        S2 = S[:]; S2[i] = rng.randint(1, 60)
        S2 = tuple(sorted(set(S2)))
        if len(S2) != 13: continue
        if 14 not in S2 and 28 not in S2: continue
        g=0
        for x in S2: g=gcd(g,x)
        if g!=1: continue
        m,_ = min_escape_m(S2)
        if m and m > cur_m:
            cur_m, cur = m, S2
    print(f"hill-climb best: minimal escaping m={cur_m}  S={cur}")
    print(f"=> ceiling 2n-1=27 {'HOLDS' if cur_m<=27 and not over27 else 'VIOLATED'} on this search")
