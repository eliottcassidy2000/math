#!/usr/bin/env python3
"""
S625 — the PRIME-SHELL DODGE: a uniform-in-n window lemma for LRC, and the cyclotomic
tower it lives in.  (Extends THM-412/415, HYP-2220/2230; the disproof port.)

Convention (THM-398): n runners, n-1 speeds, gap 1/n, M(S)=max_t min_i ||v_i t||.

THEOREM (prime-shell dodge, elementary).  Let p be prime, p ∤ v_i for all i.  For a unit
multiplier a∈(ℤ/p)*, the witness t=a/p has min_i||v_i a/p|| = (1/p) min_i dist(a v_i mod p).
A runner v_i is 'banded' (dist ≤1, i.e. a v_i ∈{0,±1}) for exactly the 2 multipliers a=±v_i^{-1}
(a v_i=0 impossible).  So at most 2(n-1) multipliers are bad; if  p-1 > 2(n-1)  a GOOD multiplier
exists, giving  M(S) ≥ 2/p.
The multipliers (ℤ/p)* = Gal(ℚ(ζ_p)/ℚ); the ±-pairing = complex conjugation; the count
'p-1 > 2(n-1)' is the finite WINDOW LEMMA (enough modulus-1 rotations to dodge the width-2 band).

Taking q = smallest prime > 2n-1 (≤ 4n-3 by Bertrand) gives the UNIFORM bound M ≥ 2/q > 1/(2n-1)
for every config with all speeds < q (in particular the near-tight regime, THM-411).
The OPTIMAL value is 2/(2n-1) (THM-415); reaching it needs the COMPOSITE shell 2n-1 itself —
the ramified/cyclotomic-tower case (e.g. n=14: 2n-1=27=3³).
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, sys

def isprime(k):
    if k < 2: return False
    i = 2
    while i*i <= k:
        if k % i == 0: return False
        i += 1
    return True

def nextprime(k):
    k += 1
    while not isprime(k): k += 1
    return k

sys.path.insert(0, '04-computation')
from lrc_n14_flowshells_s622 import gap_and_argmax

def good_multipliers(S, p):
    """a in (Z/p)* with all a*v_i off the band {0,±1} (dist>=2). returns list."""
    out = []
    for a in range(1, p):
        if all(min((a*v) % p, p-((a*v) % p)) >= 2 for v in S):
            out.append(a)
    return out

def verify_count(S, p):
    g = good_multipliers(S, p)
    bound = (p-1) - 2*len(S)
    return len(g), bound

if __name__ == "__main__":
    print("PRIME-SHELL DODGE — counting bound and uniform M>=2/q\n" + "="*70)
    # (1) verify the counting bound #good >= (p-1)-2(n-1) on random configs
    import random; rng = random.Random(3); bad = 0
    for _ in range(3000):
        n = rng.randint(4, 18); p = nextprime(2*n-1)
        S = rng.sample(range(1, p), n-1)            # speeds < p so p ∤ v_i
        cnt, bnd = verify_count(S, p)
        if cnt < max(bnd, 0): bad += 1
        if cnt < 1: print("  !! no good multiplier", n, S, p)
    print(f"(1) counting bound #good >= (p-1)-2(n-1) held on 3000 random configs: violations={bad}")

    # (2) the table: optimal 2/(2n-1) vs provable 2/q, prime vs composite 2n-1
    print("\n(2)  n | 2n-1 (prime?) | q=nextprime(2n-1) | optimal 2/(2n-1) | provable 2/q | ratio")
    for n in range(4, 31):
        s = 2*n-1; q = nextprime(s)
        print(f"  {n:2d} | {s:3d} {'prime' if isprime(s) else 'comp ':5s} | q={q:3d} |"
              f" 2/{s}={float(Fr(2,s)):.5f} | 2/{q}={float(Fr(2,q)):.5f} | {float(Fr(s,q)):.3f}")

    # (3) exhaustive sanity: M(S) >= 2/q for all bounded configs (small n)
    print("\n(3) exhaustive: min M over primitive configs (speeds<q) vs provable 2/q and optimal 2/(2n-1)")
    for n in (4, 5, 6, 7):
        q = nextprime(2*n-1); worst = None
        for S in itertools.combinations(range(1, q), n-1):
            if reduce(gcd, S) != 1: continue
            M, _ = gap_and_argmax(S)
            if worst is None or M < worst[0]: worst = (M, S)
        print(f"  n={n}: min M={worst[0]}={float(worst[0]):.4f}  provable 2/q=2/{q}={float(Fr(2,q)):.4f}"
              f"  optimal 2/(2n-1)=2/{2*n-1}  (2n-1 {'prime' if isprime(2*n-1) else 'composite'})  argmin={worst[1]}")

    # (4) composite shell (n=14): does the 27-shell carry a good multiplier (the tower)?
    print("\n(4) n=14, 2n-1=27=3^3 (composite/ramified): prime shell gives 2/29; the OPTIMAL 2/27")
    n = 14; ext = tuple([1]+list(range(3, n+1)))   # AP-minus-2 style extremal probe
    g27 = good_multipliers(ext, 27); g29 = good_multipliers(ext, 29)
    print(f"   probe S={ext}")
    print(f"   good multipliers on shell 27 (->2/27): {len(g27)}  e.g. {g27[:5]}")
    print(f"   good multipliers on shell 29 (->2/29): {len(g29)}  e.g. {g29[:5]}")
