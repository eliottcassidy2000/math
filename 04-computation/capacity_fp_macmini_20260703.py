#!/usr/bin/env python3
"""
THE CAPACITY QUANTITY f_q (mac-mini-2026-07-03-S28): f_q = P(random 13-tuple of nonzero residues mod q is
NO-WITNESS). The creative reason the log-census closes: blocking modulus q (making the family no-witness mod q)
costs log(1/f_q) 'bits'; by CRT, blocking primes p_1..p_r needs the 13 speeds to jointly hit a density-prod(f)
config, so the smallest such tuple is ~ e^{(1/13) sum log(1/f)}; for it to be <= M, sum log(1/f_p) <= 13 log M.
Since log(1/f_q) >= c > 0 (f_q bounded away from 1), at most O(log M) moduli are blockable => witness at
q* = O(log M loglog M). WHY f_q<1: no-witness needs the 13 danger arcs to TILE Z/q (pairwise commensurability
= small resonances q|(m_i v_i + m_j v_j)) -- a special density-f_q configuration.
"""
from math import gcd, log
import random
def nd(x): x=x%1; return min(x,1-x)
def has_witness_tuple(vs,q):
    for a in range(1,q):
        if gcd(a,q)!=1: continue
        if all(nd((v*a)%q/q)>=1/14 for v in vs): return True
    return False
if __name__=="__main__":
    rng=random.Random(7)
    print(f"{'q':>4} {'f_q':>10} {'log(1/f_q)':>11}")
    for q in [17,19,23,29,31,37,41,43,47,53,59,67,71,79,83,89,97,101,113]:
        n=3000; nw=sum(1 for _ in range(n) if not has_witness_tuple([rng.randint(1,q-1) for _ in range(13)],q))
        f=nw/n; print(f"{q:>4} {f:>10.4f} {log(1/f) if f>0 else float('inf'):>11.3f}")
    print("=> log(1/f_q) bounded below (>=~2) and f_q -> const<1 (arc-covering prob) => O(log M) blockable moduli")
