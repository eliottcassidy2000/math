#!/usr/bin/env python3
"""
lrc_nowhere_zero_flow_s521.py   claudebox-2026-06-01-S521

The nowhere-zero-flow view of LRC (reflection:
07-reflections/lrc-nowhere-zero-flow-attack-s521.md).

(1) CROSSING-FLOW DIVERGENCE LAW: on the directed n-cycle of sector boundaries,
    occupancy change = net crossing flow (conservation of runners). Verified as
    o_k(t2)-o_k(t1) = (crossings into S_k) - (crossings out of S_k). LRC(strict)
    <=> div=0 at the two observer sectors (o_0=o_{n-1}=0).
(2) OBSERVER AS FLOW-SOURCE: in the marked tournament indeg(observer)=N(t)=
    #runners within 1/n; observer is a nowhere-zero-flow source iff N(t)=0 iff lonely.
(3) RELATION LATTICE: nowhere-zero relations sum a_i v_i = 0 (the flow lattice of
    the speed system); the minimal one is the shortest nowhere-zero flow.
"""
from fractions import Fraction as F
from math import gcd
from itertools import product

def fr(x): return x % 1
def dist(x):
    x = x % 1; return min(x, 1 - x)
def sector(p, n): return int(fr(p) * n)
def occ(sp, t, n):
    o = [0]*n
    for v in sp: o[sector(F(v)*t, n)] += 1
    return o
def crossings_in_interval(sp, t1, t2, k, n):
    """net runners crossing boundary b_k=k/n (clockwise) during (t1,t2]."""
    c = 0
    for v in sp:
        # v s = k/n + j, s in (t1,t2]  ->  j in (v t1 - k/n, v t2 - k/n]
        import math
        lo = v*float(t1) - k/n; hi = v*float(t2) - k/n
        c += math.floor(hi + 1e-9) - math.floor(lo + 1e-9)
    return c

def main():
    sp = [1, 2, 4, 7]; n = 5
    t1, t2 = F(1, 13), F(3, 11)
    o1 = occ(sp, t1, n); o2 = occ(sp, t2, n)
    print(f"(1) Divergence law (sp={sp}): occupancy change = net crossing flow")
    ok = True
    for k in range(n):
        do = o2[k] - o1[k]
        net = crossings_in_interval(sp, t1, t2, k, n) - crossings_in_interval(sp, t1, t2, (k+1) % n, n)
        if do != net: ok = False
        print(f"   S_{k}: d(occ)={do:+d}  net flow(in b_{k} - out b_{(k+1)%n})={net:+d}  {'OK' if do==net else 'MISMATCH'}")
    print(f"   conservation law holds: {ok}")

    print(f"\n(2) observer is a flow-SOURCE  <=>  indeg(observer)=N(t)=0  <=>  lonely:")
    for t in [F(1, 5), F(3, 11), F(1, 13)]:
        N = sum(1 for v in sp if dist(F(v)*t) < F(1, n))
        print(f"   t={t}: N(t)=indeg(observer)={N}  -> observer source? {N == 0}")

    print(f"\n(3) Relation lattice (nowhere-zero flow of the speed system):")
    for s in [(1, 2, 3), (1, 2, 4), (2, 3, 5, 7)]:
        best = None
        for a in product(range(-4, 5), repeat=len(s)):
            if all(x != 0 for x in a) and sum(x*v for x, v in zip(a, s)) == 0:
                L1 = sum(abs(x) for x in a)
                if best is None or L1 < best[0]: best = (L1, a)
        print(f"   {s}: minimal nowhere-zero relation {best[1]}  (L1={best[0]})")

if __name__ == "__main__":
    main()
