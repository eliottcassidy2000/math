#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_threadC_unify_points_kpswf4.py  (kind-pasteur 2026-06-21, THREAD C)

Are the pq midpoints of all q arcs together an EQUALLY SPACED set on the length-7 circle?

Row 0 sample points (7pv values, mod 7):
   x_{k,t} = A_k + (2t+1)/(2q),   A_k = 7 frac(pk/q),  k=0..q-1, t=0..p-1.
Equivalently  x = 7 frac(pk/q) + (2t+1)/(2q).

CLAIM to test: { x_{k,t} mod 7 } = { (2m+1)/(2q) : m = 0 .. pq-1 } mod 7  (i.e. the pq points
equally spaced with gap 1/q, offset 1/(2q)).  If TRUE, row 0 is EXACTLY a histogram of pq
equally spaced points (gap 1/q) into 7 unit bins -- a pure 1D three-distance / equidistribution
object, and the sharp bound follows from counting points per unit bin.
"""
from fractions import Fraction as Fr
from math import gcd

P = 7

def row0_points(p, q):
    pts = []
    for k in range(q):
        Ak = (7 * Fr(p * k, q)) % 7
        for t in range(p):
            x = (Ak + Fr(2 * t + 1, 2 * q)) % 7
            pts.append(x)
    return sorted(pts)

def equally_spaced_ref(p, q):
    # pq points gap 1/q offset 1/(2q): (2m+1)/(2q), m=0..pq-1, mod 7
    return sorted(((Fr(2 * m + 1, 2 * q)) % 7 for m in range(p * q)))

def main():
    print("THREAD C: are the pq row-0 sample points equally spaced (gap 1/q)?")
    print("=" * 70)
    ok = True
    examples = []
    for q in range(1, 30):
        for p in range(q + 1, int(Fr(43, 20) * q) + 1):
            if gcd(p, q) != 1 or not (Fr(1) < Fr(p, q) <= Fr(43, 20)):
                continue
            a = row0_points(p, q)
            b = equally_spaced_ref(p, q)
            if a != b:
                ok = False
                if len(examples) < 5:
                    examples.append((p, q, a, b))
    print("pq row-0 points == equally spaced gap 1/q :", "YES" if ok else "NO")
    if not ok:
        for p, q, a, b in examples:
            print(f"  MISMATCH p/q={p}/{q}")
            print(f"    points : {[float(x) for x in a]}")
            print(f"    ref    : {[float(x) for x in b]}")
    else:
        print("\n=> ROW 0 = histogram of pq equally spaced points (gap 1/q, offset 1/(2q))")
        print("   into 7 unit bins on the length-7 circle.")
        print("   r_j = #{ m : (2m+1)/(2q) mod 7 in [j, j+1) }, j=0..6.")
        # show that this is just: pq points, gap 1/q, into 7 bins of width 1
        # => each bin gets either floor(7pq/7 /? ) ... compute the per-bin counts directly
        print("\n   Per-bin counts r_j for sample ratios (each in {floor, ceil} of pq/7):")
        for (p, q) in [(3,2),(2,1),(4,3),(5,3),(5,4),(9,5),(11,10),(17,10),(43,20)]:
            if gcd(p,q)!=1: continue
            r = [0]*P
            for m in range(p*q):
                x = (Fr(2*m+1,2*q)) % 7
                r[int(x)] += 1
            lo = (p*q)//7; hi = -(-(p*q)//7)
            allinrange = all(x in (lo,hi) for x in r)
            print(f"   p/q={p}/{q:<3d} pq={p*q:<4d} pq/7={float(Fr(p*q,7)):.3f} r={r}  "
                  f"each in {{{lo},{hi}}}:{allinrange}")

if __name__ == "__main__":
    main()
