#!/usr/bin/env python3
"""
mac-mini-2026-07-02-S5 -- HYP-3862: the d=3 DEVIATION TABLE (queue item).
For each of the 10 primitive d=3 simplex patterns (S2 census), the EXACT
max/min over the shift 2-torus of the triple overlap
  ov(t2,t3) = |{x : ||P x|| < r, ||Q x - t2|| < r, ||S x - t3|| < r}|,
r = 1/14.  Piecewise bilinear in (t2,t3) => extrema at breakpoint vertices
(endpoint-coincidence rationals).  Mean = (2r)^3 (Fubini).  These rows are
the THM-602(C) partial-cycle accounting data for depth-3 resonances.
"""
from fractions import Fraction as F
from math import gcd

r = F(1, 14)
PATTERNS = [(1,1,1),(1,1,2),(1,1,3),(1,1,4),(1,1,5),(1,2,2),(1,2,3),(1,2,4),(1,3,3),(2,2,3)]

def arcs(v, shift):
    out = []
    for k in range(v):
        a, b = (F(k) + shift - r) / v, (F(k) + shift + r) / v
        out.append((a, b))
    return out

def inter_len(iv1, iv2):
    tot = F(0)
    for a1, b1 in iv1:
        for a2, b2 in iv2:
            for sh in (-1, 0, 1):
                lo, hi = max(a1, a2 + sh), min(b1, b2 + sh)
                if hi > lo: tot += hi - lo
    return tot

def ov3(P, Q, S, t2, t3):
    A = arcs(P, F(0)); B = arcs(Q, t2); C = arcs(S, t3)
    # intersect A with B, then with C
    AB = []
    for a1, b1 in A:
        for a2, b2 in B:
            for sh in (-1, 0, 1):
                lo, hi = max(a1, a2 + sh), min(b1, b2 + sh)
                if hi > lo: AB.append((lo, hi))
    return inter_len(AB, C)

def bps(P, other, v):
    """breakpoints of the shift of comb v against comb P (both endpoints coincide)."""
    out = set()
    for k in range(P):
        for m in range(v):
            for s1 in (-1, 1):
                for s2 in (-1, 1):
                    for j in (-1, 0, 1):
                        th = (F(k) + s1 * r) * v / P - F(m) - s2 * r + j
                        th -= int(th)
                        if th < 0: th += 1
                        out.add(th)
    return sorted(out)

mean = (2 * r) ** 3
print(f"d=3 deviation table, r=1/14, mean=(2r)^3={mean}={float(mean):.6f}")
print(f"{'pattern':>10} {'exact min':>12} {'exact max':>12} {'max-mean':>12} {'mean-min':>12}")
for (P, Q, S) in PATTERNS:
    # vertex grid: breakpoints of t2 vs P and vs S-comb interplay -- take the union
    # of pairwise coincidence breakpoints for each shift variable (sufficient: the
    # overlap is piecewise bilinear with cell walls at pairwise coincidences)
    T2 = sorted(set(bps(P, Q, Q)) | set())
    T3 = sorted(set(bps(P, S, S)) | set(bps(Q, S, S)))
    # subsample if huge (patterns are tiny so fine); evaluate at vertices
    mn, mx = None, None
    for t2 in T2:
        for t3 in T3:
            v = ov3(P, Q, S, t2, t3)
            if mn is None or v < mn: mn = v
            if mx is None or v > mx: mx = v
    print(f"  ({P},{Q},{S})  {str(mn):>12} {str(mx):>12} {str(mx-mean):>12} {str(mean-mn):>12}")
print("\nAll min = 0 expected (sum|m| <= 7: danger-capable, THM-601(iii)); the max-mean")
print("column is the partial-cycle cost THM-602(C) charges per frozen depth-3 resonance.")
