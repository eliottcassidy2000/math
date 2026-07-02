#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE BOUNDARY MINIMAX PROBE for HYP-3835: hunt the interpolating region between the census (height 1,
tight 1.16x) and the deep-cluster limit (height inf, loose 3.8x). If any moderate-height cluster config
dips below the pentagon 313/9702 (or worse, below 1/36), the fixed-point conjecture is WRONG.

Families probed (all 11-cores, band 1/14, target 1/36, pentagon benchmark 313/9702 = 0.032261):
  A. consecutive 7-cluster {N..N+6} over ALL primitive compact 4-cores (max<=13), N on a height ladder;
  B. pentagon-difference 7-clusters (drop patterns of {0..8}) over the pentagon's own compact parts;
  C. mixed splits: compact k + consecutive (11-k)-cluster for k=0..3 on the height ladder;
  D. barely-gap-free two-scale clusters {N..N+m} u {qN..} (ratios 2..5, below any Lambda=10^4 gap).

opus-2026-07-01-S32.
"""
import sys, itertools
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

BAND = Fr(1, 14); TARGET = Fr(1, 36); PENT = Fr(313, 9702)

def safe_arcs(v):
    return [((Fr(k) + BAND) / v, (Fr(k + 1) - BAND) / v) for k in range(v)]

def inter(A, B):
    r = []; i = j = 0
    while i < len(A) and j < len(B):
        lo = A[i][0] if A[i][0] > B[j][0] else B[j][0]
        hi = A[i][1] if A[i][1] < B[j][1] else B[j][1]
        if lo < hi: r.append((lo, hi))
        if A[i][1] < B[j][1]: i += 1
        else: j += 1
    return r

def Lmeas(S):
    a = safe_arcs(S[0])
    for v in S[1:]:
        a = inter(a, safe_arcs(v))
        if not a: return Fr(0)
    return sum(h - l for l, h in a)

print("=" * 100)
print(f" BOUNDARY MINIMAX PROBE.  benchmarks: 1/36 = {float(TARGET):.6f} ; pentagon = {float(PENT):.6f}")
print("=" * 100)

glob_min = None; glob_arg = None
def track(m, C, fam):
    global glob_min, glob_arg
    if m > 0 and (glob_min is None or m < glob_min):
        glob_min, glob_arg = m, (C, fam)

print("\n A. consecutive 7-cluster {N..N+6} x ALL primitive compact 4-cores (max<=13):")
fours = [C for C in itertools.combinations(range(1, 14), 4)]
HEIGHTS = [14, 16, 19, 23, 28, 35, 45, 60, 85, 120]
print(f"    compact 4-cores: {len(fours)}; heights: {HEIGHTS}")
worstA = []
for N in HEIGHTS:
    mn = None; arg = None
    clu = set(range(N, N + 7))
    for C4 in fours:
        if set(C4) & clu: continue
        C = tuple(sorted(set(C4) | clu))
        if reduce(gcd, C) != 1: continue
        m = Lmeas(C)
        if m > 0 and (mn is None or m < mn): mn, arg = m, C4
        track(m, C, "A")
    worstA.append((N, mn, arg))
    print(f"    N={N:>4}: min meas = {float(mn):.6f} at C_low={arg}   (/pentagon: {float(mn/PENT):.3f}x)")

print("\n B. pentagon-difference clusters: 7-subsets of {N+c: c in {0..8}} with AP-with-drops patterns,")
print("    over compact parts drawn from the pentagon's small elements:")
drop_patterns = [tuple(sorted(set(range(9)) - set(d))) for d in itertools.combinations(range(1, 8), 2)]
lows = [(1,2,3,4), (1,2,3,5), (1,2,4,5), (2,3,4,5), (1,3,4,5)]
mnB = None; argB = None
for N in (14, 20, 30, 50, 80):
    for pat in drop_patterns:
        clu = set(N + c for c in pat)
        for low in lows:
            if set(low) & clu: continue
            C = tuple(sorted(set(low) | clu))
            if len(C) != 11 or reduce(gcd, C) != 1: continue
            m = Lmeas(C)
            if m > 0 and (mnB is None or m < mnB): mnB, argB = m, (low, N, pat)
            track(m, C, "B")
print(f"    min over {len(drop_patterns)}x{len(lows)}x5 configs = {float(mnB):.6f} at low={argB[0]}, N={argB[1]}, pattern={argB[2]}")
print(f"    (/pentagon: {float(mnB/PENT):.3f}x)")

print("\n C. mixed splits: compact k + consecutive (11-k)-cluster:")
for k, lows_k in [(0, [()]), (1, [(1,), (2,), (3,)]), (2, [(1,2),(1,3),(2,3),(1,5)]), (3, [(1,2,3),(1,2,5),(2,3,5),(1,3,4)])]:
    j = 11 - k
    mn = None; arg = None
    for N in (14, 20, 30, 50, 80, 120):
        clu = set(range(N, N + j))
        for low in lows_k:
            if set(low) & clu: continue
            C = tuple(sorted(set(low) | clu))
            if len(C) != 11 or reduce(gcd, C) != 1: continue
            m = Lmeas(C)
            if m > 0 and (mn is None or m < mn): mn, arg = m, (low, N)
            track(m, C, "C")
    print(f"    k={k} compact + {j}-cluster: min = {float(mn):.6f} at low={arg[0]}, N={arg[1]}   (/pentagon: {float(mn/PENT):.3f}x)")

print("\n D. barely-gap-free two-scale clusters (ratios 2..5, no Lambda-gap):")
mnD = None; argD = None
for N in (14, 20, 30, 50):
    for q in (2, 3, 4, 5):
        for split in ((6, 5), (7, 4), (5, 6), (4, 7)):
            a, b = split
            lowc = tuple(range(N, N + a))
            high = tuple(range(q * N, q * N + b))
            C = tuple(sorted(set(lowc) | set(high)))
            if len(C) != 11 or reduce(gcd, C) != 1: continue
            m = Lmeas(C)
            if m > 0 and (mnD is None or m < mnD): mnD, argD = m, (N, q, split)
            track(m, C, "D")
print(f"    min over N x q x splits = {float(mnD):.6f} at N={argD[0]}, ratio q={argD[1]}, split={argD[2]}")
print(f"    (/pentagon: {float(mnD/PENT):.3f}x)")

print("\n" + "=" * 100)
print(f" GLOBAL MIN over all boundary families: {float(glob_min):.6f} = {glob_min}")
print(f"   at C = {glob_arg[0]} (family {glob_arg[1]})")
print(f"   vs pentagon 313/9702 = {float(PENT):.6f}: ratio {float(glob_min/PENT):.3f}x")
print(f"   vs 1/36: ratio {float(glob_min/TARGET):.3f}x ; below pentagon? {glob_min < PENT} ; below 1/36? {glob_min < TARGET}")
print(" If ratio/pentagon >= 1: the fixed-point conjecture SURVIVES the boundary probe (census stays binding).")
print("DONE.")
