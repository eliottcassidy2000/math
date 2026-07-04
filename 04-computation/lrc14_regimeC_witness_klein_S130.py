#!/usr/bin/env python3
"""
klein-2026-07-04-S130 - REGIME C (near-equal-small 7-far wall): lonely? effective bound? witness?

The c=7 trichotomy (kps-S24): the seven-danger-arc wall in LRC(14) breaks 3 ways.
  A clustered-huge (w1>=7392): CLOSED (frozen arcs, sum-combo gap-deficit).
  B spread (some consecutive D*L>=2): CLOSED (first pristine pair-floor).
  C near-equal-small (w1<7392, tight, e.g. 23..29): OPEN. Both citations fail -- the allowed
    window ~1/B is too coarse to see the arc correlation as measure.

This session: DIRECTLY probe regime C. A regime-C family = 6 bounded runners + 7 near-consecutive
"far" runners {W,...,W+6} (D=1 tight, small W). Compute exact M(family), confirm lonely (M>=1/14),
find (a) the MIN M over worst-case bounded runners per W, (b) whether the hard cases are BOUNDED in
W (opus's Q: is the tight-small corner a bounded finite residual?), (c) the WITNESS structure
(three-distance? bounded combo?) that makes them lonely -- a hint at the proof technique.

Exact (Fractions).
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def cdist_q(a, q):
    r = a % q
    return min(r, q - r)

def Mval_arg(S, Qcap):
    """Exact M and one optimal (a,Q)."""
    best = F(0); arg = (0, 1)
    for Q in range(2, Qcap + 1):
        for a in range(1, Q // 2 + 1):
            if gcd(a, Q) != 1: continue
            m = min(F(cdist_q(v * a, Q), Q) for v in S)
            if m > best: best = m; arg = (a, Q)
    return best, arg

def binders(S, a, Q, M):
    return [v for v in S if F(cdist_q(v * a, Q), Q) == M]

ONE14 = F(1, 14)
print(f"1/14 = {float(ONE14):.6f}")
print("=" * 88)
print("opus's concrete regime-C example {1..6} u {23..29}:")
S = list(range(1, 7)) + list(range(23, 30))
M, (a, Q) = Mval_arg(S, 2 * max(S) + 40)
print(f"  S={S}")
print(f"  M={M} (~{float(M):.6f})  {'>= 1/14 LONELY' if M>=ONE14 else 'BELOW 1/14 !!'}  witness t*={a}/{Q}  binders={binders(S,a,Q,M)}")
print("=" * 88)

# Sweep: for each far-block start W, worst-case 6 bounded runners from a small pool.
# Bounded pool = {1,...,Bpool}\{far block}. Far block = {W,...,W+6} (7 consecutive).
# Report min M over choices of 6 bounded runners, and the witness for the extremal.
Bpool = 12
print(f"Sweep: far block {{W..W+6}} + worst-case 6 bounded runners from 1..{Bpool}.")
print(f"{'W':>4} {'min M':>12} {'~':>9} {'lonely?':>8} {'argmin bounded':>22} {'t*':>10}")
overall_min = None; overall_at = None
for W in range(8, 60):
    far = list(range(W, W + 7))
    pool = [b for b in range(1, Bpool + 1) if b not in far]
    if len(pool) < 6: continue
    wmin = None; wat = None; warg = None
    for bnd in combinations(pool, 6):
        S = sorted(list(bnd) + far)
        if len(set(S)) != 13: continue
        M, arg = Mval_arg(S, 2 * (W + 6) + 30)
        if wmin is None or M < wmin:
            wmin, wat, warg = M, bnd, arg
    lon = "yes" if wmin >= ONE14 else "NO!!"
    print(f"{W:>4} {str(wmin):>12} {float(wmin):>9.6f} {lon:>8} {str(list(wat)):>22} {warg[0]}/{warg[1]}")
    if overall_min is None or wmin < overall_min:
        overall_min = wmin; overall_at = (W, wat)

print("=" * 88)
print(f"OVERALL worst regime-C family (W in 8..59, bounded from 1..{Bpool}):")
print(f"  min M = {overall_min} (~{float(overall_min):.6f}) at W={overall_at[0]}, bounded={list(overall_at[1])}")
print(f"  all lonely (>= 1/14)? {overall_min >= ONE14}")
print()
print("READING: if all regime-C families are comfortably >= 1/14 and the min M is BOUNDED AWAY")
print("from 1/14 (not approaching it as W grows), then regime C is TRUE (lonely) and BOUNDED --")
print("a finite residual, closeable by a finite check. The witness t* structure hints the technique.")
print("DONE")
