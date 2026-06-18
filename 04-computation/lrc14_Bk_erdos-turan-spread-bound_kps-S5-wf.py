#!/usr/bin/env python3
"""
lrc14_Bk_erdos-turan-spread-bound_kps-S5-wf.py   (kps-2026-06-18-S5-wf)

ANGLE = erdos-turan-spread-bound.  GOAL: prove (or explicitly reduce to a finite check)
the SPREAD BOUND B(k):   spread(E) > B(k)  ==>  mu(E) >= mu_min^bdd(k)
so that the residual uniform-floor lemma (THM-527-G / THM-528-E / OPEN-Q-108) becomes a
FINITE / COMPACT problem.

mu(E) = meas{ x in [0,1) : the points {frac(e_i x): e_i in E} have circular maxgap > 2/7 }.
E is a set of nonneg integers with 0 in E (co-offsets), k = |E|, spread = max E.

PART 0 (this file, foundation):  a RIGOROUS exact mu computation (the prompt's tool is NOT
rigorous: it misses the gap=2/7 crossings strictly inside an order-cell).  We add the
gap=2/7 breakpoints and verify against the canon exact values mu_3..mu_7 for consecutive E.

Subsequent PARTs (added below): structure of mu vs spread, the Erdos-Turan / Weyl deviation
bound, and the explicit B(k) calibration.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

TWO7 = F(2,7)

# ---------------------------------------------------------------------------
# RIGOROUS exact mu.
#
# On an open cell (a,b) between consecutive "order breakpoints" the cyclic ORDER of the
# points {e_i x mod 1} is fixed and each point moves AFFINELY in x: p_i(x) = e_i x - floor(e_i a_mid).
# So every gap g(x) between cyclically-consecutive points is AFFINE in x on the cell.
# maxgap(x) = max of finitely many affine functions => piecewise-linear, convex-ish but the
# MAX of affines is convex; {maxgap > 2/7} on the cell is the complement of a closed interval
# (where the convex function <= 2/7), i.e. a union of at most two sub-intervals.  We compute it
# EXACTLY by: enumerating the affine gap functions on the cell, forming maxgap as their upper
# envelope, and finding exactly where the envelope crosses 2/7.
#
# Implementation: within a cell, pick the order at the midpoint; the cyclic neighbor structure
# is constant; each gap is affine g(x)=alpha*x+beta with alpha,beta in Q (exact).  maxgap is the
# max; we find the measure of {x in (a,b): max_g g(x) > 2/7} exactly.
# ---------------------------------------------------------------------------

def order_breakpoints(E):
    """x where two points collide: (e_i - e_j) x in Z."""
    bp = {F(0), F(1)}
    Es = sorted(set(E))
    for i in range(len(Es)):
        for j in range(i+1, len(Es)):
            d = Es[j] - Es[i]
            for m in range(0, d+1):
                bp.add(F(m, d))
    return sorted(b for b in bp if F(0) <= b <= F(1))

def affine_gaps_on_cell(E, a, b):
    """On open cell (a,b): return list of (alpha,beta) so each cyclic gap is alpha*x+beta.
    We fix the integer 'level' n_i = floor(e_i * mid) for each point so frac(e_i x)=e_i x - n_i
    holds throughout the cell (no wrap inside the cell since order is constant)."""
    mid = (a + b) / 2
    Es = sorted(set(E))
    # fractional value at mid, and the integer level
    pts = []
    for e in Es:
        val = e * mid
        n = val - (val % 1)        # floor as Fraction
        # frac(e x) = e x - n on this cell.  representative position at mid:
        pts.append((e * mid - n, e, n))
    # sort by position at mid to get cyclic order
    pts.sort(key=lambda t: t[0])
    m = len(pts)
    gaps = []
    for i in range(m):
        (pos_i, e_i, n_i) = pts[i]
        (pos_j, e_j, n_j) = pts[(i+1) % m]
        if i < m-1:
            # gap = (e_j x - n_j) - (e_i x - n_i)
            alpha = e_j - e_i
            beta  = -(n_j) + (n_i)
        else:
            # wrap gap = (e_{0} x - n_0 + 1) - (e_i x - n_i)
            (pos0, e0, n0) = pts[0]
            alpha = e0 - e_i
            beta  = -(n0) + (n_i) + 1
        gaps.append((alpha, beta))
    return gaps

def measure_envelope_gt(gaps, a, b, c=TWO7):
    """Exact measure of {x in (a,b): max_g (alpha*x+beta) > c}.
    max of affines is convex PL; {>c} is (a,b) minus the closed interval where envelope<=c.
    We compute the set where envelope>c by: for the upper envelope U(x)=max_g g(x),
    {U>c} = union over g of {g>c}, intersected appropriately?  NO: U>c  <=>  exists g with g>c.
    So {U>c} = UNION_g {x: alpha*x+beta>c}.  That's exact and simple (no envelope needed)."""
    # union of half-lines intersected with (a,b)
    intervals = []
    for (alpha, beta) in gaps:
        # alpha x + beta > c
        if alpha == 0:
            if beta > c:
                intervals.append((a, b))
            # else empty
        else:
            xstar = (c - beta) / alpha
            if alpha > 0:
                lo = max(a, xstar); hi = b
            else:
                lo = a; hi = min(b, xstar)
            if lo < hi:
                intervals.append((lo, hi))
    # union measure
    if not intervals:
        return F(0)
    intervals.sort()
    tot = F(0); cur_lo, cur_hi = intervals[0]
    for lo, hi in intervals[1:]:
        if lo <= cur_hi:
            cur_hi = max(cur_hi, hi)
        else:
            tot += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
    tot += cur_hi - cur_lo
    return tot

def mu_exact(E, c=TWO7):
    """Rigorous exact mu(E)."""
    E = sorted(set(E))
    if len(E) == 1:
        return F(1)  # single point: maxgap = 1 > 2/7 always
    bps = order_breakpoints(E)
    tot = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b:
            continue
        gaps = affine_gaps_on_cell(E, a, b)
        tot += measure_envelope_gt(gaps, a, b, c)
    return tot

# ---------------------------------------------------------------------------
# VERIFICATION against canon exact consecutive values (THM-527-C/E):
#   mu_3=1, mu_4=19/21, mu_5=9/14, mu_6=4/7, mu_7=13/35, mu_13=829/4620
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("="*80)
    print("PART 0: RIGOROUS exact mu — verify vs canon consecutive values")
    print("="*80)
    canon = {3:F(1), 4:F(19,21), 5:F(9,14), 6:F(4,7), 7:F(13,35), 13:F(829,4620)}
    allok = True
    for k in range(3, 14):
        E = list(range(k))
        m = mu_exact(E)
        tag = ""
        if k in canon:
            ok = (m == canon[k])
            allok &= ok
            tag = f"  canon={canon[k]}  {'OK' if ok else 'MISMATCH!!'}"
        print(f"  k={k:2d}  consecutive  mu = {m} = {float(m):.6f}{tag}")
    print(f"\n  ALL canon checks: {'PASS' if allok else 'FAIL'}")
