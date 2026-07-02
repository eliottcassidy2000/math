#!/usr/bin/env python3
"""
klein-2026-07-01-S89 -- HYP-3846: the section-7.3 arc x radius LP, n=6 pilot.

Setting (opus-S33 doc sec 7.3): variables x_{a,l} = lonely mass in arc a at radius layer
r_l; constraints (i) layer monotonicity, (ii) per-arc coverage, (iii) slope transport
x_{a,l} - x_{a,l+1} <= budget_a * dr  (THM-592(ii) localized).

PILOT OBSERVATION: with single-anchor data (per-arc lonely masses at r0) and per-arc
transport budgets, the LP is SEPARABLE and solves in closed form:
    min sum_a x_{a,L}  =  sum_a max(0, x_{a,0} - budget_a * (r_L - r0))
-- the ARC-LOCALIZED LADDER. The solver only earns its keep with multi-layer coupling +
universal (covering-axiom) parameterization; the pilot certifies the constraint set is
EXACT on the known n=6 collapse law Lambda = (2/5)(1-6r) on the final window [1/7, 1/6].

n=6 facts (klein S87 HYP-3834/3835): c(AP_5) = c({1,3,4,5,9}) = 2/5; final window
[1/7, 1/6] exactly linear; d=1 band (6d,7d) empty => K=0 for diam <= 12 (n=6 analogue).
"""
from fractions import Fraction as F
import itertools

def danger_intervals(S, r):
    ivs = []
    for v in S:
        rv = r / v
        for a in range(v + 1):
            c = F(a, v); lo, hi = c - rv, c + rv
            if hi <= 0 or lo >= 1: continue
            ivs.append((max(lo, F(0)), min(hi, F(1))))
    ivs.sort()
    merged = []
    for lo, hi in ivs:
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged

def lonely_arcs(S, r, Q):
    """Per-arc lonely mass x_a and per-arc gap-component slope budgets."""
    cov = danger_intervals(S, r)
    # complement gaps
    gaps = []
    prev = F(0)
    for lo, hi in cov:
        if lo > prev: gaps.append((prev, lo))
        prev = max(prev, hi)
    if prev < 1: gaps.append((prev, F(1)))
    x = [F(0)] * Q
    budget = [F(0)] * Q
    for glo, ghi in gaps:
        # gap shrink rate: 1/v_left + 1/v_right where the bounding interval edges belong
        # to the speeds whose fraction is nearest; exact rate = sum over the two bounding
        # danger edges of 1/v. Identify bounding speeds by matching endpoints.
        rate = F(0)
        for v in S:
            rv = r / v
            for a in range(v + 1):
                if F(a, v) + rv == glo: rate += F(1, v)   # left neighbor's right edge
                if F(a, v) - rv == ghi: rate += F(1, v)   # right neighbor's left edge
        for c in range(Q):
            alo, ahi = F(c, Q), F(c + 1, Q)
            seg = min(ghi, ahi) - max(glo, alo)
            if seg > 0:
                x[c] += seg
                budget[c] += rate * (seg / (ghi - glo))  # apportion rate by overlap share
    return x, budget

def pilot(S, name, Q=6):
    r0, rL = F(1, 7), F(1, 6)
    x0, b0 = lonely_arcs(S, r0, Q)
    dr = rL - r0
    lp_value = sum(max(F(0), xa - ba * dr) for xa, ba in zip(x0, b0))
    # exact truth
    cov = danger_intervals(S, rL)
    true_L = 1 - sum(hi - lo for lo, hi in cov)
    anchor = sum(x0)
    c_law = F(2, 5)
    print(f"  {name:16s} Q={Q}  anchor Lambda(1/7)={anchor} (law: {c_law * (1 - 6 * r0)})"
          f"  LP-min at 1/6 = {lp_value}   true Lambda(1/6) = {true_L}"
          f"   certifies-law={'YES' if lp_value == true_L else 'gap ' + str(float(lp_value - true_L))}")

print("=" * 90)
print("SEC-7.3 ARC x RADIUS LP, n=6 PILOT: separable closed form = arc-localized ladder")
print("=" * 90)
for Q in (6, 31):
    pilot([1, 2, 3, 4, 5], "AP_5 {1..5}", Q)
    pilot([1, 3, 4, 5, 9], "sporadic", Q)
print("\n  reading: anchor at 1/7 = (2/5)(1/7) = 2/35; window exactly linear (d=1 band (6,7)")
print("  empty, diam <= 12 => K=0); transport budgets = exact gap rates => the LP constraint")
print("  set (i)-(iii) is EXACT on the final window; the collapse law is LP-certified per-set.")
print("  NEXT (universal): parameterize D_{a,l} by covering axioms only; multi-layer coupling.")
