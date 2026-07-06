#!/usr/bin/env python3
"""
klein-2026-07-05-S142 (HYP-4187) - DRAG(12) VIA MULTI-TOP SPLITS: coverage measurement.

For a middle-strip 12-family (v_max/v_min > 11.5, v_max/2ndmax <= 11.5), search all
(base, tops) suffix-splits l = 1..6 for one satisfying the FORMALIZED multi-window
criterion (margin_of_window_of_fees at rho = 2/25, base of 12-l runners cited at
beta = 1/(13-l) > 2/25):

    delta = (beta - rho)/B',   B' = max(base)
    per-top fee = 2*rho*(2*delta) + 6*rho/v   [density, teethR_mass]
              or  2*rho/v  if  2*delta*v <= 1 - 2*rho   [sparse, teethR_mass_far]
    criterion:  sum fees < 2*delta

(Fee formulas EXACTLY mirror LRCTeethR; a split that passes here is a Lean-ready
instance.) Measure coverage over strip families; characterize the residual.
"""
from fractions import Fraction as F
from math import gcd
import random

rho = F(2, 25)

def split_closes(V):
    """V sorted ascending, 12 entries. Try l = 1..6 suffix splits."""
    for l in range(1, 7):
        base, tops = V[:12-l], V[12-l:]
        beta = F(1, 13 - l)
        if beta <= rho: continue
        B = base[-1]
        delta = (beta - rho) / B
        total = F(0)
        for v in tops:
            if 2 * delta * v <= 1 - 2 * rho:
                fee = 2 * rho / v                       # sparse
            else:
                fee = 2 * rho * (2 * delta) + 6 * rho / v   # density
            total += fee
        if total < 2 * delta:
            return l
    return None

random.seed(17)
cands = []
for lo_top in range(1, 7):
    nb = 12 - lo_top
    base = list(range(1, nb + 1))
    for T in range(int(11.5 * max(base)) + 2, 30 * max(base), max(base)):
        cands.append(sorted(base + [T + 2*j for j in range(lo_top)]))
for _ in range(600):
    k = random.randint(1, 6)
    base = sorted(random.sample(range(1, 28), 12 - k))
    T = random.randint(int(11.5 * min(base)) + 1, 40 * min(base))
    tops = sorted(random.sample(range(T, T + max(3, 3 * T)), k))
    V = sorted(set(base + tops))
    if len(V) == 12: cands.append(V)

tested = closed = 0
residual = []
for V in cands:
    if len(set(V)) != 12 or gcd(*V) != 1: continue
    if not (V[-1] * 2 > 23 * V[0] and V[-1] * 2 <= 23 * V[-2]): continue
    tested += 1
    l = split_closes(V)
    if l is not None:
        closed += 1
    else:
        residual.append(V)

print(f"strip families tested: {tested}")
print(f"closed by a multi-window split: {closed} ({100*closed/max(tested,1):.1f}%)")
print(f"residual (no split closes): {len(residual)}")
for V in residual[:12]:
    print("  ", V, " ratios: full", round(V[-1]/V[0],1), "top", round(V[-1]/V[-2],2))
