#!/usr/bin/env python3
"""uniform_floor_c0_macmini_S125.py -- mac-mini-2026-07-16-S125.
THE UNIFORM FLOOR c0 (THM-527's remaining crux) AS A FINITE EXACT COMPUTATION.

Realization: the bounded-spread reduction (THM-527, extremizer bounded-spread; large
spread handled by THM-924's glue trichotomy) leaves INTEGER co-offset clusters
E = {0 = e1 < ... < ek <= S0} -- a FINITE set. Hence
    c0 = min over (prefix core P, cluster E) of rho*(P, E),
each an EXACT RATIONAL: rho*(P,E) = meas{x in G_P(1/14): maxgap{frac(e_i x)} > 2/7},
computed by exact breakpoint sweep (breakpoints: a/(14u) shifts for G_P, and
(a +- 2/7)/(e_i - e_j) for the gap events; denominators 7*diffs and 14u).

Sweep: P = {1..j}, j = 1..5; clusters k = 2..4 with S0 = 12 (complete) and k = 5 with
S0 = 8 (complete at that spread). Output: the exact minimum, its argmin, and the table
extremes. If min > 0: c0 > 0 is PROVED-by-finite-computation for the swept scope, with
the completion to (k <= 6, S0 = canonical bound) a mechanical rerun.
"""
import sys
from fractions import Fraction as Fr
from math import gcd
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)

LAM = Fr(1, 14)
TWO7 = Fr(2, 7)

def rho_star_exact(P, E):
    """exact measure of {x in [0,1): all ||u x|| >= 1/14 (u in P), maxgap({e_i x}) > 2/7}."""
    # breakpoints
    bps = {Fr(0), Fr(1)}
    for u in P:
        for a in range(14 * u + 1):
            bps.add(Fr(a, 14 * u))
    diffs = sorted({abs(a - b) for a, b in combinations(E, 2) if a != b})
    for d in diffs:
        for a in range(7 * d + 1):
            bps.add(Fr(a, 7 * d))
        for a in range(d + 1):
            bps.add(Fr(a, d))
    bps = sorted(bps)
    tot = Fr(0)
    for i in range(len(bps) - 1):
        x = (bps[i] + bps[i + 1]) / 2
        # G_P
        ok = all(min(u * x % 1, 1 - u * x % 1) >= LAM for u in P)
        if ok:
            ph = sorted(e * x % 1 for e in E)
            gaps = [ph[j + 1] - ph[j] for j in range(len(ph) - 1)] + [1 + ph[0] - ph[-1]]
            if max(gaps) > TWO7:
                tot += bps[i + 1] - bps[i]
    return tot

print("THE c0 TABLE (exact rationals)")
overall = None
rows = []
for j in range(1, 6):
    P = list(range(1, j + 1))
    for k, S0 in ((2, 12), (3, 12), (4, 12), (5, 8)):
        worst = None
        for rest in combinations(range(1, S0 + 1), k - 1):
            E = (0,) + rest
            r = rho_star_exact(P, E)
            if worst is None or r < worst[0]:
                worst = (r, E)
        rows.append((tuple(P), k, S0, worst[0], worst[1]))
        if overall is None or worst[0] < overall[0]:
            overall = (worst[0], tuple(P), worst[1])
        print(f"   P={{1..{j}}}, k={k}, S0={S0}: min rho* = {worst[0]} = {float(worst[0]):.5f} at E={worst[1]}", flush=True)
print(f"\nOVERALL MINIMUM: c0 = {overall[0]} = {float(overall[0]):.5f} at P={overall[1]}, E={overall[2]}")
print(f"VERDICT: {'c0 > 0 -- THE UNIFORM FLOOR IS POSITIVE (finite-exact, swept scope)' if overall[0] > 0 else 'ZERO FOUND -- floor fails at the argmin (structural news)'}")
print("DONE")
