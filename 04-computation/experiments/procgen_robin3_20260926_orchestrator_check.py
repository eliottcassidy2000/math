#!/usr/bin/env python3
"""Orchestrator audit of lane `robin3` (Robin inequality at shifts 2 and 3). Uses the orchestrator's own exact
counters N_m(L), A_M(L) (from procgen_robin2_20260926_orchestrator_check.py, written from the robin note's
definitions); the robin3 lane's code was not read.

  1. Theorem R_2 on a finite range: N_m(L) <= 2.1285 A_(m+2)(L), and the sharper per-range constant
     1.0473 for 2 <= m <= 24 (L <= 400 for m <= 14, L <= 250 for 15 <= m <= 24).
  2. Theorem R_3 on the same range: N_m(L) <= 1.1632 A_(m+3)(L).
  3. Shift 0 is false: N_3(L)/A_3(L) grows without bound (exact values at L = 50, 100, 200, 400).
"""
import importlib.util, os, sys
from fractions import Fraction
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "procgen_robin2_20260926_orchestrator_check.py")).read()
defs = src[:src.index("# 1. N_1 = 0")]
ns = {}
exec(defs, ns)
N_barrier, A_hard, check = ns["N_barrier"], ns["A_hard"], ns["check"]

bad2 = bad2s = bad3 = 0
best2 = best3 = Fraction(0)
for m in range(2, 25):
    L = 400 if m <= 14 else 250
    N = N_barrier(m, L)
    A2 = A_hard(m + 2, L)
    A3 = A_hard(m + 3, L)
    for l in range(1, L + 1):
        if 10000 * N[l] > 21285 * A2[l]:
            bad2 += 1
        if 10000 * N[l] > 10473 * A2[l]:
            bad2s += 1
        if 10000 * N[l] > 11632 * A3[l]:
            bad3 += 1
        if A2[l]:
            best2 = max(best2, Fraction(N[l], A2[l]))
        if A3[l]:
            best3 = max(best3, Fraction(N[l], A3[l]))
check(bad2 == 0 and bad2s == 0, f"Theorem R_2 finite range (2 <= m <= 24): N_m <= 2.1285 A_(m+2) and even <= 1.0473 A_(m+2); observed max ratio {float(best2):.7f}")
check(bad3 == 0, f"Theorem R_3 finite range (2 <= m <= 24): N_m <= 1.1632 A_(m+3); observed max ratio {float(best3):.7f}")
N3 = N_barrier(3, 400)
A3 = A_hard(3, 400)
rat = {l: float(Fraction(N3[l], A3[l])) for l in (50, 100, 200, 400)}
check(rat[400] > rat[200] > rat[100] > rat[50] > 1, "shift 0 fails: N_3(L)/A_3(L) = " + ", ".join(f"{v:.3g} (L={l})" for l, v in rat.items()) + " grows with L")
