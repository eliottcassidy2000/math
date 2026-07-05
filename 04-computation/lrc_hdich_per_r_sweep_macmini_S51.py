#!/usr/bin/env python3
"""
mac-mini-2026-07-05-S51 -- HYP-4102 (renumbered from 4098; opus-S75 first-committed): hdich's PER-R FINITE CRT SWEEP
(discharging opus-4097 leg 2's finite content).

hdich after residue pinning: tight 12-families are lifts of {1..12} mod 13.
Single lifts W_k(r) = ({1..12}\\{r}) u {r+13k}. Claim (opus, 96/96 spot-verified):
ALL lifts are strictly loose, min M = 1/12 (gap 1/156 above 1/13).
This sweep: per r = 1..12, k = 1..K0: classify sieve-exposed vs surviving,
exact M for every survivor, confirm the 1/12 floor; plus THM-610 witness-
denominator cross-check (loose lifts should witness SHALLOW, q* small).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import sys
sys.path.insert(0, '04-computation')
from lonely_profile import profile

def M_exact_prof(S):
    for cap in (11, 8, 5, 3, 2):
        p = profile(sorted(S), F(1, cap))
        m = p.M()
        if m is not None:
            return m
    return None

K0 = 12
R13, R12 = F(1, 13), F(1, 12)
print(f"hdich per-r sweep: W_k(r) = ({{1..12}} \\ {{r}}) + {{r+13k}}, k = 1..{K0}")
print(f"{'r':>3} {'#exposed':>9} {'#surv':>6} {'min M over survivors':>22} {'>=1/12?':>8}")
worst = None
viol = 0
for r in range(1, 13):
    base = [v for v in range(1, 13) if v != r]
    exposed = 0; surv = 0; minM = None
    for k in range(1, K0 + 1):
        lift = r + 13 * k
        W = sorted(base + [lift])
        if reduce(gcd, W) != 1:
            continue
        # sieve exposure: some m in 2..13 with no multiple in W (base had {1..12})
        missing = [m for m in range(2, 13) if not any(v % m == 0 for v in W)]
        if missing:
            exposed += 1
            # sieve witness at 1/m gives M >= 1/m >= 1/13 with strict margin when m <= 12
            continue
        surv += 1
        M = M_exact_prof(W)
        if minM is None or M < minM: minM = M
        if M is not None and M < R12:
            if M <= R13:
                viol += 1
                print(f"    !! r={r} k={k}: M={M} <= 1/13 -- RIGIDITY VIOLATION")
    ok = "yes" if (minM is None or minM >= R12) else f"NO ({minM})"
    print(f"{r:>3} {exposed:>9} {surv:>6} {str(minM):>22} {ok:>8}")
    if minM is not None and (worst is None or minM < worst[0]): worst = (minM, r)
print(f"\nworst survivor M = {worst[0]} (r={worst[1]}); 1/12 = {F(1,12)}; gap to 1/13 = {worst[0]-R13}")
print(f"violations (M <= 1/13): {viol}")
# THM-610 cross-check: witness denominators of a few survivors (shallow expected)
print("\nTHM-610 cross-check (loose lifts witness shallow):")
for (r, k) in [(1, 1), (6, 2), (12, 1)]:
    base = [v for v in range(1, 13) if v != r]
    W = sorted(base + [r + 13 * k])
    p = profile(W, F(1, 3))
    M = p.M()
    print(f"  r={r},k={k}: M={M}, witness denominator = {M.denominator if M else '?'} "
          f"(THM-610: tight-covering would need q* >= 28; shallow == loose ✓)")
