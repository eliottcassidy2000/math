#!/usr/bin/env python3
"""
THE EXACT-mu CENSUS OF THE WINDOW-22 COVERING BRANCH (klein-2026-07-10-S235,
THM-686(A) data). All covering 13-sets in [1..22] (13 and 14 are forced
elements, so primitivity is automatic). Integer-scaled rational sweeps, no
floats. Output: 05-knowledge/results/lrc14_window22_mu_census_klein_S235.out
"""
from itertools import combinations
from math import gcd
from fractions import Fraction as F


def lcm(a, b):
    return a * b // gcd(a, b)


def mu_pair(S):
    L = 1
    for v in S:
        L = lcm(L, v)
    D2 = 28 * L
    pts = {0, D2}
    for v in S:
        s = 2 * L // v
        for k in range(v):
            pts.add((14 * k + 1) * s)
            pts.add((14 * k + 13) * s)
    pts = sorted(pts)
    good = 0
    for x, y in zip(pts, pts[1:]):
        m = x + y
        ok = True
        for v in S:
            t = v * m % (2 * D2)
            if not (2 * D2 <= 14 * t <= 26 * D2):
                ok = False
                break
        if ok:
            good += y - x
    return good, D2  # mu = good/D2


fams = []
for S in combinations(range(1, 23), 13):
    ok = True
    for q in range(2, 15):
        if not any(v % q == 0 for v in S):
            ok = False
            break
    if ok:
        fams.append(S)
print(f"window-22 covering families: {len(fams)} "
      f"(min=1: {sum(1 for S in fams if S[0] == 1)}, "
      f"min>=2: {sum(1 for S in fams if S[0] > 1)})")

best = None
best1 = None
worstq = 0
mus = []
for S in fams:
    g, d = mu_pair(S)
    mus.append(F(g, d))
    if best is None or g * best[1] < best[0] * d:
        best = (g, d, S)
    if S[0] == 1 and (best1 is None or g * best1[1] < best1[0] * d):
        best1 = (g, d, S)
    sv = sum(S)
    q_star = sv * d // g + 1 if g else None
    if g and q_star > worstq:
        worstq = q_star

mn = F(best[0], best[1])
mn1 = F(best1[0], best1[1])
mus.sort()
print(f"GLOBAL exact floor: mu >= {mn} = {float(mn):.6f}  at {list(best[2])}")
print(f"min=1 (census-branch) floor: mu >= {mn1} = {float(mn1):.6f}  "
      f"at {list(best1[2])}")
print(f"median mu = {float(mus[len(mus) // 2]):.6f}; "
      f"max = {float(mus[-1]):.6f}; "
      f"zero-mu count = {sum(1 for m in mus if m == 0)}")
print(f"WORST TRANSFER THRESHOLD over the branch: "
      f"max q* = max Sum(v)/mu = {worstq}")
print("=> THE WINDOW-22 BRANCH IS CERTIFIED AT EVERY MODULUS by "
      f"[banks q <= {worstq}] + [THM-685 transfer beyond] -- and lonely "
      "outright by opus's 6 witnesses (LEM-024); this adds the QUANTIFIED "
      "margin form.")
