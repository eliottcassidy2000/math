#!/usr/bin/env python3
"""The amplitude law of the multiplier families u 2^j (S22): family maximum / pure-family maximum against u^(theta*/ln 2).

The corridor of the family u 2^j is narrower by log2 u bits (phase e(u 2^(kappa_j - j log2 3)) in the corridor), so the
Chernoff mass law with delta -> delta + log2 u predicts max_j |mu_hat_h(u 2^j)| ~ M(h) e^(theta* log2 u) = M(h) u^(-0.438),
theta* = ln(2(1 - log_3 2)).  Each family is closed under the frequency recursion (u 2^j -> u 2^(j-a)); computed exactly
(truncation a <= 40) at the levels h = 20, 30, 40, 50 for u = 1, 5, 7, ..., 49, 55, 65, 85, 127.
Run: python 04-computation/experiments/collatz_five_mirrors_multiplier_law_20260929.py
"""
from __future__ import annotations

import cmath
import math
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "04-computation/experiments")
from collatz_five_mirrors_coherent_20260929 import LOG23  # noqa: E402
from collatz_five_mirrors_renewal_bound_20260929 import THETA  # noqa: E402

AMAX = 40


def family_max(HMAX: int, u: int, levels):
    """max_j |mu_hat_n(u 2^j mod 3^n)| over j in [-20, n log2 3 + 10] at the given levels n."""
    j_hi = int(HMAX * LOG23) + 12
    j_lo = -20
    weights = [2.0 ** (-a) for a in range(1, AMAX + 1)]
    prev = {j: 1.0 + 0j for j in range(j_lo - AMAX * HMAX, j_hi + 1)}
    out = {}
    for n in range(1, HMAX + 1):
        mod = 3 ** n
        lo = j_lo - AMAX * (HMAX - n)
        ph = {}
        for j in range(lo - AMAX, j_hi + 1):
            r = (u * pow(2, j, mod)) % mod
            ph[j] = cmath.exp(2j * math.pi * float(Fraction(r, mod)))
        cur = {}
        for j in range(lo, j_hi + 1):
            s = 0j
            for a in range(1, AMAX + 1):
                s += weights[a - 1] * ph[j - a] * prev.get(j - a, 0j)
            cur[j] = s
        prev = cur
        if n in levels:
            jm = max(range(j_lo, int(n * LOG23) + 12), key=lambda j: abs(cur[j]))
            out[n] = (abs(cur[jm]), jm)
    return out


if __name__ == "__main__":
    levels = (20, 30, 40, 50)
    us = [1, 5, 7, 11, 13, 17, 19, 23, 25, 29, 31, 35, 37, 41, 43, 47, 49, 55, 65, 85, 127]
    res = {u: family_max(50, u, levels) for u in us}
    print(f"predicted exponent theta*/ln 2 = {THETA / math.log(2):.4f}")
    for n in levels:
        pure = res[1][n][0]
        print(f"h={n}: pure-family max {pure:.4e} at j = {res[1][n][1]} (= h log2 3 - {n * LOG23 - res[1][n][1]:.2f})")
        print("   u: family max / pure, predicted u^-0.438, argmax j offset (j_u - j_1 + log2 u)")
        for u in us[1:]:
            fm, jm = res[u][n]
            print(f"   {u:3d}: {fm / pure:.3f}  {u ** (THETA / math.log(2)):.3f}  {jm - res[1][n][1] + math.log2(u):+.2f}")
        xs = np.log([u for u in us[1:]])
        ys = np.log([res[u][n][0] / pure for u in us[1:]])
        slope, icpt = np.polyfit(xs, ys, 1)
        print(f"   least-squares exponent over u = 5..127: {slope:.3f} (prefactor {math.exp(icpt):.3f})")
    print("DONE")
