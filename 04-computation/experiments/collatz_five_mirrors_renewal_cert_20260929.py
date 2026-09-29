#!/usr/bin/env python3
"""Certification and extension of the renewal picture (S22).

(i)   Ntilde_n (n <= 80) recomputed with the valuation truncation a <= 60 instead of a <= 40 (relative agreement).
(ii)  J-decomposition at h = 40, 80, 120: remainder |full - cum_J|/|full|, J_eff (remainder below 10%), phases of c_J,
      per-mass |c_J|/mass_h(J) for J <= 12 (the 1/h test of the top-part weights).
(iii) exact mass_h(J) e^(hI) e^(-theta delta) for h = 40, 80, 120, 300, 1000, J <= 24: growth (log2 3 - 1)^-J and the
      Gaussian window (O(1)-per-step recursion p_(N+1)(t) = (p_N(t-1) + p_(N+1)(t-1))/2).
(iv)  the bound sum_J mass_h(J) Ntilde_(J-1) for h up to 2000 with Ntilde_n = 1.3 * 3^(-n/2) beyond n = 80 (hypothesis H
      with the observed constant), reported as bound / e^(-hI).
Run: python 04-computation/experiments/collatz_five_mirrors_renewal_cert_20260929.py
"""
from __future__ import annotations

import cmath
import math
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, "04-computation/experiments")
from collatz_five_mirrors_coherent_20260929 import LOG23, closed_family  # noqa: E402
from collatz_five_mirrors_renewal_20260929 import closed_family_all, decompose  # noqa: E402
from collatz_five_mirrors_renewal_bound_20260929 import IRATE, THETA  # noqa: E402


def closed_family_all_amax(HMAX: int, j_lo: int, j_hi: int, AMAX: int):
    weights = [2.0 ** (-a) for a in range(1, AMAX + 1)]
    prev = {j: 1.0 + 0j for j in range(j_lo - AMAX * HMAX, j_hi + 1)}
    out = {}
    for n in range(1, HMAX + 1):
        mod = 3 ** n
        lo = j_lo - AMAX * (HMAX - n)
        pw = {j: pow(2, j, mod) for j in range(lo - AMAX, j_hi + 1)}
        ph = {j: cmath.exp(2j * math.pi * float(Fraction(pw[j], mod))) for j in pw}
        cur = {}
        for j in range(lo, j_hi + 1):
            s = 0j
            for a in range(1, AMAX + 1):
                s += weights[a - 1] * ph[j - a] * prev.get(j - a, 0j)
            cur[j] = s
        prev = cur
        out[n] = {j: cur[j] for j in range(j_lo, j_hi + 1)}
    return out


def cost_laws(h: int, smax: int):
    """p[N][t] = P(S_N = t), t <= smax, N <= h, by the geometric recursion (O(smax) per step)."""
    P = np.zeros((h + 1, smax + 1))
    P[0, 0] = 1.0
    for N in range(1, h + 1):
        prev = P[N - 1]
        cur = P[N]
        for t in range(1, smax + 1):
            cur[t] = 0.5 * (prev[t - 1] + cur[t - 1])
    return P


def mass_law_fast(h: int, s: int, JMAX: int, P):
    out = []
    ts = np.arange(0, s + 1)
    w = 2.0 ** (-(s - ts))
    for J in range(0, JMAX + 1):
        N = h - J
        if N < 0:
            out.append(0.0)
        elif J == 0:
            out.append(float(P[N, : s + 1].sum()))
        else:
            out.append(float((P[N, : s + 1] * w).sum()))
    return out


if __name__ == "__main__":
    t0 = time.time()
    print("== (i) truncation certification of the negative-power family: a <= 60 against a <= 40, n <= 80 ==")
    f40 = closed_family_all(80, -60, 0)
    f60 = closed_family_all_amax(80, -60, 0, 60)
    worst = 0.0
    for n in range(1, 81):
        N40 = sum(2.0 ** (-m) * abs(f40[n][-m]) for m in range(1, 61))
        N60 = sum(2.0 ** (-m) * abs(f60[n][-m]) for m in range(1, 61))
        rel = abs(N40 - N60) / N60
        worst = max(worst, rel)
        if n in (20, 40, 60, 80):
            e = max(abs(f40[n][-m] - f60[n][-m]) / abs(f60[n][-m]) for m in range(1, 61))
            print(f"   n={n}: Ntilde (a<=40) {N40:.6e}, (a<=60) {N60:.6e}, rel diff {rel:.1e}; worst single-coefficient rel diff over m <= 60: {e:.1e}")
    print(f"   worst relative difference of Ntilde_n over n <= 80: {worst:.1e}   ({time.time() - t0:.0f}s)")
    print("== (ii) J-decomposition at h = 40, 80, 120: remainders, phases, per-mass weights ==")
    permass = {}
    for h in (40, 80, 120):
        s0 = int(round(h * LOG23)) - 6
        famh = closed_family(h, s0 - 8, s0 + 8)
        s_star = max(famh, key=lambda j: abs(famh[j]))
        full = famh[s_star]
        mass, contrib = decompose(h, s_star, "J")
        cum = np.cumsum(contrib)
        rem = [abs(full - cum[J]) / abs(full) for J in range(31)]
        jeff = next((J for J in range(31) if rem[J] < 0.1), None)
        jeff5 = next((J for J in range(31) if all(r < 0.05 for r in rem[J:])), None)
        print(f"   h={h}: s* = {s_star}, |full| = {abs(full):.4e}, mass covered {mass.sum():.4f}, J_eff(10%) = {jeff}, J_eff(5%, stable) = {jeff5}   ({time.time() - t0:.0f}s)")
        print("      J: mass, |c_J|/|full|, arg c_J, per-mass |c_J|/mass, remainder")
        for J in range(0, 25):
            print(f"      {J:2d}: {mass[J]:.3e}  {abs(contrib[J]) / abs(full):.3f}  {cmath.phase(contrib[J]):+.2f}  {abs(contrib[J]) / mass[J] if mass[J] > 0 else 0:.3e}  {rem[J]:.3f}")
        permass[h] = [abs(contrib[J]) / mass[J] if mass[J] > 0 else 0 for J in range(13)]
        sys.stdout.flush()
    print("   per-mass weights |c_J|/mass_h(J) times h (the 1/h test), J = 2..10:")
    for h in (40, 80, 120):
        print(f"      h={h}: " + " ".join(f"{permass[h][J] * h:.2f}" for J in range(2, 11)))
    print("== (iii) exact mass_h(J) e^(hI) e^(-theta delta) (Chernoff-normalised), growth (log2 3 - 1)^-J = 1.7095^J and the window ==")
    for h in (40, 80, 120, 300, 1000):
        s = int(round(h * LOG23)) - 6
        delta = h * LOG23 - s
        P = cost_laws(h, s + 1)
        em = mass_law_fast(h, s, 24, P)
        norm = math.exp(h * IRATE) * math.exp(-THETA * delta)
        vals = [em[J] * norm for J in range(25)]
        ratios = [vals[J] / vals[J - 1] if vals[J - 1] > 0 else 0 for J in range(1, 25)]
        print(f"   h={h}: delta = {delta:.2f}; normalised mass J=0..24: " + " ".join(f"{v:.2e}" for v in vals))
        print(f"          ratios J/J-1: " + " ".join(f"{r:.2f}" for r in ratios))
    print("== (iv) bound sum_J mass_h(J) Ntilde_(J-1) / e^(-hI), Ntilde_n = 1.3 * 3^(-n/2) for n > 80 (hypothesis H, observed constant) ==")
    Nt = {0: 1.0}
    for n in range(1, 81):
        Nt[n] = sum(2.0 ** (-m) * abs(f60[n][-m]) for m in range(1, 61)) + 2.0 ** (-60)
    for h in (100, 200, 300, 500, 1000, 2000):
        s = int(round(h * LOG23)) - 6
        P = cost_laws(h, s + 1)
        JMAX = min(h, 400)
        em = mass_law_fast(h, s, JMAX, P)
        b = 0.0
        for J in range(0, JMAX + 1):
            if J <= 1:
                Ntm = 1.0
            elif J - 1 <= 80:
                Ntm = Nt[J - 1]
            else:
                Ntm = 1.3 * 3 ** (-(J - 1) / 2)
            b += em[J] * Ntm
        print(f"   h={h:4d}: bound / e^(-hI) = {b * math.exp(h * IRATE):.4f}   (delta = {h * LOG23 - s:.2f}, {time.time() - t0:.0f}s)")
    print("DONE")
