#!/usr/bin/env python3
"""Travelling ridges in the negative-power family of the 3-adic Syracuse law (S22).

The family m_n(k) = mu_hat_n(2^k mod 3^n) is periodic in k with period L_n = 2 3^(n-1).  At the lowest levels the period
is small, so the negative exponents -m contain wrapped copies of the positive resonance (n = 4: L = 54, resonance at
k = 4..6, copies at m = 48..50, 102..104, ...; n = 5: L = 162, copies at m = 154..157; n = 6: L = 486, m = 476..482).
Under the recursion a copy propagates to the next level at exponents shifted by +a (weights 2^-a), so it travels toward
m = 0 at about 1.3 per level.  This script follows the ridges of v_n(m) = 3^(n/2) |m_n(-m)| over m in [1, WIDTH]
to level NMAX: local maxima above 1.5, the window's Fourier mass, and the amplitude along the predicted ridge paths.
Run: python 04-computation/experiments/collatz_five_mirrors_ridges_20260929.py [NMAX] [WIDTH]
"""
from __future__ import annotations

import math
import sys
import time

sys.path.insert(0, "04-computation/experiments")
from collatz_five_mirrors_renewal_20260929 import closed_family_all  # noqa: E402

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    WIDTH = int(sys.argv[2]) if len(sys.argv) > 2 else 600
    t0 = time.time()
    fam = closed_family_all(NMAX, -WIDTH, 0)
    print(f"(recursion to level {NMAX}, window [-{WIDTH}, 0], {time.time() - t0:.0f}s)")
    print("periods: " + ", ".join(f"L_{n} = {2 * 3 ** (n - 1)}" for n in range(1, 8)))
    print("local maxima of v_n(m) = 3^(n/2) |m_n(-m)| above 1.5 (m, value), up to six per level, and the window's Fourier mass 3^n sum_m |m_n(-m)|^2 / WIDTH")
    for n in range(3, NMAX + 1):
        row = fam[n]
        v = [0.0] + [abs(row[-m]) * 3 ** (n / 2) for m in range(1, WIDTH + 1)]
        peaks = []
        for m in range(1, WIDTH + 1):
            if v[m] >= 1.5 and all(v[m] >= v[mm] for mm in range(max(1, m - 3), min(WIDTH, m + 3) + 1)):
                peaks.append((v[m], m))
        peaks.sort(reverse=True)
        mass = sum(abs(row[-m]) ** 2 for m in range(1, WIDTH + 1)) * 3 ** n / WIDTH
        if n <= 12 or n % 2 == 0 or peaks:
            print(f"   n={n:3d}: mass {mass:.3f}; peaks " + " ".join(f"(m={m}, {val:.2f})" for val, m in peaks[:6]))
    print("amplitude along the predicted ridge paths (max of v over m within +-12 of the predicted position; '-' if outside the window)")
    seeds = [(4, 49), (4, 103), (4, 157), (4, 211), (5, 155), (5, 317), (6, 479)]
    print("   n     " + "  ".join(f"seed({n0},{m0})" for n0, m0 in seeds))
    for n in range(6, NMAX + 1, 2):
        row = fam[n]
        cells = []
        for n0, m0 in seeds:
            mp = m0 - 1.3 * (n - n0)
            if mp < 1 or mp > WIDTH:
                cells.append("      -     ")
                continue
            lo, hi = max(1, int(mp) - 12), min(WIDTH, int(mp) + 12)
            best = max(range(lo, hi + 1), key=lambda m: abs(row[-m]))
            cells.append(f"m={best:3d} {abs(row[-best]) * 3 ** (n / 2):6.2f}")
        print(f"   {n:3d}  " + "  ".join(cells))
    print("DONE")
