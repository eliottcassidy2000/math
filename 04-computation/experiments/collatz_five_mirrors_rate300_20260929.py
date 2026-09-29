#!/usr/bin/env python3
"""Rate analysis of the resonant coefficients to level 300 (S21): reads collatz_five_mirrors_powers_of_two300_20260929.out
(the closed recursion on the powers of two to h = 300), recomputes the no-descent probability P_h to 300 by dynamic
programming, and fits both to C h^-beta r^h over 100..300 against r0 = 3^(h* - 1), h* = h(log_3 2) = 0.94996.
Run: python 04-computation/experiments/collatz_five_mirrors_rate300_20260929.py
"""
from __future__ import annotations

import math
import re

import numpy as np

if __name__ == "__main__":
    vals = {}
    for line in open("05-knowledge/results/collatz_five_mirrors_powers_of_two300_20260929.out", encoding="utf-8"):
        m = re.match(r"\s+h=\s*(\d+): ([0-9.e+-]+)\s+argmax j=(-?\d+)", line)
        if m:
            vals[int(m.group(1))] = (float(m.group(2)), int(m.group(3)))
    L23 = math.log2(3)
    dist = {0: 1.0}
    P = {}
    for j in range(1, 301):
        nd = {}
        for p, pr in dist.items():
            for a in range(1, 70):
                q = p + a
                if q < j * L23:
                    nd[q] = nd.get(q, 0.0) + pr * 2.0 ** (-a)
        dist = nd
        P[j] = sum(dist.values())
    hstar = -(1 / L23) * math.log2(1 / L23) - (1 - 1 / L23) * math.log2(1 - 1 / L23)
    r0 = 3 ** (hstar - 1)
    print(f"h* = {hstar:.5f}, r0 = 3^(h*-1) = {r0:.5f}")
    print("h, M(h) = max_j |mu_hat_h(2^j)|, s - h log2 3, M/P_h, M/r0^h, P_h/r0^h:")
    for h in (20, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300):
        M, j = vals[h]
        print(f"   {h:3d}: {M:.4e}  {j - h * L23:+.2f}  {M / P[h]:.4f}  {M / r0 ** h:.4e}  {P[h] / r0 ** h:.4e}")
    loc = lambda h1, h2: math.log(vals[h1][0] / vals[h2][0], 2) / math.log2(h2 / h1)  # noqa: E731
    print("local exponents on doublings: " + ", ".join(f"{h}->{2*h}: {loc(h, 2*h):.2f}" for h in (25, 50, 75, 100, 125, 150)))
    print("   (they grow linearly; slope per level " + f"{(loc(150, 300) - loc(50, 100)) / 100:.4f} against |log2 r0| = {abs(math.log2(r0)):.4f})")
    print("geometric-mean level ratio of M: " + ", ".join(f"{a}-{b}: {(vals[b][0]/vals[a][0])**(1/(b-a)):.5f}" for a, b in ((100, 200), (150, 300), (200, 300), (250, 300))))
    print("geometric-mean level ratio of P_h: " + ", ".join(f"{a}-{b}: {(P[b]/P[a])**(1/(b-a)):.5f}" for a, b in ((100, 200), (150, 300), (200, 300))))
    for name, series in (("M", {h: v[0] for h, v in vals.items()}), ("P_h", P)):
        hs = np.array(sorted(h for h in series if h >= 100))
        ys = np.array([series[h] for h in hs])
        A = np.vstack([np.ones_like(hs, dtype=float), -np.log(hs), hs]).T
        coef, *_ = np.linalg.lstsq(A, np.log(ys), rcond=None)
        print(f"fit over 100..300: {name} ~ C h^-beta r^h with beta = {coef[1]:.3f}, r = {math.exp(coef[2]):.5f}; max log-residual {np.abs(A @ coef - np.log(ys)).max():.4f}")
        # with the rate pinned to r0: fit the prefactor only
        A0 = np.vstack([np.ones_like(hs, dtype=float), -np.log(hs)]).T
        c0, *_ = np.linalg.lstsq(A0, np.log(ys) - hs * math.log(r0), rcond=None)
        print(f"   with r pinned to r0: beta = {c0[1]:.3f}; max log-residual {np.abs(A0 @ c0 + hs * math.log(r0) - np.log(ys)).max():.4f}")
    print("ratio M/P_h over 20..300 (every 20 levels): " + ", ".join(f"{h}: {vals[h][0] / P[h]:.3f}" for h in range(20, 301, 20)))
    print("DONE")
