#!/usr/bin/env python3
"""The J-terms of the resonant coefficient at the scale h^(-3/2) e^(-hI): do they converge? (S22 obligation 2)

c_J = contribution of the words with exactly J bottom levels at negative exponent (exact DP, Awindow = 80), at the
argmax exponent s* of the powers-of-two family; gamma_J(h) := c_J / (h^(-3/2) e^(-hI)); partial sums; the full
coefficient at the same scale and against the no-descent probability P_h.  Levels h = 40, 80, 120, 160, 200.
Run: python 04-computation/experiments/collatz_five_mirrors_gamma_20260929.py [h ...]
"""
from __future__ import annotations

import cmath
import math
import sys
import time

import numpy as np

sys.path.insert(0, "04-computation/experiments")
from collatz_five_mirrors_coherent_20260929 import LOG23, closed_family  # noqa: E402
from collatz_five_mirrors_renewal_20260929 import decompose  # noqa: E402
from collatz_five_mirrors_renewal_bound_20260929 import IRATE  # noqa: E402


def no_descent(h: int) -> float:
    lim = [0] + [int(math.floor(j * LOG23 - 1e-12)) for j in range(1, h + 1)]
    dp = np.zeros(lim[h] + 1)
    dp[0] = 1.0
    for j in range(1, h + 1):
        new = np.zeros_like(dp)
        for a in range(1, 41):
            hi = min(lim[j - 1] + a, lim[j])
            if hi < a:
                continue
            new[a: hi + 1] += 2.0 ** (-a) * dp[: hi - a + 1]
        dp = new
    return float(dp.sum())


if __name__ == "__main__":
    t0 = time.time()
    hs = [int(x) for x in sys.argv[1:]] or [40, 80, 120, 160, 200]
    table = {}
    for h in hs:
        s0 = int(round(h * LOG23)) - 6
        fam = closed_family(h, s0 - 8, s0 + 8)
        s_star = max(fam, key=lambda j: abs(fam[j]))
        full = fam[s_star]
        mass, contrib = decompose(h, s_star, "J", Awindow=80)
        scale = h ** (-1.5) * math.exp(-h * IRATE)
        Ph = no_descent(h)
        cum = np.cumsum(contrib)
        print(f"h={h}: s* = {s_star} (delta = {h * LOG23 - s_star:.2f}), |full| = {abs(full):.4e}, |full|/scale = {abs(full) / scale:.4f}, arg full = {cmath.phase(full):+.3f}, P_h = {Ph:.4e}, |full|/P_h = {abs(full) / Ph:.4f}, P_h/scale = {Ph / scale:.4f}  ({time.time() - t0:.0f}s)")
        print("   J: gamma_J = c_J/scale (modulus, arg); partial sum modulus; remainder")
        for J in range(0, 25):
            g = contrib[J] / scale
            print(f"   {J:2d}: {abs(g):.4f} {cmath.phase(g):+.2f}   {abs(cum[J]) / scale:.4f}   {abs(full - cum[J]) / abs(full):.3f}")
        table[h] = [contrib[J] / scale for J in range(25)]
        sys.stdout.flush()
    print("gamma_J(h) moduli across levels (rows J, columns h):")
    for J in range(0, 21):
        print(f"   J={J:2d}: " + "  ".join(f"{abs(table[h][J]):.4f}" for h in hs))
    print("gamma_J(h) arguments across levels:")
    for J in range(0, 21):
        print(f"   J={J:2d}: " + "  ".join(f"{cmath.phase(table[h][J]):+.2f}" for h in hs))
    print("DONE")
