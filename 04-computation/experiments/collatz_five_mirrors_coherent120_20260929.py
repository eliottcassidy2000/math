#!/usr/bin/env python3
"""Bounded-excursion families at h = 120 (S21): the words with P_j < j log2 3 + c for all j, their mass and their
share of the resonant coefficient, for c = 0, 2, 4, 6, 8, 12, 16, 24.  Uses the DP of collatz_five_mirrors_coherent_20260929.py.
Run: python 04-computation/experiments/collatz_five_mirrors_coherent120_20260929.py
"""
from __future__ import annotations

import math
import sys

sys.path.insert(0, "04-computation/experiments")
from collatz_five_mirrors_coherent_20260929 import LOG23, closed_family, coherent_split  # noqa: E402

if __name__ == "__main__":
    h = 120
    s0 = int(round(h * LOG23)) - 6
    fam = closed_family(h, s0 - 8, s0 + 8)
    s_star = max(fam, key=lambda j: abs(fam[j]))
    full = fam[s_star]
    print(f"h={h}: s* = {s_star} (h log2 3 - s* = {h * LOG23 - s_star:.2f}); |full| = {abs(full):.4e}")
    for c in (0, 2, 4, 6, 8, 12, 16, 24):
        coh, coh_abs, mass = coherent_split(h, s_star, Awindow=45 + c, c=c)
        print(f"   c={c:2d}: family mass {mass:.4e}; |coh_c| = {abs(coh):.4e} = {abs(coh) / abs(full):.4f} |full|; |full - coh_c|/|full| = {abs(full - coh) / abs(full):.4f}; coherent fraction {abs(coh) / coh_abs:.4f}")
        sys.stdout.flush()
    print("DONE")
