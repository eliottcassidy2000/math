#!/usr/bin/env python3
"""Sticky reflections: a reflection coupling j -> b - j (a = -1) whose residue pattern is symmetric, mbar_(b-j) = mbar_j for
all j, maps to a reflection under every digit (a' = a * mbar_i / mbar_j = -1).  List them with the rank of their covariance."""
import numpy as np
from block_balance import int_coords, coupling_D, ratio_group
MAPS = [(5, [1, 2, 3, 7, 1]), (5, [1, 1, 2, 3, 7]), (5, [1, 2, 3, 7, 11]), (7, [1, 1, 1, 1, 2, 3, 5]), (7, [1, 1, 1, 2, 3, 5, 11]),
        (5, [1, 4, 1, 11, 34]), (5, [1, 6, 11, 11, 4]), (5, [1, 1, 6, 39, 11]), (5, [1, 14, 4, 1, 29]), (5, [1, 4, 31, 1, 6]), (5, [1, 6, 4, 31, 1])]
for d, m in MAPS:
    rho, v = int_coords(m); mb = [x % d for x in m]; G = ratio_group(m, d)
    out = []
    for b in range(d):
        a = d - 1
        if a not in G: continue
        if all(mb[(b - j) % d] == mb[j] for j in range(d)):
            D = np.array(coupling_D(v, a, b, d, rho)); out.append((f"j -> {b} - j", int(np.linalg.matrix_rank(D)) if D.any() else 0))
    print(f"Z_{d} m = {m}: residues {mb}; sticky reflections (rank): {out}")
