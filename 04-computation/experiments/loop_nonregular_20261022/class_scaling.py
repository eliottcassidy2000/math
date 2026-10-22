"""Power counting of the leaf classes X_(e,k)(1) on random almost-regular tournaments of growing even n.
Naive counting predicts X_(e,k) = O(n^-(e+k)): each excess costs 1/n and each leaf brings delta = O(1/n).
Session opus-2026-10-22-S54.  Usage: python class_scaling.py"""
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import leaf_engine as LE
from cluster_brute_classes import almost_regular

CLASSES = [(-1, 2), (-1, 3), (-1, 4), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (2, 0)]
rng = np.random.default_rng(7)
for n in (10, 14, 20, 30, 44, 64):
    vals = {c: [] for c in CLASSES}
    for rep in range(6 if n <= 30 else 3):
        A = almost_regular(n, rng, steps=20 * n * n)
        Sh = (A - A.T).astype(float) / n
        P = LE.LEF.stack(Sh.astype(complex), 1.0)
        for c in CLASSES:
            if c == (2, 0) and n > 44:
                continue
            vals[c].append(LE.X(Sh, c[0], c[1], 1.0, P).real)
    print(f'n={n}: ' + '; '.join(f'{c}: n^{c[0] + c[1]} X = {np.mean(np.array(v) * n ** (c[0] + c[1])):+.4f} '
                                 f'(sd {np.std(np.array(v) * n ** (c[0] + c[1])):.4f})' for c, v in vals.items() if v), flush=True)
