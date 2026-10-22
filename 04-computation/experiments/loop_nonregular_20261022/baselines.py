"""Cheap 'before' baselines for the prediction tables (closed forms only).  Session opus-2026-10-22-S54.

y_mm   = -n X_(1,0) - E0[a]               (the regular two-loop formula, HYP-9365, applied blindly)
y_mm_t = y_mm - n X_(-1,2)                (plus the path trees; mac-mini's almost_regular_tree.py)
Writes pred{n}.txt-compatible rows (code H y y_mm y_mm_t) for files whose expensive predictions were not run.
Usage: python baselines.py n FILE OUT"""
import sys, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from predict_nonregular import cumulants, X1_closed, tree_closed, load

n = int(sys.argv[1]); path = sys.argv[2]; out = sys.argv[3]
with open(out, 'w') as fh:
    fh.write('code H y y_mm y_mm_t\n')
    for code, H, S in load(n, path):
        Sh = S / n
        x = np.linalg.eigvals(Sh).imag
        rho = float(np.sqrt(np.prod((1 + x ** 2) / (1 - x ** 2))))
        k = cumulants(lambda z: 0.5 * np.sum(np.log((1 + (z * x) ** 2) / (1 - (z * x) ** 2))), kmax=2)
        Ea = (k[2] + k[1] ** 2 - k[1]) / 2
        y_mm = -n * X1_closed(Sh, 1.0).real - Ea
        y_mm_t = y_mm - n * tree_closed(Sh, 1.0)
        kap = H * 2 ** (n - 1) / math.factorial(n)
        fh.write(f'{code} {H} {n * (1 - kap / rho):.10g} {y_mm:.10g} {y_mm_t:.10g}\n')
