# Validation of the three-loop engine: z-series coefficients of X_2 (tensor engine, Cauchy sums on |z| = r) against
# the brute-force excess-2 cluster sums of cluster_check.py (R_15 and P_19, L = 4, 6).  mac-mini-2026-10-22-sierp.
import sys
import numpy as np
sys.path.insert(0, '.')
import loop_engine_fast as LEF
from cluster_check import tournament
def series(Sh, excess, Lmax, K=24, r=0.5):
    vals = np.array([LEF.XB(Sh, excess, complex(r * np.exp(2j * np.pi * k / K))) for k in range(K)])
    return [(np.mean(vals * np.exp(-2j * np.pi * np.arange(K) * L / K)) / r ** L).real for L in range(Lmax + 1)]
brute = {('rot', 15): {1: {4: -1.991e-01, 6: -6.025e-02}, 2: {4: 1.244e-02, 6: 1.005e-01}, 3: {6: -1.766e-02}},
         ('paley', 19): {1: {4: -1.627e-01, 6: -6.814e-03}, 2: {4: 7.873e-03, 6: 6.033e-02}, 3: {6: -9.000e-03}}}
for (kind, n), bf in brute.items():
    A = tournament(kind, n); Sh = (A - A.T) / n
    for e in (1, 2):
        c = series(Sh, e, 8)
        print(f"{kind} n={n} excess {e}: z^4 = {c[4]:+.4e} (brute {bf[e].get(4, float('nan')):+.4e}), z^6 = {c[6]:+.4e} (brute {bf[e].get(6, float('nan')):+.4e}), z^8 = {c[8]:+.4e}; odd z^5 = {c[5]:+.1e}", flush=True)
