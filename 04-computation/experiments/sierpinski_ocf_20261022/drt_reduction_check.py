# Check of the doubly regular reduction (loop_engine4.drt_coeffs / reduce_W / drt_gamma2): for the doubly regular
# P_7 and P_11 the {P, Shat} two-matrix evaluation of X_1, X_2, X_3 equals the general five-matrix engine (itself
# checked against exact blow-up scaling in excess3_validate.out), as z-series and at z = 1.  mac-mini-2026-10-22-loop4.
import sys
import numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'audit2')
import loop_engine4 as LE4
from common import paley_A, skew
K, R = 32, 0.5
for q in (7, 11):
    S = skew(paley_A(q)); n = q; Sh = S / n
    assert LE4.check_drt(S)
    for e in (1, 2, 3):
        cores = LE4.load_cores(e)
        G2 = [LE4.drt_gamma2(S, e, i) for i in range(len(cores))]
        def val2(z):
            C = LE4.drt_coeffs(n, z)
            return sum(np.sum(LE4.reduce_W(LE4.weight_tensor(e, i), C) * G2[i]) / cores[i][2] for i in range(len(cores)))
        v2, vg = val2(1.0), LE4.X_e(Sh, e, 1.0)
        zz = R * np.exp(2j * np.pi * np.arange(K) / K)
        vals = np.array([val2(z) for z in zz])
        c2 = [(np.mean(vals * np.exp(-2j * np.pi * np.arange(K) * L / K)) / R ** L).real for L in range(13)]
        cg = LE4.series(Sh, e, 12, K=K, r=R)
        print(f"P_{q} excess {e}: z = 1: two-matrix {v2:+.13e}, general {vg:+.13e}; max series difference (L <= 12) {np.max(np.abs(np.array(c2) - cg)):.1e}", flush=True)
