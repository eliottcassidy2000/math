# Audit (HYP-9371 section 8), power counting: q^4 X_4(P_q) at z = 1 recomputed in the engine's MATRIX (general)
# mode, which contracts the sign network with the actual n x n propagators over [n]^V (no doubly regular two-label
# reduction, no {P, Shat} homomorphism tensors), for q = 19 and 23, against paley_x4_scaling.out (DRT mode):
# 6340.7731 and 6055.0577.  Usage: python3 a5_task6_paley_general.py
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.dirname(HERE)
sys.path.insert(0, D); sys.path.insert(0, os.path.join(D, 'audit2'))
import loop_engine5 as LE5
from common import paley_A, skew
ref = {19: 6340.7731, 23: 6055.0577}
for q in (19, 23):
    t0 = time.time()
    S = skew(paley_A(q))
    x4 = LE5.X_general(S / q, 4, np.array([1.0]))[0].real
    print(f"P_{q}: matrix mode q^4 X_4(1) = {q ** 4 * x4:.4f}; DRT mode (paley_x4_scaling.out) {ref[q]:.4f}; "
          f"diff {q ** 4 * x4 - ref[q]:+.1e} ({time.time() - t0:.0f}s)", flush=True)
