# Four-loop power counting on Paley tournaments: q X_1, q^2 X_2, q^3 X_3 at z = 1 for P_q (doubly regular, two-matrix
# reduction of loop_engine4), q = 19 .. 83.  HYP-9367 found q X_1 -> -7/2 and q^2 X_2 -> 29.  mac-mini-2026-10-22-loop4.
import sys, time
import numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'audit2')
import loop_engine4 as LE4
from common import paley_A, skew
for q in (19, 23, 31, 43, 47, 59, 67, 71, 79, 83):
    t = time.time()
    S = skew(paley_A(q)); assert LE4.check_drt(S)
    C = LE4.drt_coeffs(q, 1.0)
    X = {}
    for e in (1, 2, 3):
        X[e] = sum(np.sum(LE4.reduce_W(LE4.weight_tensor(e, i), C) * LE4.drt_gamma2(S, e, i)) / aut
                   for i, (V, edges, aut) in enumerate(LE4.load_cores(e)))
    print(f"P_{q}: q X_1 = {q * X[1]:+.6f}, q^2 X_2 = {q ** 2 * X[2]:+.6f}, q^3 X_3 = {q ** 3 * X[3]:+.4f}  ({time.time() - t:.0f}s)", flush=True)
