#!/usr/bin/env python3
"""audit H: negative controls for the exact balance test (it must reject what floats reject, and accept what floats accept
with a clear margin), on Z_5 (1,2,3,7,1) at k = 2, 3, 4 with the THM-4611 form, Q = I, and a perturbed form."""
import numpy as np
from hdp import Tables, balance_exact, margins_float, exact_rank_vec
tb = Tables(5, [1, 2, 3, 7, 1])
for Q in ([[20, -4, -8], [-4, 15, -3], [-8, -3, 20]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[20, -4, -8], [-4, 16, -3], [-8, -3, 20]]):
    for k in (2, 3, 4):
        T = tb.table(k).reshape(-1, 3, 3)
        nz = np.any(T.reshape(len(T), -1) != 0, axis=1)
        A = np.unique(T[nz].reshape(-1, 9), axis=0).reshape(-1, 3, 3)
        ok = balance_exact(Q, A)
        mg, al = margins_float(Q, A)
        rk = exact_rank_vec(A)
        agree = np.all(ok[mg > 1e-9]) and not np.any(ok[mg < -1e-9])
        print(f"Q = {Q}, k = {k}: {len(A)} blocks, exact-balanced {int(ok.sum())}, float margin min {mg.min():+.5f}, "
              f"exact/float agree away from 0: {agree}, min rank {rk.min()}, #rank<3 {int((rk < 3).sum())}")
