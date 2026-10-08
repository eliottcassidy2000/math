#!/usr/bin/env python3
"""audit H: recompute and DUMP the exact non-existence certificates for fixed-length block forms to dual_certificates.json:
for each (map, k): a list of (block matrix A_l (integer), rational z_l, rational weight c_l) such that
  N = sum_l c_l [A_l - 2 (A_l z_l)(A_l z_l)^T / (z_l^T A_l z_l)]   is negative definite.
Then l_verify_certs.py re-checks each certificate from the JSON alone plus an independent recomputation of the block set."""
import json, time
from fractions import Fraction as Fr
import numpy as np
from f_dual_certificates import blocks, psd_sqrt, cutting_plane, frac_vec

def cert_data(U, cuts_A, cuts_w, duals):
    out = []
    for i, w, c in zip(cuts_A, cuts_w, duals):
        if c <= 1e-12: continue
        A = U[i]; lam, V = np.linalg.eigh(A.astype(float))
        z = V @ ((V.T @ w) / np.sqrt(lam)); z = frac_vec(z / np.abs(z).max())
        cc = Fr(int(round(c * 10 ** 12)), 10 ** 12) / Fr(int(np.trace(A)))
        out.append({'A': A.tolist(), 'z': [str(x) for x in z], 'c': str(cc)})
    return out

jobs = [([1, 2, 3, 7, 1], 2), ([1, 2, 3, 7, 1], 3)]
for m in ([1, 4, 1, 11, 34], [1, 6, 11, 11, 4], [1, 1, 6, 39, 11], [1, 14, 4, 1, 29], [1, 4, 31, 1, 6], [1, 6, 4, 31, 1]):
    for k in (2, 3, 4):
        jobs.append((m, k))
res = []
for m, k in jobs:
    t0 = time.time()
    U = blocks(m, k)
    if np.linalg.matrix_rank(U.astype(float)).min() < 3:
        res.append({'m': m, 'k': k, 'rank_obstruction': True}); print(m, k, 'rank obstruction', flush=True); continue
    Uf = U.astype(float); An = Uf / np.trace(Uf, axis1=1, axis2=2)[:, None, None]
    Ah, Aih = psd_sqrt(An)
    s, best, duals, cuts_A, cuts_w = cutting_plane(An, Ah, Aih)
    res.append({'m': m, 'k': k, 'cuts': cert_data(U, cuts_A, cuts_w, duals)})
    print(m, k, f'{len(res[-1]["cuts"])} cuts  [{time.time() - t0:.1f}s]', flush=True)
json.dump(res, open('dual_certificates.json', 'w'), indent=1)
