#!/usr/bin/env python3
"""audit H: exact non-existence certificates (as in f_dual_certificates.py) for the sixth {+-1} map (1,6,4,31,1), k = 1..4."""
import time, numpy as np
from f_dual_certificates import blocks, psd_sqrt, cutting_plane, exact_cert, margins_at
from hdp import exact_rank_vec
m = [1, 6, 4, 31, 1]
for k in (1, 2, 3, 4):
    t0 = time.time()
    U = blocks(m, k)
    rk = exact_rank_vec(U)
    if rk.min() < 3:
        print(f"Z_5 m = {m}, k = {k}: {len(U)} blocks; a block of exact rank {rk.min()} -> no form (exact by rank)", flush=True); continue
    Uf = U.astype(float); An = Uf / np.trace(Uf, axis1=1, axis2=2)[:, None, None]
    Ah, Aih = psd_sqrt(An)
    s, best, duals, cuts_A, cuts_w = cutting_plane(An, Ah, Aih)
    ok, used, eig = exact_cert(U, cuts_A, cuts_w, duals)
    print(f"Z_5 m = {m}, k = {k}: {len(U)} blocks; LP max-min s = {s:+.6f}; dual certificate with {used} cuts: negative definite EXACTLY: {ok}  [{time.time() - t0:.1f}s]", flush=True)
