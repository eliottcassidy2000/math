#!/usr/bin/env python3
"""d = 4 (multipliers 1, 3, 5, 7; r_i = -m_i i mod 4): no fixed-length block certificate exists.  For every k the state
M = 1, e = 2*4^(k-1) (mod 4^k) runs k-1 identity couplings and then the translation j -> j + 2, whose covariance has rank 2,
so its block matrix has rank 2 and no form can balance it (tr <= 2 lambda_max).  Checked exactly for k = 1..5; also lists
the minimal block rank over all states."""
import numpy as np
from block_balance import int_coords, coupling_D, std_r
from block_balance2 import block_arrays
d, m = 4, [1, 3, 5, 7]; r = std_r(m)
rho, v = int_coords(m)
D12 = np.array(coupling_D(v, 1, 2, d, rho))
print(f"m = {m}, r = {r}, rank {rho}; D(translation by 2) = {D12.tolist()} (rank {np.linalg.matrix_rank(D12)})")
for k in range(1, 6):
    H, A = block_arrays(d, m, r, v, rho, k)
    iM = list(H).index(1); e = 2 * 4 ** (k - 1)
    Ak = A[iM, e]
    flat = A.reshape(-1, rho, rho); nz = np.any(flat.reshape(len(flat), -1) != 0, axis=1)
    minrank = np.linalg.matrix_rank(flat[nz].astype(float)).min()
    print(f"k = {k}: block at (M, e) = (1, 2*4^{k-1}) = {Ak.tolist()}  equals 4^(k-1) D(1,2): {np.array_equal(Ak, 4 ** (k - 1) * D12)}; "
          f"min block rank over all {nz.sum()} nonzero states: {minrank}")
