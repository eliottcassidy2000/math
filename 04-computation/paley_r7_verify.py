#!/usr/bin/env python3
"""
paley_r7_verify.py — Verify r_7 = 390 for Paley T_11

PREDICTION: If anomaly depth D = 6 (HYP-550), then r_7^(k) = 390 for ALL k.
This would confirm:
  beta_5^(k=0) = 5, beta_5^(k>=1) = 0  =>  beta_5 = 5
  beta_6^(k=0) = 5, beta_6^(k>=1) = 1   =>  beta_6 = 5 + 10 = 15

Strategy: Compute r_7^(k=0) WITHOUT building full BD matrix.
BD @ omega_m_basis.T computed column-by-column.

Author: opus-2026-03-12-S68
"""
import sys, time, gc
import numpy as np
sys.path.insert(0, '04-computation')
from circulant_homology import PaleyHomology, find_nth_root_of_unity, _gauss_rank

p = 11
h = PaleyHomology(p=p)
h._ensure_enumerated(p)
prime = h.prime
omega_p = find_nth_root_of_unity(p, prime)

m = 7

# Get allowed paths
A_m = h._diff_seqs.get(m, [])
A_m1 = h._diff_seqs.get(m - 1, [])
print(f"|A_{m}| = {len(A_m)}, |A_{m-1}| = {len(A_m1)}")
print(f"Omega_{m} = {h.omega_dims(max_degree=m)[-1]}")
print(f"Omega_{m-1} = {h.omega_dims(max_degree=m-1)[-1]}")

A_m1_idx = {d: i for i, d in enumerate(A_m1)}
n_Am1 = len(A_m1)

for k in [0, 1]:
    omega_k = pow(omega_p, k, prime) if k > 0 else 1
    label = f"k={k}"
    print(f"\n{'='*60}")
    print(f"Computing r_{m}^({label})...")
    print(f"{'='*60}")

    t0 = time.time()

    # Step 1: omega_m basis
    print(f"  Building omega_{m} basis...", flush=True)
    omega_m_basis, face_data_m = h._omega_basis_k(m, omega_k)
    dim_Om = omega_m_basis.shape[0]
    print(f"  dim(Omega_{m}) = {dim_Om}, basis shape = {omega_m_basis.shape}")
    print(f"  omega_m_basis memory: {omega_m_basis.nbytes / 1e6:.0f} MB")
    sys.stdout.flush()

    # Step 2: omega_{m-1} basis
    print(f"  Building omega_{m-1} basis...", flush=True)
    omega_m1_basis, _ = h._omega_basis_k(m - 1, omega_k)
    dim_Om1 = omega_m1_basis.shape[0]
    print(f"  dim(Omega_{m-1}) = {dim_Om1}, basis shape = {omega_m1_basis.shape}")
    print(f"  omega_m1_basis memory: {omega_m1_basis.nbytes / 1e6:.0f} MB")
    sys.stdout.flush()

    # Step 3: Compute BD @ omega_m_basis.T WITHOUT storing BD
    # Result: bd_omega has shape (n_Am1, dim_Om)
    # Process in chunks to avoid OOM
    print(f"  Computing BD @ omega_basis.T column-by-column...", flush=True)
    bd_omega = np.zeros((n_Am1, dim_Om), dtype=np.int64)

    for j, (D, faces) in enumerate(zip(A_m, face_data_m)):
        if j % 10000 == 0 and j > 0:
            print(f"    processed {j}/{len(A_m)} paths ({j*100//len(A_m)}%)", flush=True)
        omega_j = omega_m_basis[:, j] if omega_m_basis.ndim == 2 else omega_m_basis[j]
        # Wait, omega_m_basis is (dim_Om, |A_m|), so column j is omega_m_basis[:, j]
        # Actually no — need to check: _omega_basis_k returns basis with rows=Omega vectors
        # So omega_m_basis shape is (dim_Om, |A_m|)
        # Column j = omega_m_basis[:, j] = coordinates of path j in each Omega basis vector
        # But we need BD @ omega_m_basis.T = BD @ (|A_m|, dim_Om) = (|A_m1|, dim_Om)
        # For path j in A_m, BD[:,j] has entries at the face rows
        # contribution to result: BD[:,j] * omega_m_basis.T[j,:] = BD[:,j] outer omega_j
        # But omega_j = omega_m_basis[:, j] (vector of length dim_Om)
        col_j = omega_m_basis[:, j]  # shape (dim_Om,)
        for fd, sign, offset, is_allowed in faces:
            if is_allowed and fd in A_m1_idx:
                row = A_m1_idx[fd]
                w = pow(omega_k, offset, prime) if offset != 0 else 1
                entry = (sign * w) % prime
                bd_omega[row] = (bd_omega[row] + entry * col_j) % prime

    print(f"  bd_omega memory: {bd_omega.nbytes / 1e6:.0f} MB")

    # Free omega_m_basis
    del omega_m_basis
    gc.collect()

    # Step 4: d_restricted = omega_m1_basis @ bd_omega
    print(f"  Computing d_restricted = omega_{m-1}_basis @ bd_omega...", flush=True)
    d_restricted = omega_m1_basis @ bd_omega % prime % prime
    print(f"  d_restricted shape = {d_restricted.shape}")

    del bd_omega, omega_m1_basis
    gc.collect()

    # Step 5: Gauss rank
    print(f"  Computing rank...", flush=True)
    rk = _gauss_rank(d_restricted, prime)
    dt = time.time() - t0

    print(f"\n  r_{m}^({label}) = {rk}  ({dt:.1f}s)")
    print(f"  PREDICTED: 390")
    print(f"  {'CONFIRMED!' if rk == 390 else 'MISMATCH — anomaly extends to m=7!'}")
    sys.stdout.flush()

    del d_restricted
    gc.collect()

print("\nDONE.")
