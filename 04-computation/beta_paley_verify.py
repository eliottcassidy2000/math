#!/usr/bin/env python3
"""
beta_paley_verify.py - Verify Paley T_7 Betti numbers carefully

Check for numerical issues and verify chi independently.

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis, path_betti_numbers
)
sys.stdout = _saved


# Paley T_7
n = 7
qr = {1, 2, 4}
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j and ((j - i) % n) in qr:
            A[i][j] = 1

# Verify tournament
for i in range(n):
    for j in range(i+1, n):
        assert A[i][j] + A[j][i] == 1, f"Not tournament at {i},{j}"
    assert A[i][i] == 0

scores = [sum(row) for row in A]
print(f"Paley T_7 scores: {scores} (all {scores[0]}? {all(s == scores[0] for s in scores)})")

# Compute full Betti numbers up to dimension n-1 = 6
print(f"\nComputing full Betti numbers (max_dim={n-1})...")
betti = path_betti_numbers(A, n, max_dim=n-1)
print(f"  beta = {betti}")

# Independent chi computation from Omega dims
print(f"\nIndependent Euler characteristic computation:")
allowed = {}
for p in range(-1, n + 1):
    if p < 0:
        allowed[p] = []
    else:
        allowed[p] = enumerate_allowed_paths(A, n, p)

omega_dims = []
A_dims = []
for p in range(n):
    omega_p = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
    dim_omega = omega_p.shape[1] if omega_p.ndim == 2 else 0
    omega_dims.append(dim_omega)
    A_dims.append(len(allowed[p]))

print(f"  dim(A_p):     {A_dims}")
print(f"  dim(Omega_p): {omega_dims}")

chi_omega = sum((-1)**p * d for p, d in enumerate(omega_dims))
chi_A = sum((-1)**p * d for p, d in enumerate(A_dims))
chi_betti = sum((-1)**p * b for p, b in enumerate(betti))

print(f"  chi(Omega) = {chi_omega}")
print(f"  chi(A)     = {chi_A}")
print(f"  chi(betti) = {chi_betti}")
print(f"  1 - beta_1 = {1 - (betti[1] if len(betti) > 1 else 0)}")

if chi_omega != chi_betti:
    print(f"\n  WARNING: chi(Omega) != chi(betti) -- possible numerical error!")

# Detailed rank analysis
print(f"\nDetailed rank analysis:")
for p in range(n):
    omega_p = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
    dim_omega = omega_p.shape[1] if omega_p.ndim == 2 else 0

    # rank of d_p
    if dim_omega > 0 and len(allowed[p-1]) > 0:
        D = build_full_boundary_matrix(allowed[p], allowed[p-1])
        D_omega = D @ omega_p
        S = np.linalg.svd(D_omega, compute_uv=False)
        rank = sum(s > 1e-8 for s in S)
        # Also check with tighter threshold
        rank_tight = sum(s > 1e-6 for s in S)
        min_sv = min(S) if len(S) > 0 else 0
        print(f"  p={p}: dim(O)={dim_omega}, rank(d_p)={rank} (tight:{rank_tight}), "
              f"min_sv={min_sv:.2e}, dim(A)={A_dims[p]}")
        if rank != rank_tight:
            print(f"    *** THRESHOLD SENSITIVE! Some singular values near 1e-8 ***")
            # Show suspicious values
            suspicious = [s for s in S if 1e-10 < s < 1e-5]
            if suspicious:
                print(f"    Suspicious SVs: {suspicious}")
    else:
        print(f"  p={p}: dim(O)={dim_omega}, rank(d_p)=0, dim(A)={A_dims[p]}")


# Try with tighter threshold (1e-6)
print(f"\n\nRecompute betti with threshold 1e-6:")
# Manual computation with adjustable threshold
for threshold in [1e-8, 1e-6, 1e-10]:
    omega = {}
    for p in range(n + 1):
        if p < n:
            omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        else:
            omega[p] = np.zeros((0, 0))

    betti_t = []
    for p in range(n):
        dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_omega_p == 0:
            betti_t.append(0)
            continue

        bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_omega = bd_p @ omega[p]
        if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
            S_p = np.linalg.svd(bd_p_omega, compute_uv=False)
            rank_p = sum(s > threshold for s in S_p)
        else:
            rank_p = 0
        ker_dim = dim_omega_p - rank_p

        dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 and p+1 < n else 0
        if dim_omega_p1 > 0:
            bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            bd_p1_omega = bd_p1 @ omega[p+1]
            S_p1 = np.linalg.svd(bd_p1_omega, compute_uv=False)
            im_dim = sum(s > threshold for s in S_p1)
        else:
            im_dim = 0

        betti_t.append(ker_dim - im_dim)

    chi_t = sum((-1)**p * b for p, b in enumerate(betti_t))
    print(f"  threshold={threshold}: beta={betti_t}, chi={chi_t}")


print("\n\nDone.")
