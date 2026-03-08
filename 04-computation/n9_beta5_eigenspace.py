#!/usr/bin/env python3
"""
n9_beta5_eigenspace.py - Focused: eigenspace decomposition of beta_5=10 only

Only needs Om_4, Om_5, Om_6 (and boundaries d_5, d_6).
Skips the expensive Om_7, Om_8 computation.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_circulant(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S:
                A[i][j] = 1
    return A

def shift_perm_matrix(paths, n):
    idx = {tuple(p): i for i, p in enumerate(paths)}
    m = len(paths)
    P = np.zeros((m, m), dtype=float)
    for j, p in enumerate(paths):
        shifted = tuple((v + 1) % n for v in p)
        if shifted in idx:
            P[idx[shifted], j] = 1.0
    return P

print("=" * 70)
print("FOCUSED: beta_5 EIGENSPACE DECOMPOSITION (Z_9 maximizer)")
print("=" * 70)

n = 9
S = {1, 5, 6, 7}
A = build_circulant(n, S)
t0 = time.time()

# Only need dims 4, 5, 6
dims_needed = [4, 5, 6]
allowed = {}
for p in range(-1, 7):
    if p < 0:
        allowed[p] = []
    else:
        t1 = time.time()
        allowed[p] = enumerate_allowed_paths(A, n, p)
        print(f"  A_{p}: {len(allowed[p])} ({time.time()-t1:.1f}s)")

# Omega bases for dims 4, 5, 6
omega_basis = {}
for p in dims_needed:
    t1 = time.time()
    basis = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
    omega_basis[p] = basis
    dim = basis.shape[1] if basis.ndim == 2 else 0
    print(f"  Om_{p}: dim={dim} ({time.time()-t1:.1f}s)")

# Boundary matrices: d_5 (Om_5 -> A_4) and d_6 (Om_6 -> A_5)
print("\n--- Boundary matrices ---")
bd_5 = build_full_boundary_matrix(allowed[5], allowed[4]) @ omega_basis[5]
rk_5 = np.linalg.matrix_rank(bd_5, tol=1e-8)
print(f"  d_5: rank={rk_5}")

bd_6 = build_full_boundary_matrix(allowed[6], allowed[5]) @ omega_basis[6]
rk_6 = np.linalg.matrix_rank(bd_6, tol=1e-8)
print(f"  d_6: rank={rk_6}")

dim_5 = omega_basis[5].shape[1]
ker_5 = dim_5 - rk_5
beta_5 = ker_5 - rk_6
print(f"\n  beta_5 = ker(d_5) - im(d_6) = {ker_5} - {rk_6} = {beta_5}")

# ===== Eigenspace decomposition =====
print("\n" + "=" * 70)
print("EIGENSPACE DECOMPOSITION OF beta_5")
print("=" * 70)

omega_9 = np.exp(2j * np.pi / 9)

# Sigma on Omega_5 and Omega_6
sigma_Om = {}
for p in [5, 6]:
    sigma_A = shift_perm_matrix(allowed[p], n)
    om = omega_basis[p]
    sigma_Om[p] = om.T @ sigma_A @ om

print(f"\n  Om_5 dim = {dim_5}, Om_5/9 = {dim_5/9:.1f}")
dim_6 = omega_basis[6].shape[1]
print(f"  Om_6 dim = {dim_6}, Om_6/9 = {dim_6/9:.1f}")

beta5_per_k = []
for k in range(n):
    m5 = dim_5
    sig5 = sigma_Om[5]
    P5_k = np.zeros((m5, m5), dtype=complex)
    sig_j = np.eye(m5, dtype=complex)
    for j in range(n):
        P5_k += (omega_9**(-k*j)) * sig_j
        sig_j = sig_j @ sig5
    P5_k /= n
    proj_dim_5 = round(np.real(np.trace(P5_k)))

    m6 = dim_6
    sig6 = sigma_Om[6]
    P6_k = np.zeros((m6, m6), dtype=complex)
    sig_j = np.eye(m6, dtype=complex)
    for j in range(n):
        P6_k += (omega_9**(-k*j)) * sig_j
        sig_j = sig_j @ sig6
    P6_k /= n
    proj_dim_6 = round(np.real(np.trace(P6_k)))

    # ker(d_5) in eigenspace k
    bd5_proj = bd_5.astype(complex) @ P5_k
    S_vals = np.linalg.svd(bd5_proj, compute_uv=False)
    rk5_k = sum(s > 1e-6 for s in S_vals)
    ker5_k = proj_dim_5 - rk5_k

    # im(d_6) in eigenspace k
    bd6_proj = bd_6.astype(complex) @ P6_k
    S_vals = np.linalg.svd(bd6_proj, compute_uv=False)
    rk6_k = sum(s > 1e-6 for s in S_vals)

    b5_k = max(0, int(round(ker5_k - rk6_k)))
    beta5_per_k.append(b5_k)
    print(f"  k={k}: Om5_dim={proj_dim_5}, Om6_dim={proj_dim_6}, "
          f"ker(d5)={ker5_k}, im(d6)={rk6_k}, beta5={b5_k}")

print(f"\n  SUM of per-eigenspace beta_5 = {sum(beta5_per_k)}")
print(f"  Expected: 10")

# Check conjugation symmetry: beta(k) should equal beta(9-k)
print(f"\n  Conjugation symmetry check:")
for k in range(1, 5):
    print(f"    b5({k}) = {beta5_per_k[k]}, b5({9-k}) = {beta5_per_k[9-k]}, "
          f"equal={beta5_per_k[k]==beta5_per_k[9-k]}")

# Identify order-3 vs order-9 characters
print(f"\n  By character order:")
print(f"    trivial (k=0): beta5 = {beta5_per_k[0]}")
print(f"    order 3 (k=3,6): beta5 = {beta5_per_k[3]}, {beta5_per_k[6]}")
print(f"    order 9 (k=1,2,4,5,7,8): beta5 = {[beta5_per_k[k] for k in [1,2,4,5,7,8]]}")

print(f"\nTotal time: {time.time()-t0:.1f}s")
