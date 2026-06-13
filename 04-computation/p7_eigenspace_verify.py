#!/usr/bin/env python3
"""
p7_eigenspace_verify.py - Verify eigenspace decomposition code on P_7

Known result (opus): P_7 has beta = [1,0,0,0,6,0,0].
Each non-trivial eigenspace (k=1..6) contributes beta_4=1.
Trivial (k=0) contributes beta_0=1 only.

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
print("EIGENSPACE VERIFY: P_7 (Paley tournament)")
print("=" * 70)

n = 7
S = {1, 2, 4}  # QR of GF(7)
A = build_circulant(n, S)
t0 = time.time()

max_dim = 6

# Allowed paths
allowed = {}
for p in range(-1, max_dim + 2):
    if p < 0:
        allowed[p] = []
    else:
        allowed[p] = enumerate_allowed_paths(A, n, p)
        print(f"  A_{p}: {len(allowed[p])}")
        if not allowed[p]:
            max_dim = min(max_dim, p - 1)
            break

# Omega bases
omega_basis = {}
for p in range(max_dim + 2):
    if p not in allowed or not allowed[p]:
        omega_basis[p] = np.zeros((0, 0))
        continue
    basis = compute_omega_basis(A, n, p, allowed[p], allowed[p-1] if p-1 in allowed else [])
    omega_basis[p] = basis
    dim = basis.shape[1] if basis.ndim == 2 else 0
    print(f"  Om_{p}: dim={dim}")

# Boundary matrices in Omega coords
bd_omega = {}
for p in range(1, max_dim + 2):
    dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
    if dim_p == 0:
        continue
    bd = build_full_boundary_matrix(allowed[p], allowed[p-1] if p-1 in allowed else [])
    bd_om = bd @ omega_basis[p]
    bd_omega[p] = bd_om

# Global Betti
betti_global = []
for p in range(max_dim + 1):
    dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
    if dim_p == 0:
        betti_global.append(0)
        continue
    if p not in bd_omega:
        ker = dim_p
    else:
        S_vals = np.linalg.svd(bd_omega[p], compute_uv=False)
        ker = dim_p - sum(s > 1e-8 for s in S_vals)
    if p+1 not in bd_omega:
        im = 0
    else:
        S_vals = np.linalg.svd(bd_omega[p+1], compute_uv=False)
        im = sum(s > 1e-8 for s in S_vals)
    betti_global.append(ker - im)

print(f"\nGlobal b = {betti_global}")
print(f"Expected: [1, 0, 0, 0, 6, 0, 0]")

# Eigenspace decomposition
print("\n--- Eigenspace decomposition ---")
omega_n = np.exp(2j * np.pi / n)

# Sigma on Omega
sigma_Om = {}
for p in range(max_dim + 2):
    if p not in allowed or not allowed[p]:
        continue
    dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
    if dim_p == 0:
        continue
    sigma_A = shift_perm_matrix(allowed[p], n)
    om = omega_basis[p]
    sigma_Om[p] = om.T @ sigma_A @ om

for k in range(n):
    betti_k = []
    for p in range(max_dim + 1):
        dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
        if dim_p == 0:
            betti_k.append(0)
            continue

        sig = sigma_Om[p]
        m = dim_p
        P_k = np.zeros((m, m), dtype=complex)
        sig_j = np.eye(m, dtype=complex)
        for j in range(n):
            P_k += (omega_n**(-k*j)) * sig_j
            sig_j = sig_j @ sig
        P_k /= n
        proj_dim = round(np.real(np.trace(P_k)))

        if proj_dim == 0:
            betti_k.append(0)
            continue

        if p not in bd_omega:
            ker_k = proj_dim
        else:
            bd_proj = bd_omega[p].astype(complex) @ P_k
            S_vals = np.linalg.svd(bd_proj, compute_uv=False)
            rank_k = sum(s > 1e-6 for s in S_vals)
            ker_k = proj_dim - rank_k

        if p+1 not in bd_omega or p+1 not in sigma_Om:
            im_k = 0
        else:
            dim_p1 = omega_basis[p+1].shape[1]
            sig_p1 = sigma_Om[p+1]
            P_k_p1 = np.zeros((dim_p1, dim_p1), dtype=complex)
            sig_j = np.eye(dim_p1, dtype=complex)
            for j in range(n):
                P_k_p1 += (omega_n**(-k*j)) * sig_j
                sig_j = sig_j @ sig_p1
            P_k_p1 /= n
            bd_next_proj = bd_omega[p+1].astype(complex) @ P_k_p1
            S_vals = np.linalg.svd(bd_next_proj, compute_uv=False)
            im_k = sum(s > 1e-6 for s in S_vals)

        betti_k.append(max(0, int(round(ker_k - im_k))))

    print(f"  k={k}: b={betti_k}")

print(f"\nExpected: k=0 -> [1,0,0,0,0,0,0], k=1..6 -> [0,0,0,0,1,0,0]")
print(f"Total time: {time.time()-t0:.1f}s")
