#!/usr/bin/env python3
"""
n9_eigenspace_v3.py - Correct eigenspace decomposition using path_homology_v2

Uses the proper GLMY Omega computation (null space of junk map).
Decomposes b5=10 across Z/9Z eigenspaces for the circulant maximizer.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

# Suppress path_homology_v2's validation output
_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, boundary_coeffs
)
sys.stdout = _saved

def build_circulant(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S:
                A[i][j] = 1
    return A

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])

def shift_perm_matrix(paths, n):
    """Permutation matrix of sigma: i -> i+1 mod n on a list of paths."""
    idx = {tuple(p): i for i, p in enumerate(paths)}
    m = len(paths)
    P = np.zeros((m, m), dtype=float)
    for j, p in enumerate(paths):
        shifted = tuple((v + 1) % n for v in p)
        if shifted in idx:
            P[idx[shifted], j] = 1.0
    return P

# ===== Main =====
print("=" * 70)
print("EIGENSPACE DECOMPOSITION v3: n=9 CIRCULANT MAXIMIZER")
print("=" * 70)

n = 9
S = {1, 5, 6, 7}
A = build_circulant(n, S)
H = H_tournament(A, n)
print(f"Z_{n}, S={sorted(S)}, H={H}")
assert H == 3357

t0 = time.time()
max_dim = 7  # need up to dim 6 for b5 (since we need d6 for im at dim 5)

# Step 1: Enumerate allowed paths
print("\n--- Step 1: Allowed paths ---")
allowed = {}
for p in range(-1, max_dim + 2):
    if p < 0:
        allowed[p] = []
    else:
        t1 = time.time()
        allowed[p] = enumerate_allowed_paths(A, n, p)
        print(f"  A_{p}: {len(allowed[p])} paths ({time.time()-t1:.1f}s)")
        if not allowed[p]:
            max_dim = min(max_dim, p - 1)
            break

# Step 2: Compute Omega bases
print("\n--- Step 2: Omega bases ---")
omega_basis = {}
for p in range(max_dim + 2):
    if p not in allowed or not allowed[p]:
        omega_basis[p] = np.zeros((0, 0))
        continue
    t1 = time.time()
    basis = compute_omega_basis(A, n, p, allowed[p], allowed[p-1] if p-1 in allowed else [])
    omega_basis[p] = basis
    dim = basis.shape[1] if basis.ndim == 2 else 0
    print(f"  Om_{p}: dim={dim} (from {len(allowed[p])} paths) ({time.time()-t1:.1f}s)")

# Step 3: Build boundary matrices (in Omega coordinates)
print("\n--- Step 3: Boundary matrices ---")
bd_omega = {}  # bd_omega[p] maps Om_p -> R^{A_{p-1}} (or equivalently into Om_{p-1})
for p in range(1, max_dim + 2):
    dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
    if dim_p == 0:
        continue
    bd = build_full_boundary_matrix(allowed[p], allowed[p-1] if p-1 in allowed else [])
    bd_om = bd @ omega_basis[p]  # map from Om_p coords to A_{p-1} coords
    S_vals = np.linalg.svd(bd_om, compute_uv=False)
    rank = sum(s > 1e-8 for s in S_vals)
    bd_omega[p] = bd_om
    print(f"  d_{p}: rank={rank}")

# Step 4: Compute global Betti
print("\n--- Step 4: Global Betti ---")
betti_global = []
for p in range(max_dim + 1):
    dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
    if dim_p == 0:
        betti_global.append(0)
        continue

    # ker(d_p)
    if p not in bd_omega:
        ker = dim_p
    else:
        S_vals = np.linalg.svd(bd_omega[p], compute_uv=False)
        ker = dim_p - sum(s > 1e-8 for s in S_vals)

    # im(d_{p+1})
    if p+1 not in bd_omega:
        im = 0
    else:
        S_vals = np.linalg.svd(bd_omega[p+1], compute_uv=False)
        im = sum(s > 1e-8 for s in S_vals)

    betti_global.append(ker - im)

print(f"  b = {betti_global}")

# Step 5: Eigenspace decomposition
print("\n" + "=" * 70)
print("EIGENSPACE DECOMPOSITION (Z/9Z)")
print("=" * 70)

omega_9 = np.exp(2j * np.pi / 9)

# For each p, build sigma on A_p and project onto Om_p
sigma_Ap = {}
sigma_Om = {}

for p in range(max_dim + 2):
    if p not in allowed or not allowed[p]:
        continue
    dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
    if dim_p == 0:
        continue

    # sigma on A_p (permutation matrix)
    sigma_A = shift_perm_matrix(allowed[p], n)
    sigma_Ap[p] = sigma_A

    # sigma on Om_p: omega_basis[p]^T @ sigma_A @ omega_basis[p]
    # (valid since omega_basis columns are orthonormal from SVD)
    om = omega_basis[p]
    sig_om = om.T @ sigma_A @ om
    sigma_Om[p] = sig_om

# For each eigenspace k, compute Betti
for k in range(n):
    betti_k = []
    om_dims = []

    for p in range(max_dim + 1):
        dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
        if dim_p == 0:
            betti_k.append(0)
            om_dims.append(0)
            continue

        # Projector onto eigenspace k in Om_p coords
        sig = sigma_Om[p]
        m = dim_p
        P_k = np.zeros((m, m), dtype=complex)
        sig_j = np.eye(m, dtype=complex)
        for j in range(n):
            P_k += (omega_9**(-k*j)) * sig_j
            sig_j = sig_j @ sig
        P_k /= n

        proj_dim = round(np.real(np.trace(P_k)))
        om_dims.append(proj_dim)

        if proj_dim == 0:
            betti_k.append(0)
            continue

        # ker(d_p) in eigenspace k
        if p not in bd_omega:
            ker_k = proj_dim
        else:
            # d_p maps Om_p -> A_{p-1}. Restrict to eigenspace k:
            # bd_omega[p] @ P_k (in Om_p coords)
            bd_proj = bd_omega[p].astype(complex) @ P_k
            S_vals = np.linalg.svd(bd_proj, compute_uv=False)
            rank_k = sum(s > 1e-6 for s in S_vals)
            ker_k = proj_dim - rank_k

        # im(d_{p+1}) in eigenspace k
        if p+1 not in bd_omega or p+1 not in sigma_Om:
            im_k = 0
        else:
            # d_{p+1} maps Om_{p+1} -> A_p.
            # Project source to eigenspace k of Om_{p+1}
            dim_p1 = omega_basis[p+1].shape[1]
            sig_p1 = sigma_Om[p+1]
            P_k_p1 = np.zeros((dim_p1, dim_p1), dtype=complex)
            sig_j = np.eye(dim_p1, dtype=complex)
            for j in range(n):
                P_k_p1 += (omega_9**(-k*j)) * sig_j
                sig_j = sig_j @ sig_p1
            P_k_p1 /= n

            bd_next_proj = bd_omega[p+1].astype(complex) @ P_k_p1
            S_vals = np.linalg.svd(bd_next_proj, compute_uv=False)
            im_k = sum(s > 1e-6 for s in S_vals)

        betti_k.append(max(0, int(round(ker_k - im_k))))

    if any(b > 0 for b in betti_k):
        print(f"  k={k}: Om_dims={om_dims}, b={betti_k}")

# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"Expected: b = [1,0,0,0,0,10,...] (from path_homology_v2)")
print(f"Got:      b = {betti_global}")
print(f"\nFor P_7: each non-trivial eigenspace (k=1..6) contributes b4=1.")
print(f"         trivial (k=0) contributes b0=1 only.")
print(f"\nFor Z_9 maximizer: how does b5=10 distribute across k=0..8?")
print(f"  9 eigenspaces. If uniform: 10/9 not integer. Must be non-uniform.")

print(f"\nTotal time: {time.time()-t0:.1f}s")
