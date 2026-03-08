#!/usr/bin/env python3
"""
p11_beta8_eigenspace.py - Compute beta_8 of P_11 via eigenspace decomposition

Strategy: project FIRST into eigenspaces, THEN compute Omega and boundary.
This reduces matrix sizes by factor ~11 compared to global computation.

For beta_8 we need Om_7, Om_8, Om_9 (to compute ker d_8 - im d_9).

Prediction (HYP-212): beta_8 = 10 = (p-1) + 0 (delta=0 for prime p=11).
Per eigenspace: trivial (k=0) gives 0, each non-trivial (k=1..10) gives 1.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

# Import path homology utilities
_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, boundary_coeffs
)
sys.stdout = _saved

def build_circulant(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i + s) % n] = 1
    return A

def shift_path(path, n):
    """Shift path by +1 mod n."""
    return tuple((v + 1) % n for v in path)

def build_shift_matrix(paths, n):
    """Build the permutation matrix for sigma: i -> i+1 mod n."""
    idx = {tuple(p): i for i, p in enumerate(paths)}
    m = len(paths)
    P = np.zeros((m, m), dtype=float)
    for j, p in enumerate(paths):
        sp = shift_path(p, n)
        if sp in idx:
            P[idx[sp], j] = 1.0
    return P

def build_projector(sigma_mat, n, k):
    """Build projector P_k = (1/n) sum_{j=0}^{n-1} omega^{-kj} sigma^j."""
    m = sigma_mat.shape[0]
    omega = np.exp(2j * np.pi / n)
    P_k = np.zeros((m, m), dtype=complex)
    sig_j = np.eye(m, dtype=complex)
    for j in range(n):
        P_k += (omega**(-k*j)) * sig_j
        sig_j = sig_j @ sigma_mat
    P_k /= n
    return P_k

def project_basis(P_k, tol=1e-8):
    """Get orthonormal basis for range of projector P_k."""
    U, S, Vh = np.linalg.svd(P_k, full_matrices=False)
    rank = sum(s > tol for s in S)
    if rank == 0:
        return None, 0
    return U[:, :rank], rank

def build_boundary_matrix(paths_p, paths_pm1):
    """Build boundary matrix d_p: R^{A_p} -> R^{A_{p-1}}.
    Uses GLMY boundary: d(v0,...,vp) = sum_{i=0}^{p} (-1)^i (v0,...,vi-hat,...,vp)
    """
    idx_pm1 = {tuple(p): i for i, p in enumerate(paths_pm1)}
    m_p = len(paths_p)
    m_pm1 = len(paths_pm1)
    bd = np.zeros((m_pm1, m_p), dtype=float)
    for j, path in enumerate(paths_p):
        p = len(path) - 1  # path length
        for i in range(p + 1):
            face = tuple(path[:i] + path[i+1:])
            if face in idx_pm1:
                bd[idx_pm1[face], j] += (-1)**i
    return bd

def compute_omega_in_eigenspace(paths_p, paths_pm1, proj_basis_p, proj_basis_pm1, bd_global):
    """Compute Omega_p within an eigenspace.

    Omega_p = {u in A_p : d(u) in span(A_{p-1})}
    In the eigenspace: project both d(u) and A_{p-1} into eigenspace,
    then find u in eigenspace such that d(u) lies in A_{p-1}'s eigenspace projection.

    Returns the Omega basis (in eigenspace coords) and its dimension.
    """
    if proj_basis_p is None or proj_basis_p[1] == 0:
        return None, 0

    V_p, dim_p = proj_basis_p
    V_pm1, dim_pm1 = proj_basis_pm1 if proj_basis_pm1 is not None else (None, 0)

    # Boundary in eigenspace: V_pm1^H @ bd_global @ V_p
    # This maps eigenspace_p -> eigenspace_{p-1}
    # The "junk" is the component of bd(u) outside eigenspace_{p-1}

    # Full boundary of eigenspace vectors: bd_global @ V_p (in full A_{p-1} space)
    bd_proj = bd_global @ V_p  # shape: (|A_{p-1}|, dim_p)

    if dim_pm1 == 0:
        # All boundary is "junk" -> Omega = {u : bd(u) = 0}
        S_vals = np.linalg.svd(bd_proj, compute_uv=False)
        rank = sum(s > 1e-8 for s in S_vals)
        omega_dim = dim_p - rank
        if omega_dim == 0:
            return None, 0
        # Get kernel
        U, S_v, Vh = np.linalg.svd(bd_proj, full_matrices=False)
        kernel_basis = Vh[rank:].conj().T  # shape: (dim_p, omega_dim)
        return kernel_basis, omega_dim

    # Project boundary onto A_{p-1} eigenspace and its complement
    # Component in eigenspace: V_pm1 @ V_pm1^H @ bd_proj
    # Component outside (junk): (I - V_pm1 @ V_pm1^H) @ bd_proj
    junk = bd_proj - V_pm1 @ (V_pm1.conj().T @ bd_proj)  # shape: (|A_{p-1}|, dim_p)

    # Omega = kernel of junk map
    S_vals = np.linalg.svd(junk, compute_uv=False)
    rank = sum(s > 1e-8 for s in S_vals)
    omega_dim = dim_p - rank

    if omega_dim == 0:
        return None, 0

    U, S_v, Vh = np.linalg.svd(junk, full_matrices=False)
    kernel_basis = Vh[rank:].conj().T  # shape: (dim_p, omega_dim)
    return kernel_basis, omega_dim


print("=" * 70)
print("P_11 EIGENSPACE DECOMPOSITION: beta_8")
print("=" * 70)

n = 11
S = {1, 3, 4, 5, 9}  # QR(11)
A = build_circulant(n, S)

t0 = time.time()

# Enumerate allowed paths for dims 6-9 (need 6 for Om_7 junk)
print("\n--- Allowed paths ---")
allowed = {}
for p in range(10):
    t1 = time.time()
    allowed[p] = enumerate_allowed_paths(A, n, p)
    print(f"  A_{p}: {len(allowed[p])} ({time.time()-t1:.1f}s)")
    if not allowed[p]:
        break

# Build shift matrices for relevant dimensions
print("\n--- Shift matrices ---")
sigma = {}
for p in [6, 7, 8, 9]:
    if p in allowed and allowed[p]:
        t1 = time.time()
        sigma[p] = build_shift_matrix(allowed[p], n)
        print(f"  sigma_{p}: {sigma[p].shape} ({time.time()-t1:.1f}s)")

# Build global boundary matrices
print("\n--- Boundary matrices ---")
bd = {}
for p in [7, 8, 9]:
    if p in allowed and allowed[p] and (p-1) in allowed and allowed[p-1]:
        t1 = time.time()
        bd[p] = build_boundary_matrix(allowed[p], allowed[p-1])
        print(f"  d_{p}: {bd[p].shape} ({time.time()-t1:.1f}s)")

# For each eigenspace k, compute beta_8
print("\n" + "=" * 70)
print("PER-EIGENSPACE COMPUTATION")
print("=" * 70)

beta8_per_k = []
for k in range(n):
    t_k = time.time()

    # Build projectors
    proj_basis = {}
    for p in [6, 7, 8, 9]:
        if p not in sigma:
            proj_basis[p] = (None, 0)
            continue
        P_k = build_projector(sigma[p], n, k)
        basis, dim = project_basis(P_k)
        proj_basis[p] = (basis, dim)

    dims = [proj_basis[p][1] for p in [6, 7, 8, 9]]

    # Compute Omega_7, Omega_8, Omega_9 in eigenspace k
    omega = {}
    for p in [7, 8, 9]:
        if p not in bd or proj_basis[p][1] == 0:
            omega[p] = (None, 0)
            continue
        om_basis, om_dim = compute_omega_in_eigenspace(
            allowed[p], allowed[p-1],
            proj_basis[p], proj_basis[p-1],
            bd[p]
        )
        omega[p] = (om_basis, om_dim)

    om_dims = [omega[p][1] for p in [7, 8, 9]]

    # Compute beta_8: ker(d_8 on Om_8) - im(d_9 from Om_9)
    if omega[8][1] == 0:
        beta8_k = 0
    else:
        V_p8 = proj_basis[8][0]  # eigenspace basis for A_8
        om8_basis = omega[8][0]   # Omega_8 basis in eigenspace coords

        # d_8 restricted to Omega_8 in eigenspace:
        # d_8 maps A_8 -> A_7. In eigenspace: V_7^H @ bd[8] @ V_8 @ om8_basis
        V_p7 = proj_basis[7][0]
        if V_p7 is not None and V_p8 is not None:
            # Boundary of Omega_8 elements (in eigenspace 7 coords)
            bd8_om = V_p7.conj().T @ bd[8] @ V_p8 @ om8_basis
            S_vals = np.linalg.svd(bd8_om, compute_uv=False)
            rk8 = sum(s > 1e-6 for s in S_vals)
        else:
            rk8 = 0
        ker8 = omega[8][1] - rk8

        # im(d_9) from Omega_9 into eigenspace 8:
        if omega[9][1] == 0:
            im9 = 0
        else:
            V_p9 = proj_basis[9][0]
            om9_basis = omega[9][0]
            if V_p8 is not None and V_p9 is not None:
                bd9_om = V_p8.conj().T @ bd[9] @ V_p9 @ om9_basis
                S_vals = np.linalg.svd(bd9_om, compute_uv=False)
                im9 = sum(s > 1e-6 for s in S_vals)
            else:
                im9 = 0

        beta8_k = max(0, int(round(ker8 - im9)))

    beta8_per_k.append(beta8_k)
    dt = time.time() - t_k
    print(f"  k={k}: eig_dims={dims}, Om_dims={om_dims}, beta_8={beta8_k} ({dt:.1f}s)")

print(f"\n  SUM beta_8 = {sum(beta8_per_k)}")
print(f"  Prediction (HYP-212, prime p=11): beta_8 = 10 = (p-1) + 0")
print(f"  Per-eigenspace prediction: k=0 gives 0, each k=1..10 gives 1")

# Check conjugation symmetry
print(f"\n  Conjugation symmetry:")
for k in range(1, 6):
    print(f"    beta_8({k}) = {beta8_per_k[k]}, beta_8({11-k}) = {beta8_per_k[11-k]}, "
          f"equal={beta8_per_k[k]==beta8_per_k[11-k]}")

print(f"\nTotal time: {time.time()-t0:.1f}s")
