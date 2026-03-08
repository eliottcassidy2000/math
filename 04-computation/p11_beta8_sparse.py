#!/usr/bin/env python3
"""
p11_beta8_sparse.py - P_11 beta_8 via orbit-based eigenspace decomposition

Since n=11 is prime, every path orbit under sigma has size 11.
Per-eigenspace dimension = |A_p|/11.
Use sparse boundary matrices and orbit-based Fourier projection.

For beta_8: need Om_7^(k), Om_8^(k), Om_9^(k) in each eigenspace k.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import svds
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import enumerate_allowed_paths
sys.stdout = _saved

def build_circulant(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i + s) % n] = 1
    return A

def shift_path(path, n):
    return tuple((v + 1) % n for v in path)

def find_orbits(paths, n):
    """Find orbits under cyclic shift. Return list of orbit representatives."""
    idx = {tuple(p): i for i, p in enumerate(paths)}
    visited = set()
    orbits = []
    for i, p in enumerate(paths):
        if i in visited:
            continue
        orbit = [i]
        visited.add(i)
        q = tuple(p)
        for _ in range(n - 1):
            q = shift_path(q, n)
            j = idx.get(q)
            if j is not None and j not in visited:
                orbit.append(j)
                visited.add(j)
        orbits.append(orbit)
    return orbits

def build_eigenspace_basis(orbits, n_paths, n, k):
    """Build the eigenspace basis for eigenvalue omega^k.

    For orbit [i0, i1, ..., i_{n-1}] where i_j = sigma^j(i_0),
    the eigenvector is (1/sqrt(n)) * sum_j omega^{kj} e_{i_j}.
    """
    omega = np.exp(2j * np.pi / n)
    n_orbits = len(orbits)
    # Build basis as dense matrix: n_paths x n_orbits
    V = np.zeros((n_paths, n_orbits), dtype=complex)
    for r, orbit in enumerate(orbits):
        for j, idx in enumerate(orbit):
            V[idx, r] = omega**(k * j) / np.sqrt(n)
    return V

def build_sparse_boundary(paths_p, paths_pm1):
    """Build sparse boundary matrix d_p: R^{|A_p|} -> R^{|A_{p-1}|}."""
    idx_pm1 = {tuple(p): i for i, p in enumerate(paths_pm1)}
    rows, cols, vals = [], [], []
    for j, path in enumerate(paths_p):
        p = len(path) - 1
        for i in range(p + 1):
            face = tuple(path[:i] + path[i+1:])
            if face in idx_pm1:
                rows.append(idx_pm1[face])
                cols.append(j)
                vals.append((-1.0)**i)
    return sparse.csr_matrix((vals, (rows, cols)),
                              shape=(len(paths_pm1), len(paths_p)))

def compute_omega_and_betti_in_eigenspace(V_p, V_pm1, V_pp1, bd_p, bd_pp1):
    """Compute beta_p in eigenspace.

    V_p: eigenspace basis for A_p (|A_p| x dim_k)
    V_pm1: eigenspace basis for A_{p-1}
    V_pp1: eigenspace basis for A_{p+1}
    bd_p: sparse boundary d_p: A_p -> A_{p-1}
    bd_pp1: sparse boundary d_{p+1}: A_{p+1} -> A_p

    Returns: (Om_p_dim, beta_p)
    """
    dim_k_p = V_p.shape[1] if V_p is not None else 0
    if dim_k_p == 0:
        return 0, 0

    dim_k_pm1 = V_pm1.shape[1] if V_pm1 is not None else 0

    # Boundary of eigenspace vectors: bd_p @ V_p  (|A_{p-1}| x dim_k_p)
    bd_V = bd_p @ V_p if bd_p is not None else np.zeros((0, dim_k_p))

    # Project boundary onto complement of A_{p-1} eigenspace = "junk"
    if dim_k_pm1 > 0:
        # Component in eigenspace: V_pm1 @ V_pm1^H @ bd_V
        proj = V_pm1 @ (V_pm1.conj().T @ bd_V)
        junk = bd_V - proj
    else:
        junk = bd_V

    # Omega = kernel of junk map
    if junk.shape[0] == 0:
        # No boundary constraint -> all of eigenspace is Omega
        om_dim = dim_k_p
        om_basis = np.eye(dim_k_p, dtype=complex)
    else:
        S_vals = np.linalg.svd(junk, compute_uv=False)
        junk_rank = sum(s > 1e-8 for s in S_vals)
        om_dim = dim_k_p - junk_rank
        if om_dim == 0:
            return 0, 0
        _, _, Vh = np.linalg.svd(junk, full_matrices=True)
        # Last om_dim rows of Vh are the kernel
        om_basis = Vh[junk_rank:].conj().T  # (dim_k_p, om_dim)

    # d_p restricted to Om_p in eigenspace
    # Boundary of Om_p vectors in eigenspace_{p-1} coords
    if dim_k_pm1 > 0:
        bd_om_proj = V_pm1.conj().T @ bd_p @ V_p @ om_basis
        S_vals = np.linalg.svd(bd_om_proj, compute_uv=False)
        rk_dp = sum(s > 1e-6 for s in S_vals)
    else:
        rk_dp = 0
    ker_dp = om_dim - rk_dp

    # im(d_{p+1}) from Om_{p+1}
    dim_k_pp1 = V_pp1.shape[1] if V_pp1 is not None else 0
    if dim_k_pp1 == 0 or bd_pp1 is None:
        im_dp1 = 0
    else:
        # Compute Om_{p+1} too (need junk map for p+1)
        bd_V_pp1 = bd_pp1 @ V_pp1
        if dim_k_p > 0:
            proj_pp1 = V_p @ (V_p.conj().T @ bd_V_pp1)
            junk_pp1 = bd_V_pp1 - proj_pp1
        else:
            junk_pp1 = bd_V_pp1

        S_vals_pp1 = np.linalg.svd(junk_pp1, compute_uv=False)
        junk_rank_pp1 = sum(s > 1e-8 for s in S_vals_pp1)
        om_dim_pp1 = dim_k_pp1 - junk_rank_pp1

        if om_dim_pp1 == 0:
            im_dp1 = 0
        else:
            _, _, Vh_pp1 = np.linalg.svd(junk_pp1, full_matrices=True)
            om_basis_pp1 = Vh_pp1[junk_rank_pp1:].conj().T

            # Image of Om_{p+1} in Om_p via d_{p+1}
            bd_om_pp1 = V_p.conj().T @ bd_pp1 @ V_pp1 @ om_basis_pp1
            S_vals_im = np.linalg.svd(bd_om_pp1, compute_uv=False)
            im_dp1 = sum(s > 1e-6 for s in S_vals_im)

    beta_p = max(0, int(round(ker_dp - im_dp1)))
    return om_dim, beta_p


print("=" * 70)
print("P_11 BETA_8 (SPARSE ORBIT-BASED EIGENSPACE)")
print("=" * 70)

n = 11
S = {1, 3, 4, 5, 9}
A = build_circulant(n, S)
t0 = time.time()

# Step 1: Enumerate allowed paths
print("\n--- Allowed paths ---")
allowed = {}
for p in range(10):
    t1 = time.time()
    allowed[p] = enumerate_allowed_paths(A, n, p)
    print(f"  A_{p}: {len(allowed[p])} ({time.time()-t1:.1f}s)")
    if not allowed[p]:
        break

# Step 2: Find orbits
print("\n--- Orbits ---")
orbits = {}
for p in [6, 7, 8, 9]:
    if p in allowed and allowed[p]:
        t1 = time.time()
        orbits[p] = find_orbits(allowed[p], n)
        print(f"  dim {p}: {len(orbits[p])} orbits (expect {len(allowed[p])}/{n} = {len(allowed[p])//n})")
        assert len(orbits[p]) == len(allowed[p]) // n, f"Not all orbits have size {n}!"

# Step 3: Build sparse boundary matrices
print("\n--- Sparse boundary matrices ---")
bd = {}
for p in [7, 8, 9]:
    if p in allowed and allowed[p] and (p-1) in allowed:
        t1 = time.time()
        bd[p] = build_sparse_boundary(allowed[p], allowed[p-1])
        nnz = bd[p].nnz
        print(f"  d_{p}: {bd[p].shape}, nnz={nnz} ({time.time()-t1:.1f}s)")

# Step 4: Per-eigenspace beta_8
print("\n" + "=" * 70)
print("PER-EIGENSPACE beta_8")
print("=" * 70)

beta8_per_k = []
for k in range(n):
    t_k = time.time()

    # Build eigenspace bases
    V = {}
    for p in [6, 7, 8, 9]:
        if p in orbits:
            V[p] = build_eigenspace_basis(orbits[p], len(allowed[p]), n, k)
        else:
            V[p] = None

    dims = {p: V[p].shape[1] if V[p] is not None else 0 for p in [6, 7, 8, 9]}

    # Compute beta_8
    om_dim, beta8_k = compute_omega_and_betti_in_eigenspace(
        V.get(8), V.get(7), V.get(9),
        bd.get(8), bd.get(9)
    )

    beta8_per_k.append(beta8_k)
    dt = time.time() - t_k
    print(f"  k={k}: dims={[dims[p] for p in [6,7,8,9]]}, Om_8={om_dim}, "
          f"beta_8={beta8_k} ({dt:.1f}s)")

print(f"\n  SUM beta_8 = {sum(beta8_per_k)}")
print(f"  Prediction (HYP-212, prime p=11): beta_8 = 10 = (p-1) + 0")
print(f"  Trivial (k=0): beta_8 = {beta8_per_k[0]}")
print(f"  Non-trivial (k=1..10): beta_8 = {beta8_per_k[1:]}")

# Conjugation symmetry
print(f"\n  Conjugation symmetry:")
for k in range(1, 6):
    print(f"    beta_8({k}) = {beta8_per_k[k]}, beta_8({11-k}) = {beta8_per_k[11-k]}, "
          f"equal={beta8_per_k[k]==beta8_per_k[11-k]}")

print(f"\nTotal time: {time.time()-t0:.1f}s")
