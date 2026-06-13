#!/usr/bin/env python3
"""
paley11_omega.py — opus-2026-03-13-S71

Compute GLMY Ω dimensions for Paley P_11.
Uses sparse matrices and the JJ^T/J^TJ trick for rank computation.

Key conjecture to test: β_8(P_11) = 10 = p-1
"""

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
import time
import sys

p = 11
QR = {a*a % p for a in range(1, p)}
print(f"P_{p}: QR = {sorted(QR)}")

# Adjacency
adj = [set() for _ in range(p)]
for i in range(p):
    for s in QR:
        adj[i].add((i + s) % p)

def enumerate_paths(max_depth):
    """Enumerate ALL directed paths up to given depth."""
    all_paths = {0: [(v,) for v in range(p)]}

    for m in range(1, max_depth + 1):
        t0 = time.time()
        paths = []
        for path in all_paths[m-1]:
            last = path[-1]
            visited = set(path)
            for v in adj[last]:
                if v not in visited:
                    paths.append(path + (v,))
        all_paths[m] = paths
        t1 = time.time()
        print(f"  m={m}: {len(paths)} paths ({t1-t0:.1f}s)")

    return all_paths

def compute_omega_dim(all_paths, m):
    """Compute dim(Ω_m) using the junk map."""
    if m <= 1:
        return len(all_paths[m])

    paths_m = all_paths[m]
    paths_m1_set = set(all_paths[m-1])
    dim_Am = len(paths_m)

    # Find junk faces
    junk = {}
    junk_count = 0
    for path in paths_m:
        for i in range(m+1):
            face = path[:i] + path[i+1:]
            if face not in paths_m1_set and face not in junk:
                junk[face] = junk_count
                junk_count += 1

    if junk_count == 0:
        return dim_Am

    print(f"    Junk faces: {junk_count}, building sparse J matrix...")
    t0 = time.time()

    # Build sparse junk map
    J = lil_matrix((junk_count, dim_Am))
    for j, path in enumerate(paths_m):
        for i in range(m+1):
            face = path[:i] + path[i+1:]
            if face in junk:
                J[junk[face], j] += (-1)**i

    J_csr = J.tocsr()
    t1 = time.time()
    print(f"    J matrix built ({t1-t0:.1f}s)")

    # Rank computation via the smaller of JJ^T or J^TJ
    if junk_count <= dim_Am:
        print(f"    Computing JJ^T ({junk_count}×{junk_count})...")
        M = (J_csr @ J_csr.T).toarray().astype(float)
    else:
        print(f"    Computing J^TJ ({dim_Am}×{dim_Am})...")
        M = (J_csr.T @ J_csr).toarray().astype(float)

    t2 = time.time()
    rank_J = np.linalg.matrix_rank(M, tol=1e-8)
    t3 = time.time()
    print(f"    rank(J) = {rank_J} ({t3-t2:.1f}s)")

    return dim_Am - rank_J

def compute_boundary_rank(all_paths, omega_bases, m):
    """Compute rank of ∂_m: Ω_m → Ω_{m-1}."""
    if omega_bases[m].shape[1] == 0 or omega_bases[m-1].shape[1] == 0:
        return 0

    paths_m = all_paths[m]
    paths_m1 = all_paths[m-1]
    idx = {path: i for i, path in enumerate(paths_m1)}

    B = lil_matrix((len(paths_m1), len(paths_m)))
    for j, path in enumerate(paths_m):
        for i in range(m+1):
            face = path[:i] + path[i+1:]
            if face in idx:
                B[idx[face], j] += (-1)**i

    B_csr = B.tocsr()
    # B_omega = omega_{m-1}^T @ B @ omega_m
    B_omega = omega_bases[m-1].T @ B_csr.toarray() @ omega_bases[m]
    return np.linalg.matrix_rank(B_omega, tol=1e-8)

# ============================================================
print(f"\n{'='*70}")
print("COMPUTING Ω DIMENSIONS FOR P_11")
print("="*70)

# Phase 1: Enumerate paths
max_m = 8  # Try up to m=8 (where we expect β_8=10)
print(f"\nEnumerating paths up to m={max_m}+1:")
all_paths = enumerate_paths(max_m + 1)

# Phase 2: Compute Ω dims
print(f"\nComputing Ω dimensions:")
omega_dims = {}
for m in range(max_m + 2):
    if m not in all_paths:
        break
    t0 = time.time()
    dim = compute_omega_dim(all_paths, m)
    t1 = time.time()
    omega_dims[m] = dim
    div_11 = f"(={dim//11}×11)" if dim % 11 == 0 else ""
    print(f"  Ω_{m} = {dim} {div_11} ({t1-t0:.1f}s)")
    sys.stdout.flush()

print(f"\nSummary:")
om_list = [omega_dims.get(m, '?') for m in range(max_m + 2)]
print(f"  Ω = {om_list}")
om_div = [omega_dims[m]//11 if omega_dims.get(m, 0) % 11 == 0 else '?' for m in range(max_m + 2)]
print(f"  Ω/11 = {om_div}")

# Compare patterns
print(f"\n  P_7 Ω/7 = [1, 3, 6, 9, 9, 6, 3]")
print(f"  P_11 Ω/11 = {om_div}")

# Phase 3: If we have enough, compute boundary ranks and Betti
if len(omega_dims) >= 3:
    print(f"\n{'='*70}")
    print("COMPUTING BOUNDARY RANKS AND BETTI NUMBERS")
    print("="*70)

    # We need omega bases (not just dims) for the boundary computation
    # This is more expensive. Let's do it for the degrees we have.

    # Build omega bases
    omega_bases = {}
    for m in range(max_m + 2):
        if m not in all_paths:
            break
        paths_m = all_paths[m]
        dim_Am = len(paths_m)

        if m <= 1:
            omega_bases[m] = np.eye(dim_Am)
            continue

        paths_m1_set = set(all_paths[m-1])
        junk = {}
        junk_count = 0
        for path in paths_m:
            for i in range(m+1):
                face = path[:i] + path[i+1:]
                if face not in paths_m1_set and face not in junk:
                    junk[face] = junk_count
                    junk_count += 1

        if junk_count == 0:
            omega_bases[m] = np.eye(dim_Am)
            continue

        J = lil_matrix((junk_count, dim_Am))
        for j, path in enumerate(paths_m):
            for i in range(m+1):
                face = path[:i] + path[i+1:]
                if face in junk:
                    J[junk[face], j] += (-1)**i

        J_dense = J.toarray()
        U, s, Vh = np.linalg.svd(J_dense, full_matrices=True)
        rank_J = int(np.sum(s > 1e-10))
        if rank_J < dim_Am:
            omega_bases[m] = Vh[rank_J:].T
        else:
            omega_bases[m] = np.zeros((dim_Am, 0))

    # Compute boundary ranks
    print(f"\nBoundary ranks:")
    ranks = {}
    for m in range(1, max_m + 1):
        if m not in omega_bases or m-1 not in omega_bases:
            break
        if omega_bases[m].shape[1] == 0 or omega_bases[m-1].shape[1] == 0:
            ranks[m] = 0
        else:
            t0 = time.time()
            rk = compute_boundary_rank(all_paths, omega_bases, m)
            t1 = time.time()
            ranks[m] = rk
            print(f"  rk(∂_{m}) = {rk} ({t1-t0:.1f}s)")

    # Betti numbers
    print(f"\nBetti numbers:")
    for m in range(max_m + 1):
        if m not in omega_dims:
            break
        dim_m = omega_dims[m]
        rk_m = ranks.get(m, 0)
        rk_m1 = ranks.get(m+1, 0)
        beta = dim_m - rk_m - rk_m1
        print(f"  β_{m} = Ω_{m} - rk(∂_{m}) - rk(∂_{m+1}) = {dim_m} - {rk_m} - {rk_m1} = {beta}")

print("\nDONE.")
