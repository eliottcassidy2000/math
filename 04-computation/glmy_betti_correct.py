#!/usr/bin/env python3
"""
glmy_betti_correct.py — opus-2026-03-13-S71

Compute CORRECT GLMY path homology Betti numbers.

The GLMY chain complex:
  1. A_m = regular m-paths (our basis)
  2. Ω_m = {u ∈ span(A_m) : ∂u ∈ span(A_{m-1})}
     This is a SUBSPACE — linear combinations of regular paths
     whose boundary lies entirely in the span of regular (m-1)-paths.
  3. ∂ restricted to Ω_m → Ω_{m-1} is the GLMY boundary.
  4. Betti = dim(ker ∂_m) - dim(im ∂_{m+1}).

Key: Ω_m is NOT just the set of individual paths with all faces regular.
It includes linear combinations that CANCEL non-regular face components.
"""

import numpy as np

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def transitive_tournament(n):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def get_regular_paths(A, m):
    n = A.shape[0]
    paths = []
    def dfs(path, depth, prev):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v in path: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path.append(v)
            dfs(path, depth+1, last)
            path.pop()
    for start in range(n):
        dfs([start], 0, -1)
    return paths

def compute_glmy_betti(A, max_dim=None):
    """Compute correct GLMY path homology Betti numbers."""
    n = A.shape[0]
    if max_dim is None:
        max_dim = n - 1

    # Step 1: Regular paths
    A_paths = {}
    for m in range(max_dim + 2):  # need m+1 for boundary
        A_paths[m] = get_regular_paths(A, m)

    # Step 2: Build Ω_m as kernel of the "junk map"
    # For each m-path p, ∂p might have faces not in A_{m-1}.
    # Let J_{m-1} = set of non-regular (m-1)-paths that appear as faces.
    # The "junk map" sends p ↦ projection of ∂p onto J_{m-1}.
    # Ω_m = kernel of this junk map.

    omega_basis = {}  # omega_basis[m] = matrix whose columns are basis for Ω_m (in A_m coords)

    for m in range(max_dim + 2):
        dim_Am = len(A_paths[m])
        if m <= 1:
            # Ω_0 = A_0 (vertices, no faces to worry about)
            # Ω_1 = A_1 (faces of 1-paths are vertices, always in A_0)
            omega_basis[m] = np.eye(dim_Am)
            continue

        A_set_prev = set(A_paths[m-1])

        # Find all non-regular faces
        junk_faces = {}  # junk_face -> index
        junk_count = 0

        for path in A_paths[m]:
            path_len = len(path) - 1  # = m
            for i in range(path_len + 1):
                face = path[:i] + path[i+1:]
                if face not in A_set_prev and face not in junk_faces:
                    junk_faces[face] = junk_count
                    junk_count += 1

        if junk_count == 0:
            # All faces are regular → Ω_m = A_m
            omega_basis[m] = np.eye(dim_Am)
            continue

        # Build junk map: dim_Am → junk_count
        J = np.zeros((junk_count, dim_Am))
        for j, path in enumerate(A_paths[m]):
            path_len = len(path) - 1
            for i in range(path_len + 1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in junk_faces:
                    J[junk_faces[face], j] += sign

        # Ω_m = kernel of J
        # Use SVD to find null space
        U, s, Vh = np.linalg.svd(J)
        tol = 1e-10
        null_mask = np.abs(s) < tol
        # If all singular values are nonzero, extend check
        rank_J = np.sum(np.abs(s) > tol)
        null_dim = dim_Am - rank_J
        if null_dim > 0:
            omega_basis[m] = Vh[rank_J:].T  # columns = null space basis
        else:
            omega_basis[m] = np.zeros((dim_Am, 0))

    # Step 3: Build boundary maps on Ω
    omega_dims = []
    for m in range(max_dim + 2):
        omega_dims.append(omega_basis[m].shape[1])

    ranks = {}
    for m in range(1, max_dim + 2):
        dim_omega_m = omega_dims[m]
        dim_omega_m1 = omega_dims[m-1]
        if dim_omega_m == 0 or dim_omega_m1 == 0:
            ranks[m] = 0
            continue

        # Full boundary matrix B: A_{m-1} → A_m (columns of B correspond to A_m paths)
        A_m_paths = A_paths[m]
        A_m1_paths = A_paths[m-1]
        path_to_idx = {p: i for i, p in enumerate(A_m1_paths)}
        dim_Am = len(A_m_paths)
        dim_Am1 = len(A_m1_paths)

        B_full = np.zeros((dim_Am1, dim_Am))
        for j, path in enumerate(A_m_paths):
            path_len = len(path) - 1
            for i in range(path_len + 1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in path_to_idx:
                    B_full[path_to_idx[face], j] += sign

        # Restrict: ∂_m: Ω_m → A_{m-1}, then project to Ω_{m-1} coordinates
        # B_omega = (omega_basis[m-1])^T @ B_full @ omega_basis[m]
        # This gives the boundary map in Ω coordinates
        B_omega = omega_basis[m-1].T @ B_full @ omega_basis[m]

        ranks[m] = np.linalg.matrix_rank(B_omega, tol=1e-8)

    # Betti numbers
    betti = []
    for m in range(max_dim + 1):
        om = omega_dims[m]
        rk_dm = ranks.get(m, 0)
        rk_dm1 = ranks.get(m+1, 0)
        betti.append(om - rk_dm - rk_dm1)

    return betti, omega_dims, [len(A_paths[m]) for m in range(max_dim + 1)]

def compute_interior_betti(A, max_dim=None):
    """Interior-only boundary on A_m."""
    n = A.shape[0]
    if max_dim is None:
        max_dim = n - 1
    all_paths = {}
    for m in range(max_dim + 1):
        all_paths[m] = get_regular_paths(A, m)
    ranks = {}
    for m in range(1, max_dim + 1):
        if not all_paths[m] or not all_paths[m-1]:
            ranks[m] = 0
            continue
        path_to_idx = {p: i for i, p in enumerate(all_paths[m-1])}
        path_m = len(all_paths[m][0]) - 1
        B = np.zeros((len(all_paths[m-1]), len(all_paths[m])), dtype=int)
        for j, path in enumerate(all_paths[m]):
            for i in range(1, path_m):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in path_to_idx:
                    B[path_to_idx[face], j] += sign
        ranks[m] = np.linalg.matrix_rank(B)
    betti = []
    for m in range(max_dim + 1):
        omega_m = len(all_paths[m])
        rank_dm = ranks.get(m, 0)
        rank_dm_plus_1 = ranks.get(m+1, 0)
        betti.append(omega_m - rank_dm - rank_dm_plus_1)
    return betti, [len(all_paths[m]) for m in range(max_dim + 1)]

# ============================================================
print("="*70)
print("GLMY PATH HOMOLOGY: CORRECT (subspace Ω_m) vs INTERIOR-ONLY")
print("="*70)

test_cases = [
    ("Transitive T_3", transitive_tournament(3)),
    ("Cyclic C_3", circulant_tournament(3, {1})),
    ("Transitive T_5", transitive_tournament(5)),
    ("Regular n=5", circulant_tournament(5, {1,2})),
    ("Paley P_7", circulant_tournament(7, {1,2,4})),
    ("Interval n=7", circulant_tournament(7, {1,2,3})),
]

for name, A in test_cases:
    n = A.shape[0]
    print(f"\n{'='*50}")
    print(f"  {name}")
    print(f"{'='*50}")

    betti_g, omega_dims, a_dims = compute_glmy_betti(A)
    betti_i, a_dims_i = compute_interior_betti(A)

    print(f"  A_m:       {a_dims}")
    print(f"  Ω_m (dim): {omega_dims[:len(a_dims)]}")

    gap = [a_dims[m] - omega_dims[m] for m in range(len(a_dims))]
    if any(g > 0 for g in gap):
        print(f"  A_m - Ω_m: {gap}")

    print(f"\n  GLMY Betti:     {betti_g}")
    print(f"  Interior Betti: {betti_i}")

    chi_g = sum((-1)**m * b for m, b in enumerate(betti_g))
    chi_i = sum((-1)**m * b for m, b in enumerate(betti_i))
    print(f"  chi (GLMY)={chi_g}, chi (INT)={chi_i}")

    if any(b < 0 for b in betti_g):
        print(f"  *** NEGATIVE GLMY BETTI — BUG! ***")

    # Check divisibility
    if name.startswith(("Paley", "Interval", "Regular")):
        div_g = all(b % n == 0 for b in betti_g)
        div_i = all(b % n == 0 for b in betti_i)
        if div_g: print(f"  GLMY β/n: {[b//n for b in betti_g]}")
        if div_i: print(f"  INT  β/n: {[b//n for b in betti_i]}")

    # β_2 check
    if n >= 5:
        print(f"  β_2: GLMY={betti_g[2]}, INT={betti_i[2]}")

# ============================================================
# Verify d^2 = 0 for GLMY boundary on Ω_m
print(f"\n{'='*70}")
print("d² = 0 VERIFICATION ON Ω_m")
print("="*70)

for name, A in [("Paley P_7", circulant_tournament(7, {1,2,4})),
                ("Interval n=7", circulant_tournament(7, {1,2,3}))]:
    n = A.shape[0]
    print(f"\n  {name}:")

    A_paths = {}
    for m in range(n+1):
        A_paths[m] = get_regular_paths(A, m)

    # Build full boundary and omega for each degree
    omega_bases = {}
    for m in range(n+1):
        dim_Am = len(A_paths[m])
        if m <= 1:
            omega_bases[m] = np.eye(dim_Am)
            continue
        A_set_prev = set(A_paths[m-1])
        junk_faces = {}
        junk_count = 0
        for path in A_paths[m]:
            plen = len(path)-1
            for i in range(plen+1):
                face = path[:i] + path[i+1:]
                if face not in A_set_prev and face not in junk_faces:
                    junk_faces[face] = junk_count
                    junk_count += 1
        if junk_count == 0:
            omega_bases[m] = np.eye(dim_Am)
            continue
        J = np.zeros((junk_count, dim_Am))
        for j, path in enumerate(A_paths[m]):
            plen = len(path)-1
            for i in range(plen+1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in junk_faces:
                    J[junk_faces[face], j] += sign
        U, s, Vh = np.linalg.svd(J)
        rank_J = np.sum(np.abs(s) > 1e-10)
        null_dim = dim_Am - rank_J
        if null_dim > 0:
            omega_bases[m] = Vh[rank_J:].T
        else:
            omega_bases[m] = np.zeros((dim_Am, 0))

    # Build boundary on Omega and check d^2
    def boundary_on_omega(m):
        if omega_bases[m].shape[1] == 0 or omega_bases[m-1].shape[1] == 0:
            return np.zeros((omega_bases[m-1].shape[1], omega_bases[m].shape[1]))
        path_to_idx = {p: i for i, p in enumerate(A_paths[m-1])}
        dim_Am = len(A_paths[m])
        dim_Am1 = len(A_paths[m-1])
        B = np.zeros((dim_Am1, dim_Am))
        for j, path in enumerate(A_paths[m]):
            plen = len(path)-1
            for i in range(plen+1):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in path_to_idx:
                    B[path_to_idx[face], j] += sign
        return omega_bases[m-1].T @ B @ omega_bases[m]

    for m in range(2, n):
        B_m = boundary_on_omega(m)
        B_m1 = boundary_on_omega(m-1)
        if B_m.shape[1] > 0 and B_m1.shape[1] > 0:
            d2 = B_m1 @ B_m
            maxabs = np.max(np.abs(d2)) if d2.size > 0 else 0
            print(f"    m={m}: d^2 max|entry| = {maxabs:.2e}, zero? {maxabs < 1e-8}")

print("\nDONE.")
