#!/usr/bin/env python3
"""
boundary_convention_check.py — opus-2026-03-13-S71

CRITICAL CHECK: Two different boundary conventions have been used in this repo.
1. FULL boundary: ∂(v_0,...,v_m) = Σ_{i=0}^{m} (-1)^i (v_0,...,v̂_i,...,v_m)
   Used in: path_homology_v2.py, circulant_homology.py
2. INTERIOR boundary: ∂(v_0,...,v_m) = Σ_{i=1}^{m-1} (-1)^i (v_0,...,v̂_i,...,v_m)
   Used in: betti_omega_connection.py, per_eigenspace_betti.py, betti_divisibility.py

Which is the correct GLMY path homology?
Answer: The FULL boundary is the standard GLMY convention.

Let's compute Betti with BOTH conventions for small cases and compare.
"""

import numpy as np

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
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

def boundary_matrix_FULL(paths_m, paths_m_minus_1):
    """Full GLMY boundary: delete ALL vertices."""
    if not paths_m or not paths_m_minus_1:
        return np.zeros((len(paths_m_minus_1) if paths_m_minus_1 else 0,
                         len(paths_m) if paths_m else 0), dtype=int)
    path_to_idx = {p: i for i, p in enumerate(paths_m_minus_1)}
    m = len(paths_m[0]) - 1
    matrix = np.zeros((len(paths_m_minus_1), len(paths_m)), dtype=int)
    for j, path in enumerate(paths_m):
        for i in range(m + 1):  # ALL vertices: 0 to m
            face = path[:i] + path[i+1:]
            sign = (-1)**i
            if face in path_to_idx:
                matrix[path_to_idx[face], j] += sign
    return matrix

def boundary_matrix_INTERIOR(paths_m, paths_m_minus_1):
    """Interior-only boundary: delete vertices 1 to m-1."""
    if not paths_m or not paths_m_minus_1:
        return np.zeros((len(paths_m_minus_1) if paths_m_minus_1 else 0,
                         len(paths_m) if paths_m else 0), dtype=int)
    path_to_idx = {p: i for i, p in enumerate(paths_m_minus_1)}
    m = len(paths_m[0]) - 1
    matrix = np.zeros((len(paths_m_minus_1), len(paths_m)), dtype=int)
    for j, path in enumerate(paths_m):
        for i in range(1, m):  # Interior only: 1 to m-1
            face = path[:i] + path[i+1:]
            sign = (-1)**i
            if face in path_to_idx:
                matrix[path_to_idx[face], j] += sign
    return matrix

def compute_betti(A, boundary_fn, label, max_dim=None):
    n = A.shape[0]
    if max_dim is None:
        max_dim = n - 1

    all_paths = {}
    for m in range(max_dim + 1):
        all_paths[m] = get_regular_paths(A, m)

    ranks = {}
    for m in range(1, max_dim + 1):
        if all_paths[m] and all_paths[m-1]:
            B = boundary_fn(all_paths[m], all_paths[m-1])
            # Verify d^2 = 0
            if m >= 2 and all_paths[m-2]:
                B_prev = boundary_fn(all_paths[m-1], all_paths[m-2])
                d2 = B_prev @ B
                if np.any(np.abs(d2) > 1e-10):
                    print(f"    WARNING: d^2 != 0 at m={m} for {label}!")
            ranks[m] = np.linalg.matrix_rank(B)
        else:
            ranks[m] = 0

    betti = []
    for m in range(max_dim + 1):
        omega_m = len(all_paths[m])
        rank_dm = ranks.get(m, 0)
        rank_dm_plus_1 = ranks.get(m+1, 0)
        beta_m = omega_m - rank_dm - rank_dm_plus_1
        betti.append(beta_m)

    omegas = [len(all_paths[m]) for m in range(max_dim + 1)]
    return betti, omegas, ranks

# ============================================================
print("="*70)
print("BOUNDARY CONVENTION COMPARISON")
print("="*70)

test_cases = [
    ("Transitive T_5", 5, None),
    ("Paley P_7", 7, {1,2,4}),
    ("Interval n=7", 7, {1,2,3}),
    ("Regular n=5", 5, {1,2}),
]

for name, n, S in test_cases:
    if S is None:
        # Transitive tournament
        A = np.zeros((n,n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                A[i][j] = 1
    else:
        A = circulant_tournament(n, S)

    print(f"\n{'='*50}")
    print(f"  {name}")
    print(f"{'='*50}")

    betti_full, omegas, ranks_full = compute_betti(A, boundary_matrix_FULL, "FULL")
    betti_int, _, ranks_int = compute_betti(A, boundary_matrix_INTERIOR, "INTERIOR")

    print(f"  Omega:  {omegas}")
    print(f"  FULL boundary:     β = {betti_full}   ranks = {[ranks_full.get(m,0) for m in range(1,n)]}")
    print(f"  INTERIOR boundary: β = {betti_int}   ranks = {[ranks_int.get(m,0) for m in range(1,n)]}")

    chi_full = sum((-1)**m * betti_full[m] for m in range(len(betti_full)))
    chi_int = sum((-1)**m * betti_int[m] for m in range(len(betti_int)))
    print(f"  chi (FULL) = {chi_full}, chi (INTERIOR) = {chi_int}")

    if betti_full == betti_int:
        print(f"  *** SAME BETTI ***")
    else:
        print(f"  *** DIFFERENT BETTI! ***")
        # Show which dims differ
        for m in range(len(betti_full)):
            if betti_full[m] != betti_int[m]:
                print(f"    β_{m}: FULL={betti_full[m]}, INT={betti_int[m]}")

# ============================================================
print(f"\n{'='*70}")
print("CHECKING d^2 = 0 FOR BOTH CONVENTIONS")
print("="*70)

# Test d^2 = 0 explicitly for Paley P_7
A = circulant_tournament(7, {1,2,4})
all_paths = {}
for m in range(7):
    all_paths[m] = get_regular_paths(A, m)

print("\n  FULL boundary:")
for m in range(2, 7):
    if all_paths[m] and all_paths[m-1] and all_paths[m-2]:
        B_m = boundary_matrix_FULL(all_paths[m], all_paths[m-1])
        B_m1 = boundary_matrix_FULL(all_paths[m-1], all_paths[m-2])
        d2 = B_m1 @ B_m
        print(f"    m={m}: d^2 zero? {np.allclose(d2, 0)}, max|d^2|={np.max(np.abs(d2))}")

print("\n  INTERIOR boundary:")
for m in range(2, 7):
    if all_paths[m] and all_paths[m-1] and all_paths[m-2]:
        B_m = boundary_matrix_INTERIOR(all_paths[m], all_paths[m-1])
        B_m1 = boundary_matrix_INTERIOR(all_paths[m-1], all_paths[m-2])
        d2 = B_m1 @ B_m
        print(f"    m={m}: d^2 zero? {np.allclose(d2, 0)}, max|d^2|={np.max(np.abs(d2))}")

# ============================================================
print(f"\n{'='*70}")
print("WHICH CONVENTION MATCHES path_homology_v2.py?")
print("="*70)

# path_homology_v2.py uses full boundary AND the Omega subspace construction.
# Let's check if that gives different results from direct FULL boundary on
# regular paths. For tournaments, Omega_p = A_p (all regular path faces are
# regular), so the Omega construction shouldn't change anything.

# But WAIT: path_homology_v2.py computes Ω_p as the subspace of A_p whose
# boundary lies in A_{p-1}. Even though all faces of regular tournament paths
# are regular, some ENDPOINT deletions might give non-regular faces!

# Check: for (v_0, v_1, v_2), deleting v_0 gives (v_1, v_2). Is this regular?
# Need v_1 → v_2 (which we have). For a 1-path, there's no regularity condition
# beyond the edge, so yes.

# For (v_0,v_1,v_2,v_3), deleting v_0 gives (v_1,v_2,v_3). Need v_1→v_2, v_2→v_3,
# AND v_1→v_3. We have the first two from regularity. v_1→v_3 is from the skip-one
# condition on the original path. So yes, the face is regular.

# Actually, the regularity of the original path (v_0,...,v_m) requires v_i → v_{i+2}
# for all valid i. When we delete v_0, we need v_1 → v_3 (which was v_{1} → v_{0+3}
# in the original, guaranteed by regularity). So ALL endpoint deletions give regular
# faces for tournament regular paths.

# Therefore Ω_p = A_p for tournaments, and path_homology_v2.py with full boundary
# should give the same result as our boundary_matrix_FULL.

print("  For tournaments: Ω_p = A_p (all faces of regular paths are regular)")
print("  Therefore: FULL boundary on regular paths = GLMY path homology")

# Let's also check: do the S70 results match FULL or INTERIOR?
# S70 session log says: Transitive T_n: β=(n, n-1, 0, ..., 0)
# FULL gives T_5: β=[5,4,0,0,0]? Let me check above.

print(f"\n  S70 reported: Transitive T_5 β = (5, 4, 0, 0, 0)")
print(f"  Our FULL result for T_5: see above")

# Also compare with known GLMY results from literature
# For the complete tournament T_3 (transitive on 3 vertices):
A3 = np.zeros((3,3), dtype=int)
A3[0][1] = A3[0][2] = A3[1][2] = 1
betti3_f, omega3, _ = compute_betti(A3, boundary_matrix_FULL, "FULL")
betti3_i, _, _ = compute_betti(A3, boundary_matrix_INTERIOR, "INTERIOR")
print(f"\n  Transitive T_3:")
print(f"    Omega = {omega3}")
print(f"    FULL: β = {betti3_f}")
print(f"    INT:  β = {betti3_i}")
print(f"    Expected (GLMY literature): β = [1, 0, 0] for acyclic")

# For the cyclic tournament C_3:
A3c = np.zeros((3,3), dtype=int)
A3c[0][1] = A3c[1][2] = A3c[2][0] = 1
betti3c_f, omega3c, _ = compute_betti(A3c, boundary_matrix_FULL, "FULL")
betti3c_i, _, _ = compute_betti(A3c, boundary_matrix_INTERIOR, "INTERIOR")
print(f"\n  Cyclic C_3:")
print(f"    Omega = {omega3c}")
print(f"    FULL: β = {betti3c_f}")
print(f"    INT:  β = {betti3c_i}")
print(f"    Expected (GLMY): β = [1, 1, 0] (one 1-cycle)")

print("\nDONE.")
