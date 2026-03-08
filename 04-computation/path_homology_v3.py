#!/usr/bin/env python3
"""
GLMY PATH HOMOLOGY v3 — Optimized implementation.

Optimizations over v2:
1. DFS-based path enumeration instead of permutation filtering.
   For sparse digraphs, this is dramatically faster.
2. Bitset-based visited tracking for DFS.
3. Same mathematical definitions as v2 (drop-in replacement).

DEFINITIONS (GLMY):
- Elementary p-path: sequence (v_0, ..., v_p) of DISTINCT vertices
- ALLOWED p-path: one where v_i -> v_{i+1} is a directed edge for all i
- Boundary: ∂(v_0...v_p) = sum_{i=0}^p (-1)^i (v_0...v̂_i...v_p)
- A_p = vector space spanned by allowed p-paths
- Ω_p = {u ∈ A_p : ∂u ∈ A_{p-1}}  (for p ≥ 1; Ω_0 = A_0)
- H_p^{path}(G) = ker(∂_p: Ω_p → Ω_{p-1}) / im(∂_{p+1}: Ω_{p+1} → Ω_p)

Author: opus-2026-03-08-S48
"""
import numpy as np
from collections import defaultdict


def enumerate_allowed_paths(A, n, p):
    """All sequences of p+1 distinct vertices following directed edges.

    Uses DFS instead of permutation enumeration — much faster for
    sparse or moderately dense digraphs like tournaments.
    """
    if p < 0:
        return []
    if p == 0:
        return [(v,) for v in range(n)]

    # Precompute adjacency lists
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                adj[i].append(j)

    paths = []
    # DFS from each starting vertex
    stack = []  # (current_path_as_list, visited_bitmask)

    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))

    return paths


def boundary_coeffs(path):
    """Returns list of (sign, face_tuple) for ∂(path)."""
    p = len(path) - 1
    result = []
    for i in range(p + 1):
        face = path[:i] + path[i+1:]
        result.append(((-1)**i, face))
    return result


def build_full_boundary_matrix(allowed_p, allowed_pm1):
    """Build ∂: A_p → A_{p-1} as a matrix."""
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))

    idx_pm1 = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))

    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx_pm1:
                M[idx_pm1[face], j] += sign
    return M


def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    """Compute a basis for Ω_p = {u ∈ A_p : ∂u ∈ A_{p-1}}."""
    dim_Ap = len(allowed_p)
    if dim_Ap == 0:
        return np.zeros((0, 0))

    if p == 0:
        return np.eye(dim_Ap)

    allowed_pm1_set = set(allowed_pm1)

    non_allowed_faces = {}
    na_count = 0

    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed_faces:
                    non_allowed_faces[face] = na_count
                    na_count += 1

    if na_count == 0:
        return np.eye(dim_Ap)

    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed_faces:
                P[non_allowed_faces[face], j] += sign

    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    null_space = Vt[rank:].T

    if null_space.shape[1] == 0:
        return np.zeros((dim_Ap, 0))
    return null_space


def path_betti_numbers(A, n, max_dim=None):
    """Compute GLMY path Betti numbers β_0, β_1, ..., β_{max_dim}."""
    if max_dim is None:
        max_dim = n - 1

    allowed = {}
    for p in range(-1, max_dim + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)

    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])

    betti = []
    for p in range(max_dim + 1):
        dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0

        if dim_omega_p == 0:
            betti.append(0)
            continue

        bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_omega = bd_p @ omega[p]

        if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
            S_p = np.linalg.svd(bd_p_omega, compute_uv=False)
            rank_p = sum(s > 1e-8 for s in S_p)
        else:
            rank_p = 0

        ker_dim = dim_omega_p - rank_p

        dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_omega_p1 > 0:
            bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            bd_p1_omega = bd_p1 @ omega[p+1]
            S_p1 = np.linalg.svd(bd_p1_omega, compute_uv=False)
            im_dim = sum(s > 1e-8 for s in S_p1)
        else:
            im_dim = 0

        beta_p = ker_dim - im_dim
        betti.append(max(0, beta_p))

    return betti


# ===== VALIDATION & BENCHMARK =====
if __name__ == "__main__":
    from time import perf_counter
    from itertools import permutations

    print("=" * 70)
    print("GLMY PATH HOMOLOGY v3 — VALIDATION & BENCHMARK")
    print("=" * 70)

    # Test 1: Directed 3-cycle
    print("\n--- Directed 3-cycle C_3 ---")
    A = [[0,1,0],[0,0,1],[1,0,0]]
    print(f"  β = {path_betti_numbers(A, 3)}  (expected [1, 1, 0])")

    # Test 2: Transitive tournament T_3
    A = [[0,1,1],[0,0,1],[0,0,0]]
    print(f"  T_3: β = {path_betti_numbers(A, 3)}  (expected [1, 0, 0])")

    # Benchmark: enumerate_allowed_paths v2 vs v3
    print("\n--- Benchmark: path enumeration v2 vs v3 ---")

    def enumerate_v2(A, n, p):
        """Old v2 method."""
        paths = []
        for perm in permutations(range(n), p + 1):
            ok = True
            for i in range(p):
                if A[perm[i]][perm[i+1]] != 1:
                    ok = False
                    break
            if ok:
                paths.append(perm)
        return paths

    # Random tournament
    import random
    random.seed(42)
    for n_test in [7, 8, 9]:
        A = [[0]*n_test for _ in range(n_test)]
        for i in range(n_test):
            for j in range(i+1, n_test):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        for p_test in [2, 3, 4]:
            t0 = perf_counter()
            paths_v2 = enumerate_v2(A, n_test, p_test)
            dt_v2 = perf_counter() - t0

            t0 = perf_counter()
            paths_v3 = enumerate_allowed_paths(A, n_test, p_test)
            dt_v3 = perf_counter() - t0

            match = set(paths_v2) == set(paths_v3)
            speedup = dt_v2 / dt_v3 if dt_v3 > 0 else float('inf')
            print(f"  n={n_test} p={p_test}: v2={dt_v2:.4f}s v3={dt_v3:.4f}s "
                  f"speedup={speedup:.1f}x  {len(paths_v3)} paths  {'OK' if match else 'MISMATCH'}")

    # Full betti benchmark
    print("\n--- Full betti benchmark ---")
    for n_test in [5, 6, 7]:
        A = [[0]*n_test for _ in range(n_test)]
        for i in range(n_test):
            for j in range(i+1, n_test):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        t0 = perf_counter()
        betti = path_betti_numbers(A, n_test)
        dt = perf_counter() - t0
        print(f"  n={n_test}: β = {betti}  [{dt:.4f}s]")

    print("\nDone.")
