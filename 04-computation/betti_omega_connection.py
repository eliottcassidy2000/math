#!/usr/bin/env python3
"""
Betti-Omega Connection for Tournaments
=======================================
Computes GLMY path homology Betti numbers and Omega profiles for:
1. All 12 tournament isomorphism types on n=5
2. Paley tournament on p=7 (QR = {1,2,4})
3. Interval tournament on n=7 (S = {1,2,3})
4. Interval tournament on n=9 (S = {1,2,3,4})

Uses the GLMY "regular path" formulation:
- Regular m-path (v_0,...,v_m): all v_i distinct, v_i -> v_{i+1} (edges),
  and v_{i-1} -> v_{i+1} for 1 <= i <= m-1 (skip-one / regularity condition)
- Boundary: d_m(v_0,...,v_m) = sum_{i=1}^{m-1} (-1)^i (v_0,...,v_hat_i,...,v_m)
  (only interior vertices deleted)
- Betti_m = dim(ker(d_m)) - dim(im(d_{m+1}))
"""

import numpy as np
from itertools import permutations
from collections import defaultdict

# ============================================================
# Core: Regular path enumeration
# ============================================================

def enumerate_regular_paths(A, n, m):
    """
    Enumerate all regular m-paths in tournament with adjacency matrix A.

    A regular m-path (v_0, ..., v_m) satisfies:
    1. All vertices distinct
    2. A[v_i][v_{i+1}] = 1 for all 0 <= i < m  (path edges)
    3. A[v_{i-1}][v_{i+1}] = 1 for all 1 <= i <= m-1  (skip-one edges)

    Returns list of tuples.
    """
    if m < 0:
        return []
    if m == 0:
        return [(v,) for v in range(n)]

    # Build adjacency lists
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                adj[i].append(j)

    # DFS with regularity check
    paths = []
    # Stack: (current_path, visited_bitmask)
    for start in range(n):
        stack = [([start], 1 << start)]
        while stack:
            path, visited = stack.pop()
            if len(path) == m + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if visited & (1 << u):
                    continue
                # Check regularity: if path has at least 2 vertices,
                # the second-to-last must have edge to u
                if len(path) >= 2:
                    prev = path[-2]
                    if A[prev][u] != 1:
                        continue
                stack.append((path + [u], visited | (1 << u)))

    return paths


def boundary_interior(path):
    """
    Boundary map deleting only interior vertices.
    d(v_0,...,v_m) = sum_{i=1}^{m-1} (-1)^i (v_0,...,v_hat_i,...,v_m)

    Returns list of (sign, face_tuple).
    """
    m = len(path) - 1
    if m <= 0:
        return []
    result = []
    for i in range(1, m):  # interior indices only: 1, ..., m-1
        face = path[:i] + path[i+1:]
        result.append(((-1)**i, face))
    return result


def build_boundary_matrix(paths_m, paths_m1):
    """
    Build the boundary matrix d_m: Omega_m -> Omega_{m-1}.
    Columns indexed by regular m-paths, rows by regular (m-1)-paths.
    Uses interior-only deletion.
    """
    if not paths_m or not paths_m1:
        return np.zeros((max(len(paths_m1), 0), max(len(paths_m), 0)))

    idx = {p: i for i, p in enumerate(paths_m1)}
    M = np.zeros((len(paths_m1), len(paths_m)))

    for j, path in enumerate(paths_m):
        for sign, face in boundary_interior(path):
            if face in idx:
                M[idx[face], j] += sign

    return M


def compute_betti_and_omega(A, n, max_dim=None):
    """
    Compute Omega profile and Betti numbers for tournament A on n vertices.

    Returns (omega_dims, betti_nums).
    """
    if max_dim is None:
        max_dim = n - 1

    # Enumerate regular paths for each dimension
    regular = {}
    for m in range(max_dim + 2):
        regular[m] = enumerate_regular_paths(A, n, m)

    omega_dims = [len(regular[m]) for m in range(max_dim + 1)]

    # Build boundary matrices
    bd = {}
    for m in range(max_dim + 2):
        bd[m] = build_boundary_matrix(regular[m], regular.get(m-1, []))

    betti = []
    for m in range(max_dim + 1):
        # dim(ker(d_m))
        if bd[m].shape[1] == 0:
            ker_dim = 0
        elif bd[m].shape[0] == 0:
            ker_dim = bd[m].shape[1]  # d_m maps to 0-space, everything is kernel
        else:
            sv = np.linalg.svd(bd[m], compute_uv=False)
            rank_dm = int(np.sum(sv > 1e-8))
            ker_dim = bd[m].shape[1] - rank_dm

        # dim(im(d_{m+1}))
        if m + 1 > max_dim + 1 or bd[m+1].shape[1] == 0:
            im_dim = 0
        else:
            sv1 = np.linalg.svd(bd[m+1], compute_uv=False)
            im_dim = int(np.sum(sv1 > 1e-8))

        beta_m = ker_dim - im_dim
        betti.append(max(0, beta_m))  # numerical safety

    return omega_dims, betti


# ============================================================
# Tournament constructors
# ============================================================

def transitive_tournament(n):
    """Transitive tournament: i -> j iff i < j."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def circulant_tournament(n, S):
    """Circulant tournament on Z_n with connection set S.
    A[i][j] = 1 iff (j - i) mod n in S.
    """
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S:
                A[i][j] = 1
    return A

def all_tournaments_by_iso(n):
    """Enumerate all isomorphism classes of tournaments on n vertices."""
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    seen = set()
    reps = []
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(edges):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        # Canonical form
        canon = None
        for perm in permutations(range(n)):
            B = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    B[perm[i]][perm[j]] = A[i][j]
            key = tuple(tuple(row) for row in B)
            if canon is None or key < canon:
                canon = key
        if canon not in seen:
            seen.add(canon)
            reps.append(A)
    return reps

def score_sequence(A, n):
    """Score sequence of tournament (sorted)."""
    scores = [sum(A[i]) for i in range(n)]
    return tuple(sorted(scores))

def count_3cycles(A, n):
    """Count directed 3-cycles."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[i][k] and A[k][j] and A[j][i]):
                    count += 1
    return count


# ============================================================
# Verification
# ============================================================

def verify_boundary_squared_zero(A, n, max_dim=4):
    """Verify d_{m-1} . d_m = 0 for the chain complex."""
    regular = {}
    for m in range(max_dim + 2):
        regular[m] = enumerate_regular_paths(A, n, m)

    for m in range(2, max_dim + 1):
        bd_m = build_boundary_matrix(regular[m], regular[m-1])
        bd_m1 = build_boundary_matrix(regular[m-1], regular[m-2])
        if bd_m.shape[1] > 0 and bd_m1.shape[1] > 0:
            comp = bd_m1 @ bd_m
            if np.max(np.abs(comp)) > 1e-10:
                return False, m
    return True, None


# ============================================================
# Main computation
# ============================================================

def print_tournament_info(name, A, n):
    """Compute and print Omega + Betti for a tournament."""
    omega, betti = compute_betti_and_omega(A, n)
    chi = sum((-1)**m * b for m, b in enumerate(betti))

    print(f"\n  {name}")
    print(f"    Score seq: {score_sequence(A, n)},  3-cycles: {count_3cycles(A, n)}")
    print(f"    Omega: {omega}")
    print(f"    Betti: {betti}")
    print(f"    chi  : {chi}")
    return omega, betti, chi


def main():
    print("=" * 72)
    print("BETTI-OMEGA CONNECTION FOR TOURNAMENTS")
    print("=" * 72)

    # ----------------------------------------------------------
    # Verify chain complex property
    # ----------------------------------------------------------
    print("\n--- Chain complex verification (d^2 = 0) ---")
    T5_trans = transitive_tournament(5)
    ok, dim = verify_boundary_squared_zero(T5_trans, 5)
    print(f"  Transitive T_5: d^2=0 ? {ok}")

    P7 = circulant_tournament(7, {1, 2, 4})
    ok, dim = verify_boundary_squared_zero(P7, 7, max_dim=6)
    print(f"  Paley P_7:      d^2=0 ? {ok}")

    # ----------------------------------------------------------
    # Part 1: All 12 tournament types on n=5
    # ----------------------------------------------------------
    print("\n" + "=" * 72)
    print("PART 1: ALL 12 TOURNAMENT ISOMORPHISM TYPES ON n=5")
    print("=" * 72)

    reps = all_tournaments_by_iso(5)
    print(f"  Found {len(reps)} isomorphism classes")

    results_n5 = []
    for idx, A in enumerate(reps):
        omega, betti, chi = print_tournament_info(
            f"Type {idx+1}", A, 5
        )
        results_n5.append((score_sequence(A, 5), count_3cycles(A, 5), omega, betti, chi))

    # Summary table
    print("\n  --- Summary ---")
    print(f"  {'Type':>5} {'Score':>15} {'t3':>4} {'Omega':>30} {'Betti':>20} {'chi':>5}")
    for idx, (ss, t3, om, bt, chi) in enumerate(results_n5):
        print(f"  {idx+1:>5} {str(ss):>15} {t3:>4} {str(om):>30} {str(bt):>20} {chi:>5}")

    # Check beta_2 = 0
    all_b2_zero = all(bt[2] == 0 for _, _, _, bt, _ in results_n5)
    print(f"\n  beta_2 = 0 for ALL types? {all_b2_zero}")

    # ----------------------------------------------------------
    # Part 2: Paley tournament on p=7
    # ----------------------------------------------------------
    print("\n" + "=" * 72)
    print("PART 2: PALEY TOURNAMENT P_7 (QR = {1,2,4})")
    print("=" * 72)

    P7 = circulant_tournament(7, {1, 2, 4})
    print_tournament_info("Paley P_7", P7, 7)

    # ----------------------------------------------------------
    # Part 3: Interval tournament on n=7, S={1,2,3}
    # ----------------------------------------------------------
    print("\n" + "=" * 72)
    print("PART 3: INTERVAL TOURNAMENT n=7, S={1,2,3}")
    print("=" * 72)

    I7 = circulant_tournament(7, {1, 2, 3})
    print_tournament_info("Interval C_7^{1,2,3}", I7, 7)

    # ----------------------------------------------------------
    # Part 4: Interval tournament on n=9, S={1,2,3,4}
    # ----------------------------------------------------------
    print("\n" + "=" * 72)
    print("PART 4: INTERVAL TOURNAMENT n=9, S={1,2,3,4}")
    print("=" * 72)

    I9 = circulant_tournament(9, {1, 2, 3, 4})
    print_tournament_info("Interval C_9^{1,2,3,4}", I9, 9)

    # ----------------------------------------------------------
    # Part 5: Comparison of Paley vs Interval at n=7
    # ----------------------------------------------------------
    print("\n" + "=" * 72)
    print("COMPARISON: PALEY vs INTERVAL at n=7")
    print("=" * 72)

    om_p, bt_p, chi_p = compute_betti_and_omega(P7, 7)[0], *compute_betti_and_omega(P7, 7)
    om_i, bt_i, chi_i = compute_betti_and_omega(I7, 7)[0], *compute_betti_and_omega(I7, 7)

    # Fix: recompute properly
    om_p, bt_p = compute_betti_and_omega(P7, 7)
    om_i, bt_i = compute_betti_and_omega(I7, 7)
    chi_p = sum((-1)**m * b for m, b in enumerate(bt_p))
    chi_i = sum((-1)**m * b for m, b in enumerate(bt_i))

    print(f"  {'':>20} {'Paley':>15} {'Interval':>15}")
    for m in range(7):
        print(f"  dim {m}: Omega={om_p[m]:>5} beta={bt_p[m]:>3}    Omega={om_i[m]:>5} beta={bt_i[m]:>3}")
    print(f"  {'chi':>20} {chi_p:>15} {chi_i:>15}")

    # ----------------------------------------------------------
    # Bonus: Transitive tournaments comparison
    # ----------------------------------------------------------
    print("\n" + "=" * 72)
    print("BONUS: TRANSITIVE TOURNAMENTS T_n for n=3..8")
    print("=" * 72)

    for n in range(3, 9):
        T = transitive_tournament(n)
        max_d = min(n-1, 7)
        omega, betti = compute_betti_and_omega(T, n, max_dim=max_d)
        chi = sum((-1)**m * b for m, b in enumerate(betti))
        print(f"  T_{n}: Omega={omega},  Betti={betti},  chi={chi}")

    print("\n" + "=" * 72)
    print("DONE")
    print("=" * 72)


if __name__ == "__main__":
    main()
