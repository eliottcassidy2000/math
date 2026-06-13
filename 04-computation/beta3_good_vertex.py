# ⚠️ WARNING: This script assumes beta_3 ≤ 1 for all tournaments, which is
# FALSE at n=8 (beta_3=2 found in 0.08% of tournaments). See MISTAKE-018.
# Results are valid only for n ≤ 7.

#!/usr/bin/env python3
"""
beta3_good_vertex.py — opus-2026-03-09-S52

The LES gives: beta_3(T) <= beta_3(T\v) + dim H_3(T,T\v) <= 1 + 1 = 2.
To get beta_3 <= 1, we need ONE of:
  (a) "Good vertex": exists v with beta_3(T\v) = 0 or H_3(T,T\v) = 0
  (b) Sharper bound: when both are 1, rank(i_*) + rank(j_*) < 2

Strategy: For each tournament T with beta_3(T)=1, check that there EXISTS
a vertex v where the LES contribution is ≤ 1.

Also check the "simultaneous" condition: when beta_3(T\v)=1 AND H_3(T,T\v)=1,
does i_* = 0 necessarily?
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    return T

def is_allowed_path(T, path):
    for i in range(len(path)-1):
        if not T[path[i]][path[i+1]]:
            return False
    return len(path) == len(set(path))

def compute_betti_full(T, n, max_p=4):
    """Compute all Betti numbers for tournament T."""
    tol = 1e-8

    all_paths = {}
    for p in range(0, max_p + 2):
        paths = []
        if p + 1 <= n:
            for verts in combinations(range(n), p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths.append(perm)
        all_paths[p] = paths

    def build_omega(paths_dict, mp):
        omega = {}
        for p in range(0, mp + 1):
            a_p = paths_dict[p]
            if not a_p:
                omega[p] = np.zeros((0, 0))
                continue
            a_pm1_set = set(paths_dict[p-1]) if p > 0 else set()
            na_faces = {}
            for sigma in a_p:
                for i in range(1, len(sigma)-1):
                    face = sigma[:i] + sigma[i+1:]
                    if p > 0 and face not in a_pm1_set:
                        if face not in na_faces:
                            na_faces[face] = len(na_faces)
            if not na_faces:
                omega[p] = np.eye(len(a_p))
            else:
                mat = np.zeros((len(na_faces), len(a_p)))
                for j, sigma in enumerate(a_p):
                    for i in range(1, len(sigma)-1):
                        face = sigma[:i] + sigma[i+1:]
                        if face in na_faces:
                            mat[na_faces[face], j] += (-1)**i
                U, S, Vt = np.linalg.svd(mat, full_matrices=True)
                rank = int(np.sum(S > tol))
                null_dim = len(a_p) - rank
                if null_dim == 0:
                    omega[p] = np.zeros((len(a_p), 0))
                else:
                    omega[p] = Vt[rank:].T
        return omega

    omega = build_omega(all_paths, max_p)

    # Boundary maps in Omega basis
    boundary = {}
    for p in range(1, max_p + 1):
        a_p = all_paths[p]
        a_pm1 = all_paths[p-1]
        if not a_p or not a_pm1:
            boundary[p] = np.zeros((0, 0))
            continue
        idx = {path: i for i, path in enumerate(a_pm1)}
        raw = np.zeros((len(a_pm1), len(a_p)))
        for j, sigma in enumerate(a_p):
            for i in range(len(sigma)):
                face = sigma[:i] + sigma[i+1:]
                if face in idx:
                    raw[idx[face], j] += (-1)**i
        boundary[p] = raw

    # Omega boundary
    omega_bd = {}
    for p in range(1, max_p + 1):
        Om_p = omega[p]
        Om_pm1 = omega[p-1]
        if Om_p.ndim < 2 or Om_p.shape[1] == 0 or Om_pm1.ndim < 2 or Om_pm1.shape[1] == 0:
            omega_bd[p] = np.zeros((0, 0))
            continue
        omega_bd[p] = Om_pm1.T @ boundary[p] @ Om_p

    # Ranks and Betti
    ranks = {}
    dims = {}
    for p in range(max_p + 1):
        dims[p] = omega[p].shape[1] if omega[p].ndim == 2 else 0
    for p in range(1, max_p + 1):
        if omega_bd[p].size == 0:
            ranks[p] = 0
        else:
            S = np.linalg.svd(omega_bd[p], compute_uv=False)
            ranks[p] = int(np.sum(S > tol))

    betti = {}
    for p in range(max_p + 1):
        ker = dims[p] - ranks.get(p, 0)
        im_next = ranks.get(p+1, 0)
        betti[p] = ker - im_next

    return betti, omega, omega_bd, ranks, dims, all_paths

def compute_relative_betti(T, n, v, max_p=4):
    """Compute relative Betti H_p(T, T\\v)."""
    tol = 1e-8

    all_T = {}
    all_Tv = {}
    for p in range(0, max_p + 2):
        paths_T = []
        if p + 1 <= n:
            for verts in combinations(range(n), p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths_T.append(perm)
        all_T[p] = paths_T
        all_Tv[p] = [path for path in paths_T if v not in path]

    def build_omega(paths_dict, mp):
        omega = {}
        for p in range(0, mp + 1):
            a_p = paths_dict[p]
            if not a_p:
                omega[p] = np.zeros((0, 0))
                continue
            a_pm1_set = set(paths_dict[p-1]) if p > 0 else set()
            na_faces = {}
            for sigma in a_p:
                for i in range(1, len(sigma)-1):
                    face = sigma[:i] + sigma[i+1:]
                    if p > 0 and face not in a_pm1_set:
                        if face not in na_faces:
                            na_faces[face] = len(na_faces)
            if not na_faces:
                omega[p] = np.eye(len(a_p))
            else:
                mat = np.zeros((len(na_faces), len(a_p)))
                for j, sigma in enumerate(a_p):
                    for i in range(1, len(sigma)-1):
                        face = sigma[:i] + sigma[i+1:]
                        if face in na_faces:
                            mat[na_faces[face], j] += (-1)**i
                U, S, Vt = np.linalg.svd(mat, full_matrices=True)
                rank = int(np.sum(S > tol))
                null_dim = len(a_p) - rank
                if null_dim == 0:
                    omega[p] = np.zeros((len(a_p), 0))
                else:
                    omega[p] = Vt[rank:].T
        return omega

    omega_T = build_omega(all_T, max_p)
    omega_Tv = build_omega(all_Tv, max_p)

    boundary_T = {}
    for p in range(1, max_p + 1):
        a_p = all_T[p]
        a_pm1 = all_T[p-1]
        if not a_p or not a_pm1:
            boundary_T[p] = np.zeros((len(a_pm1) if a_pm1 else 0, len(a_p) if a_p else 0))
            continue
        idx = {path: i for i, path in enumerate(a_pm1)}
        mat = np.zeros((len(a_pm1), len(a_p)))
        for j, sigma in enumerate(a_p):
            for i in range(len(sigma)):
                face = sigma[:i] + sigma[i+1:]
                if face in idx:
                    mat[idx[face], j] += (-1)**i
        boundary_T[p] = mat

    omega_Tv_in_T = {}
    for p in range(0, max_p + 1):
        a_p_T = all_T[p]
        a_p_Tv = all_Tv[p]
        om_Tv = omega_Tv[p]
        if om_Tv.ndim < 2 or om_Tv.shape[1] == 0 or not a_p_T:
            omega_Tv_in_T[p] = np.zeros((len(a_p_T) if a_p_T else 0, 0))
            continue
        idx_T = {path: i for i, path in enumerate(a_p_T)}
        mapping = np.zeros((len(a_p_T), len(a_p_Tv)))
        for j, path in enumerate(a_p_Tv):
            if path in idx_T:
                mapping[idx_T[path], j] = 1.0
        omega_Tv_in_T[p] = mapping @ om_Tv

    complement = {}
    rel_dims = {}
    for p in range(0, max_p + 1):
        om_T = omega_T[p]
        om_Tv = omega_Tv_in_T[p]
        if om_T.ndim < 2 or om_T.shape[1] == 0:
            complement[p] = np.zeros((0, 0))
            rel_dims[p] = 0
            continue
        if om_Tv.shape[1] == 0:
            complement[p] = om_T
            rel_dims[p] = om_T.shape[1]
            continue
        proj = om_T.T @ om_Tv
        U, S, Vt = np.linalg.svd(proj, full_matrices=True)
        sub_rank = int(np.sum(S > tol))
        if sub_rank < om_T.shape[1]:
            comp_coords = U[:, sub_rank:]
            complement[p] = om_T @ comp_coords
            rel_dims[p] = comp_coords.shape[1]
        else:
            complement[p] = np.zeros((om_T.shape[0], 0))
            rel_dims[p] = 0

    rel_rank = {}
    for p in range(1, max_p + 1):
        comp_p = complement[p]
        om_T_pm1 = omega_T[p-1]
        comp_pm1 = complement[p-1]
        if comp_p.ndim < 2 or comp_p.shape[1] == 0 or comp_pm1.ndim < 2 or comp_pm1.shape[1] == 0:
            rel_rank[p] = 0
            continue
        if p not in boundary_T or boundary_T[p].shape[0] == 0:
            rel_rank[p] = 0
            continue
        dp_comp = boundary_T[p] @ comp_p
        dp_in_omega = om_T_pm1.T @ dp_comp
        comp_pm1_in_omega = om_T_pm1.T @ comp_pm1
        if comp_pm1_in_omega.shape[1] == 0:
            rel_rank[p] = 0
            continue
        Q, R = np.linalg.qr(comp_pm1_in_omega, mode='reduced')
        proj_dp = Q.T @ dp_in_omega
        S = np.linalg.svd(proj_dp, compute_uv=False)
        rel_rank[p] = int(np.sum(S > tol))

    rel_betti = {}
    for p in range(0, max_p + 1):
        ker = rel_dims[p] - rel_rank.get(p, 0)
        im_next = rel_rank.get(p+1, 0)
        rel_betti[p] = ker - im_next

    return rel_betti

def main():
    print("=" * 70)
    print("GOOD VERTEX ANALYSIS: Can we prove beta_3 <= 1?")
    print("=" * 70)

    # =====================================================================
    # Part 1: At n=7, for every T with beta_3=1, check LES contributions
    # =====================================================================
    print("\n--- Part 1: n=7 sampled ---")
    n = 7
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs
    rng = np.random.RandomState(42)
    sample_size = 200

    # Track the critical question: for T with beta_3(T)=1,
    # does there exist v with beta_3(T\v) + H_3(T,T\v) <= 1?
    total_b3_eq_1 = 0
    has_good_vertex = 0
    worst_min_sum = 0

    # Also track: when beta_3(T\v)=1 AND H_3(T,T\v)=1, what happens?
    simultaneous_cases = 0

    for trial in range(sample_size):
        if trial % 50 == 0:
            print(f"  ... trial {trial}/{sample_size}")
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)

        betti_T, _, _, _, _, _ = compute_betti_full(T, n, max_p=4)
        b3_T = betti_T.get(3, 0)

        if b3_T == 0:
            continue

        total_b3_eq_1 += 1
        min_sum = float('inf')

        for v in range(n):
            # Compute beta_3(T\v)
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            betti_sub, _, _, _, _, _ = compute_betti_full(T_sub, n-1, max_p=4)
            b3_sub = betti_sub.get(3, 0)

            # Compute H_3(T, T\v)
            rel_betti = compute_relative_betti(T, n, v, max_p=4)
            h3_rel = rel_betti.get(3, 0)

            this_sum = b3_sub + h3_rel
            min_sum = min(min_sum, this_sum)

            if b3_sub == 1 and h3_rel == 1:
                simultaneous_cases += 1

        if min_sum <= 1:
            has_good_vertex += 1
        worst_min_sum = max(worst_min_sum, min_sum)

    print(f"\n  Tournaments with beta_3 = 1: {total_b3_eq_1}")
    print(f"  Has good vertex (min(b3_sub + h3_rel) <= 1): {has_good_vertex}/{total_b3_eq_1}")
    print(f"  Worst min sum: {worst_min_sum}")
    print(f"  Simultaneous beta_3(T\\v)=1 AND H_3(T,T\\v)=1: {simultaneous_cases} vertex-pairs")

    # =====================================================================
    # Part 2: n=6 exhaustive — detailed LES analysis for beta_3=1 cases
    # =====================================================================
    print(f"\n{'='*70}")
    print("--- Part 2: n=6 exhaustive ---")
    n = 6
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    total_b3_eq_1 = 0
    has_good_vertex = 0
    worst_min_sum = 0
    simultaneous_cases = 0

    # Detailed: track what happens for each beta_3=1 tournament
    # For each vertex v: (beta_3(T\v), H_3(T,T\v))
    vertex_pair_dist = Counter()

    for bits in range(n_total):
        if bits % 5000 == 0 and bits > 0:
            print(f"  ... {bits}/{n_total}")

        T = tournament_from_bits(n, bits)
        betti_T, _, _, _, _, _ = compute_betti_full(T, n, max_p=4)
        b3_T = betti_T.get(3, 0)

        if b3_T == 0:
            continue

        total_b3_eq_1 += 1
        min_sum = float('inf')

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            betti_sub, _, _, _, _, _ = compute_betti_full(T_sub, n-1, max_p=4)
            b3_sub = betti_sub.get(3, 0)

            rel_betti = compute_relative_betti(T, n, v, max_p=4)
            h3_rel = rel_betti.get(3, 0)

            vertex_pair_dist[(b3_sub, h3_rel)] += 1

            this_sum = b3_sub + h3_rel
            min_sum = min(min_sum, this_sum)

            if b3_sub == 1 and h3_rel == 1:
                simultaneous_cases += 1

        if min_sum <= 1:
            has_good_vertex += 1
        worst_min_sum = max(worst_min_sum, min_sum)

    print(f"\n  Tournaments with beta_3 = 1: {total_b3_eq_1}")
    print(f"  Has good vertex: {has_good_vertex}/{total_b3_eq_1}")
    print(f"  Worst min sum: {worst_min_sum}")
    print(f"  Simultaneous beta_3(T\\v)=1 AND H_3(T,T\\v)=1: {simultaneous_cases}")
    print(f"\n  (beta_3(T\\v), H_3(T,T\\v)) distribution:")
    for key in sorted(vertex_pair_dist.keys()):
        print(f"    {key}: {vertex_pair_dist[key]}")

    # =====================================================================
    # Part 3: n=8 sampled — does the good vertex condition hold?
    # =====================================================================
    print(f"\n{'='*70}")
    print("--- Part 3: n=8 sampled ---")
    n = 8
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs
    rng = np.random.RandomState(123)
    sample_size = 50

    total_b3_eq_1 = 0
    has_good_vertex = 0
    worst_min_sum = 0
    simultaneous_cases = 0

    for trial in range(sample_size):
        if trial % 10 == 0:
            print(f"  ... trial {trial}/{sample_size}")
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)

        try:
            betti_T, _, _, _, _, _ = compute_betti_full(T, n, max_p=4)
        except np.linalg.LinAlgError:
            continue

        b3_T = betti_T.get(3, 0)
        if b3_T == 0:
            continue

        total_b3_eq_1 += 1
        min_sum = float('inf')

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            try:
                betti_sub, _, _, _, _, _ = compute_betti_full(T_sub, n-1, max_p=4)
                b3_sub = betti_sub.get(3, 0)
                rel_betti = compute_relative_betti(T, n, v, max_p=4)
                h3_rel = rel_betti.get(3, 0)
            except np.linalg.LinAlgError:
                continue

            this_sum = b3_sub + h3_rel
            min_sum = min(min_sum, this_sum)

            if b3_sub == 1 and h3_rel == 1:
                simultaneous_cases += 1

        if min_sum <= 1:
            has_good_vertex += 1
        worst_min_sum = max(worst_min_sum, min_sum)

    print(f"\n  Tournaments with beta_3 = 1: {total_b3_eq_1}")
    print(f"  Has good vertex: {has_good_vertex}/{total_b3_eq_1}")
    print(f"  Worst min sum: {worst_min_sum}")
    print(f"  Simultaneous cases: {simultaneous_cases}")

    # =====================================================================
    # Part 4: Check an even stronger condition — adjacent-odd seesaw
    # =====================================================================
    print(f"\n{'='*70}")
    print("SEESAW INTERACTION WITH GOOD VERTEX")
    print("="*70)
    print("\nBy adjacent-odd seesaw: beta_1 * beta_3 = 0.")
    print("So if beta_3(T) = 1, then beta_1(T) = 0.")
    print("Question: does beta_1(T\\v) = 0 also hold for some/all v?")
    print("If so, the LES at level 1 simplifies.\n")

    n = 6
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    b1_when_b3 = Counter()
    for bits in range(n_total):
        T = tournament_from_bits(n, bits)
        betti_T, _, _, _, _, _ = compute_betti_full(T, n, max_p=4)
        if betti_T.get(3, 0) != 1:
            continue
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            betti_sub, _, _, _, _, _ = compute_betti_full(T_sub, n-1, max_p=4)
            b1_when_b3[(betti_sub.get(1,0), betti_sub.get(3,0))] += 1

    print(f"  When beta_3(T)=1 at n=6, (beta_1(T\\v), beta_3(T\\v)) distribution:")
    for key in sorted(b1_when_b3.keys()):
        print(f"    {key}: {b1_when_b3[key]}")

if __name__ == '__main__':
    main()
