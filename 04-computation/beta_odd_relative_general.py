#!/usr/bin/env python3
"""
beta_odd_relative_general.py — opus-2026-03-09-S52

Check H_{2k+1}(T, T\v) <= 1 for all k.
If true for all k, gives Boolean odd Betti for all tournaments.

Focus on H_1, H_3, H_5 at n=7,8.
We know H_1(T,T\v) <= 1 (from beta_1 <= 1 proof).
We verified H_3(T,T\v) <= 1.
Now check H_5(T,T\v).
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter
import sys

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

def compute_all_relative_betti(T, n, v, max_p=6):
    """Compute all relative Betti numbers H_p(T, T\v) for p up to max_p."""
    tol = 1e-8

    remaining = [i for i in range(n) if i != v]

    # Build allowed paths for T and T\v
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

    # Build Omega projectors
    def build_omega(all_allowed, max_p_val):
        omega_proj = {}
        for p in range(0, max_p_val + 1):
            a_p = all_allowed[p]
            if not a_p:
                omega_proj[p] = np.zeros((0, 0))
                continue

            a_pm1_set = set(all_allowed[p-1]) if p > 0 else set()
            na_faces_map = {}
            for sigma in a_p:
                for i in range(1, len(sigma)-1):
                    face = sigma[:i] + sigma[i+1:]
                    if p > 0 and face not in a_pm1_set:
                        if face not in na_faces_map:
                            na_faces_map[face] = len(na_faces_map)

            if not na_faces_map:
                omega_proj[p] = np.eye(len(a_p))
            else:
                na_mat = np.zeros((len(na_faces_map), len(a_p)))
                for j, sigma in enumerate(a_p):
                    for i in range(1, len(sigma)-1):
                        face = sigma[:i] + sigma[i+1:]
                        if face in na_faces_map:
                            na_mat[na_faces_map[face], j] += (-1)**i

                U, S, Vt = np.linalg.svd(na_mat, full_matrices=True)
                rank = int(np.sum(S > tol))
                null_dim = len(a_p) - rank
                if null_dim == 0:
                    omega_proj[p] = np.zeros((len(a_p), 0))
                else:
                    omega_proj[p] = Vt[rank:].T
        return omega_proj

    omega_T = build_omega(all_T, max_p)
    omega_Tv = build_omega(all_Tv, max_p)

    # Map Omega_Tv into A_p(T) coordinate system
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

    # Build boundary maps for T
    boundary_T = {}
    for p in range(1, max_p + 1):
        a_p = all_T[p]
        a_pm1 = all_T[p-1]
        if not a_p or not a_pm1:
            boundary_T[p] = np.zeros((len(a_pm1) if a_pm1 else 0, len(a_p) if a_p else 0))
            continue
        idx_pm1 = {path: i for i, path in enumerate(a_pm1)}
        mat = np.zeros((len(a_pm1), len(a_p)))
        for j, sigma in enumerate(a_p):
            for i in range(len(sigma)):
                face = sigma[:i] + sigma[i+1:]
                if face in idx_pm1:
                    mat[idx_pm1[face], j] += (-1)**i
        boundary_T[p] = mat

    # Compute complement (quotient) bases and relative boundary maps
    complement = {}
    rel_dims = {}
    for p in range(0, max_p + 1):
        om_T_p = omega_T[p]
        om_Tv_p = omega_Tv_in_T[p]

        if om_T_p.ndim < 2 or om_T_p.shape[1] == 0:
            complement[p] = np.zeros((0, 0))
            rel_dims[p] = 0
            continue

        if om_Tv_p.shape[1] == 0:
            complement[p] = om_T_p
            rel_dims[p] = om_T_p.shape[1]
            continue

        proj = om_T_p.T @ om_Tv_p
        U, S, Vt = np.linalg.svd(proj, full_matrices=True)
        sub_rank = int(np.sum(S > tol))

        if sub_rank < om_T_p.shape[1]:
            comp_coords = U[:, sub_rank:]
            complement[p] = om_T_p @ comp_coords
            rel_dims[p] = comp_coords.shape[1]
        else:
            complement[p] = np.zeros((om_T_p.shape[0], 0))
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
    print("GENERAL RELATIVE ODD BETTI: H_{2k+1}(T, T\\v) <= 1?")
    print("=" * 70)

    for n in [7, 8]:
        print(f"\n{'='*70}")
        print(f"n = {n}")
        print(f"{'='*70}")

        num_arcs = n*(n-1)//2
        n_total = 1 << num_arcs
        rng = np.random.RandomState(42)
        sample_size = 100 if n == 8 else 200

        max_p = min(n-1, 6)

        # Track max values for each odd p
        max_rel = {p: 0 for p in range(1, max_p+1, 2)}
        counts = {p: Counter() for p in range(1, max_p+1, 2)}
        violations = {p: 0 for p in range(1, max_p+1, 2)}

        for trial in range(sample_size):
            if trial % 50 == 0:
                print(f"  ... {trial}/{sample_size}")

            bits = rng.randint(0, n_total)
            T = tournament_from_bits(n, bits)
            for v in range(n):
                rel_betti = compute_all_relative_betti(T, n, v, max_p=max_p)
                for p in range(1, max_p+1, 2):
                    h = rel_betti.get(p, 0)
                    counts[p][h] += 1
                    if h > max_rel[p]:
                        max_rel[p] = h
                    if h > 1:
                        violations[p] += 1

        for p in range(1, max_p+1, 2):
            print(f"\n  H_{p}(T, T\\v) distribution: {dict(sorted(counts[p].items()))}")
            print(f"  Max: {max_rel[p]}, Violations (>1): {violations[p]}")

        # Also check even relative Betti for comparison
        print(f"\n  --- Even relative Betti (for comparison) ---")
        max_rel_even = {p: 0 for p in range(2, max_p+1, 2)}
        counts_even = {p: Counter() for p in range(2, max_p+1, 2)}

        for trial in range(min(50, sample_size)):
            bits = rng.randint(0, n_total)
            T = tournament_from_bits(n, bits)
            for v in range(n):
                rel_betti = compute_all_relative_betti(T, n, v, max_p=max_p)
                for p in range(2, max_p+1, 2):
                    h = rel_betti.get(p, 0)
                    counts_even[p][h] += 1
                    if h > max_rel_even[p]:
                        max_rel_even[p] = h

        for p in range(2, max_p+1, 2):
            print(f"  H_{p}(T, T\\v) distribution: {dict(sorted(counts_even[p].items()))}")
            print(f"  Max: {max_rel_even[p]}")

if __name__ == '__main__':
    main()
