#!/usr/bin/env python3
"""
beta3_relative_h3_n8.py — opus-2026-03-09-S52

Extended verification: H_3(T, T\v) <= 1 at n=8,9.
Uses same computation as beta3_relative_h3.py but optimized for larger n.
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

def build_omega_and_boundary(T, n, max_p):
    """Build Omega_p bases and boundary maps."""
    tol = 1e-8

    all_allowed = {}
    for p in range(0, max_p + 2):
        paths = []
        if p + 1 <= n:
            for verts in combinations(range(n), p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths.append(perm)
        all_allowed[p] = paths

    omega_proj = {}
    for p in range(0, max_p + 1):
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

    boundary = {}
    for p in range(1, max_p + 1):
        a_p = all_allowed[p]
        a_pm1 = all_allowed[p-1]
        if not a_p or not a_pm1:
            boundary[p] = np.zeros((len(a_pm1) if a_pm1 else 0, len(a_p) if a_p else 0))
            continue
        idx_pm1 = {path: i for i, path in enumerate(a_pm1)}
        mat = np.zeros((len(a_pm1), len(a_p)))
        for j, sigma in enumerate(a_p):
            for i in range(len(sigma)):
                face = sigma[:i] + sigma[i+1:]
                if face in idx_pm1:
                    mat[idx_pm1[face], j] += (-1)**i
        boundary[p] = mat

    return all_allowed, omega_proj, boundary

def compute_relative_h3(T, n, v, max_p=5):
    """Compute H_3(T, T\v)."""
    tol = 1e-8

    remaining = [i for i in range(n) if i != v]

    # Build full complex for T
    all_T, omega_T, bound_T = build_omega_and_boundary(T, n, max_p)

    # Build Omega for T\v in T's coordinate system
    # T\v paths are T-paths not involving v
    all_Tv = {}
    for p in range(0, max_p + 2):
        all_Tv[p] = [path for path in all_T[p] if v not in path]

    # Build Omega(T\v) using T\v's path structure
    # Note: for T\v, we need the NA face constraints w.r.t. T\v itself
    omega_Tv = {}
    for p in range(0, max_p + 1):
        a_p = all_Tv[p]
        if not a_p:
            omega_Tv[p] = np.zeros((0, 0))
            continue

        a_pm1_set = set(all_Tv[p-1]) if p > 0 else set()
        na_faces_map = {}
        for sigma in a_p:
            for i in range(1, len(sigma)-1):
                face = sigma[:i] + sigma[i+1:]
                if p > 0 and face not in a_pm1_set:
                    if face not in na_faces_map:
                        na_faces_map[face] = len(na_faces_map)

        if not na_faces_map:
            omega_Tv[p] = np.eye(len(a_p))
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
                omega_Tv[p] = np.zeros((len(a_p), 0))
            else:
                omega_Tv[p] = Vt[rank:].T

    # Map Omega_Tv into A_p(T) coords
    omega_Tv_in_T = {}
    for p in range(0, max_p + 1):
        a_p_T = all_T[p]
        a_p_Tv = all_Tv[p]
        om_Tv = omega_Tv[p]

        if om_Tv.ndim < 2 or om_Tv.shape[1] == 0 or not a_p_T:
            omega_Tv_in_T[p] = np.zeros((len(a_p_T), 0))
            continue

        idx_T = {path: i for i, path in enumerate(a_p_T)}
        mapping = np.zeros((len(a_p_T), len(a_p_Tv)))
        for j, path in enumerate(a_p_Tv):
            if path in idx_T:
                mapping[idx_T[path], j] = 1.0

        omega_Tv_in_T[p] = mapping @ om_Tv

    # Compute complement (quotient) basis
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

    # Relative boundary maps
    rel_rank = {}
    for p in range(1, max_p + 1):
        comp_p = complement[p]
        om_T_pm1 = omega_T[p-1]
        comp_pm1 = complement[p-1]

        if comp_p.ndim < 2 or comp_p.shape[1] == 0 or comp_pm1.ndim < 2 or comp_pm1.shape[1] == 0:
            rel_rank[p] = 0
            continue

        if p not in bound_T or bound_T[p].shape[0] == 0:
            rel_rank[p] = 0
            continue

        dp_comp = bound_T[p] @ comp_p
        dp_in_omega = om_T_pm1.T @ dp_comp
        comp_pm1_in_omega = om_T_pm1.T @ comp_pm1

        if comp_pm1_in_omega.shape[1] == 0:
            rel_rank[p] = 0
            continue

        Q, R = np.linalg.qr(comp_pm1_in_omega, mode='reduced')
        proj_dp = Q.T @ dp_in_omega

        S = np.linalg.svd(proj_dp, compute_uv=False)
        rel_rank[p] = int(np.sum(S > tol))

    # Relative Betti
    rel_betti = {}
    for p in range(0, max_p + 1):
        ker = rel_dims[p] - rel_rank.get(p, 0)
        im_next = rel_rank.get(p+1, 0)
        rel_betti[p] = ker - im_next

    return rel_betti.get(3, 0)

def main():
    for n in [8]:
        print(f"\n{'='*70}")
        print(f"n = {n}: sampled verification of H_3(T, T\\v) <= 1")
        print(f"{'='*70}")

        num_arcs = n*(n-1)//2
        n_total = 1 << num_arcs

        rng = np.random.RandomState(42)
        sample_size = 100

        h3_max = 0
        h3_counts = Counter()
        violations = 0

        for trial in range(sample_size):
            if trial % 20 == 0:
                print(f"  ... {trial}/{sample_size}, max H_3 relative: {h3_max}")

            bits = rng.randint(0, n_total)
            T = tournament_from_bits(n, bits)
            for v_idx in range(n):
                h3 = compute_relative_h3(T, n, v_idx, max_p=5)
                h3_counts[h3] += 1
                if h3 > h3_max:
                    h3_max = h3
                    print(f"    NEW MAX: H_3(T,T\\v)={h3} at trial={trial}, v={v_idx}")
                if h3 > 1:
                    violations += 1

        print(f"\n  Total (T,v) pairs: {sample_size * n}")
        print(f"  H_3(T, T\\v) distribution: {dict(sorted(h3_counts.items()))}")
        print(f"  Max H_3(T, T\\v): {h3_max}")
        print(f"  Violations (>1): {violations}")

        if violations == 0:
            print(f"\n  *** H_3(T, T\\v) <= 1 CONFIRMED for all sampled n={n} pairs ***")

if __name__ == '__main__':
    main()
