#!/usr/bin/env python3
"""
beta3_relative_h3.py — opus-2026-03-09-S52

CRITICAL COMPUTATION: Is dim H_3(T, T\v) <= 1 for all (T,v)?

If yes, combined with the LES and beta_2=0, this gives:
  beta_3(T) = rank(i_*) + dim H_3(T,T\v) <= beta_3(T\v) + 1

With base case beta_3=0 at n<=5, induction gives beta_3 <= 1 for all n.

Method: Compute H_3(T,T\v) via the quotient complex Omega_*(T)/Omega_*(T\v).
The relative complex R_p = Omega_p(T) / Omega_p(T\v) consists of
elements of Omega_p(T) not lying in Omega_p(T\v).

Actually more precisely: we embed Omega_p(T\v) into Omega_p(T) via inclusion
of allowed paths, and take the quotient. The relative boundary is induced.
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

def compute_relative_h3(T, n, v, max_p=5):
    """Compute H_3(T, T\v) using the relative chain complex.

    H_p(T, T\v) = ker(d_p^rel) / im(d_{p+1}^rel)
    where d_p^rel: R_p -> R_{p-1} is the induced boundary on R_p = Omega_p(T)/Omega_p(T\v).
    """
    tol = 1e-8

    # Get Omega_p(T) and Omega_p(T\v) for p = 2,3,4,5
    # T\v is the (n-1)-tournament on V\{v}

    remaining = [i for i in range(n) if i != v]
    n_sub = n - 1

    # T\v adjacency
    T_sub = [[T[remaining[i]][remaining[j]] for j in range(n_sub)] for i in range(n_sub)]

    # For each level p, compute:
    # 1. Omega_p(T) basis (as linear combinations of allowed p-paths of T)
    # 2. Omega_p(T\v) basis (mapped into the same coordinate system)
    # 3. Quotient R_p = Omega_p(T) / Omega_p(T\v)

    # First, get all allowed p-paths for T and T\v
    all_allowed_T = {}
    all_allowed_Tv = {}  # in terms of original vertex labels
    for p in range(0, max_p + 2):
        paths_T = []
        if p + 1 <= n:
            for verts in combinations(range(n), p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths_T.append(perm)
        all_allowed_T[p] = paths_T

        # T\v: paths on remaining vertices
        paths_Tv = []
        if p + 1 <= n_sub:
            for verts in combinations(remaining, p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):  # same check since T restricted
                        paths_Tv.append(perm)
        all_allowed_Tv[p] = paths_Tv

    # Build Omega_p(T) for each level
    def build_omega(all_allowed, n_val, max_p_val):
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

    omega_T = build_omega(all_allowed_T, n, max_p)

    # For T\v, we need Omega in the SAME coordinate system as T
    # The allowed paths of T\v are a subset of allowed paths of T
    # So Omega_p(T\v) lives inside the space spanned by A_p(T)
    omega_Tv = build_omega(all_allowed_Tv, n, max_p)

    # Map Omega_p(T\v) into the A_p(T) coordinate system
    omega_Tv_in_T = {}
    for p in range(0, max_p + 1):
        a_p_T = all_allowed_T[p]
        a_p_Tv = all_allowed_Tv[p]
        om_Tv = omega_Tv[p]

        if om_Tv.ndim < 2 or om_Tv.shape[1] == 0 or not a_p_T:
            omega_Tv_in_T[p] = np.zeros((len(a_p_T), 0))
            continue

        # Map: column j of om_Tv is a linear combination of A_p(T\v)
        # We need to express this in terms of A_p(T)
        idx_T = {path: i for i, path in enumerate(a_p_T)}
        mapping = np.zeros((len(a_p_T), len(a_p_Tv)))
        for j, path in enumerate(a_p_Tv):
            if path in idx_T:
                mapping[idx_T[path], j] = 1.0

        omega_Tv_in_T[p] = mapping @ om_Tv

    # Build the relative complex R_p
    # R_p = Omega_p(T) / Omega_p(T\v)
    # We compute R_p as the quotient of the column space of omega_T[p]
    # by the subspace omega_Tv_in_T[p] projected onto omega_T[p]

    # Step 1: Project omega_Tv_in_T[p] onto omega_T[p] to get the subcomplex
    # Step 2: Compute quotient dimensions and maps

    rel_dims = {}
    rel_dp = {}  # relative boundary maps

    # Boundary of T
    boundary_T = {}
    for p in range(1, max_p + 1):
        a_p = all_allowed_T[p]
        a_pm1 = all_allowed_T[p-1]
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

    # For each p, compute relative complex
    for p in range(0, max_p + 1):
        om_T_p = omega_T[p]  # columns span Omega_p(T)
        om_Tv_p = omega_Tv_in_T[p]  # columns span Omega_p(T\v) inside A_p(T)

        if om_T_p.ndim < 2 or om_T_p.shape[1] == 0:
            rel_dims[p] = 0
            continue

        dim_T = om_T_p.shape[1]

        if om_Tv_p.shape[1] == 0:
            # Omega_p(T\v) = 0 in Omega_p(T), so R_p = Omega_p(T)
            rel_dims[p] = dim_T
            continue

        # Project Omega_Tv into Omega_T basis
        # om_Tv_p columns are in A_p(T) coords
        # Project to Omega_T: om_T_p.T @ om_Tv_p gives coords in Omega_T basis
        proj = om_T_p.T @ om_Tv_p  # (dim_T x dim_Tv)

        # Dimension of the image (subspace of Omega_T)
        S = np.linalg.svd(proj, compute_uv=False)
        sub_dim = int(np.sum(S > tol))

        rel_dims[p] = dim_T - sub_dim

    # Compute relative boundary maps
    # d_p^rel: R_p -> R_{p-1}
    # These are induced by the boundary of T, modulo the subcomplex

    # For the actual relative homology computation, we work with:
    # Omega_p(T) in a basis that separates the sub and quotient parts

    # Use the approach: compute d_p on Omega_p(T), then take quotient

    # For each p, find complement basis of Omega_Tv in Omega_T
    complement = {}
    for p in range(0, max_p + 1):
        om_T_p = omega_T[p]
        om_Tv_p = omega_Tv_in_T[p]

        if om_T_p.ndim < 2 or om_T_p.shape[1] == 0:
            complement[p] = np.zeros((0, 0))
            continue

        if om_Tv_p.shape[1] == 0:
            complement[p] = om_T_p
            continue

        proj = om_T_p.T @ om_Tv_p
        U, S, Vt = np.linalg.svd(proj, full_matrices=True)
        sub_rank = int(np.sum(S > tol))
        # Complement = last (dim_T - sub_rank) columns of U
        if sub_rank < om_T_p.shape[1]:
            comp_coords = U[:, sub_rank:]  # coords in Omega_T basis
            complement[p] = om_T_p @ comp_coords  # back to A_p(T) coords
        else:
            complement[p] = np.zeros((om_T_p.shape[0], 0))

    # Now compute relative boundary d_p^rel using complement basis
    # d_p^rel is the boundary restricted to complement(p) projected to complement(p-1)
    # But we need to be careful: d_p maps to Omega_{p-1}(T), and we project
    # to the quotient R_{p-1} = complement part.

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

        # d_p on complement: boundary_T[p] @ comp_p -> A_{p-1}(T)
        dp_comp = boundary_T[p] @ comp_p  # in A_{p-1}(T) coords

        # Project to Omega_{p-1}(T): om_T_pm1.T @ dp_comp
        dp_in_omega = om_T_pm1.T @ dp_comp  # in Omega_{p-1}(T) coords

        # Project further to complement_{p-1} (quotient part)
        # complement_{p-1} = om_T_{p-1} @ comp_pm1_coords
        # We need the coordinates of dp_in_omega in the split basis:
        # Omega_{p-1} = sub + comp
        # To project to comp part, we need the projection matrix

        # comp_pm1 is in A_{p-1} coords
        # comp_pm1_in_omega = om_T_pm1.T @ comp_pm1
        comp_pm1_in_omega = om_T_pm1.T @ comp_pm1

        # dp_in_omega projected onto comp_pm1_in_omega
        # Since comp_pm1_in_omega columns may not be orthonormal:
        Q, R = np.linalg.qr(comp_pm1_in_omega, mode='reduced')
        # Project: Q @ Q.T @ dp_in_omega
        proj_dp = Q.T @ dp_in_omega  # relative boundary matrix (rel_dim_{p-1} x rel_dim_p)

        S = np.linalg.svd(proj_dp, compute_uv=False)
        rel_rank[p] = int(np.sum(S > tol))

    # Relative Betti numbers
    rel_betti = {}
    for p in range(0, max_p + 1):
        ker = rel_dims[p] - rel_rank.get(p, 0)
        im_next = rel_rank.get(p+1, 0)
        rel_betti[p] = ker - im_next

    return rel_betti, rel_dims, rel_rank

def main():
    print("=" * 70)
    print("RELATIVE H_3(T, T\\v) COMPUTATION")
    print("Is dim H_3(T, T\\v) <= 1 for all (T,v)?")
    print("=" * 70)

    # n=6 exhaustive
    n = 6
    print(f"\nn = {n}: exhaustive")
    num_arcs = n*(n-1)//2
    total = 1 << num_arcs

    h3_rel_max = 0
    h3_rel_counts = Counter()
    violations = 0

    for bits in range(total):
        if bits % 5000 == 0 and bits > 0:
            print(f"  ... {bits}/{total}, max H_3 relative so far: {h3_rel_max}")

        T = tournament_from_bits(n, bits)
        for v in range(n):
            rel_betti, rel_dims, rel_rank = compute_relative_h3(T, n, v)
            h3 = rel_betti.get(3, 0)
            h3_rel_counts[h3] += 1
            if h3 > h3_rel_max:
                h3_rel_max = h3
            if h3 > 1:
                violations += 1
                print(f"    VIOLATION: bits={bits}, v={v}, H_3(T,T\\v)={h3}")
                print(f"      rel_dims: {rel_dims}")
                print(f"      rel_rank: {rel_rank}")
                print(f"      rel_betti: {rel_betti}")

    print(f"\n  H_3(T, T\\v) distribution: {dict(sorted(h3_rel_counts.items()))}")
    print(f"  Max H_3(T, T\\v): {h3_rel_max}")
    print(f"  Violations (>1): {violations}")

    if violations == 0:
        print(f"\n  *** H_3(T, T\\v) <= 1 CONFIRMED for all n={n} (T,v) pairs ***")
        print(f"  This gives beta_3(T) <= beta_3(T\\v) + 1 by LES")
        print(f"  Combined with beta_3 = 0 at n<=5 (base), gives beta_3 <= 1 for all n")
        print(f"  (modulo verification at n>=7)")

    # n=7 sampled
    print(f"\n{'='*70}")
    n = 7
    print(f"n = {n}: sampled")
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    rng = np.random.RandomState(42)
    sample_size = 200

    h3_rel_max7 = 0
    h3_rel_counts7 = Counter()
    violations7 = 0

    for trial in range(sample_size):
        if trial % 50 == 0 and trial > 0:
            print(f"  ... {trial}/{sample_size}, max H_3 relative: {h3_rel_max7}")

        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)
        for v in range(n):
            rel_betti, rel_dims, rel_rank = compute_relative_h3(T, n, v)
            h3 = rel_betti.get(3, 0)
            h3_rel_counts7[h3] += 1
            if h3 > h3_rel_max7:
                h3_rel_max7 = h3
            if h3 > 1:
                violations7 += 1
                print(f"    VIOLATION: bits={bits}, v={v}, H_3(T,T\\v)={h3}")

    print(f"\n  H_3(T, T\\v) distribution: {dict(sorted(h3_rel_counts7.items()))}")
    print(f"  Max H_3(T, T\\v): {h3_rel_max7}")
    print(f"  Violations (>1): {violations7}")

    if violations7 == 0:
        print(f"\n  *** H_3(T, T\\v) <= 1 CONFIRMED for all sampled n={n} (T,v) pairs ***")

if __name__ == '__main__':
    main()
