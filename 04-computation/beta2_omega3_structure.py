#!/usr/bin/env python3
"""
Ω_3 structure analysis for β_2 = 0 proof.

KEY QUESTION: Is Ω_3 = span(DT paths) = span(4-clique paths)?
If so, the proof reduces to showing DT boundaries span ker(∂_2|Ω_2).

DT 4-path (a,b,c,d): a→b→c→d, a→c, b→d.
4-clique path: DT + a→d.

Also check: are there Ω_3 elements that are NOT individual A_3 paths?
(I.e., linear combinations of 3-paths needed for cancellation?)
"""
import numpy as np
from itertools import combinations, permutations
import sys, time
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def classify_3paths(A, n):
    """Classify 3-paths into types."""
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a2_set = set(tuple(p) for p in a2)

    clique_paths = []  # 4-clique (transitive 4-set)
    dt_nonclique = []  # DT but not 4-clique (d→a)
    other_a3 = []      # allowed 3-path but not DT

    for p in a3:
        a, b, c, d = tuple(p)
        # Already have a→b→c→d (allowed 3-path)
        has_ac = A[a][c] == 1
        has_bd = A[b][d] == 1
        has_ad = A[a][d] == 1

        if has_ac and has_bd:
            if has_ad:
                clique_paths.append(tuple(p))
            else:
                dt_nonclique.append(tuple(p))
        else:
            other_a3.append(tuple(p))

    return clique_paths, dt_nonclique, other_a3

def analyze_omega3(A, n):
    """Compare Ω_3 with DT and 4-clique spans."""
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]
    a3_idx = {p: i for i, p in enumerate(a3_list)}

    clique, dt_nc, other = classify_3paths(A, n)

    # Ω_3 = compute_omega_basis
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    # Which A_3 paths are individually in Ω_3?
    # A path p is in Ω_3 iff ALL its faces are in A_2.
    # Face d_i of (a,b,c,d): omit i-th vertex.
    # d_0 = (b,c,d), d_1 = (a,c,d), d_2 = (a,b,d), d_3 = (a,b,c)
    a2_set = set(a2_list)

    indiv_in_om3 = []
    for p in a3_list:
        a, b, c, d = p
        faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
        all_in_a2 = all(f in a2_set for f in faces)
        if all_in_a2:
            indiv_in_om3.append(p)

    # Compute ranks of subspaces
    # Embed each set as columns in A_3 coordinate space
    def embed(paths, a3_idx, dim):
        if not paths:
            return np.zeros((dim, 0))
        M = np.zeros((dim, len(paths)))
        for j, p in enumerate(paths):
            if p in a3_idx:
                M[a3_idx[p], j] = 1.0
        return M

    dim_a3 = len(a3_list)

    M_clique = embed(clique, a3_idx, dim_a3)
    M_dt = embed(clique + dt_nc, a3_idx, dim_a3)
    M_indiv = embed(indiv_in_om3, a3_idx, dim_a3)

    rank_clique = np.linalg.matrix_rank(M_clique) if M_clique.shape[1] > 0 else 0
    rank_dt = np.linalg.matrix_rank(M_dt) if M_dt.shape[1] > 0 else 0
    rank_indiv = np.linalg.matrix_rank(M_indiv) if M_indiv.shape[1] > 0 else 0

    # Check if indiv_in_om3 == DT paths
    indiv_set = set(indiv_in_om3)
    dt_set = set(clique + dt_nc)

    return {
        'dim_om3': dim_om3,
        'n_clique': len(clique),
        'n_dt_nc': len(dt_nc),
        'n_other': len(other),
        'n_indiv': len(indiv_in_om3),
        'rank_clique': rank_clique,
        'rank_dt': rank_dt,
        'rank_indiv': rank_indiv,
        'indiv_eq_dt': indiv_set == dt_set,
        'indiv_subset_dt': indiv_set <= dt_set,
    }

# ===== n=5 exhaustive =====
print("=" * 70)
print("Ω_3 STRUCTURE ANALYSIS")
print("=" * 70)

for n in [4, 5]:
    print(f"\n--- n={n} ---")
    type_dist = Counter()
    mismatch_count = 0

    for A in all_tournaments_gen(n):
        info = analyze_omega3(A, n)

        key = (info['dim_om3'], info['n_clique'], info['n_dt_nc'],
               info['n_indiv'], info['rank_clique'], info['rank_dt'],
               info['rank_indiv'], info['indiv_eq_dt'])
        type_dist[key] += 1

        # Check: is Ω_3 = span(individual paths)?
        if info['dim_om3'] != info['rank_indiv']:
            mismatch_count += 1

    print(f"  Total: {sum(type_dist.values())} tournaments")
    print(f"  Ω_3 ≠ span(individual paths): {mismatch_count}")
    print(f"\n  (dim_Ω3, #clique, #dt_nc, #indiv, rk_cliq, rk_dt, rk_indiv, indiv==dt): count")
    for key in sorted(type_dist.keys()):
        dim_om3, nc, ndt, ni, rc, rd, ri, eq = key
        print(f"    ({dim_om3}, {nc}, {ndt}, {ni}, {rc}, {rd}, {ri}, {'Y' if eq else 'N'}): {type_dist[key]}")

# ===== KEY TEST: Is DT always = individually-in-Ω_3? =====
print(f"\n\n{'='*70}")
print("KEY QUESTION: Is individually-in-Ω_3 = DT paths?")
print("="*70)

for n in [4, 5]:
    always_eq = True
    for A in all_tournaments_gen(n):
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)
        a2_set = set(tuple(p) for p in a2)
        a3_list = [tuple(p) for p in a3]

        clique, dt_nc, other = classify_3paths(A, n)
        dt_set = set(clique + dt_nc)

        indiv = set()
        for p in a3_list:
            a, b, c, d = p
            faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
            if all(f in a2_set for f in faces):
                indiv.add(p)

        if indiv != dt_set:
            always_eq = False
            # Show the difference
            extra_indiv = indiv - dt_set
            missing_indiv = dt_set - indiv
            if extra_indiv or missing_indiv:
                print(f"  n={n}: indiv\\DT = {extra_indiv}, DT\\indiv = {missing_indiv}")
                # Check what the extra paths look like
                for p in extra_indiv:
                    a, b, c, d = p
                    print(f"    {p}: a→c? {A[a][c]}, b→d? {A[b][d]}, a→d? {A[a][d]}")
                break

    if always_eq:
        print(f"  n={n}: YES, always equal! DT = {{p ∈ A_3 : all faces in A_2}}")

# ===== Does Ω_3 ever have cancellation chains (non-individual elements)? =====
print(f"\n\n{'='*70}")
print("Does Ω_3 have non-individual elements (cancellation 3-chains)?")
print("="*70)

for n in [4, 5]:
    has_cancel = 0
    total = 0
    for A in all_tournaments_gen(n):
        total += 1
        info = analyze_omega3(A, n)
        if info['dim_om3'] > info['rank_indiv']:
            has_cancel += 1

    print(f"  n={n}: {has_cancel}/{total} tournaments have Ω_3 > span(individual paths)")

# ===== MOST IMPORTANT: Does span(DT boundaries) = ker(∂_2|Ω_2)? =====
print(f"\n\n{'='*70}")
print("Does im(∂_3|DT) = ker(∂_2|Ω_2)?")
print("="*70)

for n in [4, 5]:
    failures = 0
    total_with_ker = 0
    for A in all_tournaments_gen(n):
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)
        a2_list = [tuple(p) for p in a2]
        a3_list = [tuple(p) for p in a3]

        om2 = compute_omega_basis(A, n, 2, a2, a1)
        dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
        if dim_om2 == 0:
            continue

        bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
        bd2_om = bd2 @ om2
        rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
        ker_dim = dim_om2 - rank2
        if ker_dim == 0:
            continue

        total_with_ker += 1

        # DT paths
        clique, dt_nc, other = classify_3paths(A, n)
        dt_all = clique + dt_nc

        if not dt_all:
            failures += 1
            continue

        # ∂_3 on DT paths, projected to Ω_2
        bd3 = build_full_boundary_matrix(a3_list, a2_list)
        a3_idx = {p: i for i, p in enumerate(a3_list)}

        # Build DT boundary in A_2 coords
        dt_bds = np.zeros((len(a2_list), len(dt_all)))
        for j, p in enumerate(dt_all):
            if p in a3_idx:
                dt_bds[:, j] = bd3[:, a3_idx[p]]

        # Project to Ω_2 coords
        dt_bds_om2, _, _, _ = np.linalg.lstsq(om2, dt_bds, rcond=None)
        err = np.max(np.abs(om2 @ dt_bds_om2 - dt_bds))
        if err > 1e-6:
            print(f"  WARNING: DT bd not in Ω_2, err={err:.2e}")

        im_dt_rank = np.linalg.matrix_rank(dt_bds_om2, tol=1e-8)
        if im_dt_rank < ker_dim:
            failures += 1
            print(f"  FAILURE at n={n}: ker={ker_dim}, im_DT={im_dt_rank}")

    print(f"  n={n}: {total_with_ker} tournaments with ker>0, {failures} failures")

print("\nDone.")
