#!/usr/bin/env python3
"""Analyze WHY H_2(T, T\v) = 0 for all tournaments.

The relative chain complex Ω_p(T)/Ω_p(T\v) can be understood as:
  Ω_p^{rel}(v) = chains in Ω_p(T) that essentially use vertex v.

For H_2^{rel}(v) = 0, we need every 2-cycle in the relative complex
to be a relative boundary.

Key questions:
1. What does Ω_2^{rel}(v) look like? (dimension, generators)
2. What does the relative boundary ∂_2^{rel} look like?
3. What fills relative 2-cycles?
"""
import numpy as np
from itertools import combinations, permutations
import sys
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
    boundary_coeffs
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

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def compute_rel_dims(A, n, v):
    """Compute dimensions of relative chain complex at vertex v."""
    # Compute Ω_p(T) for p=0,...,4
    a = {}
    for p in range(5):
        a[p] = [tuple(x) for x in enumerate_allowed_paths(A, n, p)]

    om = {}
    for p in range(5):
        if p == 0:
            om[p] = np.eye(n)
        else:
            om[p] = compute_omega_basis(A, n, p, a[p], a[p-1])

    results = {}

    for p in [1, 2, 3]:
        dim_om = om[p].shape[1] if om[p].ndim == 2 and om[p].shape[1] > 0 else 0
        if dim_om == 0:
            results[p] = {'dim_om': 0, 'dim_V': 0, 'dim_rel': 0}
            continue

        # v-path indices
        v_idx = [i for i in range(len(a[p])) if v in a[p][i]]
        nov_idx = [i for i in range(len(a[p])) if v not in a[p][i]]

        # Ω_p ∩ V_p (no-v subspace of Ω_p)
        if nov_idx:
            R_nov = om[p][nov_idx, :]  # project onto no-v rows
        else:
            R_nov = np.zeros((0, dim_om))

        if v_idx:
            R_v = om[p][v_idx, :]  # project onto v rows
        else:
            R_v = np.zeros((0, dim_om))

        # dim(Ω_p ∩ V_p) = dim(ker R_v)
        S_v = np.linalg.svd(R_v, compute_uv=False) if R_v.shape[0] > 0 else []
        rank_v = sum(s > 1e-8 for s in S_v)
        dim_V = dim_om - rank_v
        dim_rel = rank_v  # = dim_om - dim_V

        # Count: how many paths use v?
        n_v_paths = len(v_idx)
        n_nov_paths = len(nov_idx)

        # Among v-paths, classify by position of v
        pos_counts = Counter()
        for i in v_idx:
            path = a[p][i]
            pos = path.index(v)
            pos_counts[pos] += 1

        results[p] = {
            'dim_om': dim_om,
            'dim_V': dim_V,
            'dim_rel': dim_rel,
            'n_v_paths': n_v_paths,
            'n_nov_paths': n_nov_paths,
            'v_positions': dict(pos_counts)
        }

    return results

# ===== Part 1: Relative complex dimensions =====
print("=" * 70)
print("RELATIVE CHAIN COMPLEX DIMENSIONS")
print("=" * 70)

n = 5
rel_dim_patterns = Counter()

for A in all_tournaments_gen(n):
    for v in range(n):
        res = compute_rel_dims(A, n, v)
        pattern = tuple(res[p]['dim_rel'] for p in [1, 2, 3])
        rel_dim_patterns[pattern] += 1

print(f"n={n}: Relative dimension patterns (dim_1^rel, dim_2^rel, dim_3^rel):")
for pattern in sorted(rel_dim_patterns):
    count = rel_dim_patterns[pattern]
    print(f"  {pattern}: {count}")

# ===== Part 2: Is the relative complex always exact at dim 2? =====
print(f"\n\n{'='*70}")
print("RELATIVE COMPLEX EXACTNESS AT DIM 2")
print("=" * 70)

n = 5
exact_count = 0
nonexact_count = 0

for A in all_tournaments_gen(n):
    a = {}
    for p in range(5):
        a[p] = [tuple(x) for x in enumerate_allowed_paths(A, n, p)]

    om = {}
    for p in range(5):
        if p == 0:
            om[p] = np.eye(n)
        else:
            om[p] = compute_omega_basis(A, n, p, a[p], a[p-1])

    for v in range(n):
        dim_om2 = om[2].shape[1] if om[2].ndim == 2 else 0
        dim_om3 = om[3].shape[1] if om[3].ndim == 2 else 0

        if dim_om2 == 0:
            exact_count += 1
            continue

        # v-path filter for Ω_2 and Ω_3
        v_idx_2 = [i for i in range(len(a[2])) if v in a[2][i]]
        v_idx_3 = [i for i in range(len(a[3])) if v in a[3][i]]

        # R_v matrices (project onto v-rows)
        R_v2 = om[2][v_idx_2, :] if v_idx_2 else np.zeros((0, dim_om2))
        Sv2 = np.linalg.svd(R_v2, compute_uv=False) if R_v2.shape[0] > 0 else []
        rank_v2 = sum(s > 1e-8 for s in Sv2)
        dim_rel2 = rank_v2

        if dim_rel2 == 0:
            exact_count += 1
            continue

        # ∂_2 on Ω_2
        bd2 = build_full_boundary_matrix(a[2], a[1])
        bd2_om = bd2 @ om[2]
        S2 = np.linalg.svd(bd2_om, compute_uv=False)
        rank_bd2 = sum(s > 1e-8 for s in S2)
        ker_dim2 = dim_om2 - rank_bd2

        # V_2 = ker(R_v2) in Ω_2
        _, _, Vt_v2 = np.linalg.svd(R_v2)
        V2_basis = Vt_v2[rank_v2:].T if rank_v2 < dim_om2 else np.zeros((dim_om2, 0))

        # ker(∂_2) basis
        _, _, Vt_bd2 = np.linalg.svd(bd2_om)
        ker2_basis = Vt_bd2[rank_bd2:].T if rank_bd2 < dim_om2 else np.zeros((dim_om2, 0))

        # ker(∂_2) ∩ V_2
        if ker2_basis.shape[1] > 0 and V2_basis.shape[1] > 0:
            combined = np.hstack([ker2_basis, V2_basis])
            rc = np.linalg.matrix_rank(combined, tol=1e-8)
            intersect_dim = ker2_basis.shape[1] + V2_basis.shape[1] - rc
        else:
            intersect_dim = 0

        ker_rel2 = ker_dim2 - intersect_dim

        # im(∂_3) in relative sense
        if dim_om3 > 0:
            bd3 = build_full_boundary_matrix(a[3], a[2])
            bd3_om = bd3 @ om[3]
            im3_om2, _, _, _ = np.linalg.lstsq(om[2], bd3_om, rcond=None)

            if V2_basis.shape[1] > 0:
                combined_im = np.hstack([im3_om2, V2_basis])
                rc_im = np.linalg.matrix_rank(combined_im, tol=1e-8)
                im3_rel = rc_im - V2_basis.shape[1]
            else:
                im3_rel = np.linalg.matrix_rank(im3_om2, tol=1e-8)
        else:
            im3_rel = 0

        h2_rel = max(0, ker_rel2 - im3_rel)

        if h2_rel == 0:
            exact_count += 1
        else:
            nonexact_count += 1
            out_v = sum(A[v])
            in_v = n - 1 - out_v
            print(f"  NON-EXACT at v={v}, out_deg={out_v}: H_2^rel={h2_rel}")
            print(f"    dim_rel2={dim_rel2}, ker_rel2={ker_rel2}, im3_rel={im3_rel}")

print(f"\n  Exact: {exact_count}, Non-exact: {nonexact_count}")

# ===== Part 3: What makes H_2^rel = 0? =====
print(f"\n\n{'='*70}")
print("WHY IS H_2^rel = 0? DIMENSION ANALYSIS")
print("=" * 70)

# Two cases that would give H_2^rel = 0:
# (a) ker_rel2 = 0 (no relative 2-cycles)
# (b) im3_rel = ker_rel2 (all relative cycles are boundaries)

n = 5
case_a = 0
case_b = 0
case_c = 0  # something else

for A in all_tournaments_gen(n):
    a = {}
    for p in range(5):
        a[p] = [tuple(x) for x in enumerate_allowed_paths(A, n, p)]

    om = {}
    for p in range(5):
        if p == 0:
            om[p] = np.eye(n)
        else:
            om[p] = compute_omega_basis(A, n, p, a[p], a[p-1])

    for v in range(n):
        dim_om2 = om[2].shape[1] if om[2].ndim == 2 else 0
        if dim_om2 == 0:
            case_a += 1
            continue

        v_idx_2 = [i for i in range(len(a[2])) if v in a[2][i]]
        R_v2 = om[2][v_idx_2, :] if v_idx_2 else np.zeros((0, dim_om2))
        Sv2 = np.linalg.svd(R_v2, compute_uv=False) if R_v2.shape[0] > 0 else []
        rank_v2 = sum(s > 1e-8 for s in Sv2)

        if rank_v2 == 0:
            case_a += 1
            continue

        _, _, Vt_v2 = np.linalg.svd(R_v2)
        V2_basis = Vt_v2[rank_v2:].T if rank_v2 < dim_om2 else np.zeros((dim_om2, 0))

        bd2 = build_full_boundary_matrix(a[2], a[1])
        bd2_om = bd2 @ om[2]
        _, S2, Vt_bd2 = np.linalg.svd(bd2_om)
        rank_bd2 = sum(s > 1e-8 for s in S2)
        ker2_basis = Vt_bd2[rank_bd2:].T if rank_bd2 < dim_om2 else np.zeros((dim_om2, 0))

        if ker2_basis.shape[1] > 0 and V2_basis.shape[1] > 0:
            combined = np.hstack([ker2_basis, V2_basis])
            rc = np.linalg.matrix_rank(combined, tol=1e-8)
            intersect_dim = ker2_basis.shape[1] + V2_basis.shape[1] - rc
        else:
            intersect_dim = 0

        ker_rel2 = ker2_basis.shape[1] - intersect_dim

        if ker_rel2 == 0:
            case_a += 1  # No relative 2-cycles (trivially exact)
        else:
            case_b += 1  # Has relative cycles but all are boundaries

print(f"n={n}:")
print(f"  Case (a) ker_rel2=0 (no relative cycles): {case_a}")
print(f"  Case (b) ker_rel2>0 (cycles exist, all filled): {case_b}")

# Does case (a) always hold? That would be the simplest!
print(f"\n  Is ker_rel2 ALWAYS 0? {'YES' if case_b == 0 else 'NO'}")

# ===== Part 4: Check ker_rel2 at n=6 =====
print(f"\n\n{'='*70}")
print("CHECKING ker_rel2 AT n=6")
print("=" * 70)

n = 6
case_a_6 = 0
case_b_6 = 0

for tidx, A in enumerate(all_tournaments_gen(n)):
    if tidx % 5000 == 0:
        print(f"  ... {tidx}", flush=True)

    a = {}
    for p in range(5):
        a[p] = [tuple(x) for x in enumerate_allowed_paths(A, n, p)]

    om = {}
    for p in range(5):
        if p == 0:
            om[p] = np.eye(n)
        else:
            om[p] = compute_omega_basis(A, n, p, a[p], a[p-1])

    for v in range(n):
        dim_om2 = om[2].shape[1] if om[2].ndim == 2 else 0
        if dim_om2 == 0:
            case_a_6 += 1
            continue

        v_idx_2 = [i for i in range(len(a[2])) if v in a[2][i]]
        R_v2 = om[2][v_idx_2, :] if v_idx_2 else np.zeros((0, dim_om2))
        Sv2 = np.linalg.svd(R_v2, compute_uv=False) if R_v2.shape[0] > 0 else []
        rank_v2 = sum(s > 1e-8 for s in Sv2)

        if rank_v2 == 0:
            case_a_6 += 1
            continue

        _, _, Vt_v2 = np.linalg.svd(R_v2)
        V2_basis = Vt_v2[rank_v2:].T if rank_v2 < dim_om2 else np.zeros((dim_om2, 0))

        bd2 = build_full_boundary_matrix(a[2], a[1])
        bd2_om = bd2 @ om[2]
        _, S2, Vt_bd2 = np.linalg.svd(bd2_om)
        rank_bd2 = sum(s > 1e-8 for s in S2)
        ker2_basis = Vt_bd2[rank_bd2:].T if rank_bd2 < dim_om2 else np.zeros((dim_om2, 0))

        if ker2_basis.shape[1] > 0 and V2_basis.shape[1] > 0:
            combined = np.hstack([ker2_basis, V2_basis])
            rc = np.linalg.matrix_rank(combined, tol=1e-8)
            intersect_dim = ker2_basis.shape[1] + V2_basis.shape[1] - rc
        else:
            intersect_dim = 0

        ker_rel2 = ker2_basis.shape[1] - intersect_dim

        if ker_rel2 == 0:
            case_a_6 += 1
        else:
            case_b_6 += 1

print(f"\nn=6:")
print(f"  Case (a) ker_rel2=0 (no relative cycles): {case_a_6}")
print(f"  Case (b) ker_rel2>0 (cycles exist, all filled): {case_b_6}")
print(f"  Is ker_rel2 ALWAYS 0? {'YES' if case_b_6 == 0 else 'NO'}")

print("\nDone.")
