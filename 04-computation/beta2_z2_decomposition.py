#!/usr/bin/env python3
"""beta2_z2_decomposition.py - Analyze Z_2 decomposition into swap cycles

KEY QUESTIONS:
1. Is Z_2 spanned by swap cycles from all vertices?
2. Does the multi-vertex cone (union of B_fill^v for all v) always cover Z_2?
3. What is the structure of Z_2 relative to swap cycles?

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


def all_tournaments_gen(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def compute_z2_and_cones(A, n):
    """Compute Z_2, and check if multi-vertex cone covers it."""
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths3 = enumerate_allowed_paths(A, n, 3)

    if not paths2:
        return None

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    if omega2.ndim < 2 or omega2.shape[1] == 0:
        return None

    D2 = build_full_boundary_matrix([tuple(p) for p in paths2], [tuple(p) for p in paths1])
    D2_om = D2 @ omega2
    svals = np.linalg.svd(D2_om, compute_uv=False)
    rk_d2 = int(sum(s > 1e-8 for s in svals))
    z2_dim = omega2.shape[1] - rk_d2

    if z2_dim == 0:
        return {'z2_dim': 0}

    # Get Z_2 basis (columns in A_2 coordinates)
    _, _, Vt = np.linalg.svd(D2_om, full_matrices=True)
    Z2_basis = omega2 @ Vt[rk_d2:].T  # dim_A2 x z2_dim

    path2_idx = {tuple(p): i for i, p in enumerate(paths2)}

    # For each vertex v, build the cone B matrix
    all_cone_cols = []
    swap_cycle_vecs = []
    per_vertex_swap_dims = []

    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

        if len(arcs_PQ) < 2:
            per_vertex_swap_dims.append(0)
            continue

        m = len(arcs_PQ)
        arc_idx = {arc: i for i, arc in enumerate(arcs_PQ)}

        # Build C
        rows = []
        for a in P:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if a2 == a:
                    row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        for b in Q:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if b2 == b:
                    row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)

        if not rows:
            per_vertex_swap_dims.append(0)
            continue

        C_mat = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C_mat, compute_uv=False)
        rank_C = int(sum(s > 1e-8 for s in Sc))
        ker_dim = m - rank_C
        per_vertex_swap_dims.append(ker_dim)

        if ker_dim == 0:
            continue

        # Get swap cycle vectors in A_2 coordinates
        _, _, Vt_C = np.linalg.svd(C_mat, full_matrices=True)
        ker_basis = Vt_C[rank_C:]

        for ki in range(ker_basis.shape[0]):
            M_vec = ker_basis[ki]
            z = np.zeros(len(paths2))
            for j, (a, b) in enumerate(arcs_PQ):
                coeff = M_vec[j]
                if abs(coeff) < 1e-12:
                    continue
                if (a, b, v) in path2_idx and (v, a, b) in path2_idx:
                    z[path2_idx[(a, b, v)]] += coeff
                    z[path2_idx[(v, a, b)]] -= coeff
            if np.max(np.abs(z)) > 1e-12:
                swap_cycle_vecs.append(z)

        # Build cone B_fill for vertex v
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
        Tp_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]

        for (a, b, c) in Tp_orig:
            col = np.zeros(len(paths2))
            terms = [
                ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
                ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
            ]
            for path, coeff in terms:
                if path in path2_idx:
                    col[path2_idx[path]] += coeff
            if np.max(np.abs(col)) > 1e-12:
                all_cone_cols.append(col)

    # Check: do swap cycles from all vertices span Z_2?
    if swap_cycle_vecs:
        swap_mat = np.column_stack(swap_cycle_vecs)
        rank_swap_all = np.linalg.matrix_rank(swap_mat, tol=1e-8)

        # Does Z_2 lie in span of swap cycles?
        combined = np.hstack([swap_mat, Z2_basis])
        rank_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        z2_in_swap = (rank_combined == rank_swap_all)
    else:
        rank_swap_all = 0
        z2_in_swap = (z2_dim == 0)

    # Check: does multi-vertex cone cover Z_2?
    if all_cone_cols:
        cone_mat = np.column_stack(all_cone_cols)
        rank_cone = np.linalg.matrix_rank(cone_mat, tol=1e-8)

        combined2 = np.hstack([cone_mat, Z2_basis])
        rank_combined2 = np.linalg.matrix_rank(combined2, tol=1e-8)
        z2_in_cone = (rank_combined2 == rank_cone)
    else:
        rank_cone = 0
        z2_in_cone = (z2_dim == 0)

    return {
        'z2_dim': z2_dim,
        'total_swap_dim': rank_swap_all,
        'z2_in_swap_span': z2_in_swap,
        'num_cone_cols': len(all_cone_cols),
        'rank_cone': rank_cone,
        'z2_in_cone_span': z2_in_cone,
        'per_vertex_swap': per_vertex_swap_dims,
    }


# ============================================================
# Part 1: n=5 exhaustive
# ============================================================
print("=" * 70)
print("PART 1: Z_2 DECOMPOSITION AT n=5 (exhaustive)")
print("=" * 70)

n = 5
z2_in_swap_count = 0
z2_in_cone_count = 0
z2_total = 0
swap_deficits = []

for A in all_tournaments_gen(n):
    result = compute_z2_and_cones(A, n)
    if result is None or result['z2_dim'] == 0:
        continue

    z2_total += 1
    if result['z2_in_swap_span']:
        z2_in_swap_count += 1
    else:
        deficit = result['z2_dim'] - result['total_swap_dim']
        swap_deficits.append(deficit)
        if len(swap_deficits) <= 3:
            print(f"  SWAP DEFICIT: z2_dim={result['z2_dim']}, "
                  f"swap_dim={result['total_swap_dim']}, "
                  f"per_v={result['per_vertex_swap']}")

    if result['z2_in_cone_span']:
        z2_in_cone_count += 1
    else:
        print(f"  CONE FAILS: z2_dim={result['z2_dim']}, "
              f"rank_cone={result['rank_cone']}")

print(f"\nn=5: {z2_total} tournaments with z2_dim > 0")
print(f"  Z_2 in span(swap cycles from all v): {z2_in_swap_count}/{z2_total}")
print(f"  Z_2 in span(cone cols from all v): {z2_in_cone_count}/{z2_total}")
if swap_deficits:
    print(f"  Swap deficit range: {min(swap_deficits)} to {max(swap_deficits)}")


# ============================================================
# Part 2: n=6 exhaustive
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: Z_2 DECOMPOSITION AT n=6 (exhaustive)")
print("=" * 70)

import time
n = 6
z2_in_swap_count = 0
z2_in_cone_count = 0
z2_total = 0
swap_deficits = []
t0 = time.time()

for tidx, A in enumerate(all_tournaments_gen(n)):
    result = compute_z2_and_cones(A, n)
    if result is None or result['z2_dim'] == 0:
        continue

    z2_total += 1
    if result['z2_in_swap_span']:
        z2_in_swap_count += 1
    else:
        deficit = result['z2_dim'] - result['total_swap_dim']
        swap_deficits.append(deficit)
        if len(swap_deficits) <= 3:
            print(f"  SWAP DEFICIT T#{tidx}: z2_dim={result['z2_dim']}, "
                  f"swap_dim={result['total_swap_dim']}")

    if result['z2_in_cone_span']:
        z2_in_cone_count += 1
    else:
        print(f"  CONE FAILS T#{tidx}: z2_dim={result['z2_dim']}, "
              f"rank_cone={result['rank_cone']}")

    if (tidx + 1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {tidx + 1}/32768 ({elapsed:.0f}s) swap_ok={z2_in_swap_count} cone_ok={z2_in_cone_count}")

elapsed = time.time() - t0
print(f"\nn=6 ({elapsed:.0f}s): {z2_total} tournaments with z2_dim > 0")
print(f"  Z_2 in span(swap cycles from all v): {z2_in_swap_count}/{z2_total}")
print(f"  Z_2 in span(cone cols from all v): {z2_in_cone_count}/{z2_total}")
if swap_deficits:
    print(f"  Swap deficit range: {min(swap_deficits)} to {max(swap_deficits)}")


# ============================================================
# Part 3: n=7,8 sampled
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: Z_2 DECOMPOSITION AT n=7,8 (sampled)")
print("=" * 70)

for n in [7, 8]:
    random.seed(42)
    num_trials = 300 if n == 7 else 100
    z2_in_swap_count = 0
    z2_in_cone_count = 0
    z2_total = 0
    swap_deficits = []
    t0 = time.time()

    for trial in range(num_trials):
        A = random_tournament(n)
        result = compute_z2_and_cones(A, n)
        if result is None or result['z2_dim'] == 0:
            continue

        z2_total += 1
        if result['z2_in_swap_span']:
            z2_in_swap_count += 1
        else:
            deficit = result['z2_dim'] - result['total_swap_dim']
            swap_deficits.append(deficit)

        if result['z2_in_cone_span']:
            z2_in_cone_count += 1
        else:
            print(f"  CONE FAILS n={n} trial {trial}: z2_dim={result['z2_dim']}, "
                  f"rank_cone={result['rank_cone']}")

        if (trial + 1) % (num_trials // 4) == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial+1}/{num_trials} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): {z2_total}/{num_trials} with z2_dim > 0")
    print(f"  Z_2 in span(swap from all v): {z2_in_swap_count}/{z2_total}")
    print(f"  Z_2 in span(cone from all v): {z2_in_cone_count}/{z2_total}")
    if swap_deficits:
        print(f"  Swap deficits: {len(swap_deficits)}, "
              f"range {min(swap_deficits)} to {max(swap_deficits)}")


print("\n\nDone.")
