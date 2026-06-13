#!/usr/bin/env python3
"""
beta2_relative_n4_check.py - Check H_2(T, T\\v) at n=4

Key discovery: H_2(T, T\\v) might be NONZERO at n=4!
Specifically, when v is a sink/source in a tournament with a 3-cycle,
removing v "exposes" a cycle that was filled by paths through v.

The inductive proof of beta_2 = 0 still works:
  Base: beta_2 = 0 for n <= 4 (direct)
  Step: H_2(T, T\\v) = 0 for n >= 5 => beta_2(T) = 0

This script verifies the n=4 situation exhaustively.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, path_betti_numbers
)
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def compute_relative_h2(A, n, v):
    """Correct relative H_2 computation from beta2_relative_homology.py."""
    paths = {}
    omega = {}
    for p in range(5):
        paths[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif len(paths[p]) > 0 and len(paths[p-1]) > 0:
            omega[p] = compute_omega_basis(A, n, p, paths[p], paths[p-1])
        else:
            omega[p] = np.zeros((max(1, len(paths[p])), 0))

    dim_om2 = omega[2].shape[1] if omega[2].ndim == 2 else 0
    dim_om3 = omega[3].shape[1] if omega[3].ndim == 2 else 0

    if dim_om2 == 0:
        return 0, {}

    bd2 = build_full_boundary_matrix(paths[2], paths[1])
    bd2_om = bd2 @ omega[2]
    U2, S2, Vt2 = np.linalg.svd(bd2_om, full_matrices=True)
    rk_d2 = int(np.sum(np.abs(S2) > 1e-8))
    ker2_dim = dim_om2 - rk_d2

    def get_V_basis(p):
        dim = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim == 0: return np.zeros((dim, 0))
        v_idx = [i for i in range(len(paths[p])) if v in paths[p][i]]
        if not v_idx: return np.eye(dim)
        R = omega[p][v_idx, :]
        U, S, Vt = np.linalg.svd(R, full_matrices=True)
        rk = int(np.sum(np.abs(S) > 1e-8))
        if rk >= dim: return np.zeros((dim, 0))
        return Vt[rk:].T

    V2 = get_V_basis(2)
    dim_V2 = V2.shape[1]

    if ker2_dim == 0:
        return 0, {'ker2_rel': 0, 'im3_rel': 0, 'dim_V2': dim_V2}

    ker2_basis = Vt2[rk_d2:].T

    if dim_V2 > 0:
        combined = np.hstack([ker2_basis, V2])
        rk_comb = np.linalg.matrix_rank(combined, tol=1e-8)
        inter_dim = ker2_dim + dim_V2 - rk_comb
    else:
        inter_dim = 0

    ker2_rel = ker2_dim - inter_dim

    if ker2_rel == 0:
        return 0, {'ker2_rel': 0, 'im3_rel': 0, 'dim_V2': dim_V2}

    if dim_om3 > 0:
        bd3 = build_full_boundary_matrix(paths[3], paths[2])
        bd3_om = bd3 @ omega[3]
        im3_coords, _, _, _ = np.linalg.lstsq(omega[2], bd3_om, rcond=None)

        if dim_V2 > 0:
            combined_im = np.hstack([im3_coords, V2])
            rk_im_comb = np.linalg.matrix_rank(combined_im, tol=1e-8)
            im3_rel = rk_im_comb - dim_V2
        else:
            im3_rel = np.linalg.matrix_rank(im3_coords, tol=1e-8)
    else:
        im3_rel = 0

    h2_rel = max(0, ker2_rel - im3_rel)
    info = {
        'ker2_rel': ker2_rel, 'im3_rel': im3_rel,
        'dim_V2': dim_V2, 'dim_om2': dim_om2, 'dim_om3': dim_om3,
    }
    return h2_rel, info


# ===== n=3 =====
print("=" * 70)
print("RELATIVE HOMOLOGY AT n=3 and n=4")
print("=" * 70)

for n in [3, 4]:
    print(f"\n--- n={n} exhaustive ---")
    n_arcs = n*(n-1)//2
    total = 1 << n_arcs

    violations = 0
    violation_details = []

    for bits in range(total):
        A = build_adj(n, bits)
        betti = path_betti_numbers(A, n, max_dim=3)

        for v_test in range(n):
            h2_rel, info = compute_relative_h2(A, n, v_test)

            if h2_rel > 0:
                violations += 1
                d_v = sum(A[v_test])
                scores = tuple(sorted(sum(row) for row in A))

                # Count 3-cycles
                t3 = 0
                from itertools import combinations
                for triple in combinations(range(n), 3):
                    a, b, c = triple
                    if A[a][b] and A[b][c] and A[c][a]: t3 += 1
                    if A[a][c] and A[c][b] and A[b][a]: t3 += 1

                # Betti of T\v
                B = [[A[i][j] for j in range(n) if j != v_test] for i in range(n) if i != v_test]
                betti_sub = path_betti_numbers(B, n-1, max_dim=2)

                violation_details.append({
                    'bits': bits, 'v': v_test, 'd_v': d_v,
                    'scores': scores, 't3': t3,
                    'betti': betti, 'betti_sub': betti_sub,
                    'h2_rel': h2_rel, 'info': info
                })

    print(f"  Total pairs: {total * n}")
    print(f"  H_2^rel > 0 violations: {violations}")

    if violations > 0:
        print(f"\n  Details of H_2^rel > 0 cases:")
        # Group by (scores, d_v)
        groups = Counter()
        for d in violation_details:
            key = (d['scores'], d['d_v'], d['t3'], tuple(d['betti']),
                   tuple(d['betti_sub']), d['h2_rel'])
            groups[key] += 1

        for key, count in sorted(groups.items()):
            scores, dv, t3, betti, betti_sub, h2r = key
            print(f"    scores={scores}, d_v={dv}, t3={t3}: "
                  f"beta(T)={betti}, beta(T\\v)={betti_sub}, "
                  f"H_2^rel={h2r}, count={count}")

# ===== n=5 quick verify =====
print(f"\n--- n=5 exhaustive ---")
n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs
violations = 0

for bits in range(total):
    A = build_adj(n, bits)
    for v_test in range(n):
        h2_rel, _ = compute_relative_h2(A, n, v_test)
        if h2_rel > 0:
            violations += 1

print(f"  Total pairs: {total * n}")
print(f"  H_2^rel > 0 violations: {violations}")

print("\nDone.")
