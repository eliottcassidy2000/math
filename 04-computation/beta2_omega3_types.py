#!/usr/bin/env python3
"""
beta2_omega3_types.py - Complete Omega_3 type decomposition

3-path types by bad face structure:
  DD: a->c AND b->d (0 bad faces) -> individually in Omega_3
  02: c->a AND b->d (1 bad face at d_1)
  13: a->c AND d->b (1 bad face at d_2)
  NN: c->a AND d->b (2 bad faces at d_1 and d_2)

dim(Omega_3) = DD + cancel_02 + cancel_13 + NN_contribution

The NN contribution is the number of NN linear combinations where
BOTH bad faces simultaneously cancel. This requires solving a
system of linear constraints.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
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


print("=" * 70)
print("OMEGA_3 TYPE DECOMPOSITION")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

# For each tournament, decompose dim(Omega_3) into components
gap_dist = Counter()
nn_contrib_data = []

for bits in range(total):
    A = build_adj(n, bits)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    dd_count = sum(1 for p in a3 if A[p[0]][p[2]] and A[p[1]][p[3]])

    # 02 cancellation
    t02_idx = [i for i, p in enumerate(a3) if not A[p[0]][p[2]] and A[p[1]][p[3]]]
    bad02 = defaultdict(list)
    for i in t02_idx:
        p = a3[i]
        bad02[(p[0], p[2], p[3])].append(i)
    cancel_02 = sum(max(0, len(g)-1) for g in bad02.values())

    # 13 cancellation
    t13_idx = [i for i, p in enumerate(a3) if A[p[0]][p[2]] and not A[p[1]][p[3]]]
    bad13 = defaultdict(list)
    for i in t13_idx:
        p = a3[i]
        bad13[(p[0], p[1], p[3])].append(i)
    cancel_13 = sum(max(0, len(g)-1) for g in bad13.values())

    # NN paths
    nn_idx = [i for i, p in enumerate(a3) if not A[p[0]][p[2]] and not A[p[1]][p[3]]]
    n_nn = len(nn_idx)

    # NN d1 and d2 face groups
    nn_d1_groups = defaultdict(list)
    nn_d2_groups = defaultdict(list)
    for i in nn_idx:
        p = a3[i]
        nn_d1_groups[(p[0], p[2], p[3])].append(i)
        nn_d2_groups[(p[0], p[1], p[3])].append(i)

    # NN contribution = dim(Omega_3) - DD - cancel_02 - cancel_13
    nn_contrib = d_om3 - dd_count - cancel_02 - cancel_13
    gap_dist[nn_contrib] += 1

    # Compute NN contribution directly via constraint counting
    # Each NN path has 2 constraints: d1 face and d2 face must cancel
    # Constraint matrix: rows = distinct bad faces, columns = NN paths
    # For d1 face f: constraint = sum of coefficients of NN paths with d1=f is 0
    # For d2 face f: constraint = sum of coefficients of NN paths with d2=f is 0
    if n_nn > 0:
        d1_faces = sorted(set(nn_d1_groups.keys()))
        d2_faces = sorted(set(nn_d2_groups.keys()))
        n_constraints = len(d1_faces) + len(d2_faces)

        C = np.zeros((n_constraints, n_nn))
        for row, face in enumerate(d1_faces):
            for j, idx in enumerate(nn_idx):
                p = a3[idx]
                if (p[0], p[2], p[3]) == face:
                    C[row, j] = 1
        for row_offset, face in enumerate(d2_faces):
            row = len(d1_faces) + row_offset
            for j, idx in enumerate(nn_idx):
                p = a3[idx]
                if (p[0], p[1], p[3]) == face:
                    C[row, j] = 1

        rk_C = np.linalg.matrix_rank(C, tol=1e-8)
        nn_predicted = n_nn - rk_C  # freedom = |NN| - rank(constraints)

        nn_contrib_data.append({
            'bits': bits, 'n_nn': n_nn, 'd1': len(d1_faces),
            'd2': len(d2_faces), 'rk_C': rk_C, 'predicted': nn_predicted,
            'actual': nn_contrib
        })

print(f"NN contribution distribution: {dict(sorted(gap_dist.items()))}")

# Check if predicted = actual
mismatches = 0
for d in nn_contrib_data:
    if d['predicted'] != d['actual']:
        mismatches += 1
        if mismatches <= 5:
            print(f"  MISMATCH: bits={d['bits']}, NN={d['n_nn']}, d1={d['d1']}, "
                  f"d2={d['d2']}, rk_C={d['rk_C']}, pred={d['predicted']}, actual={d['actual']}")

print(f"\nMismatches (predicted vs actual NN contrib): {mismatches}/{len(nn_contrib_data)}")
if mismatches == 0:
    print("  ** FORMULA CONFIRMED: dim(Omega_3) = DD + cancel_02 + cancel_13 + (|NN| - rk(C)) **")
    print("  where C is the constraint matrix for NN bad face cancellation")


# ============================================================
# ANALYSIS 2: Check formula at n=6
# ============================================================
print(f"\n{'='*70}")
print("FORMULA CHECK at n=6")
print("=" * 70)

n = 6
total = 1 << (n*(n-1)//2)

mismatches_6 = 0
import time
t0 = time.time()
for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s), mismatches={mismatches_6}")

    A = build_adj(n, bits)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    dd_count = sum(1 for p in a3 if A[p[0]][p[2]] and A[p[1]][p[3]])

    t02_idx = [i for i, p in enumerate(a3) if not A[p[0]][p[2]] and A[p[1]][p[3]]]
    bad02 = defaultdict(list)
    for i in t02_idx:
        p = a3[i]
        bad02[(p[0], p[2], p[3])].append(i)
    cancel_02 = sum(max(0, len(g)-1) for g in bad02.values())

    t13_idx = [i for i, p in enumerate(a3) if A[p[0]][p[2]] and not A[p[1]][p[3]]]
    bad13 = defaultdict(list)
    for i in t13_idx:
        p = a3[i]
        bad13[(p[0], p[1], p[3])].append(i)
    cancel_13 = sum(max(0, len(g)-1) for g in bad13.values())

    nn_idx = [i for i, p in enumerate(a3) if not A[p[0]][p[2]] and not A[p[1]][p[3]]]
    n_nn = len(nn_idx)

    if n_nn > 0:
        nn_d1_groups = defaultdict(list)
        nn_d2_groups = defaultdict(list)
        for i in nn_idx:
            p = a3[i]
            nn_d1_groups[(p[0], p[2], p[3])].append(i)
            nn_d2_groups[(p[0], p[1], p[3])].append(i)

        d1_faces = sorted(set(nn_d1_groups.keys()))
        d2_faces = sorted(set(nn_d2_groups.keys()))
        n_constraints = len(d1_faces) + len(d2_faces)

        C = np.zeros((n_constraints, n_nn))
        for row, face in enumerate(d1_faces):
            for j, idx in enumerate(nn_idx):
                p = a3[idx]
                if (p[0], p[2], p[3]) == face:
                    C[row, j] = 1
        for row_offset, face in enumerate(d2_faces):
            row = len(d1_faces) + row_offset
            for j, idx in enumerate(nn_idx):
                p = a3[idx]
                if (p[0], p[1], p[3]) == face:
                    C[row, j] = 1

        rk_C = np.linalg.matrix_rank(C, tol=1e-8)
        nn_contrib = n_nn - rk_C
    else:
        nn_contrib = 0

    predicted = dd_count + cancel_02 + cancel_13 + nn_contrib
    if predicted != d_om3:
        mismatches_6 += 1
        if mismatches_6 <= 3:
            scores = tuple(sorted([sum(row) for row in A]))
            print(f"  MISMATCH at bits={bits}, scores={scores}: "
                  f"DD={dd_count}, c02={cancel_02}, c13={cancel_13}, NN_c={nn_contrib}, "
                  f"pred={predicted}, actual={d_om3}")

dt = time.time() - t0
print(f"\nDone in {dt:.0f}s")
print(f"Mismatches at n=6: {mismatches_6}/{total}")
if mismatches_6 == 0:
    print("  ** FORMULA CONFIRMED at n=6: dim(Omega_3) = DD + cancel_02 + cancel_13 + (|NN| - rk(C)) **")


print("\nDone.")
