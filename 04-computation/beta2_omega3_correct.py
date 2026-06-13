#!/usr/bin/env python3
"""
beta2_omega3_correct.py - Correct Omega_3 dimension formula

The FULL constraint matrix for Omega_3:
- Rows: ALL invalid 2-sequences appearing as d_1 or d_2 faces
- Columns: ALL non-DD 3-paths (02, 13, NN types)
- An invalid face (a,c,d) [d_1 type] constrains: sum of alpha for (a,b,c,d) = 0
- An invalid face (a,b,d) [d_2 type] constrains: sum of alpha for (a,b,c,d) = 0

dim(Omega_3) = |DD| + |non-DD| - rk(full constraint matrix)

The cross-links from NN paths (which appear in BOTH d_1 and d_2 constraints)
can reduce the rank of the constraint matrix, giving MORE Omega_3 dimensions.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
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


def compute_omega3_formula(A, n):
    """Compute dim(Omega_3) using the full constraint matrix."""
    a3 = enumerate_allowed_paths(A, n, 3)
    if not a3:
        return 0, 0, 0

    # Classify paths
    dd_idx = []
    nondd_idx = []

    for i, p in enumerate(a3):
        if A[p[0]][p[2]] and A[p[1]][p[3]]:
            dd_idx.append(i)
        else:
            nondd_idx.append(i)

    if not nondd_idx:
        return len(dd_idx), 0, len(dd_idx)

    # Build constraint matrix
    # Find all invalid faces
    d1_faces = {}  # (a,c,d) -> row index
    d2_faces = {}  # (a,b,d) -> row index

    for i in nondd_idx:
        p = a3[i]
        a_, b_, c_, d_ = p
        if not A[a_][c_]:  # d_1 face (a,c,d) is invalid
            face = (a_, c_, d_)
            if face not in d1_faces:
                d1_faces[face] = len(d1_faces)
        if not A[b_][d_]:  # d_2 face (a,b,d) is invalid
            face = (a_, b_, d_)
            if face not in d2_faces:
                d2_faces[face] = len(d2_faces)

    n_rows = len(d1_faces) + len(d2_faces)
    n_cols = len(nondd_idx)

    C = np.zeros((n_rows, n_cols))

    for j, idx in enumerate(nondd_idx):
        p = a3[idx]
        a_, b_, c_, d_ = p
        if not A[a_][c_]:
            row = d1_faces[(a_, c_, d_)]
            C[row, j] = 1  # d_1 face has coefficient -1 in boundary, constraint is sum = 0
        if not A[b_][d_]:
            row = len(d1_faces) + d2_faces[(a_, b_, d_)]
            C[row, j] = 1

    rk_C = np.linalg.matrix_rank(C, tol=1e-8)
    predicted = len(dd_idx) + n_cols - rk_C

    return len(dd_idx), rk_C, predicted


# ============================================================
# TEST at n=5
# ============================================================
print("=" * 70)
print("CORRECT OMEGA_3 FORMULA TEST")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

mismatches = 0
for bits in range(total):
    A = build_adj(n, bits)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    dd, rk_C, predicted = compute_omega3_formula(A, n)

    if predicted != d_om3:
        mismatches += 1
        if mismatches <= 5:
            scores = tuple(sorted([sum(row) for row in A]))
            print(f"  MISMATCH: bits={bits}, scores={scores}, "
                  f"DD={dd}, rk_C={rk_C}, pred={predicted}, actual={d_om3}")

print(f"\nn=5: Mismatches: {mismatches}/{total}")
if mismatches == 0:
    print("  ** CORRECT FORMULA VERIFIED at n=5! **")
    print("  dim(Omega_3) = |DD| + |non-DD| - rk(full constraint matrix)")


# ============================================================
# TEST at n=6
# ============================================================
print(f"\n{'='*70}")
print("FORMULA TEST at n=6")
print("=" * 70)

n = 6
total = 1 << (n*(n-1)//2)

mismatches_6 = 0
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

    dd, rk_C, predicted = compute_omega3_formula(A, n)

    if predicted != d_om3:
        mismatches_6 += 1
        if mismatches_6 <= 3:
            scores = tuple(sorted([sum(row) for row in A]))
            print(f"  MISMATCH: bits={bits}, scores={scores}, "
                  f"DD={dd}, rk_C={rk_C}, pred={predicted}, actual={d_om3}")

dt = time.time() - t0
print(f"\nn=6: Done in {dt:.0f}s")
print(f"Mismatches: {mismatches_6}/{total}")
if mismatches_6 == 0:
    print("  ** CORRECT FORMULA VERIFIED at n=6! **")
    print("  dim(Omega_3) = |DD| + |non-DD| - rk(constraint matrix)")
    print("  = |A_3| - rk(constraint matrix)")
    print("  where constraint matrix has one row per invalid d_1/d_2 face")
    print("  and one column per non-DD 3-path")


# ============================================================
# ANALYSIS: constraint matrix structure
# ============================================================
print(f"\n{'='*70}")
print("CONSTRAINT MATRIX STRUCTURE")
print("=" * 70)

# The constraint matrix has a bipartite structure:
# d_1 constraints (top rows): each involves non-DD paths with fixed (a,c,d)
# d_2 constraints (bottom rows): each involves non-DD paths with fixed (a,b,d)
# NN paths appear in BOTH sections, creating cross-links.
#
# If there were NO NN paths, the matrix would be block-diagonal:
# rk = rk(d1 block) + rk(d2 block) = |d1 faces| + |d2 faces|
# (since each face has at least one path, and the constraints are just
# "sum of coefficients = 0" which is rank 1 per face)
#
# With NN paths, the cross-links can create rank deficiency.
# Key question: by how much does the cross-link reduce rank?

print("\nRank reduction from cross-links at n=5:")
n = 5
total = 1 << (n*(n-1)//2)

reduction_dist = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    a3 = enumerate_allowed_paths(A, n, 3)

    # Compute constraint matrix components
    nondd_idx = [i for i, p in enumerate(a3) if not (A[p[0]][p[2]] and A[p[1]][p[3]])]
    if not nondd_idx:
        continue

    d1_faces = {}
    d2_faces = {}
    for i in nondd_idx:
        p = a3[i]
        if not A[p[0]][p[2]]:
            d1_faces.setdefault((p[0], p[2], p[3]), len(d1_faces))
        if not A[p[1]][p[3]]:
            d2_faces.setdefault((p[0], p[1], p[3]), len(d2_faces))

    n_d1 = len(d1_faces)
    n_d2 = len(d2_faces)

    # Full matrix
    C = np.zeros((n_d1 + n_d2, len(nondd_idx)))
    for j, idx in enumerate(nondd_idx):
        p = a3[idx]
        if not A[p[0]][p[2]]:
            C[d1_faces[(p[0], p[2], p[3])], j] = 1
        if not A[p[1]][p[3]]:
            C[n_d1 + d2_faces[(p[0], p[1], p[3])], j] = 1

    rk_full = np.linalg.matrix_rank(C, tol=1e-8)
    reduction = (n_d1 + n_d2) - rk_full
    reduction_dist[reduction] += 1

print(f"Rank reduction distribution: {dict(sorted(reduction_dist.items()))}")
print(f"  reduction=0: block-diagonal (no cross-link effect)")
print(f"  reduction>0: cross-links reduce constraints")


print("\nDone.")
