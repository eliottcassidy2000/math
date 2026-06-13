#!/usr/bin/env python3
"""
paley11_glmy_estimate.py — opus-2026-03-13-S71

Estimate GLMY Betti numbers for Paley P_11.

Strategy: Use the circulant structure to work with paths starting at 0.
For a circulant tournament on Z_p, vertex-transitivity means:
  |directed m-paths| = p * |directed m-paths starting at 0|
  dim(Ω_m) = p * dim(Ω_m restricted to paths starting at 0)

This reduces the computation by a factor of p.

Even better: for paths starting at 0, we can enumerate efficiently
and compute the junk map restricted to these paths.

Conjecture to test: β_8(P_11) = 10 = p-1
"""

import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import svds
import time

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def enumerate_paths_from_start(A, n, m, start):
    """Directed m-paths starting at 'start'."""
    if m == 0: return [(start,)]
    paths = []
    def dfs(path, depth):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        visited = set(path)
        for v in range(n):
            if v not in visited and A[last][v]:
                path.append(v)
                dfs(path, depth+1)
                path.pop()
    dfs([start], 0)
    return paths

def enumerate_all_directed_paths(A, n, m):
    """All directed m-paths."""
    paths = []
    for start in range(n):
        paths.extend(enumerate_paths_from_start(A, n, m, start))
    return paths

# ============================================================
# Phase 1: Count directed paths for P_11
# ============================================================
print("="*70)
print("PALEY P_11 DIRECTED PATH COUNTS")
print("="*70)

QR11 = {a*a % 11 for a in range(1, 11)}
print(f"  QR mod 11 = {sorted(QR11)}")
A = circulant_tournament(11, QR11)

for m in range(11):
    t0 = time.time()
    paths = enumerate_paths_from_start(A, 11, m, 0)
    t1 = time.time()
    total = len(paths) * 11  # by vertex transitivity
    print(f"  m={m:2d}: from_0={len(paths):8d}, total={total:10d}, time={t1-t0:.2f}s")

# ============================================================
# Phase 2: Compute Ω dims using "paths from 0" reduction
# ============================================================
print(f"\n{'='*70}")
print("PALEY P_11 OMEGA DIMS (VERTEX-TRANSITIVE REDUCTION)")
print("="*70)

# For GLMY: Ω_m = {u ∈ span(A_m) : ∂u ∈ span(A_{m-1})}
# For a circulant tournament, the Z_p action on paths gives:
#   dim(Ω_m) = p * dim(Ω_m ∩ {paths starting at 0})
# This is because the Z_p action commutes with ∂ and the junk map.

# Actually, that's not quite right. The Ω subspace involves LINEAR COMBINATIONS
# of paths. The Z_p action on Ω gives a representation that decomposes into
# eigenspaces. But since we're just computing dim(Ω_m), we need to be more careful.

# CORRECT approach: compute dim(Ω_m) = |A_m| - rank(J_m)
# where J_m is the junk map from A_m to junk faces.
# By vertex transitivity, J_m has a block structure under cyclic shifts.
# But the full J_m matrix might be too large.

# Let me try computing GLMY up to m=5 or so directly

max_m = 6  # conservative limit
allowed = {}
omega_dims = []

for m in range(max_m + 2):
    t0 = time.time()
    allowed[m] = enumerate_all_directed_paths(A, 11, m)
    t1 = time.time()
    print(f"  m={m}: |A_m|={len(allowed[m])}, enum time={t1-t0:.2f}s")

    if m <= 1:
        omega_dims.append(len(allowed[m]))
        print(f"         Ω_m={len(allowed[m])}")
        continue

    # Compute junk map
    am1_set = set(allowed[m-1])
    junk = {}
    junk_count = 0
    for path in allowed[m]:
        for i in range(m+1):
            face = path[:i] + path[i+1:]
            if face not in am1_set and face not in junk:
                junk[face] = junk_count
                junk_count += 1

    if junk_count == 0:
        omega_dims.append(len(allowed[m]))
        print(f"         No junk faces → Ω_m={len(allowed[m])}")
        continue

    print(f"         {junk_count} junk faces, building J matrix...")
    t2 = time.time()

    # Build sparse junk map
    rows_list = []
    cols_list = []
    vals_list = []
    for j, path in enumerate(allowed[m]):
        for i in range(m+1):
            face = path[:i] + path[i+1:]
            if face in junk:
                rows_list.append(junk[face])
                cols_list.append(j)
                vals_list.append((-1)**i)

    J = csr_matrix((vals_list, (rows_list, cols_list)),
                    shape=(junk_count, len(allowed[m])))

    # Rank of J via dense SVD (if small enough) or sparse
    if junk_count * len(allowed[m]) < 5_000_000:
        J_dense = J.toarray()
        U, s, Vh = np.linalg.svd(J_dense, full_matrices=False)
        rank_J = int(np.sum(s > 1e-10))
    else:
        # Use J^T J approach (smaller matrix if one dimension is small)
        JTJ = (J.T @ J).toarray()
        eigvals = np.linalg.eigvalsh(JTJ)
        rank_J = int(np.sum(eigvals > 1e-10))

    omega_m = len(allowed[m]) - rank_J
    omega_dims.append(omega_m)
    t3 = time.time()
    print(f"         rank(J)={rank_J}, Ω_m={omega_m}, time={t3-t2:.2f}s")

print(f"\n  Summary: Ω = {omega_dims}")
if all(o % 11 == 0 for o in omega_dims):
    print(f"  Ω/11 = {[o//11 for o in omega_dims]}")

# ============================================================
# Phase 3: Compare with known P_7 values
# ============================================================
print(f"\n{'='*70}")
print("COMPARISON: P_7 Ω/7 vs P_11 Ω/11")
print("="*70)

p7_omega_n = [1, 3, 6, 9, 9, 6, 3]
print(f"  P_7:  Ω/7  = {p7_omega_n}")
if len(omega_dims) > 0:
    p11_omega_n = [o//11 if o % 11 == 0 else f"{o}/11" for o in omega_dims]
    print(f"  P_11: Ω/11 = {p11_omega_n} (computed up to m={max_m+1})")

print("\nDONE.")
