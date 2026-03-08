#!/usr/bin/env python3
"""
beta_2 = 0 for tournaments: algebraic proof exploration.
"""
import numpy as np
from itertools import combinations
import sys, time
from collections import Counter, defaultdict
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

# ===== Cancellation space formula =====
print("=" * 70)
print("CANCELLATION SPACE FORMULA")
print("=" * 70)

for n in [5]:
    match = 0
    total = 0
    for A in all_tournaments_gen(n):
        total += 1
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a2_list = [tuple(p) for p in a2]
        n_tt = sum(1 for p in a2_list if A[p[0]][p[2]] == 1)
        om2 = compute_omega_basis(A, n, 2, a2, a1)
        dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
        face_mult = defaultdict(int)
        for p in a2_list:
            a, b, c = p
            if A[c][a] == 1:
                face_mult[(a, c)] += 1
        predicted_gap = sum(max(0, m - 1) for m in face_mult.values())
        actual_gap = dim_om2 - n_tt
        if predicted_gap == actual_gap:
            match += 1
    print(f"  n={n}: gap = sum(mult-1) matches: {match}/{total}")

# ===== Dimension table by t3 at n=5 =====
print(f"\n{'='*70}")
print("DIMENSION TABLE n=5 (by t3)")
print("="*70)

by_t3 = defaultdict(list)
for A in all_tournaments_gen(5):
    a1 = enumerate_allowed_paths(A, 5, 1)
    a2 = enumerate_allowed_paths(A, 5, 2)
    a3 = enumerate_allowed_paths(A, 5, 3)
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]
    a2_set = set(a2_list)

    n_tt = sum(1 for p in a2_list if A[p[0]][p[2]] == 1)
    n_dt = sum(1 for p in a3_list if
               all(f in a2_set for f in [(p[1],p[2],p[3]),
                   (p[0],p[2],p[3]), (p[0],p[1],p[3]), (p[0],p[1],p[2])]))

    om2 = compute_omega_basis(A, 5, 2, a2, a1)
    om3 = compute_omega_basis(A, 5, 3, a3, a2)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2 if dim_om2 > 0 else np.zeros((len(a1), 0))
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8) if dim_om2 > 0 else 0
    ker2 = dim_om2 - rank2

    t3 = sum(1 for a, b, c in combinations(range(5), 3)
             if (A[a][b] and A[b][c] and A[c][a]) or
                (A[b][a] and A[a][c] and A[c][b]))
    by_t3[t3].append((n_tt, n_dt, dim_om2, dim_om3, rank2, ker2))

print(f"  t3 -> (|TT|, |DT|, dim_Om2, dim_Om3, rk_bd2, ker_bd2)")
for t3 in sorted(by_t3.keys()):
    vals = by_t3[t3]
    unique = Counter(tuple(v) for v in vals)
    for key in sorted(unique.keys()):
        print(f"    t3={t3}: TT={key[0]}, DT={key[1]}, Om2={key[2]}, Om3={key[3]}, "
              f"rk={key[4]}, ker={key[5]} [{unique[key]}x]")

# ===== ker(bd_2) vs |DT| =====
print(f"\n{'='*70}")
print("ker(bd_2|Om_2) vs |DT| relationship")
print("="*70)

for n in [4, 5]:
    print(f"\n--- n={n} ---")
    ker_dt = Counter()
    for A in all_tournaments_gen(n):
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)
        a2_list = [tuple(p) for p in a2]
        a3_list = [tuple(p) for p in a3]
        a2_set = set(a2_list)

        n_dt = sum(1 for p in a3_list if
                   all(f in a2_set for f in [(p[1],p[2],p[3]),
                       (p[0],p[2],p[3]), (p[0],p[1],p[3]), (p[0],p[1],p[2])]))

        om2 = compute_omega_basis(A, n, 2, a2, a1)
        dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
        if dim_om2 == 0:
            ker_dt[(0, n_dt)] += 1
            continue

        bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
        bd2_om = bd2 @ om2
        rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
        ker2 = dim_om2 - rank2
        ker_dt[(ker2, n_dt)] += 1

    print(f"  (ker, |DT|): count")
    for key in sorted(ker_dt.keys()):
        ker2, ndt = key
        # Compute ratio
        ratio = ker2 / ndt if ndt > 0 else 0
        print(f"    ker={ker2}, |DT|={ndt}, ratio={ratio:.3f}: {ker_dt[key]}")
    always_le = all(k <= d for (k, d) in ker_dt.keys())
    print(f"  ker <= |DT| always? {always_le}")

# ===== Per-edge kernel contributions =====
print(f"\n{'='*70}")
print("EDGE-BASED ANALYSIS: which edges are 'saturated'?")
print("="*70)

# For each 2-cycle, compute its "support" on each edge
# Question: does the cycle vanish on all edges incident to some vertex?

n = 5
vertex_zero_dist = Counter()
for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a2_list = [tuple(p) for p in a2]
    a1_list = [tuple(p) for p in a1]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2_list, a1_list)
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        continue

    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    ker_basis = (om2 @ Vt[rank2:, :].T).T

    for i in range(ker_dim):
        cyc = ker_basis[i]
        # For each vertex v, check if cycle has zero support on all paths through v
        for v in range(n):
            zero_v = all(abs(cyc[j]) < 1e-8 for j, p in enumerate(a2_list)
                        if v in p)
            if zero_v:
                vertex_zero_dist['has_zero_vertex'] += 1
                break
        else:
            vertex_zero_dist['all_vertices_active'] += 1

print(f"  n=5: {vertex_zero_dist}")

print("\nDone.")
