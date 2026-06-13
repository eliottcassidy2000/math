#!/usr/bin/env python3
"""beta2_cone_rank_analysis.py - Analyze rank formula for cone B matrix

Looking for closed-form relationship between:
- swap_dim(v): dim of swap cycle space at vertex v
- rank(B_v): rank of unfiltered cone B matrix at vertex v
- The surplus: rank(B_v) - swap_dim(v)

Observed pattern:
  n=5: surplus always = 2 = C(3,2)
  n=6: surplus 2-5, varies
  n=7: surplus 6-12
  n=8: surplus 11-17
  n=9: surplus 15-23

Test: is min surplus = C(n-2, 2)?
  C(3,2)=3 (n=5: observed 2, WRONG)

Let me find the right formula.

Author: kind-pasteur-2026-03-08-S42
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


def analyze_vertex(A, n, v, paths2, path2_idx):
    """Compute swap_dim, rank_B, and related quantities for vertex v."""
    d_out = sum(A[v])
    d_in = n - 1 - d_out
    P = [a for a in range(n) if a != v and A[v][a] == 1]
    Q = [b for b in range(n) if b != v and A[b][v] == 1]
    arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

    if len(arcs_PQ) < 2:
        return None

    m = len(arcs_PQ)
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
        return None
    C_mat = np.array(rows, dtype=float)
    Sc = np.linalg.svd(C_mat, compute_uv=False)
    rank_C = sum(s > 1e-8 for s in Sc)
    ker_dim = m - rank_C
    if ker_dim == 0:
        return None

    _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)
    ker_basis = Vt[rank_C:]

    swap_vecs = []
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
            swap_vecs.append(z)
    if not swap_vecs:
        return None

    swap_dim = np.linalg.matrix_rank(np.column_stack(swap_vecs), tol=1e-8)

    others = [x for x in range(n) if x != v]
    B_sub, vlist = get_induced(A, n, others)
    paths2_Tp = enumerate_allowed_paths(B_sub, len(others), 2)
    num_Tp = len(paths2_Tp)
    Tp_orig = [tuple(vlist[x] for x in p) for p in paths2_Tp]

    B_fill = np.zeros((len(paths2), num_Tp))
    for j, (a, b, c) in enumerate(Tp_orig):
        terms = [
            ((v, b, c), -1), ((v, a, c), +1), ((v, a, b), -1),
            ((b, c, v), +1), ((a, c, v), -1), ((a, b, v), +1),
        ]
        for path, coeff in terms:
            if path in path2_idx:
                B_fill[path2_idx[path], j] += coeff

    rank_B = np.linalg.matrix_rank(B_fill, tol=1e-8)

    return {
        'd_out': d_out, 'd_in': d_in,
        'edges_PQ': len(arcs_PQ), 'ker_dim': ker_dim,
        'swap_dim': swap_dim, 'num_Tp': num_Tp,
        'rank_B': rank_B, 'surplus': rank_B - swap_dim,
    }


# ============================================================
# Collect data for n=5 through n=9
# ============================================================
print("=" * 70)
print("CONE RANK ANALYSIS")
print("=" * 70)

for n in [5, 6, 7, 8, 9]:
    random.seed(42)

    if n <= 6:
        gen = list(all_tournaments_gen(n))
    else:
        gen = [random_tournament(n) for _ in range(200)]

    data = []
    for A in gen:
        paths2 = enumerate_allowed_paths(A, n, 2)
        path2_idx = {tuple(p): i for i, p in enumerate(paths2)}

        for v in range(n):
            result = analyze_vertex(A, n, v, paths2, path2_idx)
            if result:
                data.append(result)

    if not data:
        print(f"\nn={n}: no cases")
        continue

    surpluses = [d['surplus'] for d in data]
    swap_dims = [d['swap_dim'] for d in data]
    rank_Bs = [d['rank_B'] for d in data]
    num_Tps = [d['num_Tp'] for d in data]

    print(f"\nn={n}: {len(data)} cases")
    print(f"  swap_dim: min={min(swap_dims)}, max={max(swap_dims)}")
    print(f"  rank_B:   min={min(rank_Bs)}, max={max(rank_Bs)}")
    print(f"  #Tp:      min={min(num_Tps)}, max={max(num_Tps)}")
    print(f"  surplus:  min={min(surpluses)}, max={max(surpluses)}, "
          f"mean={np.mean(surpluses):.1f}")

    # Check candidate formulas for min surplus
    # n=5: min=2, n=6: min=2, n=7: min=6, n=8: min=11, n=9: min=15
    # Differences: 0, 4, 5, 4
    # C(n-3, 2): C(2,2)=1, C(3,2)=3, C(4,2)=6, C(5,2)=10, C(6,2)=15
    # n*(n-5)/2+1: 0, 3, 7, 12, 18
    # (n-3)*(n-4)/2: 1, 3, 6, 10, 15
    for formula_name, formula in [
        ("C(n-3,2)", (n-3)*(n-4)//2),
        ("C(n-2,2)", (n-2)*(n-3)//2),
        ("(n-3)^2/2", (n-3)**2 // 2),
        ("n-3", n-3),
        ("2*(n-3)", 2*(n-3)),
    ]:
        if formula == min(surpluses):
            print(f"  *** min surplus = {formula_name} = {formula} ***")

    # Check: swap_dim formula
    # At n=5, swap_dim is always 1. ker_dim = 1 always.
    # Bipartite constraint: #P + #Q - 1 active rows
    # ker_dim = #edges_PQ - (#P + #Q - 1) [for connected bipartite]
    # swap_dim <= ker_dim
    print(f"\n  Swap dim vs ker_dim:")
    for s in sorted(set(swap_dims)):
        matching_k = [d['ker_dim'] for d in data if d['swap_dim'] == s]
        print(f"    swap_dim={s}: ker_dim in {sorted(set(matching_k))}")


# ============================================================
# Summary table
# ============================================================
print(f"\n{'=' * 70}")
print("SUMMARY: Min surplus by n")
print("=" * 70)
print("n | min_surplus | C(n-3,2) | C(n-2,2) | 2*(n-3)")
for n_val, min_s in [(5, 2), (6, 2), (7, 6), (8, 11), (9, 15)]:
    cn3 = (n_val-3)*(n_val-4)//2
    cn2 = (n_val-2)*(n_val-3)//2
    n3_2 = 2*(n_val-3)
    print(f"{n_val} | {min_s:>10} | {cn3:>8} | {cn2:>8} | {n3_2:>7}")


print("\n\nDone.")
