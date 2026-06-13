#!/usr/bin/env python3
"""
beta2_h2rel_combinatorial.py — Combinatorial characterization of h₂_rel

Goal: find a purely combinatorial condition for h₂_rel(T,T\\v) = 1.

Strategy: For each (T, v) pair with h₂_rel = 1, extract the
structural features that distinguish it from h₂_rel = 0 cases.

Key candidate: the "3-cycle neighborhood" of v.
h₂_rel > 0 requires β₁(T\\v) > 0, which requires 3-cycles in T\\v.
Which 3-cycles? And how does v relate to them?

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def compute_h2_rel(A, n, v):
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
    remap = {i: others[i] for i in range(n1)}

    ap2_T = enumerate_allowed_paths(A, n, 2)
    if not ap2_T: return 0
    ap1_T = enumerate_allowed_paths(A, n, 1)
    ap3_T = enumerate_allowed_paths(A, n, 3)
    ap0_T = enumerate_allowed_paths(A, n, 0)

    om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T)
    om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
    om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))

    d2_T = dim_om(om2_T)
    if d2_T == 0: return 0

    ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
    ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
    ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
    om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))
    om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
    d2_sub = dim_om(om2_sub); d1_T = dim_om(om1_T); d1_sub = dim_om(om1_sub)
    d3_T = dim_om(om3_T)

    ap2_T_list = [tuple(p) for p in ap2_T]
    ap1_T_list = [tuple(p) for p in ap1_T]

    if ap2_sub and d2_sub > 0:
        embed = np.zeros((len(ap2_T_list), d2_sub))
        for j in range(d2_sub):
            for k, path_sub in enumerate(ap2_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap2_T_list:
                    embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
        phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]
    else:
        phi = np.zeros((d2_T, 0))

    rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
    if rk_phi > 0:
        U_phi, _, _ = np.linalg.svd(phi, full_matrices=True)
        Q = U_phi[:, rk_phi:]
    else:
        Q = np.eye(d2_T)

    d_rel = Q.shape[1]
    if d_rel == 0: return 0

    coords2_T = np.linalg.lstsq(om1_T, build_full_boundary_matrix(ap2_T, ap1_T) @ om2_T, rcond=None)[0]

    if d1_sub > 0:
        embed1 = np.zeros((len(ap1_T_list), d1_sub))
        for j in range(d1_sub):
            for k, path_sub in enumerate(ap1_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap1_T_list:
                    embed1[ap1_T_list.index(path_T), j] = om1_sub[k, j]
        psi = np.linalg.lstsq(om1_T, embed1, rcond=None)[0]
    else:
        psi = np.zeros((d1_T, 0))

    rk_psi = np.linalg.matrix_rank(psi, tol=1e-8)
    R = np.linalg.svd(psi, full_matrices=True)[0][:, rk_psi:] if rk_psi > 0 else np.eye(d1_T)

    coords2_rel = R.T @ coords2_T @ Q
    rk_d2_rel = np.linalg.matrix_rank(coords2_rel, tol=1e-8)
    z2_rel = d_rel - rk_d2_rel
    if z2_rel == 0: return 0

    if d3_T > 0:
        ap3_sub = enumerate_allowed_paths(A_sub, n1, 3)
        coords3_T = np.linalg.lstsq(om2_T, build_full_boundary_matrix(ap3_T, ap2_T) @ om3_T, rcond=None)[0]
        om3_sub = compute_omega_basis(A_sub, n1, 3, ap3_sub, ap2_sub) if ap3_sub else np.zeros((0,0))
        d3_sub = dim_om(om3_sub)
        d3_proj = Q.T @ coords3_T
        if d3_sub > 0:
            ap3_T_list = [tuple(p) for p in ap3_T]
            embed3 = np.zeros((len(ap3_T_list), d3_sub))
            for j in range(d3_sub):
                for k, path_sub in enumerate(ap3_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap3_T_list:
                        embed3[ap3_T_list.index(path_T), j] = om3_sub[k, j]
            chi = np.linalg.lstsq(om3_T, embed3, rcond=None)[0]
            rk_chi = np.linalg.matrix_rank(chi, tol=1e-8)
            if rk_chi > 0:
                U_chi, _, _ = np.linalg.svd(chi, full_matrices=True)
                d3_proj_rel = d3_proj @ U_chi[:, rk_chi:]
            else:
                d3_proj_rel = d3_proj
        else:
            d3_proj_rel = d3_proj
    else:
        d3_proj_rel = np.zeros((d_rel, 0))

    rk_d3_rel = np.linalg.matrix_rank(d3_proj_rel, tol=1e-8)
    return z2_rel - rk_d3_rel


print("=" * 70)
print("COMBINATORIAL CHARACTERIZATION OF h₂_rel > 0")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

# For each interior vertex v with h₂_rel = 1:
# - What is the relationship of v to 3-cycles in T\v?
# - Specifically: does v "see" the 3-cycle? (edges from v to all 3 cycle vertices)

print(f"\n--- Relationship of v to 3-cycles in T\\v ---")

# For each 3-cycle {a,b,c} in T\v (say a→b→c→a):
# v has 3 edges to {a,b,c}. There are 8 possible patterns.
# Which patterns occur when h₂_rel = 1 vs 0?

cycle_relation_h2_1 = Counter()
cycle_relation_h2_0 = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    for v in range(n):
        if scores[v] == 0 or scores[v] == n-1:
            continue

        h2r = compute_h2_rel(A, n, v)

        # Find 3-cycles in T\v
        others = [i for i in range(n) if i != v]
        cycles = []
        for a_idx in range(len(others)):
            for b_idx in range(a_idx+1, len(others)):
                for c_idx in range(b_idx+1, len(others)):
                    a, b, c = others[a_idx], others[b_idx], others[c_idx]
                    if A[a][b] and A[b][c] and A[c][a]:
                        cycles.append((a, b, c))
                    elif A[a][c] and A[c][b] and A[b][a]:
                        cycles.append((a, c, b))

        if not cycles:
            continue  # No 3-cycles in T\v

        # For each cycle, classify v's relationship
        for cyc in cycles:
            a, b, c = cyc  # a→b→c→a
            # v's edges to a,b,c: (v→a or a→v), etc.
            pattern = (A[v][a], A[v][b], A[v][c])  # 1 if v→x, 0 if x→v
            outdeg_to_cycle = sum(pattern)

            if h2r > 0:
                cycle_relation_h2_1[(outdeg_to_cycle, pattern)] += 1
            else:
                cycle_relation_h2_0[(outdeg_to_cycle, pattern)] += 1

print(f"h₂_rel = 1 cases — v's edges to 3-cycle vertices:")
for key in sorted(cycle_relation_h2_1):
    print(f"  d+(v→cycle)={key[0]}, pattern={key[1]}: {cycle_relation_h2_1[key]}")

print(f"\nh₂_rel = 0 cases — v's edges to 3-cycle vertices:")
for key in sorted(cycle_relation_h2_0):
    print(f"  d+(v→cycle)={key[0]}, pattern={key[1]}: {cycle_relation_h2_0[key]}")


# Part 2: When h₂_rel = 1, how does v relate to ALL 3-cycles in T\v?
# (not just individual ones)
print(f"\n--- Part 2: v's total relationship to 3-cycles ---")

# For each vertex: count (total_cycles_in_Tv, cycles_where_v_beats_2, cycles_where_v_beats_1, etc.)
relationship_h2_1 = Counter()
relationship_h2_0 = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    for v in range(n):
        if scores[v] == 0 or scores[v] == n-1:
            continue

        h2r = compute_h2_rel(A, n, v)

        others = [i for i in range(n) if i != v]
        n_cycles = 0
        beat_pattern = Counter()
        for a_idx in range(len(others)):
            for b_idx in range(a_idx+1, len(others)):
                for c_idx in range(b_idx+1, len(others)):
                    a, b, c = others[a_idx], others[b_idx], others[c_idx]
                    is_cycle = False
                    if A[a][b] and A[b][c] and A[c][a]:
                        is_cycle = True
                    elif A[a][c] and A[c][b] and A[b][a]:
                        is_cycle = True
                    if is_cycle:
                        n_cycles += 1
                        outdeg = A[v][a] + A[v][b] + A[v][c]
                        beat_pattern[outdeg] += 1

        key = (n_cycles, tuple(sorted(beat_pattern.items())))
        if h2r > 0:
            relationship_h2_1[key] += 1
        else:
            relationship_h2_0[key] += 1

print(f"h₂_rel = 1: (n_cycles, beat_pattern)")
for key in sorted(relationship_h2_1, key=lambda k: -relationship_h2_1[k]):
    print(f"  {key}: {relationship_h2_1[key]}")

print(f"\nh₂_rel = 0 (only with cycles): (n_cycles, beat_pattern)")
for key in sorted(relationship_h2_0, key=lambda k: -relationship_h2_0[k])[:10]:
    print(f"  {key}: {relationship_h2_0[key]}")


# Part 3: What is special about the Σ=3 tournaments?
print(f"\n--- Part 3: Anatomy of Σ h₂_rel = 3 tournaments ---")

for bits in range(total):
    A = build_adj(n, bits)
    h2_sum = sum(compute_h2_rel(A, n, v) for v in range(n))

    if h2_sum == 3:
        scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
        print(f"\n  bits={bits}, scores={sorted(scores)}")
        adj_str = [''.join(str(A[i][j]) for j in range(n)) for i in range(n)]
        print(f"  Adj: {adj_str}")

        for v in range(n):
            h2r = compute_h2_rel(A, n, v)
            # Count 3-cycles in T\v
            others = [i for i in range(n) if i != v]
            cyc_count = 0
            for a in range(len(others)):
                for b in range(a+1, len(others)):
                    for c in range(b+1, len(others)):
                        aa, bb, cc = others[a], others[b], others[c]
                        if (A[aa][bb] and A[bb][cc] and A[cc][aa]) or (A[aa][cc] and A[cc][bb] and A[bb][aa]):
                            cyc_count += 1
            print(f"    v={v}: d+={scores[v]}, h₂_rel={h2r}, 3-cycles_in_T\\v={cyc_count}")

        # Count total 3-cycles
        t3 = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                        t3 += 1
        print(f"  t₃ = {t3}")
        break  # just show first one

print("\nDone.")
