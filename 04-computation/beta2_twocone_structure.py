#!/usr/bin/env python3
"""
beta2_twocone_structure.py — Why do 2 cones always suffice?

For each tournament where a single cone fails:
1. Which pairs of vertices work?
2. What's the relationship between the two vertices?
3. Is there a simple criterion for choosing the pair?

Key insight: cone c_v prepends v. It works on paths NOT starting with v's
predecessors (since v→path[0] is needed). So c_v covers the "v-reachable"
part of Z₂, and a second cone covers the rest.

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from itertools import combinations
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

def compute_cone_image(A, n, v, ap2, ap3_list):
    C = np.zeros((len(ap3_list), len(ap2)))
    for j, path in enumerate(ap2):
        if v in path:
            continue
        if not A[v][path[0]]:
            continue
        cone_path = tuple([v] + list(path))
        if cone_path in ap3_list:
            C[ap3_list.index(cone_path), j] = 1.0
    return C


print("=" * 70)
print("TWO-CONE STRUCTURE ANALYSIS")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

# For each 2-cone-needed tournament: which pairs work?
pair_stats = Counter()
pair_score_stats = defaultdict(list)
single_cone_ranks = []
detail_count = 0

for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))

    d2 = dim_om(om2)
    if d2 == 0:
        continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
    U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
    rk = sum(s > 1e-8 for s in S)
    z2_dim = d2 - rk

    if z2_dim == 0:
        continue

    z2_om = Vt[rk:, :]
    z2_A2 = om2 @ z2_om.T
    z2_pinv = np.linalg.pinv(z2_A2)

    ap3_list = [tuple(p) for p in ap3]
    bd3 = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2), 0))

    projected = {}
    cone_ranks = {}
    for v in range(n):
        C_v = compute_cone_image(A, n, v, ap2, ap3_list)
        fill_v = bd3 @ C_v @ z2_A2
        projected[v] = z2_pinv @ fill_v
        cone_ranks[v] = np.linalg.matrix_rank(projected[v], tol=1e-8)

    # Check if single cone works
    single_works = any(cone_ranks[v] >= z2_dim for v in range(n))
    if single_works:
        continue

    # This tournament needs 2 cones
    # Find which pairs work
    working_pairs = []
    for v1, v2 in combinations(range(n), 2):
        cols = np.hstack([projected[v1], projected[v2]])
        if np.linalg.matrix_rank(cols, tol=1e-8) >= z2_dim:
            working_pairs.append((v1, v2))

    # Characterize the working pairs
    for v1, v2 in working_pairs:
        # Edge direction between v1, v2
        if A[v1][v2]:
            edge = "v1→v2"
        else:
            edge = "v2→v1"
        pair_stats[edge] += 1

        # Score relationship
        s1, s2 = scores[v1], scores[v2]
        pair_score_stats[(min(s1,s2), max(s1,s2))].append(bits)

    single_cone_ranks.append((bits, sorted(cone_ranks.values(), reverse=True), z2_dim,
                              len(working_pairs), scores))

    detail_count += 1

print(f"\nn={n}: {detail_count} tournaments need 2 cones")
print(f"\n  Edge direction in working pairs: {dict(pair_stats)}")
print(f"\n  Score pairs in working pairs:")
for sp, lst in sorted(pair_score_stats.items()):
    print(f"    scores ({sp[0]},{sp[1]}): {len(lst)} working pairs")

print(f"\n  Cone rank patterns for 2-cone tournaments (top 10):")
rank_patterns = Counter(tuple(r) for _, r, _, _, _ in single_cone_ranks)
for pattern, cnt in sorted(rank_patterns.items(), key=lambda x: -x[1])[:10]:
    z2 = [z for _, r, z, _, _ in single_cone_ranks if tuple(r) == pattern][0]
    print(f"    ranks={pattern}, z2_dim={z2}: {cnt} tournaments")

# Detailed analysis of a specific 2-cone tournament
print(f"\n{'='*70}")
print("DETAILED ANALYSIS OF ONE 2-CONE TOURNAMENT")
print("=" * 70)

bits, ranks, z2d, nwp, scores = single_cone_ranks[0]
A = build_adj(n, bits)
print(f"\nbits={bits}, scores={scores}, z2_dim={z2d}")
print(f"  Cone ranks per vertex: {dict(zip(range(n), [r for r in ranks]))}")

# Show which 2-paths are in Z₂
ap0 = enumerate_allowed_paths(A, n, 0)
ap1 = enumerate_allowed_paths(A, n, 1)
ap2 = enumerate_allowed_paths(A, n, 2)
ap3 = enumerate_allowed_paths(A, n, 3)
om1 = compute_omega_basis(A, n, 1, ap1, ap0)
om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))

bd2 = build_full_boundary_matrix(ap2, ap1)
d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
rk = sum(s > 1e-8 for s in S)
z2_om = Vt[rk:, :]
z2_A2 = om2 @ z2_om.T

print(f"\n  Z₂ basis vectors (in A₂ coords):")
for j in range(z2_A2.shape[1]):
    z = z2_A2[:, j]
    paths = []
    for i, p in enumerate(ap2):
        if abs(z[i]) > 1e-8:
            tt = "TT" if A[p[0]][p[2]] else "NT"
            paths.append((tuple(p), round(z[i], 3), tt))
    print(f"    z_{j}: {paths}")

# For each vertex: what does the cone cover?
ap3_list = [tuple(p) for p in ap3]
bd3 = build_full_boundary_matrix(ap3, ap2)

for v in range(n):
    print(f"\n  Cone c_{v} (d+={scores[v]}, out-neighbors={[j for j in range(n) if A[v][j]]}):")
    C_v = compute_cone_image(A, n, v, ap2, ap3_list)
    fill_v = bd3 @ C_v @ z2_A2
    z2_pinv = np.linalg.pinv(z2_A2)
    proj_v = z2_pinv @ fill_v
    cr = np.linalg.matrix_rank(proj_v, tol=1e-8)

    # Which Z₂ directions does this cone cover?
    if cr > 0:
        U_p, S_p, Vt_p = np.linalg.svd(proj_v)
        covered_dirs = U_p[:, :cr]
        print(f"    rank={cr}, covers Z₂ directions: {covered_dirs.T.round(3).tolist()}")
    else:
        print(f"    rank=0, covers nothing")

    # Which 2-paths can be coned by v?
    coneable = []
    for i, p in enumerate(ap2):
        if v not in p and A[v][p[0]]:
            cone_path = tuple([v] + list(p))
            if cone_path in ap3_list:
                tt = "TT" if A[p[0]][p[2]] else "NT"
                coneable.append((tuple(p), tt))
    print(f"    Can cone {len(coneable)} paths: {coneable[:5]}{'...' if len(coneable)>5 else ''}")

print("\nDone.")
