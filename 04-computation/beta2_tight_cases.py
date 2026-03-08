#!/usr/bin/env python3
"""Analyze the tightest β_2=0 cases: tournaments where dim(Ω_3) - dim(Z_2) = 1.

These are the "barely exact" tournaments where one extra Ω_3 element is the
only thing preventing β_2 > 0. Understanding their structure should reveal
WHY β_2=0 is forced.

At n=6: tightest cases have score (1,1,1,4,4,4), t3=2,
        Ω_3=9, Z_2=8, B_2=8, ker(∂_3)=1.
"""
import numpy as np
from itertools import combinations
import sys
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
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

# First pass: find tight cases
print("=" * 70)
print("FINDING TIGHT CASES AT n=6")
print("=" * 70)

n = 6
tight_tournaments = []

for tidx, A in enumerate(all_tournaments_gen(n)):
    if tidx % 10000 == 0:
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

    dim_om2 = om[2].shape[1] if om[2].ndim == 2 else 0
    dim_om3 = om[3].shape[1] if om[3].ndim == 2 else 0

    if dim_om2 == 0 or dim_om3 == 0:
        continue

    bd2 = build_full_boundary_matrix(a[2], a[1])
    bd2_om = bd2 @ om[2]
    S2 = np.linalg.svd(bd2_om, compute_uv=False)
    rank_d2 = sum(s > 1e-8 for s in S2)
    z2 = dim_om2 - rank_d2

    surplus = dim_om3 - z2
    if surplus <= 2:
        tight_tournaments.append((tidx, [row[:] for row in A], surplus, dim_om2, dim_om3, z2))

print(f"\nFound {len(tight_tournaments)} tight tournaments (surplus ≤ 2)")

# Analyze tight cases in detail
print(f"\n{'='*70}")
print(f"DETAILED ANALYSIS OF TIGHT CASES")
print(f"{'='*70}")

for tidx, A, surplus, dom2, dom3, z2 in tight_tournaments[:5]:
    print(f"\n--- Tournament #{tidx}, surplus={surplus} ---")
    scores = [sum(A[i]) for i in range(n)]
    print(f"  Scores: {scores}")
    print(f"  Adjacency matrix:")
    for i in range(n):
        print(f"    {A[i]}")

    # Full chain complex computation
    a = {}
    for p in range(5):
        a[p] = [tuple(x) for x in enumerate_allowed_paths(A, n, p)]

    om = {}
    for p in range(5):
        if p == 0:
            om[p] = np.eye(n)
        else:
            om[p] = compute_omega_basis(A, n, p, a[p], a[p-1])

    dim_om2 = om[2].shape[1] if om[2].ndim == 2 else 0
    dim_om3 = om[3].shape[1] if om[3].ndim == 2 else 0

    print(f"  |A_2|={len(a[2])}, |A_3|={len(a[3])}")
    print(f"  dim(Ω_2)={dim_om2}, dim(Ω_3)={dim_om3}")
    print(f"  Z_2={z2}, surplus={surplus}")

    # Identify bad pairs (a,c) where c→a
    bad_pairs = [(i,j) for i in range(n) for j in range(n) if i!=j and A[j][i]]
    print(f"  Bad pairs (c→a): {len(bad_pairs)}")

    # For each bad pair, how many 2-paths (a,b,c) have it?
    for (ai, ci) in bad_pairs[:10]:
        mids = [b for b in range(n) if b != ai and b != ci and A[ai][b] and A[b][ci]]
        if len(mids) > 0:
            print(f"    ({ai},{ci}): {len(mids)} paths through {mids}")

    # What are the DT 4-paths?
    dt_paths = []
    for path in a[3]:
        a_,b_,c_,d_ = path
        if A[a_][c_] and A[b_][d_]:
            dt_paths.append(path)
    print(f"  DT 4-paths: {len(dt_paths)}")

    # What are the non-DT Ω_3 elements?
    if dim_om3 > 0:
        # Check each Ω_3 basis element
        om3_basis = om[3]  # |A_3| × dim_om3
        print(f"  Ω_3 basis vectors ({dim_om3} total):")
        for j in range(min(dim_om3, 15)):
            col = om3_basis[:, j]
            nonzero = [(a[3][i], col[i]) for i in range(len(col)) if abs(col[i]) > 1e-10]
            # Check if this is a single DT path
            if len(nonzero) == 1:
                path, coeff = nonzero[0]
                is_dt = A[path[0]][path[2]] and A[path[1]][path[3]]
                print(f"    v{j}: {path} (DT={is_dt})")
            else:
                paths_str = ", ".join(f"{c:.0f}*{p}" for p,c in nonzero[:4])
                if len(nonzero) > 4:
                    paths_str += f" + {len(nonzero)-4} more"
                print(f"    v{j}: {paths_str}")

    # ker(∂_3): the single kernel element
    bd3 = build_full_boundary_matrix(a[3], a[2])
    bd3_om = bd3 @ om[3]
    im3_coords, _, _, _ = np.linalg.lstsq(om[2], bd3_om, rcond=None)
    _, S3, Vt3 = np.linalg.svd(im3_coords)
    rank_d3 = sum(s > 1e-8 for s in S3)
    ker_d3_dim = dim_om3 - rank_d3
    print(f"  rank(∂_3|_{{Ω_3}})={rank_d3}, ker(∂_3)={ker_d3_dim}")

    if ker_d3_dim > 0:
        ker_basis = Vt3[rank_d3:]
        for k_idx in range(ker_d3_dim):
            ker_vec = ker_basis[k_idx]
            # Express in A_3 coordinates
            ker_chain = om[3] @ ker_vec
            nonzero = [(a[3][i], ker_chain[i]) for i in range(len(ker_chain)) if abs(ker_chain[i]) > 1e-10]
            paths_str = ", ".join(f"{c:.2f}*{p}" for p,c in nonzero[:6])
            if len(nonzero) > 6:
                paths_str += f" + {len(nonzero)-6} more"
            print(f"  ker(∂_3) element: {paths_str}")

    # Check: is the ker(∂_3) element related to β_3?
    if dim_om3 > 0:
        dim_om4 = om[4].shape[1] if om[4].ndim == 2 else 0
        if dim_om4 > 0:
            bd4 = build_full_boundary_matrix(a[4], a[3])
            bd4_om = bd4 @ om[4]
            im4_coords, _, _, _ = np.linalg.lstsq(om[3], bd4_om, rcond=None)
            rank_d4 = np.linalg.matrix_rank(im4_coords, tol=1e-8)
            beta3 = ker_d3_dim - rank_d4
            print(f"  rank(∂_4)={rank_d4}, β_3={beta3}")
        else:
            print(f"  Ω_4=0, β_3={ker_d3_dim}")

# ===== Structural analysis =====
print(f"\n\n{'='*70}")
print("STRUCTURAL PATTERNS IN TIGHT CASES")
print("='*70")

# What score sequences appear in tight cases?
score_counter = Counter()
for tidx, A, surplus, dom2, dom3, z2 in tight_tournaments:
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    score_counter[scores] += 1

print(f"\nScore sequences in tight cases:")
for scores, count in sorted(score_counter.items()):
    print(f"  {scores}: {count}")

# How many 3-cycles?
t3_counter = Counter()
for tidx, A, surplus, dom2, dom3, z2 in tight_tournaments:
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: t3 += 1
                if A[j][i] and A[i][k] and A[k][j]: t3 += 1
    t3_counter[t3] += 1

print(f"\n3-cycle counts in tight cases:")
for t3, count in sorted(t3_counter.items()):
    print(f"  t3={t3}: {count}")

# ===== Compare: what makes DT insufficiency happen? =====
print(f"\n\n{'='*70}")
print("DT SUFFICIENCY vs TIGHTNESS")
print("=" * 70)

# Check: are tight cases exactly the DT-insufficient ones?
dt_insuff = 0
dt_suff = 0
for tidx, A, surplus, dom2, dom3, z2 in tight_tournaments:
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
    dt_count = sum(1 for path in a3 if A[path[0]][path[2]] and A[path[1]][path[3]])
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]

    # DT subspace
    dt_indices = [i for i, path in enumerate(a3) if A[path[0]][path[2]] and A[path[1]][path[3]]]

    # Check if DT alone fills Z_2
    bd3 = build_full_boundary_matrix(a3, a2)
    bd2 = build_full_boundary_matrix(a2, a1)

    dt_bd_matrix = bd3[np.ix_(range(len(a2)), dt_indices)]
    rank_dt_bd = np.linalg.matrix_rank(dt_bd_matrix, tol=1e-8)

    # rank of bd2
    rank_bd2 = np.linalg.matrix_rank(bd2, tol=1e-8)

    # Need rank_dt_bd = z2 for DT sufficiency
    # But z2 involves Ω_2, not A_2...
    if rank_dt_bd < z2:
        dt_insuff += 1
    else:
        dt_suff += 1

print(f"  Tight cases: DT-sufficient={dt_suff}, DT-insufficient={dt_insuff}")

print("\nDone.")
