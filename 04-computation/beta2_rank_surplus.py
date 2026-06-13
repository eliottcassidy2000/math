#!/usr/bin/env python3
"""Analyze the rank surplus of ∂_3: Ω_3 → ker(∂_2) for all tournaments.

Key question: How much "extra" capacity does Ω_3 have for filling 2-cycles?

Define:
  Z_2 = ker(∂_2) ∩ Ω_2 (space of 2-cycles)
  B_2 = im(∂_3|_{Ω_3}) (space of 2-boundaries)
  surplus = dim(Ω_3) - dim(Z_2) (excess 3-chains beyond what's needed)
  rank_surplus = dim(B_2) - dim(Z_2) (must be ≥ 0 for β_2=0)

Also analyze: does dim(Ω_3) ≥ dim(Z_2) ALWAYS hold?
If dim(Ω_3) < dim(Z_2), then β_2=0 requires ∂_3 to be injective
AND have full rank on Z_2 — a tighter constraint.
"""
import numpy as np
from itertools import combinations
import sys
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
    path_betti_numbers
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

def analyze_tournament(A, n):
    """Full rank analysis for β_2."""
    a = {}
    for p in range(5):
        a[p] = [tuple(x) for x in enumerate_allowed_paths(A, n, p)]

    om = {}
    for p in range(5):
        if p == 0:
            om[p] = np.eye(n)
        else:
            om[p] = compute_omega_basis(A, n, p, a[p], a[p-1])

    dim_om = {}
    for p in range(5):
        if om[p].ndim == 2:
            dim_om[p] = om[p].shape[1]
        else:
            dim_om[p] = 0

    # ∂_2: Ω_2 → Ω_1
    if dim_om[2] == 0:
        return None  # trivial

    bd2 = build_full_boundary_matrix(a[2], a[1])
    bd2_om = bd2 @ om[2]
    S2 = np.linalg.svd(bd2_om, compute_uv=False)
    rank_d2 = sum(s > 1e-8 for s in S2)
    dim_Z2 = dim_om[2] - rank_d2  # = dim ker(∂_2)

    # ∂_3: Ω_3 → Ω_2
    if dim_om[3] == 0:
        rank_d3 = 0
        dim_ker_d3 = 0
    else:
        bd3 = build_full_boundary_matrix(a[3], a[2])
        bd3_om = bd3 @ om[3]
        # Express im(∂_3) in Ω_2 coordinates
        im3_coords, _, _, _ = np.linalg.lstsq(om[2], bd3_om, rcond=None)
        S3 = np.linalg.svd(im3_coords, compute_uv=False)
        rank_d3 = sum(s > 1e-8 for s in S3)  # = dim(B_2)
        dim_ker_d3 = dim_om[3] - rank_d3

    beta_2 = dim_Z2 - rank_d3
    surplus_dim = dim_om[3] - dim_Z2  # raw dimension surplus
    surplus_rank = rank_d3 - dim_Z2  # boundary rank surplus (should be 0 if β_2=0)
    injectivity = dim_ker_d3  # how far ∂_3 is from injective

    t3 = sum(1 for i in range(n) for j in range(n) for k in range(n)
             if i<j<k and A[i][j] and A[j][k] and A[i][k]) + \
         sum(1 for i in range(n) for j in range(n) for k in range(n)
             if i<j<k and A[j][i] and A[k][j] and A[k][i])
    # Actually just count 3-cycles
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check all 3 possible cycle orientations
                if A[i][j] and A[j][k] and A[k][i]: t3 += 1
                if A[j][i] and A[i][k] and A[k][j]: t3 += 1

    return {
        'dim_om2': dim_om[2],
        'dim_om3': dim_om[3],
        'dim_Z2': dim_Z2,
        'dim_B2': rank_d3,
        'beta_2': beta_2,
        'surplus_dim': surplus_dim,
        'surplus_rank': surplus_rank,
        'dim_ker_d3': dim_ker_d3,
        'score': score_seq(A, n),
        't3': t3
    }

# ===== n=5 =====
print("=" * 70)
print("RANK SURPLUS ANALYSIS: β_2 = 0")
print("=" * 70)

for n in [5, 6]:
    print(f"\n--- n={n} ---")

    surplus_dist = Counter()
    dim_surplus_dist = Counter()
    pattern_dist = Counter()
    min_surplus = float('inf')
    max_surplus = float('-inf')
    min_dim_surplus = float('inf')

    total = 0
    for tidx, A in enumerate(all_tournaments_gen(n)):
        if n == 6 and tidx % 5000 == 0:
            print(f"  ... {tidx}", flush=True)

        res = analyze_tournament(A, n)
        if res is None:
            continue
        total += 1

        surplus_dist[res['surplus_rank']] += 1
        dim_surplus_dist[res['surplus_dim']] += 1

        pattern = (res['dim_om2'], res['dim_om3'], res['dim_Z2'], res['dim_B2'])
        pattern_dist[pattern] += 1

        min_surplus = min(min_surplus, res['surplus_dim'])
        max_surplus = max(max_surplus, res['surplus_dim'])
        min_dim_surplus = min(min_dim_surplus, res['surplus_dim'])

        if n == 5 or (n == 6 and res['surplus_dim'] <= 1):
            if res['surplus_dim'] <= 2:
                print(f"    Tight: Ω3={res['dim_om3']}, Z2={res['dim_Z2']}, B2={res['dim_B2']}, "
                      f"surplus_dim={res['surplus_dim']}, ker∂3={res['dim_ker_d3']}, "
                      f"score={res['score']}, t3={res['t3']}")

    print(f"\n  Total: {total}")
    print(f"  Rank surplus (B2 - Z2) distribution:")
    for k in sorted(surplus_dist):
        print(f"    {k}: {surplus_dist[k]}")

    print(f"  Dimension surplus (dim Ω_3 - dim Z_2) distribution:")
    for k in sorted(dim_surplus_dist):
        print(f"    {k}: {dim_surplus_dist[k]}")

    print(f"  Min dimension surplus: {min_dim_surplus}")

    print(f"\n  (dim_Ω2, dim_Ω3, dim_Z2, dim_B2) patterns:")
    for pat in sorted(pattern_dist):
        count = pattern_dist[pat]
        om2, om3, z2, b2 = pat
        ratio = om3/z2 if z2 > 0 else float('inf')
        print(f"    {pat}  β2={z2-b2}  Ω3/Z2={ratio:.2f}  ×{count}")

# ===== n=7 sample =====
print(f"\n\n--- n=7 (sample) ---")
import random
random.seed(42)

n = 7
sample_size = 2000
min_surplus_7 = float('inf')
pattern_dist_7 = Counter()

for _ in range(sample_size):
    A = [[0]*n for _ in range(n)]
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for i,j in edges:
        if random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

    res = analyze_tournament(A, n)
    if res is None:
        continue

    pattern = (res['dim_om2'], res['dim_om3'], res['dim_Z2'], res['dim_B2'])
    pattern_dist_7[pattern] += 1
    min_surplus_7 = min(min_surplus_7, res['surplus_dim'])

print(f"  Sample: {sample_size}")
print(f"  Min dimension surplus: {min_surplus_7}")
print(f"\n  (dim_Ω2, dim_Ω3, dim_Z2, dim_B2) patterns:")
for pat in sorted(pattern_dist_7):
    count = pattern_dist_7[pat]
    om2, om3, z2, b2 = pat
    ratio = om3/z2 if z2 > 0 else float('inf')
    print(f"    {pat}  β2={z2-b2}  Ω3/Z2={ratio:.2f}  ×{count}")

print("\nDone.")
