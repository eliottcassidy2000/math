#!/usr/bin/env python3
"""
beta4_structure_n7.py — Why is beta_4 = 6 for H-maximizers at n=7?

beta_4 = 6 = dim H_4(path complex of T_7)
This means there are 6 independent 4-dimensional "holes" in the path complex.

What does the path complex look like at dim 4?
- Allowed 4-paths: v0->v1->v2->v3->v4 (5 distinct vertices, all edges present)
- These are Hamiltonian paths of 5-vertex induced subtournaments
- The boundary map partial_4: Omega_4 -> Omega_3 takes alternating sums

Questions:
1. What is dim(Omega_4) for T_7? (Number of allowed 4-paths)
2. What is dim(ker partial_4)?
3. What is dim(im partial_5)? (From 5-paths = Ham paths on 6-vertex subtournaments)
4. beta_4 = dim(ker partial_4) / dim(im partial_5) = 6

The answer 6 might be: C(7,5) - dim(im partial_5) = 21 - 15 = 6?
Or something else entirely.

Author: kind-pasteur-2026-03-08-S40
"""
import sys, os, time
import numpy as np
from itertools import combinations, permutations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers, enumerate_allowed_paths, build_full_boundary_matrix
sys.stdout = _saved

def paley_tournament(p):
    qr = set((a*a) % p for a in range(1, p))
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                A[i][j] = 1
    return A

# Paley T_7
A7 = paley_tournament(7)
n = 7

print("=" * 70)
print("PALEY T_7: Path complex structure")
print("=" * 70)

# Count allowed paths at each dimension
for dim in range(8):
    paths = enumerate_allowed_paths(A7, n, dim)
    print(f"dim {dim}: {len(paths)} allowed paths")

# Build boundary matrices
print("\n--- Boundary matrices ---")
dims = {}
for dim in range(8):
    paths = enumerate_allowed_paths(A7, n, dim)
    dims[dim] = len(paths)

all_paths = {}
for dim in range(8):
    all_paths[dim] = enumerate_allowed_paths(A7, n, dim)

for dim in range(1, 7):
    B = build_full_boundary_matrix(all_paths[dim], all_paths[dim-1])
    if B is not None and B.size > 0:
        r = np.linalg.matrix_rank(B)
        print(f"partial_{dim}: {B.shape[0]} x {B.shape[1]}, rank = {r}")
    else:
        print(f"partial_{dim}: empty or None")

# Compute Betti manually
print("\n--- Betti numbers (manual) ---")
print("beta_p = dim(ker partial_p) - dim(im partial_{p+1})")
print("      = (dim Omega_p - rank partial_p) - rank partial_{p+1}")

ranks = {}
for dim in range(1, 7):
    B = build_full_boundary_matrix(all_paths[dim], all_paths[dim-1])
    if B is not None and B.size > 0:
        ranks[dim] = np.linalg.matrix_rank(B)
    else:
        ranks[dim] = 0

for p in range(7):
    rank_p = ranks.get(p, 0)
    rank_p1 = ranks.get(p+1, 0)
    ker_p = dims[p] - rank_p
    beta_p = ker_p - rank_p1
    print(f"  beta_{p} = {dims[p]} - {rank_p} - {rank_p1} = {beta_p}")

# Specifically for beta_4:
print("\n" + "=" * 70)
print("BETA_4 DECOMPOSITION")
print("=" * 70)

paths_4 = enumerate_allowed_paths(A7, n, 4)
paths_5 = enumerate_allowed_paths(A7, n, 5)
paths_6 = enumerate_allowed_paths(A7, n, 6)

print(f"Omega_4: {len(paths_4)} allowed 4-paths (on 5 vertices each)")
print(f"Omega_5: {len(paths_5)} allowed 5-paths (on 6 vertices each)")
print(f"Omega_6: {len(paths_6)} allowed 6-paths (on 7 vertices = Ham paths)")

# Which 5-vertex subsets have allowed 4-paths?
subsets_4 = set()
for p in paths_4:
    subsets_4.add(frozenset(p))
print(f"\nDistinct 5-vertex subsets with allowed 4-paths: {len(subsets_4)}")
print(f"Total 5-vertex subsets: {len(list(combinations(range(n), 5)))}")

# How many 4-paths per 5-vertex subset?
paths_per_subset = {}
for p in paths_4:
    s = frozenset(p)
    paths_per_subset[s] = paths_per_subset.get(s, 0) + 1

from collections import Counter
count_dist = Counter(paths_per_subset.values())
print(f"4-paths per subset: {dict(sorted(count_dist.items()))}")

# Which 6-vertex subsets have allowed 5-paths?
subsets_5 = set()
for p in paths_5:
    subsets_5.add(frozenset(p))
print(f"\nDistinct 6-vertex subsets with allowed 5-paths: {len(subsets_5)}")
print(f"5-paths per subset: ", end="")
pp5 = {}
for p in paths_5:
    s = frozenset(p)
    pp5[s] = pp5.get(s, 0) + 1
print(dict(sorted(Counter(pp5.values()).items())))

# Ham paths on 7 vertices
print(f"\n6-paths (Ham paths on 7 vertices): {len(paths_6)}")
# Should be H(T_7) = 189
# But wait, paths_6 counts ORDERED sequences v0->...->v6
# Each Ham path is a specific ordering, so this should indeed be 189
# (or each path counted once since it's a specific sequence)

# Understanding beta_4 = 6:
# beta_4 = dim(ker partial_4) - dim(im partial_5)
# = (|Omega_4| - rank partial_4) - rank partial_5
B4 = build_full_boundary_matrix(all_paths[4], all_paths[3])
B5 = build_full_boundary_matrix(all_paths[5], all_paths[4])
r4 = np.linalg.matrix_rank(B4) if B4 is not None and B4.size > 0 else 0
r5 = np.linalg.matrix_rank(B5) if B5 is not None and B5.size > 0 else 0

print(f"\nbeta_4 = ({len(paths_4)} - {r4}) - {r5} = {len(paths_4) - r4 - r5}")
print(f"      = ker partial_4 ({len(paths_4) - r4}) - im partial_5 ({r5})")

# Let's also check for OTHER regular n=7 tournaments
print("\n" + "=" * 70)
print("NON-PALEY REGULAR n=7 TOURNAMENTS")
print("=" * 70)

import random
random.seed(42)

def random_tournament(nn):
    A = [[0]*nn for _ in range(nn)]
    for i in range(nn):
        for j in range(i+1, nn):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def H_tournament(A, nn):
    dp = [[0]*nn for _ in range(1 << nn)]
    for v in range(nn):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << nn):
        for v in range(nn):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(nn):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << nn)][u] += dp[mask][v] if (1 << nn) != (mask | (1 << u)) else 0
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << nn) - 1
    return sum(dp[full])

# The 3 rigid classes at n=7 regular:
# 1. BIBD (alpha_2=7): H=189, 240 tours
# 2. alpha_2=10: H=171, 1680 tours
# 3. alpha_2=14: H=175, 720 tours

# Find representatives of each class
classes_found = {}
for trial in range(100000):
    A = random_tournament(n)
    scores = [sum(A[i]) for i in range(n)]
    if all(s == 3 for s in scores):
        H = H_tournament(A, n)
        if H not in classes_found:
            classes_found[H] = A
        if len(classes_found) >= 3:
            break

print(f"Found regular tournament classes: {sorted(classes_found.keys())}")

for H_val in sorted(classes_found.keys(), reverse=True):
    A = classes_found[H_val]
    try:
        beta = path_betti_numbers(A, n, max_dim=6)
        beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(7)]
    except:
        beta_list = "FAILED"

    # Count dim Omega_4 for this tournament
    paths_4_t = enumerate_allowed_paths(A, n, 4)

    print(f"\nH={H_val}: beta={beta_list}, |Omega_4|={len(paths_4_t)}")
