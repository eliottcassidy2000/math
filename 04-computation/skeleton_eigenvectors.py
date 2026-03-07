#!/usr/bin/env python3
"""
Skeleton eigenvectors at n=5 — what do they mean?

The skeleton has 8 SC classes with eigenvalues ±(1+√2), ±1, ±1, ±(√2-1).
The ±(1+√2) eigenspaces are t3-independent.

Questions:
1. What are the eigenvectors explicitly?
2. Do they correspond to known invariants?
3. What is the graph-theoretic meaning of the silver ratio?

Also: the skeleton is bipartite, so eigenvalues come in ±λ pairs.
The positive eigenspace for each |λ| is the "symmetric" combination,
the negative is the "antisymmetric" (across the bipartition).

kind-pasteur-2026-03-07-S26
"""
import numpy as np
from itertools import permutations, combinations
from fractions import Fraction
from collections import defaultdict

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_5cycles(A, n):
    t5 = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                t5 += 1
    return t5 // 5

def canonical_form(A, n):
    best = None
    for perm in permutations(range(n)):
        flat = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or flat < best:
            best = flat
    return best

def tiling_transpose_pairs(n):
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    edge_to_idx = {e: idx for idx, e in enumerate(edges)}
    pairs, fixed, seen = [], [], set()
    for idx, (i, j) in enumerate(edges):
        if idx in seen: continue
        ti, tj = n-1-j, n-1-i
        if ti > tj: ti, tj = tj, ti
        if (ti, tj) == (i, j):
            fixed.append(idx); seen.add(idx)
        elif (ti, tj) in edge_to_idx:
            tidx = edge_to_idx[(ti, tj)]
            pairs.append((idx, tidx)); seen.add(idx); seen.add(tidx)
    return pairs, fixed

def gen_gs_tilings(n, pairs, fixed):
    gs_dof = len(pairs) + len(fixed)
    result = []
    for free_val in range(2**gs_dof):
        bits = 0
        for k, (idx1, idx2) in enumerate(pairs):
            if (free_val >> k) & 1:
                bits |= (1 << idx1) | (1 << idx2)
        for k, fidx in enumerate(fixed):
            if (free_val >> (len(pairs) + k)) & 1:
                bits |= (1 << fidx)
        result.append(bits)
    return result

def flip_tiling(tiling_bits, m):
    return tiling_bits ^ ((1 << m) - 1)

n = 5
m = num_tiling_bits(n)
pairs, fixed = tiling_transpose_pairs(n)
gs_tilings = gen_gs_tilings(n, pairs, fixed)

# Build class data
class_map = {}
class_info = []
for idx, bits in enumerate(gs_tilings):
    A = tournament_from_tiling(n, bits)
    canon = canonical_form(A, n)
    if canon not in class_map:
        class_id = len(class_map)
        class_map[canon] = class_id
        t3 = count_3cycles(A, n)
        t5 = count_5cycles(A, n)
        H = 1 + 2*t3 + 2*t5  # OCF at n=5
        class_info.append({'t3': t3, 't5': t5, 'H': H, 'canon': canon,
                          'tilings': [idx], 'scores': tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))})
    else:
        class_info[class_map[canon]]['tilings'].append(idx)

num_classes = len(class_map)

# Build skeleton
K = np.zeros((num_classes, num_classes))
for idx, bits in enumerate(gs_tilings):
    flipped = flip_tiling(bits, m)
    A_orig = tournament_from_tiling(n, bits)
    A_flip = tournament_from_tiling(n, flipped)
    canon_orig = canonical_form(A_orig, n)
    canon_flip = canonical_form(A_flip, n)
    K[class_map[canon_orig], class_map[canon_flip]] += 1

# Eigendecomposition
evals, evecs = np.linalg.eigh(K)

print(f"n={n}: {num_classes} SC classes")
print(f"\nClass details:")
for i in range(num_classes):
    info = class_info[i]
    side = 'A' if info['t3'] % 2 == 1 else 'B'  # Bipartition side
    print(f"  Class {i}: t3={info['t3']}, t5={info['t5']}, H={info['H']}, "
          f"scores={info['scores']}, |class|={len(info['tilings'])}, side={side}")

print(f"\nSkeleton K:")
print(K)

print(f"\nEigenvalues and eigenvectors:")
sqrt2 = np.sqrt(2)
for i in range(num_classes):
    ev = evals[i]
    vec = evecs[:, i]
    # Identify eigenvalue
    if abs(ev - (1+sqrt2)) < 0.01: name = "1+sqrt2"
    elif abs(ev + (1+sqrt2)) < 0.01: name = "-(1+sqrt2)"
    elif abs(ev - 1) < 0.01: name = "1"
    elif abs(ev + 1) < 0.01: name = "-1"
    elif abs(ev - (sqrt2-1)) < 0.01: name = "sqrt2-1"
    elif abs(ev + (sqrt2-1)) < 0.01: name = "-(sqrt2-1)"
    else: name = f"{ev:.4f}"

    # Normalize for readability
    # Scale so max abs component is 1
    scale = max(abs(vec))
    vec_scaled = vec / scale

    print(f"\n  eigenvalue = {name} ({ev:.6f}):")
    print(f"    eigenvector (scaled): {np.round(vec_scaled, 4)}")

    # Express in terms of class properties
    t3_values = [class_info[j]['t3'] for j in range(num_classes)]
    print(f"    t3 values by class:   {t3_values}")
    print(f"    side (odd/even t3):   {''.join('A' if t%2==1 else 'B' for t in t3_values)}")

# Check: what functions are eigenvectors?
# The bipartition sign: +1 for side A (odd t3), -1 for side B (even t3)
bipartition = np.array([(-1)**info['t3'] for info in class_info], dtype=float)
# Note: (-1)^t3 = +1 when t3 even (side B), -1 when t3 odd (side A)
# Let's use: +1 = side B (even), -1 = side A (odd)

print(f"\nBipartition vector ((-1)^t3): {bipartition}")
Kb = K @ bipartition
print(f"K * bipartition = {Kb}")
# If this is proportional to bipartition, it's an eigenvector

# t3 - centered
t3_vec = np.array([info['t3'] for info in class_info], dtype=float)
t3c = t3_vec - np.mean(t3_vec)

# Degree vector (row sum of K = number of GS tilings in class that are flip-mapped)
degree = K.sum(axis=1)
print(f"\nDegree vector (row sums): {degree}")
print(f"Class sizes (# GS tilings): {[len(info['tilings']) for info in class_info]}")
# At n=5: degree = class size (each GS tiling maps to exactly one other class under flip)

# The uniform vector on each side
side_A = np.array([1.0 if info['t3'] % 2 == 1 else 0.0 for info in class_info])
side_B = np.array([1.0 if info['t3'] % 2 == 0 else 0.0 for info in class_info])

KA = K @ side_A
KB = K @ side_B
print(f"\nK * (side A indicator): {KA}")
print(f"K * (side B indicator): {KB}")
print(f"side_A = {side_A}")
print(f"side_B = {side_B}")

# Check: does K map side_A to some linear combination of side_B?
# Since K is bipartite, K maps A -> B and B -> A
print(f"\nSince K is bipartite: K*v_A lives in span(B), K*v_B lives in span(A)")

# The H vector
H_vec = np.array([info['H'] for info in class_info], dtype=float)
print(f"\nH vector: {H_vec}")
print(f"K * H = {K @ H_vec}")
print(f"H components on eigenvectors:")
for i in range(num_classes):
    comp = np.dot(H_vec, evecs[:, i])
    print(f"  ev={evals[i]:+.4f}: component = {comp:.6f}")
