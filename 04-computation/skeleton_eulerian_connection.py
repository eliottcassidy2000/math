#!/usr/bin/env python3
"""
Connection between skeleton eigenvalues and Eulerian structure.

At n=5, the skeleton adjacency matrix K has eigenvalues:
  {±(1+√2), ±1, ±1, ±(√2-1)}

The master polynomials at n=5 are:
  F_4(r) = 120r^4 - 30r^2 + 1   (background)
  F_2(r) = 6r^2 - 1/2            (t3 coefficient / 2)
  F_0(r) = 1                      (t5 coefficient / 2)

W(r) = F_4(r) + 2*F_2(r)*t3 + 2*F_0(r)*t5

For GS tilings, the skeleton connects tiling classes by GS flip.
The flip changes t3 by delta_t3 = t3(T) - t3(flip(T)).

The KEY observation: the skeleton K matrix operates on the space of
GS tournament classes. The W-polynomial changes under flip via:
  W(T) - W(flip(T)) = C_t3 * delta_t3 + C_t5 * delta_t5

So the W-polynomial values on skeleton vertices are EIGENVECTORS of K
up to the invariant contributions.

Let me compute the skeleton adjacency matrix at n=5 and check what
happens when we project W(r) coefficients onto eigenspaces.

kind-pasteur-2026-03-07-S26
"""
from itertools import permutations, combinations
from fractions import Fraction
from collections import defaultdict
import numpy as np

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

def compute_W_poly(A, n):
    result = defaultdict(lambda: Fraction(0))
    for perm in permutations(range(n)):
        poly = {0: Fraction(1)}
        for i in range(n-1):
            s = Fraction(A[perm[i]][perm[i+1]]) - Fraction(1, 2)
            new_poly = {}
            for power, coeff in poly.items():
                new_poly[power+1] = new_poly.get(power+1, Fraction(0)) + coeff
                new_poly[power] = new_poly.get(power, Fraction(0)) + coeff * s
            poly = new_poly
        for power, coeff in poly.items():
            result[power] += coeff
    return dict(result)

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

def canonical_form(A, n):
    """Simple canonical form: min over vertex permutations of adjacency matrix."""
    best = None
    for perm in permutations(range(n)):
        flat = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or flat < best:
            best = flat
    return best

n = 5
m = num_tiling_bits(n)
pairs, fixed = tiling_transpose_pairs(n)
gs_tilings = gen_gs_tilings(n, pairs, fixed)

# Group GS tilings by canonical form (= SC isomorphism class)
class_map = {}
class_data = {}
for idx, bits in enumerate(gs_tilings):
    A = tournament_from_tiling(n, bits)
    canon = canonical_form(A, n)
    if canon not in class_map:
        class_id = len(class_map)
        class_map[canon] = class_id
        t3 = count_3cycles(A, n)
        W = compute_W_poly(A, n)
        class_data[class_id] = {
            'canon': canon, 't3': t3, 'W': W,
            'H': sum(c * Fraction(1,2)**p for p, c in W.items()),
            'tilings': []
        }
    class_data[class_map[canon]]['tilings'].append(idx)

num_classes = len(class_map)
print(f"n={n}: {num_classes} SC classes from {len(gs_tilings)} GS tilings")

# Build skeleton adjacency matrix
K = np.zeros((num_classes, num_classes))
for idx, bits in enumerate(gs_tilings):
    flipped = flip_tiling(bits, m)
    A_orig = tournament_from_tiling(n, bits)
    A_flip = tournament_from_tiling(n, flipped)
    canon_orig = canonical_form(A_orig, n)
    canon_flip = canonical_form(A_flip, n)
    c1 = class_map[canon_orig]
    c2 = class_map[canon_flip]
    K[c1, c2] += 1

print(f"\nSkeleton adjacency matrix K ({num_classes}x{num_classes}):")
print(K)

# Eigenvalues
evals, evecs = np.linalg.eigh(K)
print(f"\nEigenvalues: {sorted(evals)[::-1]}")

# The t3 vector and H vector on classes
t3_vec = np.array([class_data[i]['t3'] for i in range(num_classes)])
H_vec = np.array([float(class_data[i]['H']) for i in range(num_classes)])
W0_vec = np.array([float(sum(c * Fraction(0)**p if p > 0 else c
                            for p, c in class_data[i]['W'].items()))
                   for i in range(num_classes)])
W2_vec = np.array([float(class_data[i]['W'].get(2, 0)) for i in range(num_classes)])

print(f"\nt3 vector: {t3_vec}")
print(f"H vector:  {H_vec}")
print(f"W(0) vector: {W0_vec}")
print(f"w_2 vector: {W2_vec}")

# Project onto eigenspaces
print(f"\nProjections onto eigenvectors:")
for i, (ev, vec) in enumerate(zip(evals, evecs.T)):
    vec_norm = vec / np.linalg.norm(vec)
    t3_proj = np.dot(t3_vec, vec_norm)
    H_proj = np.dot(H_vec, vec_norm)
    W0_proj = np.dot(W0_vec, vec_norm)
    W2_proj = np.dot(W2_vec, vec_norm)
    print(f"  lambda={ev:+.4f}: t3_proj={t3_proj:.4f}, H_proj={H_proj:.4f}, "
          f"W(0)_proj={W0_proj:.4f}, w2_proj={W2_proj:.4f}")

# Check: is t3 an eigenvector of K?
Kt3 = K @ t3_vec
print(f"\nK * t3 = {Kt3}")
print(f"t3     = {t3_vec}")
# Ratio
ratios = Kt3 / t3_vec
print(f"Kt3/t3 = {ratios}")

# Check: is (t3 - mean(t3)) an eigenvector?
t3_centered = t3_vec - np.mean(t3_vec)
Kt3c = K @ t3_centered
if np.all(t3_centered != 0):
    ratios_c = Kt3c / t3_centered
    print(f"\nK * (t3-mean) = {Kt3c}")
    print(f"t3-mean = {t3_centered}")
    print(f"ratio = {ratios_c}")
