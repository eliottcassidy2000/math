#!/usr/bin/env python3
"""
Refine rank computation for degree-4 Fourier at n=9.

Since rank >= 20 from 20 samples, the dimension is large.
Let's focus on a single graph type and determine its internal rank.

For type (5, (2,2,2,1,1)) = P4 paths on 5 vertices:
There are C(9,5) * (number of P4's on 5 vertices) such monomials.

Actually at n=7, the key insight was that monomials of the SAME graph type
had proportional coefficients. Let me check if this is still true at n=9.

If NOT proportional, the types further decompose into sub-types.

opus-2026-03-07-S37
"""
from itertools import combinations
from collections import defaultdict
import random
import time
import numpy as np

n = 9
edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
m = len(edges)  # 36
edge_idx = {e: k for k, e in enumerate(edges)}

def degree4_coefficients(A):
    coeffs = defaultdict(int)
    path_count = 0

    def backtrack(path, visited):
        nonlocal path_count
        if len(path) == n:
            path_count += 1
            path_edges = []
            path_signs = []
            for i in range(n-1):
                u, v = path[i], path[i+1]
                if u < v:
                    path_edges.append(edge_idx[(u, v)])
                    path_signs.append(1)
                else:
                    path_edges.append(edge_idx[(v, u)])
                    path_signs.append(-1)
            for combo in combinations(range(n-1), 4):
                mono = tuple(sorted(path_edges[i] for i in combo))
                sign = 1
                for i in combo:
                    sign *= path_signs[i]
                coeffs[mono] += sign
            return

        v = path[-1]
        for u in range(n):
            if u not in visited and A[v][u]:
                path.append(u)
                visited.add(u)
                backtrack(path, visited)
                path.pop()
                visited.remove(u)

    for start in range(n):
        backtrack([start], {start})
    return dict(coeffs), path_count

def classify_monomial(mono_edges):
    elist = [edges[i] for i in mono_edges]
    verts = set()
    deg = defaultdict(int)
    for u, v in elist:
        verts.add(u)
        verts.add(v)
        deg[u] += 1
        deg[v] += 1
    nv = len(verts)
    deg_seq = tuple(sorted([deg[v] for v in verts], reverse=True))
    return (nv, deg_seq)

# Focus: pick a SINGLE specific graph type and look at its coefficient structure
# For type (5, (2,2,2,1,1)): these are P4 paths (4 edges on 5 vertices forming a path)
# But also could be: C4 cycle (4 edges on 4 vertices with deg (2,2,2,2))
# Wait no — (5, (2,2,2,1,1)) has 5 vertices, so it's a P4.

# At n=7, the P4 coefficient was: c_{P4} * (product of signs based on orientation)
# And c_{P4} was the SAME for all P4's (depends only on the labeled P4 path, not T).

# The key question is whether at n=9, all P4 coefficients are proportional to
# a SINGLE function of the tournament, or whether they decompose further.

# To test this: for each labeled 5-vertex path, the coefficient should be
# some function of the tournament's edges on those 5 vertices + the remaining 4 vertices.
# At n=7, only 2 vertices remain, so the function is simple.
# At n=9, 4 vertices remain, giving more freedom.

# Let me check: for a FIXED tournament, are the P4 coefficients +/- some constant?

rng = random.Random(42)
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if rng.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

print("Computing for seed=42...")
coeffs, H = degree4_coefficients(A)
print(f"H = {H}")

# Group P4 monomials by their vertex set
p4_by_vset = defaultdict(list)
for mono, coeff in coeffs.items():
    if coeff == 0:
        continue
    t = classify_monomial(mono)
    if t == (5, (2, 2, 2, 1, 1)):
        elist = [edges[i] for i in mono]
        verts = frozenset()
        for u, v in elist:
            verts = verts | {u, v}
        p4_by_vset[verts].append((mono, coeff))

# For each 5-vertex subset, list the P4 coefficients
print(f"\nP4 monomials grouped by vertex set:")
# Sample a few
vset_list = sorted(p4_by_vset.keys(), key=lambda x: tuple(sorted(x)))
for vs in vset_list[:5]:
    print(f"  Vertices {sorted(vs)}:")
    for mono, coeff in p4_by_vset[vs]:
        print(f"    edges={[edges[i] for i in mono]}: coeff={coeff}")

# Check: within each vertex set, are all coefficients +/- the same value?
print(f"\nChecking proportionality within vertex sets:")
anomalies = 0
for vs in vset_list:
    abs_coeffs = set(abs(c) for _, c in p4_by_vset[vs])
    if len(abs_coeffs) > 1:
        anomalies += 1
        if anomalies <= 3:
            print(f"  ANOMALY: vertices {sorted(vs)}, |coeffs| = {sorted(abs_coeffs)}")

print(f"Total vertex sets with non-constant |coeff|: {anomalies} / {len(vset_list)}")

# Also check: are all P4 abs coefficients the same across ALL vertex sets?
all_abs_p4 = set()
for vs in vset_list:
    for _, c in p4_by_vset[vs]:
        all_abs_p4.add(abs(c))
print(f"Distinct |P4 coefficients| across all: {sorted(all_abs_p4)}")

# Now do the PROPER rank test: within the P4 type, use many tournaments
print(f"\n=== Computing P4 rank across 50 tournaments ===")
# First, enumerate ALL possible P4 monomials
all_p4_monos = set()
for seed in range(50):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    coeffs, _ = degree4_coefficients(A)
    for mono, coeff in coeffs.items():
        if coeff != 0 and classify_monomial(mono) == (5, (2, 2, 2, 1, 1)):
            all_p4_monos.add(mono)
    if seed < 3 or seed % 10 == 9:
        print(f"  seed={seed}: total P4 monos so far = {len(all_p4_monos)}")

p4_mono_list = sorted(all_p4_monos)
print(f"Total P4 monomials found: {len(p4_mono_list)}")

# Build matrix
M = np.zeros((50, len(p4_mono_list)))
for seed in range(50):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    coeffs, _ = degree4_coefficients(A)
    for j, mono in enumerate(p4_mono_list):
        M[seed, j] = coeffs.get(mono, 0)

U, S, Vh = np.linalg.svd(M, full_matrices=False)
threshold = S[0] * 1e-10
rank = np.sum(S > threshold)
print(f"P4 rank: {rank}")
print(f"Top singular values: {S[:rank+3]}")
