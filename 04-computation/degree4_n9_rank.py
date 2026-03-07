#!/usr/bin/env python3
"""
Determine the RANK of the degree-4 Fourier space at n=9.

Compute degree-4 Walsh coefficients for several random tournaments,
form a matrix, find its rank.

If rank = 2, the n=7 structure persists.
If rank > 2, there are new independent types.

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
    """Compute degree-4 Walsh coefficients by enumerating Hamiltonian paths."""
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

# Gather all monomials that appear across tournaments
all_monomials = set()
tournament_coeffs = []

num_tournaments = 20
print(f"Computing degree-4 coefficients for {num_tournaments} random n=9 tournaments...")

for seed in range(num_tournaments):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t0 = time.time()
    coeffs, H = degree4_coefficients(A)
    t1 = time.time()
    print(f"  seed={seed}: H={H}, nonzero monomials={len(coeffs)}, time={t1-t0:.1f}s")

    for mono in coeffs:
        if coeffs[mono] != 0:
            all_monomials.add(mono)
    tournament_coeffs.append(coeffs)

# Build matrix: rows = tournaments, columns = monomials
mono_list = sorted(all_monomials)
print(f"\nTotal unique nonzero monomials across all tournaments: {len(mono_list)}")

M = np.zeros((num_tournaments, len(mono_list)))
for i, coeffs in enumerate(tournament_coeffs):
    for j, mono in enumerate(mono_list):
        M[i, j] = coeffs.get(mono, 0)

# Compute rank
U, S, Vh = np.linalg.svd(M, full_matrices=False)
# Count significant singular values
threshold = S[0] * 1e-10
rank = np.sum(S > threshold)
print(f"\nSingular values: {S[:10]}")
print(f"Rank of degree-4 coefficient matrix: {rank}")
print(f"(Threshold: {threshold:.2e})")

# Also check by graph type
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

# Group monomials by type and check rank within each type
type_indices = defaultdict(list)
for j, mono in enumerate(mono_list):
    t = classify_monomial(mono)
    type_indices[t].append(j)

print(f"\nRank by monomial type:")
for t in sorted(type_indices.keys()):
    cols = type_indices[t]
    M_sub = M[:, cols]
    U_s, S_s, Vh_s = np.linalg.svd(M_sub, full_matrices=False)
    r = np.sum(S_s > S_s[0] * 1e-10) if S_s[0] > 0 else 0
    print(f"  {t}: {len(cols)} monomials, rank = {r}")
    if r <= 5:
        print(f"    singular values: {S_s[:r+2]}")

# Check: is each type proportional (rank 1)?
# If so, the total rank equals the number of types with nonzero contribution.
print(f"\nTotal types: {len(type_indices)}")
print(f"If each type has rank 1, total rank should be {len(type_indices)}")
