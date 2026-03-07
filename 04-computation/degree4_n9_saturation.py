#!/usr/bin/env python3
"""
Check saturation of the degree-4 Fourier rank at n=9.

Instead of using the full coefficient vector (too large),
pick a random projection and track rank growth.

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

def degree4_coefficients_sparse(A, target_monos):
    """Compute degree-4 coefficients only for specified monomials."""
    coeffs = defaultdict(int)
    target_set = set(target_monos)

    def backtrack(path, visited):
        if len(path) == n:
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
                if mono in target_set:
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
    return dict(coeffs)

# Strategy: pick ~200 random degree-4 monomials and track rank as we add tournaments
rng = random.Random(999)

# Generate random monomials (4-element subsets of edges)
num_probe_monos = 300
probe_monos = []
while len(probe_monos) < num_probe_monos:
    mono = tuple(sorted(rng.sample(range(m), 4)))
    if mono not in probe_monos:
        probe_monos.append(mono)

print(f"Using {num_probe_monos} random probe monomials")
print(f"Computing degree-4 coefficients for increasing tournament count...\n")

max_tournaments = 200
rows = []

for seed in range(max_tournaments):
    r = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if r.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    coeffs = degree4_coefficients_sparse(A, probe_monos)
    row = [coeffs.get(mono, 0) for mono in probe_monos]
    rows.append(row)

    if (seed + 1) % 10 == 0 or seed < 5:
        M = np.array(rows)
        U, S, Vh = np.linalg.svd(M, full_matrices=False)
        threshold = max(S[0] * 1e-10, 1e-6)
        rank = np.sum(S > threshold)
        print(f"  {seed+1} tournaments: rank = {rank}, smallest kept SV = {S[rank-1]:.2f}, next = {S[rank] if rank < len(S) else 0:.2f}")

M = np.array(rows)
U, S, Vh = np.linalg.svd(M, full_matrices=False)
threshold = max(S[0] * 1e-10, 1e-6)
rank = np.sum(S > threshold)

print(f"\nFinal: {max_tournaments} tournaments, rank = {rank}")
print(f"Smallest nonzero SV: {S[rank-1]:.4f}")
print(f"Next SV: {S[rank]:.6f}" if rank < len(S) else "All SVs are nonzero")
print(f"\nSingular value spectrum (first 20): {S[:20]}")
