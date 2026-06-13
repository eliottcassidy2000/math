#!/usr/bin/env python3
"""
beta2_non_tournament.py — Find digraphs with β₂ > 0

Tournaments always have β₂=0. To understand WHY, find non-tournament
digraphs where β₂ > 0. This will reveal what structural property
of tournaments prevents nonzero β₂.

Key tournament property: between every pair, exactly one arc.
What happens if we:
1. Remove an arc (making a "near-tournament")?
2. Add a reversed arc (making a bidirectional edge)?
3. Use a random sparse digraph?

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def compute_betti_2(A, n):
    """Compute β₂ for digraph A."""
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)

    if not ap2:
        return 0, {}

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    d2 = om2.shape[1] if om2.ndim == 2 and om2.shape[0] > 0 else 0
    if d2 == 0:
        return 0, {}

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2 = d2 - rk2

    if not ap3:
        return z2, {'z2': z2, 'b2': 0, 'd2': d2}

    om3 = compute_omega_basis(A, n, 3, ap3, ap2)
    d3 = om3.shape[1] if om3.ndim == 2 and om3.shape[0] > 0 else 0
    if d3 == 0:
        return z2, {'z2': z2, 'b2': 0, 'd2': d2, 'd3': 0}

    bd3 = build_full_boundary_matrix(ap3, ap2)
    bd3_om = bd3 @ om3
    bd3_coords = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
    b2 = np.linalg.matrix_rank(bd3_coords, tol=1e-8)

    return z2 - b2, {'z2': z2, 'b2': b2, 'd2': d2, 'd3': d3}


print("=" * 70)
print("FINDING DIGRAPHS WITH β₂ > 0")
print("=" * 70)

n = 5

# Test 1: Remove one arc from each tournament
print("\n--- Test 1: Remove one arc from tournament ---")
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)
found_b2 = 0

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    # Try removing each arc
    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue
            B = [row[:] for row in A]
            B[u][v] = 0
            beta2, info = compute_betti_2(B, n)
            if beta2 > 0:
                found_b2 += 1
                if found_b2 <= 3:
                    print(f"  T#{bits} remove ({u}→{v}): β₂={beta2}, {info}")
    if bits % 200 == 0 and bits > 0:
        print(f"  ... {bits}/{1 << m}, found β₂>0: {found_b2}")

print(f"Total near-tournaments with β₂>0: {found_b2}")

# Test 2: Add a reverse arc (creating bidirectional edge)
print(f"\n--- Test 2: Add reverse arc (bidirectional edge) ---")
found_b2_bi = 0

for bits in range(min(256, 1 << m)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v] or A[v][u]:
                continue
            B = [row[:] for row in A]
            B[v][u] = 1  # now both u→v and v→u
            beta2, info = compute_betti_2(B, n)
            if beta2 > 0:
                found_b2_bi += 1
                if found_b2_bi <= 3:
                    print(f"  T#{bits} add ({v}→{u}): β₂={beta2}, {info}")

print(f"Total bidirectional-edge digraphs with β₂>0: {found_b2_bi} (out of 256 tournaments)")

# Test 3: Small complete bipartite digraphs
print(f"\n--- Test 3: Known digraphs ---")

# Cycle C_3 (already a tournament)
C3 = [[0,1,0],[0,0,1],[1,0,0]]
b2, info = compute_betti_2(C3, 3)
print(f"  C₃ cycle: β₂={b2}")

# Cycle C_4
C4 = [[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,0]]
b2, info = compute_betti_2(C4, 4)
print(f"  C₄ cycle: β₂={b2}, {info}")

# Complete bipartite K_{2,2}
K22 = [[0,0,1,1],[0,0,1,1],[0,0,0,0],[0,0,0,0]]
b2, info = compute_betti_2(K22, 4)
print(f"  K₂,₂ (all left→right): β₂={b2}, {info}")

# Bidirectional complete bipartite
K22b = [[0,0,1,1],[0,0,1,1],[1,1,0,0],[1,1,0,0]]
b2, info = compute_betti_2(K22b, 4)
print(f"  K₂,₂ bidirectional: β₂={b2}, {info}")

# Complete digraph K₄
K4 = [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]]
b2, info = compute_betti_2(K4, 4)
print(f"  Complete digraph K₄: β₂={b2}, {info}")

# "Almost complete" minus one arc
K4m = [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,0,0]]
b2, info = compute_betti_2(K4m, 4)
print(f"  K₄ minus (3→2): β₂={b2}, {info}")

# Complete digraph K₅
K5 = [[1 if i != j else 0 for j in range(5)] for i in range(5)]
b2, info = compute_betti_2(K5, 5)
print(f"  Complete digraph K₅: β₂={b2}, {info}")

# Test 4: Systematic search over all 4-vertex digraphs with β₂ > 0
print(f"\n--- Test 4: All 4-vertex digraphs ---")
n4 = 4
count_4 = 0
beta2_pos = 0
beta2_dist = Counter()

# 12 possible directed edges, 2^12 = 4096 digraphs (no self-loops)
for mask in range(1 << 12):
    A = [[0]*4 for _ in range(4)]
    idx = 0
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            if mask & (1 << idx):
                A[i][j] = 1
            idx += 1
    count_4 += 1
    b2, info = compute_betti_2(A, 4)
    beta2_dist[b2] += 1
    if b2 > 0:
        beta2_pos += 1
        if beta2_pos <= 5:
            arcs = [(i,j) for i in range(4) for j in range(4) if i != j and A[i][j]]
            print(f"  mask={mask}: β₂={b2}, arcs={arcs}")

print(f"\nβ₂ distribution over {count_4} 4-vertex digraphs: {dict(beta2_dist)}")
print(f"Digraphs with β₂>0: {beta2_pos}")

# Of these, how many are tournaments?
tour_count = 0
tour_b2 = 0
for mask in range(1 << 12):
    A = [[0]*4 for _ in range(4)]
    idx = 0
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            if mask & (1 << idx):
                A[i][j] = 1
            idx += 1
    # Check if tournament
    is_tour = True
    for i in range(4):
        for j in range(i+1, 4):
            if A[i][j] + A[j][i] != 1:
                is_tour = False
                break
        if not is_tour:
            break
    if is_tour:
        tour_count += 1
        b2, _ = compute_betti_2(A, 4)
        if b2 > 0:
            tour_b2 += 1

print(f"\n4-vertex tournaments: {tour_count}, with β₂>0: {tour_b2}")

print("\nDone.")
