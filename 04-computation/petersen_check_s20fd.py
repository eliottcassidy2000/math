#!/usr/bin/env python3
"""
petersen_check_s20fd.py — Is G_5/Z_2 the Petersen graph?
kind-pasteur-2026-03-24-S20fd

The Petersen graph has:
  - 10 vertices, 15 edges
  - 3-regular
  - Girth 5 (smallest cycle = pentagon)
  - Diameter 2
  - Chromatic number 3
  - Not Hamiltonian (no Hamiltonian cycle)
  - Eigenvalues: 3 (x1), 1 (x5), -2 (x4)

The n=5 merged metagraph G_5/Z_2 has:
  - 10 vertices (verified)
  - ? edges, ? regularity

Check ALL Petersen graph invariants for G_5/Z_2.
Also check the WIGGLY-WEIGHTED version.
"""

import sys
import numpy as np
from math import factorial
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  IS G_5/Z_2 THE PETERSEN GRAPH?")
print("  kind-pasteur-2026-03-24-S20fd")
print("=" * 80)

n = 5
TILES = []
for y in range(1, n-1):
    for x in range(n, y+1, -1):
        TILES.append((x, y))
m = len(TILES)
N = n
VERTS = list(range(n, 0, -1))
all_perms = list(permutations(range(N)))

def bits_to_adj(bits):
    A = [[0]*N for _ in range(N)]
    for k in range(N-1): A[k][k+1] = 1
    for i in range(m):
        x, y = TILES[i]
        xi, yi = VERTS.index(x), VERTS.index(y)
        if bits[i] == 0: A[xi][yi] = 1
        else: A[yi][xi] = 1
    return A

def adj_complement(A):
    return [[1-A[i][j] if i!=j else 0 for j in range(N)] for i in range(N)]

def canonicalize(A):
    best = None
    for p in all_perms:
        s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
        if best is None or s < best: best = s
    return best

def count_hp(A):
    dp = [[0]*N for _ in range(1 << N)]
    for v in range(N): dp[1 << v][v] = 1
    for mask in range(1, 1 << N):
        for v in range(N):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(N):
                if mask & (1 << u): continue
                if A[v][u]: dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << N) - 1])

def count_aut(A):
    count = 0
    for p in all_perms:
        if all(A[p[i]][p[j]] == A[i][j] for i in range(N) for j in range(N)):
            count += 1
    return count

def score_seq(A):
    return tuple(sorted(sum(A[i]) for i in range(N)))

# Build all tilings
canon_map = {}; adj_map = {}
for mask in range(1 << m):
    bits = [(mask >> k) & 1 for k in range(m)]
    A = bits_to_adj(bits)
    adj_map[mask] = A; canon_map[mask] = canonicalize(A)

# Classes
class_data = {}
for cn in set(canon_map.values()):
    if cn in class_data: continue
    for mask, c in canon_map.items():
        if c == cn:
            A = adj_map[mask]
            class_data[cn] = {
                'comp': canonicalize(adj_complement(A)),
                'aut': count_aut(A), 'H': count_hp(A),
                'scores': score_seq(A),
                'size': sum(1 for m2, c2 in canon_map.items() if c2 == cn)
            }
            class_data[cn]['sc'] = class_data[cn]['comp'] == cn
            break

# Merged
merged_map = {cn: min(cn, class_data[cn]['comp']) for cn in class_data}
merged_list = sorted(set(merged_map.values()))
V = len(merged_list)
mcn_idx = {mcn: i for i, mcn in enumerate(merged_list)}

# Build adjacency (wiggly edges on merged graph)
A_mat = np.zeros((V, V), dtype=int)
W_mat = np.zeros((V, V), dtype=int)

for mask in range(1 << m):
    i = mcn_idx[merged_map[canon_map[mask]]]
    for wi in range(m):
        j = mcn_idx[merged_map[canon_map[mask ^ (1 << wi)]]]
        W_mat[i, j] += 1

np.fill_diagonal(W_mat, 0)
A_mat = (W_mat > 0).astype(int)
E = A_mat.sum() // 2

print(f"\n  G_5/Z_2: V = {V}, E = {E}")
print(f"  Petersen: V = 10, E = 15")

# Degree sequence
degrees = A_mat.sum(axis=1)
print(f"\n  Degree sequence: {sorted(degrees)}")
print(f"  Regular? {len(set(degrees)) == 1}")
print(f"  Petersen is 3-regular: {all(d == 3 for d in degrees)}")

# Eigenvalues
evals = sorted(np.linalg.eigvalsh(A_mat.astype(float)), reverse=True)
print(f"\n  Eigenvalues: {[f'{x:.4f}' for x in evals]}")
print(f"  Petersen eigenvalues: 3 (x1), 1 (x5), -2 (x4)")

# Girth (smallest cycle)
from collections import deque
def bfs_girth(adj, V):
    min_girth = float('inf')
    for start in range(V):
        dist = [-1] * V
        dist[start] = 0
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in range(V):
                if adj[u][v] == 0: continue
                if dist[v] == -1:
                    dist[v] = dist[u] + 1
                    queue.append(v)
                elif v != start and dist[v] >= dist[u]:
                    girth = dist[u] + dist[v] + 1
                    min_girth = min(min_girth, girth)
    return min_girth

girth = bfs_girth(A_mat, V)
print(f"\n  Girth: {girth}")
print(f"  Petersen girth: 5")

# Diameter
dist_mat = np.full((V, V), -1)
for start in range(V):
    dist_mat[start, start] = 0
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v in range(V):
            if A_mat[u][v] and dist_mat[start][v] == -1:
                dist_mat[start][v] = dist_mat[start][u] + 1
                queue.append(v)

diameter = dist_mat.max()
print(f"  Diameter: {diameter}")
print(f"  Petersen diameter: 2")

# Triangle count
triangles = 0
for i in range(V):
    for j in range(i+1, V):
        if A_mat[i][j]:
            for k in range(j+1, V):
                if A_mat[i][k] and A_mat[j][k]:
                    triangles += 1
print(f"\n  Triangles: {triangles}")
print(f"  Petersen triangles: 0 (triangle-free)")

# Chromatic number (brute force for small V)
def is_proper_coloring(adj, colors, V):
    for i in range(V):
        for j in range(i+1, V):
            if adj[i][j] and colors[i] == colors[j]:
                return False
    return True

def chromatic_number(adj, V, max_colors=V):
    from itertools import product
    for k in range(1, max_colors+1):
        for colors in product(range(k), repeat=V):
            if is_proper_coloring(adj, colors, V):
                return k
    return max_colors

chi = chromatic_number(A_mat, V, 5)
print(f"\n  Chromatic number: {chi}")
print(f"  Petersen chromatic number: 3")

# Independence number (max independent set)
max_indep = 0
for mask in range(1 << V):
    verts = [i for i in range(V) if mask & (1 << i)]
    is_indep = True
    for i in range(len(verts)):
        for j in range(i+1, len(verts)):
            if A_mat[verts[i]][verts[j]]:
                is_indep = False
                break
        if not is_indep: break
    if is_indep:
        max_indep = max(max_indep, len(verts))

print(f"  Independence number: {max_indep}")
print(f"  Petersen independence number: 4")

# Clique number
max_clique = 0
for mask in range(1 << V):
    verts = [i for i in range(V) if mask & (1 << i)]
    is_clique = True
    for i in range(len(verts)):
        for j in range(i+1, len(verts)):
            if not A_mat[verts[i]][verts[j]]:
                is_clique = False
                break
        if not is_clique: break
    if is_clique:
        max_clique = max(max_clique, len(verts))

print(f"  Clique number: {max_clique}")
print(f"  Petersen clique number: 2")

# Vertex connectivity
# (Petersen is 3-connected)

# Per-node data
print(f"\n  PER-NODE DATA:")
print(f"  {'idx':>4} {'H':>4} {'|Aut|':>5} {'SC':>3} {'deg':>4} {'Score':>20}")
for i, mcn in enumerate(merged_list):
    d = class_data[mcn]
    print(f"  {i:4d} {d['H']:4d} {d['aut']:5d} {'Y' if d['sc'] else 'N':>3} {degrees[i]:4d} {str(d['scores'])}")

# SC/NS decomposition
SC_nodes = [i for i, mcn in enumerate(merged_list) if class_data[mcn]['sc']]
NS_nodes = [i for i, mcn in enumerate(merged_list) if not class_data[mcn]['sc']]
print(f"\n  SC nodes: {len(SC_nodes)}: {SC_nodes}")
print(f"  NS nodes: {len(NS_nodes)}: {NS_nodes}")

# Check bipartite structure of SC-NS edges
sc_ns_edges = [(i,j) for i in range(V) for j in range(i+1,V) if A_mat[i][j] and
               (i in SC_nodes) != (j in SC_nodes)]
sc_sc_edges = [(i,j) for i in range(V) for j in range(i+1,V) if A_mat[i][j] and
               i in SC_nodes and j in SC_nodes]
ns_ns_edges = [(i,j) for i in range(V) for j in range(i+1,V) if A_mat[i][j] and
               i in NS_nodes and j in NS_nodes]

print(f"\n  SC-SC edges (spine): {len(sc_sc_edges)}")
print(f"  SC-NS edges (ribs): {len(sc_ns_edges)}")
print(f"  NS-NS edges (sea): {len(ns_ns_edges)}")

# THE VERDICT
print(f"\n{'='*60}")
print(f"  THE VERDICT: IS G_5/Z_2 THE PETERSEN GRAPH?")
print(f"{'='*60}")

petersen_match = (V == 10 and E == 15 and all(d == 3 for d in degrees) and
                  girth == 5 and diameter == 2 and triangles == 0 and chi == 3 and
                  max_indep == 4 and max_clique == 2)

if petersen_match:
    print(f"  YES! G_5/Z_2 IS the Petersen graph!")
else:
    print(f"  NO. G_5/Z_2 is NOT the Petersen graph.")
    print(f"  Mismatches:")
    if E != 15: print(f"    E = {E} (Petersen: 15)")
    if not all(d == 3 for d in degrees): print(f"    Not 3-regular (degrees: {sorted(set(degrees))})")
    if girth != 5: print(f"    Girth = {girth} (Petersen: 5)")
    if diameter != 2: print(f"    Diameter = {diameter} (Petersen: 2)")
    if triangles != 0: print(f"    {triangles} triangles (Petersen: 0)")
    if chi != 3: print(f"    chi = {chi} (Petersen: 3)")
    if max_indep != 4: print(f"    alpha = {max_indep} (Petersen: 4)")

# Adjacency list for reference
print(f"\n  ADJACENCY LIST:")
for i in range(V):
    neighbors = [j for j in range(V) if A_mat[i][j]]
    print(f"    {i}: {neighbors}")

print("\nDONE.")
print("=" * 80)
