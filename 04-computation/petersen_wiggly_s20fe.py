#!/usr/bin/env python3
"""
petersen_wiggly_s20fe.py — Petersen connections in the wiggly structure
kind-pasteur-2026-03-24-S20fe

The Petersen graph = Kneser K(5,2):
  Vertices = 2-element subsets of {1,2,3,4,5} (10 vertices)
  Edges = disjoint pairs

At n=5: 10 merged iso classes AND 6 wiggly classes (tiles).
At n=6: 34 merged classes AND 10 wiggly classes.

Is there a Petersen structure in the n=6 wiggly classes?
10 wiggly classes = vertices of Petersen = K(5,2)?

The tiles at n=6 are: (x,y) for x-y >= 2. These are pairs of
NON-ADJACENT vertices in the base path. Can we map them to
2-subsets of {1,...,5} such that the Petersen structure emerges?

Also: the 5 complement-flip overlap pairs at n=6 — do they correspond
to the 5 "spokes" of the Petersen graph?
"""

import sys
import numpy as np
from math import factorial, comb
from itertools import permutations, combinations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  PETERSEN CONNECTIONS IN WIGGLY STRUCTURE")
print("  kind-pasteur-2026-03-24-S20fe")
print("=" * 80)

# ================================================================
# PART 1: Wiggly classes at n=6 as Kneser K(5,2)?
# ================================================================
print("\n  PART 1: 10 wiggly classes at n=6")

n = 6
TILES = []
for y in range(1, n-1):
    for x in range(n, y+1, -1):
        TILES.append((x, y))
labels = [chr(65+i) for i in range(len(TILES))]

print(f"  10 tiles = 10 wiggly classes:")
for i, (x, y) in enumerate(TILES):
    print(f"    {labels[i]}: ({x},{y}), skip={x-y-1}")

# The tiles are pairs (x,y) with x > y+1 from {1,...,n}
# At n=6: pairs from {1,...,6} with gap >= 2
# Can we find a bijection to 2-subsets of {1,...,5}?

# Kneser K(5,2): vertices = C(5,2) = 10, matching!
kneser_verts = list(combinations(range(1,6), 2))
print(f"\n  Kneser K(5,2) vertices: {kneser_verts}")

# Build Kneser adjacency: disjoint pairs
K_adj = np.zeros((10, 10), dtype=int)
for i in range(10):
    for j in range(i+1, 10):
        if not set(kneser_verts[i]) & set(kneser_verts[j]):
            K_adj[i][j] = K_adj[j][i] = 1

kneser_evals = sorted(np.linalg.eigvalsh(K_adj.astype(float)), reverse=True)
print(f"  Kneser eigenvalues: {[f'{x:.1f}' for x in kneser_evals]}")
print(f"  Kneser edges: {K_adj.sum()//2}")
print(f"  Kneser degrees: {K_adj.sum(axis=1)[0]} (3-regular)")

# ================================================================
# PART 2: Build wiggly-class interaction graph
# ================================================================
# Two wiggly classes X and Y "interact" if flipping X then Y
# gives a different result than flipping Y then X on some tiling.
# More precisely: define the COMMUTATOR graph where X-Y is an edge
# if there exists a tiling where [flip X, flip Y] changes the class.

# OR: Two wiggly classes are "similar" if they produce similar
# edge sets in the metagraph. Use Jaccard similarity.

print(f"\n  PART 2: Wiggly class interaction graph at n=6")

N = n
VERTS = list(range(n, 0, -1))
all_perms = list(permutations(range(N)))
m = len(TILES)

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

canon_map = {}; adj_map = {}
for mask in range(1 << m):
    bits = [(mask >> k) & 1 for k in range(m)]
    A = bits_to_adj(bits)
    adj_map[mask] = A; canon_map[mask] = canonicalize(A)

comp_cn = {}
for cn in set(canon_map.values()):
    for mask, c in canon_map.items():
        if c == cn:
            comp_cn[cn] = canonicalize(adj_complement(adj_map[mask]))
            break

merged_map = {cn: min(cn, comp_cn.get(cn, cn)) for cn in set(canon_map.values())}

# Per wiggly class: which MERGED metagraph edges does it generate?
class_edges = [set() for _ in range(m)]
for mask in range(1 << m):
    mcn1 = merged_map[canon_map[mask]]
    for wi in range(m):
        mcn2 = merged_map[canon_map[mask ^ (1 << wi)]]
        if mcn1 != mcn2:
            class_edges[wi].add((min(mcn1,mcn2), max(mcn1,mcn2)))

# Jaccard similarity between wiggly classes
print(f"\n  JACCARD SIMILARITY between wiggly classes:")
print(f"  {'':>3}", end='')
for j in range(m):
    print(f" {labels[j]:>5}", end='')
print()

jaccard = np.zeros((m, m))
for i in range(m):
    print(f"  {labels[i]:>3}", end='')
    for j in range(m):
        if i == j:
            jaccard[i][j] = 1.0
            print(f" {'1.00':>5}", end='')
        else:
            common = len(class_edges[i] & class_edges[j])
            union = len(class_edges[i] | class_edges[j])
            jaccard[i][j] = common / union if union > 0 else 0
            print(f" {jaccard[i][j]:5.2f}", end='')
    print()

# Threshold the Jaccard to make a graph
threshold = 0.8
print(f"\n  WIGGLY CLASS GRAPH (Jaccard > {threshold}):")
wc_adj = (jaccard > threshold).astype(int)
np.fill_diagonal(wc_adj, 0)
wc_edges = wc_adj.sum() // 2
wc_degrees = wc_adj.sum(axis=1)
print(f"  Edges: {wc_edges}")
print(f"  Degrees: {list(wc_degrees)}")

# Compare with Kneser
print(f"\n  COMPARISON WITH KNESER K(5,2):")
print(f"  Kneser: 15 edges, 3-regular")
print(f"  Wiggly class (J>{threshold}): {wc_edges} edges, degrees {list(wc_degrees)}")

# ================================================================
# PART 3: Tile pairs as vertex subsets
# ================================================================
# The tiles at n=6 are: arcs (x,y) with x-y >= 2 in {1,...,6}
# The base path is 6-5-4-3-2-1 (arcs with x-y=1)
# Each tile corresponds to a "chord" of the base path.
#
# Map each tile to the SET of base-path vertices it "skips over":
# Tile (x,y) skips vertices {y+1, y+2, ..., x-1}
# e.g., tile (6,1) skips {2,3,4,5} (4 vertices)
#       tile (4,2) skips {3} (1 vertex)
#       tile (6,4) skips {5} (1 vertex)

print(f"\n  PART 3: Tiles as chord-skip sets")
for i, (x, y) in enumerate(TILES):
    skip_set = set(range(y+1, x))
    print(f"    {labels[i]}: ({x},{y}), skips = {skip_set}")

# Two tiles are "disjoint" if their skip sets don't overlap
print(f"\n  DISJOINT skip sets (= independent chords):")
disjoint_adj = np.zeros((m, m), dtype=int)
for i in range(m):
    for j in range(i+1, m):
        skip_i = set(range(TILES[i][1]+1, TILES[i][0]))
        skip_j = set(range(TILES[j][1]+1, TILES[j][0]))
        if not skip_i & skip_j:
            disjoint_adj[i][j] = disjoint_adj[j][i] = 1
            print(f"    {labels[i]} ({TILES[i]}) -- {labels[j]} ({TILES[j]})")

disjoint_evals = sorted(np.linalg.eigvalsh(disjoint_adj.astype(float)), reverse=True)
print(f"\n  Disjoint graph eigenvalues: {[f'{x:.2f}' for x in disjoint_evals]}")
print(f"  Disjoint graph edges: {disjoint_adj.sum()//2}")
print(f"  Disjoint graph degrees: {list(disjoint_adj.sum(axis=1))}")

# Is this the Petersen graph?
d_match = (disjoint_adj.sum()//2 == 15 and all(d == 3 for d in disjoint_adj.sum(axis=1)))
print(f"\n  Disjoint-chord graph = Petersen? {d_match}")

# ================================================================
# PART 4: The 5 complement-flip overlap pairs as Petersen spokes
# ================================================================
print(f"\n  PART 4: The 5 overlap pairs")

# Find overlap pairs
wiggly_edges = set()
for mask in range(1 << m):
    cn1 = canon_map[mask]
    for wi in range(m):
        cn2 = canon_map[mask ^ (1 << wi)]
        if cn1 != cn2:
            wiggly_edges.add((min(cn1,cn2), max(cn1,cn2)))

comp_pairs = set()
for cn in set(canon_map.values()):
    cc = comp_cn.get(cn, cn)
    if cc != cn:
        comp_pairs.add((min(cn, cc), max(cn, cc)))

overlap = comp_pairs & wiggly_edges
print(f"  5 overlap pairs (complement = wiggly neighbor)")
print(f"  22 total complement pairs at n=6")
print(f"  The 5 are 22.7% of complement pairs")
print(f"  Petersen has 5 spokes connecting outer pentagon to inner pentagram")
print(f"  Coincidence of 5? Or deeper?")

# At n=5 the inner/outer structure of the Petersen comes from
# 2-subsets of {1,...,5}: outer = adjacent pairs in cycle, inner = non-adjacent
# At n=6: the 5 overlap pairs might partition the 22 complement pairs similarly

print("\nDONE.")
print("=" * 80)
