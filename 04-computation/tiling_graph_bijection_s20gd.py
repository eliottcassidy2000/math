#!/usr/bin/env python3
"""
tiling_graph_bijection_s20gd.py — The tiling-to-graph map
kind-pasteur-2026-03-25-S20gd

THE INSIGHT: A tournament tiling (C(n-1,2) bits on delta_{n-2})
naturally maps to a graph on n vertices:
  - Fix the base path as "transitive direction"
  - Tile = 0 (forward arc) -> edge PRESENT
  - Tile = 1 (backward arc) -> edge ABSENT
  - Base path arcs are ALWAYS present (n-1 edges fixed)

The resulting graph has C(n-1,2) - hw + (n-1) edges total
(where hw = Hamming weight of tiling = number of "flipped" arcs).

QUESTION 1: Does this map send tournaments to even graphs?
QUESTION 2: What's the structure of the image?
QUESTION 3: The "trivial vertex" — which vertex is redundant?
QUESTION 4: The staircase delta_{n-2} vs delta_{n-1} — what's in the extra row?

Also: Redei says every tournament has a HP.
So: every tiling has a base path where ALL arcs are "forward."
In the graph interpretation: the graph contains a Hamiltonian path.
Complement: the "missing edges" graph also has structure from the tiling.
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  TILING -> GRAPH MAP")
print("  kind-pasteur-2026-03-25-S20gd")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [4, 5, 6]:
    t0 = time.time()
    TILES = get_tiles(n)
    m = len(TILES)
    N = n
    VERTS = list(range(n, 0, -1))
    all_perms = list(permutations(range(N)))
    tv = [(VERTS.index(x), VERTS.index(y)) for x, y in TILES]
    edges = [(i,j) for i in range(N) for j in range(i+1,N)]
    n_edges = len(edges)

    def b2a(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(m):
            xi, yi = tv[i]
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    def tiling_to_graph(bits):
        """Map tiling to undirected graph: forward arcs = edges present."""
        A = b2a(bits)
        # Graph: edge (i,j) present if the arc goes from lower to higher index
        # in the base path ordering. Base path: 0->1->2->...->n-1.
        # "Forward" = i->j where i < j (in the path ordering).
        # So: edge (i,j) present if A[i][j] = 1 (i < j).
        G = [[0]*N for _ in range(N)]
        for i in range(N):
            for j in range(i+1, N):
                if A[i][j]:  # forward arc = edge present
                    G[i][j] = 1
                    G[j][i] = 1
        return G

    def graph_degrees(G):
        return tuple(sorted(sum(G[i]) for i in range(N)))

    def is_even_degree(G):
        return all(sum(G[i]) % 2 == 0 for i in range(N))

    def graph_canon(G):
        best = None
        for p in all_perms:
            s = tuple(G[p[i]][p[j]] for i in range(N) for j in range(i+1, N))
            if best is None or s < best: best = s
        return best

    def tournament_canon(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m} tiles")
    print(f"{'#'*60}")

    # Map each tiling to a graph and check properties
    even_count = 0
    odd_count = 0
    graph_classes = set()
    even_graph_classes = set()
    tournament_classes = set()

    # Track which tournament classes map to which graph classes
    tourn_to_graph = defaultdict(set)

    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        G = tiling_to_graph(bits)

        t_cn = tournament_canon(A)
        g_cn = graph_canon(G)
        tournament_classes.add(t_cn)
        graph_classes.add(g_cn)
        tourn_to_graph[t_cn].add(g_cn)

        if is_even_degree(G):
            even_count += 1
            even_graph_classes.add(g_cn)
        else:
            odd_count += 1

    print(f"\n  TILING -> GRAPH MAP:")
    print(f"    Total tilings: {2**m}")
    print(f"    Produces even-degree graphs: {even_count} ({even_count/2**m*100:.1f}%)")
    print(f"    Produces odd-degree graphs: {odd_count} ({odd_count/2**m*100:.1f}%)")
    print(f"    Distinct graph classes in image: {len(graph_classes)}")
    print(f"    Distinct even graph classes: {len(even_graph_classes)}")
    print(f"    Distinct tournament classes: {len(tournament_classes)}")

    # How many graph classes per tournament class?
    graphs_per_tourn = [len(gcs) for gcs in tourn_to_graph.values()]
    print(f"\n  GRAPH CLASSES PER TOURNAMENT CLASS:")
    print(f"    min={min(graphs_per_tourn)}, max={max(graphs_per_tourn)}, avg={sum(graphs_per_tourn)/len(graphs_per_tourn):.2f}")

    # The "extra row" interpretation: base path has n-1 edges.
    # These are ALWAYS present. The graph is a subgraph of K_n
    # that includes the base path plus some non-base edges.
    # The graph has at least n-1 edges (the path) and at most
    # n-1 + C(n-1,2) = C(n,2) edges (all present).

    # Edge count distribution
    edge_counts = Counter()
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        G = tiling_to_graph(bits)
        ec = sum(G[i][j] for i in range(N) for j in range(i+1, N))
        edge_counts[ec] += 1

    print(f"\n  EDGE COUNT DISTRIBUTION:")
    print(f"    Min edges: {min(edge_counts.keys())} (base path only)")
    print(f"    Max edges: {max(edge_counts.keys())} (K_n = all forward)")
    for ec in sorted(edge_counts.keys()):
        ev = sum(1 for mask in range(1 << m)
                 if is_even_degree(tiling_to_graph([(mask >> k) & 1 for k in range(m)]))
                 and sum(tiling_to_graph([(mask >> k) & 1 for k in range(m)])[i][j]
                         for i in range(N) for j in range(i+1, N)) == ec)
        print(f"    {ec:3d} edges: {edge_counts[ec]:5d} tilings ({ev:4d} even)")

    # The key question: does EVERY even graph appear in the image?
    # Build all even graphs directly
    all_even = set()
    for mask in range(1 << n_edges):
        G = [[0]*N for _ in range(N)]
        for k, (i,j) in enumerate(edges):
            if mask & (1 << k):
                G[i][j] = G[j][i] = 1
        if is_even_degree(G):
            all_even.add(graph_canon(G))

    print(f"\n  ALL even graph classes: {len(all_even)}")
    print(f"  Even classes in tiling image: {len(even_graph_classes)}")
    print(f"  Missing from image: {len(all_even - even_graph_classes)}")
    print(f"  Extra in image: {len(even_graph_classes - all_even)}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
