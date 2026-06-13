#!/usr/bin/env python3
"""
cycle_space_bijection_s20ge.py — The tiling-to-even-graph bijection via cycle space
kind-pasteur-2026-03-25-S20ge

THE CORRECT MAP: tiling -> even graph via FUNDAMENTAL CYCLES.

The base path = spanning tree T of K_n.
Each non-base edge e creates a FUNDAMENTAL CYCLE: the path in T from
one endpoint to the other, plus the edge e.

A TILING = a binary vector on the non-base edges.
Tile = 1 means "include this fundamental cycle."

The EVEN GRAPH = XOR (symmetric difference) of all included fundamental cycles.
Since each cycle has all even degrees, any XOR of cycles also has all even degrees.
This maps 2^{C(n-1,2)} tilings to 2^{C(n-1,2)} even graphs BIJECTIVELY.

This is the GF(2) isomorphism between:
  Tiling space = {0,1}^{C(n-1,2)}
  Cycle space of K_n = {edge sets with all even degrees}

The ISO CLASS version: two tilings map to isomorphic even graphs iff...
the tournaments are related by the pullback of the graph isomorphism?
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  CYCLE SPACE BIJECTION: Tiling <-> Even Graph")
print("  kind-pasteur-2026-03-25-S20ge")
print("=" * 80)

for n in [4, 5, 6]:
    t0 = time.time()
    N = n
    edges = [(i,j) for i in range(N) for j in range(i+1,N)]
    m_full = len(edges)  # C(n,2)
    base_edges = [(i,i+1) for i in range(N-1)]  # the base path
    non_base = [e for e in edges if e not in base_edges]
    m = len(non_base)  # C(n-1,2)
    all_perms = list(permutations(range(N)))

    edge_idx = {e: i for i, e in enumerate(edges)}

    # Fundamental cycle for each non-base edge (i,j) with i < j:
    # Path in base tree: i -> i+1 -> ... -> j, plus edge (i,j)
    # Edge set = {(i,i+1), (i+1,i+2), ..., (j-1,j)} ∪ {(i,j)}
    fund_cycles = {}
    for e in non_base:
        i, j = e
        cycle_edges = set()
        for k in range(i, j):  # base path from i to j
            cycle_edges.add((k, k+1))
        cycle_edges.add((i, j))  # the non-base edge itself
        # Convert to bit vector
        cycle_bits = 0
        for ce in cycle_edges:
            cycle_bits |= (1 << edge_idx[ce])
        fund_cycles[e] = cycle_bits

    def tiling_to_even_graph(tile_bits):
        """Map tiling (m bits on non-base edges) to even graph (m_full bit vector)."""
        result = 0
        for k, e in enumerate(non_base):
            if tile_bits & (1 << k):
                result ^= fund_cycles[e]  # XOR with fundamental cycle
        return result

    def is_even_degree(graph_bits):
        deg = [0] * N
        for k, (i,j) in enumerate(edges):
            if graph_bits & (1 << k):
                deg[i] += 1
                deg[j] += 1
        return all(d % 2 == 0 for d in deg)

    def graph_canon(graph_bits):
        best = None
        for p in all_perms:
            nb = 0
            for k, (i,j) in enumerate(edges):
                pi, pj = min(p[i],p[j]), max(p[i],p[j])
                nk = edge_idx[(pi, pj)]
                if graph_bits & (1 << k):
                    nb |= (1 << nk)
            if best is None or nb < best:
                best = nb
        return best

    def tournament_canon(tile_bits):
        """Canon of the tournament encoded by this tiling."""
        # Build adjacency
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1  # base path
        for k, (i,j) in enumerate(non_base):
            xi, yi = i, j  # already in order
            if tile_bits & (1 << k):
                A[j][i] = 1  # backward = tile is 1
            else:
                A[i][j] = 1  # forward = tile is 0
        # But we also need to set A[i+1][i] = 0 for base path
        for k in range(N-1):
            A[k+1][k] = 0  # base path is always forward

        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    print(f"\n{'#'*60}")
    print(f"  n = {n}, tiles = {m}, full edges = {m_full}")
    print(f"{'#'*60}")

    # Map all tilings to even graphs
    even_images = set()
    even_image_classes = set()
    tourn_classes = set()
    tourn_to_even = defaultdict(set)  # tournament class -> even graph classes

    all_even = 0
    for tile_bits in range(1 << m):
        eg = tiling_to_even_graph(tile_bits)
        assert is_even_degree(eg), f"NOT EVEN at tile_bits={tile_bits}!"
        even_images.add(eg)
        eg_cn = graph_canon(eg)
        even_image_classes.add(eg_cn)
        t_cn = tournament_canon(tile_bits)
        tourn_classes.add(t_cn)
        tourn_to_even[t_cn].add(eg_cn)
        all_even += 1

    # Build all even graphs independently
    all_even_graphs = set()
    for bits in range(1 << m_full):
        if is_even_degree(bits):
            all_even_graphs.add(graph_canon(bits))

    print(f"\n  CYCLE SPACE MAP RESULTS:")
    print(f"    Tilings: {2**m}")
    print(f"    All produce even graphs: {all_even == 2**m}")
    print(f"    Distinct even graph images (labeled): {len(even_images)}")
    print(f"    Even graph iso classes in image: {len(even_image_classes)}")
    print(f"    All even graph iso classes: {len(all_even_graphs)}")
    print(f"    Image covers ALL even classes: {even_image_classes == all_even_graphs}")
    print(f"    Tournament iso classes: {len(tourn_classes)}")
    print(f"    Tournament classes = Even graph classes? {len(tourn_classes) == len(all_even_graphs)}")

    # Is the map INJECTIVE on labeled even graphs?
    print(f"    Map is injective (labeled): {len(even_images) == 2**m}")

    # The KEY: does each tournament class map to a SINGLE even graph class?
    multi_image = sum(1 for gcs in tourn_to_even.values() if len(gcs) > 1)
    print(f"\n  TOURNAMENT -> EVEN GRAPH CLASS MAP:")
    print(f"    1-to-1 (each tourn class -> one even class): {multi_image == 0}")
    print(f"    Multi-image tournament classes: {multi_image}/{len(tourn_classes)}")

    if multi_image > 0:
        for t_cn, gcs in sorted(tourn_to_even.items(), key=lambda x: -len(x[1])):
            if len(gcs) > 1:
                print(f"      Tournament class -> {len(gcs)} even graph classes")

    # Is the map SURJECTIVE on iso classes?
    missing = all_even_graphs - even_image_classes
    extra = even_image_classes - all_even_graphs
    print(f"    Missing even classes: {len(missing)}")
    print(f"    Extra classes: {len(extra)}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\n" + "=" * 60)
print("THE BIJECTION")
print("=" * 60)
print("""
The cycle space map: tiling -> XOR of fundamental cycles
is a LABELED BIJECTION: 2^m tilings <-> 2^m even graphs (1-to-1).

At the ISO CLASS level:
  #tournament classes = #even graph classes (Royle-Praeger)
  But the map is NOT necessarily 1-to-1 on classes!
  Different tilings of the same tournament can produce
  different even graph classes (because the base path choice
  breaks the graph symmetry differently from the tournament symmetry).

The "trivial vertex" insight: the base path pins n-1 edges.
These become part of every fundamental cycle that crosses them.
The vertex at the end of the path is the one with the most
constrained cycle structure — it's the "trivial" direction.
""")

print("DONE.")
print("=" * 80)
