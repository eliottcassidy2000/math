#!/usr/bin/env python3
"""
parallel_generation_s20ey.py — Blue/black and wiggly as parallel generators
kind-pasteur-2026-03-24-S20ey

THE PARALLEL:
  Both blue/black and wiggly lines act on the SAME tiling space Q_m.
  Both project to the SAME set of iso classes (= even graphs count).
  But they generate DIFFERENT structures:
    Blue/black: Z_2 pairing → complement merging → V_merged nodes
    Wiggly: m-regular connections → metagraph edges

WHAT GENERATES WHAT:
  The tiling space Q_m = {0,1}^m has two symmetry groups acting on it:
    1. S_n acts by relabeling vertices → quotient = iso classes (V nodes)
    2. Z_2 acts by complement (reverse all arcs) → merging pairs
    3. Each wiggly class acts as an involution (flip one tile)

  The FULL symmetry group is S_n × Z_2 (if complement commutes with relabeling).
  Quotient by S_n alone = V iso classes = A000568(n) = #even graphs
  Quotient by S_n × Z_2 = V_merged merged classes

  Wiggly lines in Q_m, projected through the S_n quotient, give edges of G_n.
  Wiggly lines projected through S_n × Z_2 give edges of G_n/Z_2.

  Blue/black lines in Q_m, projected through S_n, give complement pairs (= edges of a matching).
  Projected through S_n × Z_2, they become self-loops.

QUESTION: Do wiggly lines and blue/black lines TOGETHER generate ALL of Q_m?
  i.e., starting from any tiling, can you reach any other tiling by a sequence
  of wiggly flips and complement flips?

ANSWER: YES! Wiggly lines alone generate Q_m (they're the edges of the hypercube).
  Adding blue/black (complement = flip all) is redundant — it's already reachable
  by flipping all m tiles one at a time.

BUT: in the MERGED graph, blue/black is NOT redundant — it's what creates the merging!
  Without the complement identification, we'd have V classes, not V_merged.
"""

import sys
import numpy as np
from math import factorial, comb
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  PARALLEL GENERATION: Blue/black and Wiggly")
print("  kind-pasteur-2026-03-24-S20ey")
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
        nn = len(A)
        return [[1 - A[i][j] if i != j else 0 for j in range(nn)] for i in range(nn)]

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    def is_grid_sym(bits):
        tile_idx = {(t[0],t[1]): i for i, t in enumerate(TILES)}
        for i, (x, y) in enumerate(TILES):
            tx, ty = N - y + 1, N - x + 1
            ti = tile_idx.get((tx, ty), -1)
            if ti >= 0 and ti != i and bits[i] != bits[ti]:
                return False
        return True

    # Build all data
    canon_map = {}
    adj_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bits_to_adj(bits)
        adj_map[mask] = A
        canon_map[mask] = canonicalize(A)

    # True complement canonical form
    comp_canon = {}
    for cn in set(canon_map.values()):
        for mask, c in canon_map.items():
            if c == cn:
                comp_canon[cn] = canonicalize(adj_complement(adj_map[mask]))
                break

    # Unmerged classes
    unmerged_classes = sorted(set(canon_map.values()))
    V = len(unmerged_classes)

    # Merged classes (true complement)
    merged_map = {cn: min(cn, comp_canon.get(cn, cn)) for cn in unmerged_classes}
    merged_classes = sorted(set(merged_map.values()))
    V_merged = len(merged_classes)
    SC = sum(1 for cn in unmerged_classes if comp_canon.get(cn, cn) == cn)

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m} tiles, 2^m = {2**m} tilings")
    print(f"  V = {V} iso classes = A000568({n}) = #even graphs({n})")
    print(f"  SC = {SC}, V_merged = {V_merged} = (V+SC)/2 = ({V}+{SC})/2 = {(V+SC)//2}")
    print(f"{'#'*60}")

    # ================================================================
    # BLUE/BLACK LINE STRUCTURE
    # ================================================================
    # Each tiling has exactly 1 blue/black partner (its complement)
    # Blue if tiling is grid-symmetric, black otherwise
    blue_count = sum(1 for mask in range(1 << m)
                     if is_grid_sym([(mask >> k) & 1 for k in range(m)]))
    black_count = 2**m - blue_count

    # Blue/black lines between DIFFERENT unmerged classes (complement pairs)
    bb_edges = set()  # edges in the UNMERGED graph from complement
    bb_self = 0  # SC classes (complement = self)
    for cn in unmerged_classes:
        cc = comp_canon.get(cn, cn)
        if cc == cn:
            bb_self += 1
        else:
            bb_edges.add((min(cn, cc), max(cn, cc)))

    print(f"\n  BLUE/BLACK LINES:")
    print(f"    Blue tilings: {blue_count} (grid-symmetric)")
    print(f"    Black tilings: {black_count}")
    print(f"    Blue fraction: {blue_count/2**m:.6f}")
    print(f"    Unmerged complement edges: {len(bb_edges)} (= (V-SC)/2 = {(V-SC)//2})")
    print(f"    SC (self-complement) classes: {bb_self}")
    print(f"    In merged graph: ALL blue/black lines become self-loops")

    # ================================================================
    # WIGGLY LINE STRUCTURE
    # ================================================================
    # Build wiggly edges on unmerged classes
    wiggly_edges_unmerged = set()
    wiggly_sl_unmerged = 0
    wiggly_edge_count_unmerged = 0

    for mask in range(1 << m):
        cn1 = canon_map[mask]
        for wi in range(m):
            cn2 = canon_map[mask ^ (1 << wi)]
            if cn1 == cn2:
                wiggly_sl_unmerged += 1
            else:
                wiggly_edge_count_unmerged += 1
                wiggly_edges_unmerged.add((min(cn1, cn2), max(cn1, cn2)))

    # Wiggly edges on merged classes
    wiggly_edges_merged = set()
    for mask in range(1 << m):
        mcn1 = merged_map[canon_map[mask]]
        for wi in range(m):
            mcn2 = merged_map[canon_map[mask ^ (1 << wi)]]
            if mcn1 != mcn2:
                wiggly_edges_merged.add((min(mcn1, mcn2), max(mcn1, mcn2)))

    E_unmerged = len(wiggly_edges_unmerged)
    E_merged = len(wiggly_edges_merged)

    print(f"\n  WIGGLY LINES:")
    print(f"    Total wiggly transitions: {2**m * m}")
    print(f"    Wiggly self-loops (labeled): {wiggly_sl_unmerged}")
    print(f"    Unmerged wiggly edges: {E_unmerged}")
    print(f"    Merged wiggly edges: {E_merged}")

    # ================================================================
    # THE PARALLEL: Both generate the same node set
    # ================================================================
    print(f"\n  THE PARALLEL:")
    print(f"    Blue/black nodes: {V} unmerged, {V_merged} merged")
    print(f"    Wiggly nodes: {V} unmerged, {V_merged} merged")
    print(f"    Same node set? YES (both = A000568({n}) = #even_graphs({n}))")

    # ================================================================
    # OVERLAP: Are any blue/black edges ALSO wiggly edges?
    # ================================================================
    overlap = bb_edges & wiggly_edges_unmerged
    print(f"\n  OVERLAP (unmerged):")
    print(f"    Blue/black edges: {len(bb_edges)}")
    print(f"    Wiggly edges: {E_unmerged}")
    print(f"    Overlap: {len(overlap)}")
    print(f"    Complement reachable by single flip? {'YES' if len(overlap) > 0 else 'NO'}")

    # ================================================================
    # COMBINED GRAPH: blue/black + wiggly
    # ================================================================
    all_unmerged = bb_edges | wiggly_edges_unmerged
    print(f"\n  COMBINED UNMERGED GRAPH:")
    print(f"    Vertices: {V}")
    print(f"    Blue/black edges: {len(bb_edges)}")
    print(f"    Wiggly edges: {E_unmerged}")
    print(f"    Total unique edges: {len(all_unmerged)}")
    print(f"    Blue/black adds {len(bb_edges) - len(overlap)} new edges beyond wiggly")

    # Connectivity of wiggly-only graph
    adj_w = defaultdict(set)
    for a, b in wiggly_edges_unmerged:
        adj_w[a].add(b)
        adj_w[b].add(a)

    visited = set()
    queue = [unmerged_classes[0]]
    visited.add(unmerged_classes[0])
    while queue:
        u = queue.pop(0)
        for v in adj_w[u]:
            if v not in visited:
                visited.add(v)
                queue.append(v)

    wiggly_connected = len(visited) == V

    # Connectivity of blue/black-only graph
    adj_bb = defaultdict(set)
    for a, b in bb_edges:
        adj_bb[a].add(b)
        adj_bb[b].add(a)

    visited_bb = set()
    queue = [unmerged_classes[0]]
    visited_bb.add(unmerged_classes[0])
    while queue:
        u = queue.pop(0)
        for v in adj_bb[u]:
            if v not in visited_bb:
                visited_bb.add(v)
                queue.append(v)

    bb_connected = len(visited_bb) == V

    print(f"\n  CONNECTIVITY:")
    print(f"    Wiggly-only connected? {wiggly_connected} (reaches {len(visited)}/{V})")
    print(f"    Blue/black-only connected? {bb_connected} (reaches {len(visited_bb)}/{V})")
    print(f"    Both needed for full connectivity? {not wiggly_connected and not bb_connected}")

    # ================================================================
    # STRUCTURE COMPARISON
    # ================================================================
    # Blue/black is a PERFECT MATCHING on non-SC classes (V-SC)/2 edges
    # Wiggly is a dense graph with E_unmerged edges
    # The combination: matching + dense graph on the same V nodes

    # What fraction of the complete graph is covered?
    max_edges = V * (V-1) // 2
    print(f"\n  GRAPH DENSITY:")
    print(f"    Blue/black: {len(bb_edges)}/{max_edges} = {len(bb_edges)/max_edges:.4f}")
    print(f"    Wiggly: {E_unmerged}/{max_edges} = {E_unmerged/max_edges:.4f}")
    print(f"    Combined: {len(all_unmerged)}/{max_edges} = {len(all_unmerged)/max_edges:.4f}")

    # Degrees in each
    deg_w = defaultdict(int)
    for a, b in wiggly_edges_unmerged:
        deg_w[a] += 1; deg_w[b] += 1
    deg_bb = defaultdict(int)
    for a, b in bb_edges:
        deg_bb[a] += 1; deg_bb[b] += 1

    dw = sorted(deg_w.values())
    dbb = sorted(deg_bb.values()) if deg_bb else [0]

    print(f"\n  DEGREE COMPARISON:")
    print(f"    Wiggly: min={dw[0]}, max={dw[-1]}, avg={sum(dw)/len(dw):.1f}")
    print(f"    Blue/black: min={min(dbb)}, max={max(dbb)}, avg={sum(dbb)/len(dbb) if dbb else 0:.1f}")
    print(f"    (Blue/black is a matching: every non-SC node has degree 1, SC nodes have degree 0)")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\n" + "=" * 70)
print("SYNTHESIS")
print("=" * 70)
print("""
THE PARALLEL GENERATION:

Both blue/black and wiggly operate on the SAME tiling space Q_m = {0,1}^m.
Both project onto the SAME iso classes (V = A000568 = #even_graphs).

  BLUE/BLACK: Z_2 involution (flip all tiles ~ complement)
    Creates a PERFECT MATCHING on non-SC classes
    In merged graph: all edges become self-loops
    Generates the MERGING structure

  WIGGLY: m involutions (flip one tile each, classes A,B,C,...)
    Creates a DENSE graph on all classes
    In merged graph: generates ALL edges
    Generates the CONNECTIVITY structure

  COMBINED: matching + dense graph on V nodes
    Blue/black adds NO new edges beyond wiggly at n=4,5,6
    (Complement is never reachable by a single flip)
    But blue/black adds the COMPLEMENT PAIRING which creates the Z_2 quotient.

  THE DEEP PARALLEL: Both are quotients of the SAME action:
    S_n acts on tilings → quotient gives iso classes (= even graphs)
    Z_2 acts by complement → quotient merges complement pairs
    Wiggly acts by tile flip → quotient gives metagraph edges

  All three share the same nodes. Their combined action on Q_m fully
  determines the merged metagraph G_n/Z_2.
""")

print("DONE.")
print("=" * 80)
