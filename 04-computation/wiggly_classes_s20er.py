#!/usr/bin/env python3
"""
wiggly_classes_s20er.py — Wiggly class analysis (matches explorer exactly)
kind-pasteur-2026-03-23-S20er

A WIGGLY LINE = clicking one tile in the tournament-tiling-explorer.
Each tile position defines a WIGGLY CLASS labeled A, B, C, ...

For n=4: A=(4,1), B=(3,1), C=(4,2) — 3 classes
For n=5: A=(5,1), B=(4,1), C=(3,1), D=(5,2), E=(4,2), F=(5,3) — 6 classes
For n=6: A=(6,1)..J=(6,4) — 10 classes

The base path n->n-1->...->1 is FIXED (not flippable).
Tiling space = Q_{C(n-1,2)} (hypercube on tiles).

For each wiggly class X:
  - How many tilings have a self-loop when flipping tile X?
  - How many have an edge? To which classes?
  - What metagraph edges does class X generate?
  - Is class X "special" in any way (grid-symmetric, overlap vs boundary, etc.)?
"""

import sys
from math import factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  WIGGLY CLASS ANALYSIS")
print("  kind-pasteur-2026-03-23-S20er")
print("=" * 80)

def get_tiles(n):
    """Return tile list matching explorer's order."""
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in range(4, 7):
    t0 = time.time()
    TILES = get_tiles(n)
    num_tiles = len(TILES)
    labels = [chr(65+i) for i in range(num_tiles)]  # A, B, C, ...

    # Vertices: n, n-1, ..., 1 (explorer convention)
    VERTS = list(range(n, 0, -1))
    N = n

    print(f"\n{'#'*70}")
    print(f"  n = {n}, {num_tiles} tiles = {num_tiles} wiggly classes")
    print(f"  Base path: {' -> '.join(str(v) for v in VERTS)}")
    print(f"  Tiles: {', '.join(f'{labels[i]}=({TILES[i][0]},{TILES[i][1]})' for i in range(num_tiles))}")
    print(f"  Tiling space: Q_{num_tiles} = 2^{num_tiles} = {2**num_tiles} tilings")
    print(f"{'#'*70}")

    # Build adjacency from tiling bits (matching explorer exactly)
    def bits_to_adj(bits):
        A = [[0]*N for _ in range(N)]
        # Base path: k -> k+1 (in 0-indexed adjacency)
        for k in range(N-1):
            A[k][k+1] = 1
        # Tiles
        for i in range(num_tiles):
            x_label, y_label = TILES[i]
            xi = VERTS.index(x_label)
            yi = VERTS.index(y_label)
            if bits[i] == 0:
                A[xi][yi] = 1  # x dominates y (forward)
            else:
                A[yi][xi] = 1  # y dominates x (backward)
        return A

    def adj_to_str(A):
        return ''.join(str(A[i][j]) for i in range(N) for j in range(N))

    all_perms = list(permutations(range(N)))

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best:
                best = s
        return best

    def is_grid_sym(bits):
        """Check grid symmetry: (x,y) -> (n-y+1, n-x+1)."""
        trans_map = []
        tile_idx = {(t[0],t[1]): i for i, t in enumerate(TILES)}
        for i, (x, y) in enumerate(TILES):
            tx, ty = N - y + 1, N - x + 1
            ti = tile_idx.get((tx, ty), -1)
            trans_map.append(ti)
        for i in range(num_tiles):
            ti = trans_map[i]
            if ti != i and bits[i] != bits[ti]:
                return False
        return True

    # Enumerate all tilings
    tilings = []
    canon_map = {}
    for mask in range(1 << num_tiles):
        bits = [(mask >> k) & 1 for k in range(num_tiles)]
        A = bits_to_adj(bits)
        cn = canonicalize(A)
        canon_map[mask] = cn
        tilings.append({'mask': mask, 'bits': bits, 'canon': cn, 'grid_sym': is_grid_sym(bits)})

    # Merged classes (complement = all bits inverted)
    comp_map = {}
    for cn in set(canon_map.values()):
        # Find a mask with this canon, complement it
        for mask, c in canon_map.items():
            if c == cn:
                comp_mask = mask ^ ((1 << num_tiles) - 1)
                comp_cn = canon_map[comp_mask]
                comp_map[cn] = comp_cn
                break

    merged_map = {}
    for cn in set(canon_map.values()):
        merged_map[cn] = min(cn, comp_map.get(cn, cn))

    merged_classes = sorted(set(merged_map.values()))
    V_merged = len(merged_classes)

    # Per wiggly class analysis
    print(f"\n  {'Class':>5} {'Tile':>7} {'SL':>6} {'Edge':>6} {'SL%':>7} {'#Edges':>7} {'Blue':>5} {'Blk':>5}")

    all_metagraph_edges = set()
    class_data = []

    for wi in range(num_tiles):
        label = labels[wi]
        tile = TILES[wi]
        sl_count = 0
        edge_count = 0
        edges_from_wi = set()
        blue_sl = 0
        black_sl = 0
        blue_edge = 0
        black_edge = 0

        for mask in range(1 << num_tiles):
            cn1 = canon_map[mask]
            mcn1 = merged_map[cn1]

            flip_mask = mask ^ (1 << wi)
            cn2 = canon_map[flip_mask]
            mcn2 = merged_map[cn2]

            bits = [(mask >> k) & 1 for k in range(num_tiles)]
            gs = is_grid_sym(bits)

            if mcn1 == mcn2:
                sl_count += 1
                if gs: blue_sl += 1
                else: black_sl += 1
            else:
                edge_count += 1
                edge = (min(mcn1, mcn2), max(mcn1, mcn2))
                edges_from_wi.add(edge)
                all_metagraph_edges.add(edge)
                if gs: blue_edge += 1
                else: black_edge += 1

        total = sl_count + edge_count
        sl_pct = sl_count / total * 100 if total > 0 else 0

        print(f"  {label:>5} ({tile[0]},{tile[1]}) {sl_count:6d} {edge_count:6d} {sl_pct:6.1f}% {len(edges_from_wi):7d} {blue_sl+blue_edge:5d} {black_sl+black_edge:5d}")

        class_data.append({
            'label': label, 'tile': tile, 'sl': sl_count, 'edge': edge_count,
            'sl_pct': sl_pct, 'n_edges': len(edges_from_wi),
            'blue': blue_sl + blue_edge, 'black': black_sl + black_edge,
            'edges': edges_from_wi
        })

    E_total = len(all_metagraph_edges)
    print(f"\n  Total metagraph edges: {E_total}")

    # Check: are all wiggly classes IDENTICAL?
    sl_counts = [d['sl'] for d in class_data]
    edge_counts = [d['edge'] for d in class_data]
    n_edges_counts = [d['n_edges'] for d in class_data]

    all_same_sl = len(set(sl_counts)) == 1
    all_same_edges = len(set(n_edges_counts)) == 1
    print(f"\n  All classes have same SL count? {all_same_sl}")
    print(f"  All classes reach same # of metagraph edges? {all_same_edges}")

    if not all_same_sl:
        print(f"    SL counts: {sorted(set(sl_counts))}")
    if not all_same_edges:
        print(f"    Edge counts: {sorted(set(n_edges_counts))}")

    # Which metagraph edges are reached by ALL wiggly classes?
    common_edges = class_data[0]['edges']
    for d in class_data[1:]:
        common_edges = common_edges & d['edges']

    # Which are reached by SOME but not all?
    any_edges = set()
    for d in class_data:
        any_edges |= d['edges']

    print(f"\n  Metagraph edges reached by ALL wiggly classes: {len(common_edges)}/{E_total}")
    print(f"  Metagraph edges reached by SOME wiggly class: {len(any_edges)}/{E_total}")

    # Tile position geometry: distance from base path
    print(f"\n  TILE GEOMETRY:")
    for i in range(num_tiles):
        x, y = TILES[i]
        skip = x - y - 1  # how many vertices are skipped between x and y in the base path
        row = y  # which row of the staircase
        print(f"    {labels[i]}: ({x},{y}), skip={skip}, row={row}, SL%={class_data[i]['sl_pct']:.1f}%")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
