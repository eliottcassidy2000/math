#!/usr/bin/env python3
"""
wiggly_skip_geometry_s20es.py — Why do wiggly classes have different SL rates?
kind-pasteur-2026-03-24-S20es

The "skip" = x - y - 1 of tile (x,y) measures how far apart the vertices are
in the base path. Tiles with skip=1 (adjacent non-base arcs) have LOWER SL%.
Tiles with higher skip have HIGHER SL%.

HYPOTHESIS: SL rate depends on skip, and possibly row (y) and column (x) independently.

Also: which METAGRAPH EDGES does each wiggly class generate?
Do classes with same skip generate the same edges?
"""

import sys
from math import factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  WIGGLY SKIP GEOMETRY")
print("  kind-pasteur-2026-03-24-S20es")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in range(4, 8):
    t0 = time.time()
    TILES = get_tiles(n)
    num_tiles = len(TILES)
    labels = [chr(65+i) if i < 26 else chr(65+i-26)+str(i//26) for i in range(num_tiles)]
    N = n
    VERTS = list(range(n, 0, -1))

    print(f"\n{'#'*70}")
    print(f"  n = {n}, {num_tiles} wiggly classes")
    print(f"{'#'*70}")

    all_perms = list(permutations(range(N)))

    def bits_to_adj(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(num_tiles):
            x, y = TILES[i]
            xi, yi = VERTS.index(x), VERTS.index(y)
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    # Compute canonical forms
    canon_map = {}
    for mask in range(1 << num_tiles):
        bits = [(mask >> k) & 1 for k in range(num_tiles)]
        A = bits_to_adj(bits)
        canon_map[mask] = canonicalize(A)

    # Merged classes
    comp_map = {}
    for cn in set(canon_map.values()):
        for mask, c in canon_map.items():
            if c == cn:
                comp_mask = mask ^ ((1 << num_tiles) - 1)
                comp_map[cn] = canon_map[comp_mask]
                break

    merged_map = {cn: min(cn, comp_map.get(cn, cn)) for cn in set(canon_map.values())}

    # Per wiggly class: self-loop count and edges generated
    results = []
    for wi in range(num_tiles):
        x, y = TILES[wi]
        skip = x - y - 1
        row = y
        col = N - x  # column in the staircase grid

        sl_count = 0
        edge_count = 0
        edges = set()

        for mask in range(1 << num_tiles):
            mcn1 = merged_map[canon_map[mask]]
            flip_mask = mask ^ (1 << wi)
            mcn2 = merged_map[canon_map[flip_mask]]

            if mcn1 == mcn2:
                sl_count += 1
            else:
                edge_count += 1
                edges.add((min(mcn1, mcn2), max(mcn1, mcn2)))

        total = sl_count + edge_count
        sl_pct = sl_count / total * 100

        results.append({
            'label': labels[wi], 'tile': (x,y), 'skip': skip,
            'row': row, 'col': col,
            'sl': sl_count, 'edge': edge_count, 'sl_pct': sl_pct,
            'n_edges': len(edges), 'edges': edges
        })

    # Group by skip
    skip_groups = defaultdict(list)
    for r in results:
        skip_groups[r['skip']].append(r)

    print(f"\n  BY SKIP (distance between vertices in base path - 1):")
    print(f"  {'Skip':>4} {'Count':>5} {'Avg SL%':>8} {'SL% range':>15} {'Avg #edges':>10}")
    for skip in sorted(skip_groups.keys()):
        group = skip_groups[skip]
        sl_pcts = [r['sl_pct'] for r in group]
        n_edges_list = [r['n_edges'] for r in group]
        print(f"  {skip:4d} {len(group):5d} {sum(sl_pcts)/len(sl_pcts):8.2f} {min(sl_pcts):6.2f}-{max(sl_pcts):6.2f} {sum(n_edges_list)/len(n_edges_list):10.1f}")

    # Group by row
    row_groups = defaultdict(list)
    for r in results:
        row_groups[r['row']].append(r)

    print(f"\n  BY ROW (y coordinate in staircase):")
    for row in sorted(row_groups.keys()):
        group = row_groups[row]
        sl_pcts = [r['sl_pct'] for r in group]
        print(f"  row={row}: {len(group)} tiles, avg SL%={sum(sl_pcts)/len(sl_pcts):.2f}, range={min(sl_pcts):.2f}-{max(sl_pcts):.2f}")

    # Group by diagonal (skip = x-y-1, equivalently x-y)
    # The anti-diagonal is x+y = const
    diag_groups = defaultdict(list)
    for r in results:
        x, y = r['tile']
        diag = x + y  # anti-diagonal
        diag_groups[diag].append(r)

    print(f"\n  BY ANTI-DIAGONAL (x+y = const):")
    for diag in sorted(diag_groups.keys()):
        group = diag_groups[diag]
        sl_pcts = [r['sl_pct'] for r in group]
        tiles = [(r['tile']) for r in group]
        print(f"  x+y={diag}: tiles={tiles}, avg SL%={sum(sl_pcts)/len(sl_pcts):.2f}")

    # Edge overlap between wiggly classes with same skip
    print(f"\n  EDGE OVERLAP BETWEEN WIGGLY CLASSES:")
    for skip in sorted(skip_groups.keys()):
        group = skip_groups[skip]
        if len(group) >= 2:
            # Jaccard similarity of edge sets
            for i in range(len(group)):
                for j in range(i+1, len(group)):
                    e1 = group[i]['edges']
                    e2 = group[j]['edges']
                    common = len(e1 & e2)
                    union = len(e1 | e2)
                    jaccard = common / union if union > 0 else 0
                    print(f"    skip={skip}: {group[i]['label']} vs {group[j]['label']}: "
                          f"|common|={common}, |union|={union}, Jaccard={jaccard:.3f}")

    # Per-class detail
    print(f"\n  ALL WIGGLY CLASSES:")
    print(f"  {'Label':>5} {'Tile':>7} {'Skip':>4} {'Row':>3} {'SL%':>7} {'#Edges':>7}")
    for r in sorted(results, key=lambda r: -r['sl_pct']):
        print(f"  {r['label']:>5} ({r['tile'][0]},{r['tile'][1]}) {r['skip']:4d} {r['row']:3d} {r['sl_pct']:6.2f}% {r['n_edges']:7d}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

    if elapsed > 120:
        print("  Skipping larger n")
        break

print("\nDONE.")
print("=" * 80)
