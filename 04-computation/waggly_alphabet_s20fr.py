#!/usr/bin/env python3
"""
waggly_alphabet_s20fr.py — Waggly as letter combinations
kind-pasteur-2026-03-24-S20fr

d=1: letters A, B, C, ... (one per tile)
d=2: pairs AB, AC, AD, BC, BD, CD, ... (C(m,2) pairs)
d=3: triples ABC, ABD, ... (C(m,3) triples)

For each letter combination:
  - Which metagraph edges does it generate?
  - How many are NEW (not covered by lower d)?
  - Self-loop count
  - Geometric properties (skip sum, row pattern, adjacency)

The GRID SYMMETRY pairs letter combinations: {X,Y} <-> {grid(X), grid(Y)}.
"""

import sys
from math import factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  WAGGLY ALPHABET: Letter Combinations")
print("  kind-pasteur-2026-03-24-S20fr")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [5, 6]:
    t0 = time.time()
    TILES = get_tiles(n)
    m = len(TILES)
    N = n
    VERTS = list(range(n, 0, -1))
    labels = [chr(65+i) for i in range(m)]
    all_perms = list(permutations(range(N)))

    tv = [(VERTS.index(x), VERTS.index(y)) for x, y in TILES]

    # Grid symmetry
    tile_idx = {(t[0],t[1]): i for i, t in enumerate(TILES)}
    grid_map = []
    for i, (x, y) in enumerate(TILES):
        gx, gy = N+1-y, N+1-x
        grid_map.append(tile_idx.get((gx, gy), i))

    def b2a(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(m):
            xi, yi = tv[i]
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    canon_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        canon_map[mask] = canonicalize(A)

    classes = sorted(set(canon_map.values()))
    V = len(classes)

    # d=1 edges (wiggly)
    d1_edges = set()
    for mask in range(1 << m):
        cn1 = canon_map[mask]
        for wi in range(m):
            cn2 = canon_map[mask ^ (1 << wi)]
            if cn1 != cn2:
                d1_edges.add((min(cn1,cn2), max(cn1,cn2)))

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m} tiles ({','.join(labels)})")
    print(f"  V = {V}, d=1 edges = {len(d1_edges)}")
    print(f"{'#'*60}")

    # ================================================================
    # d=2 LETTER PAIRS
    # ================================================================
    print(f"\n  d=2: {len(list(combinations(range(m), 2)))} letter pairs")

    pair_data = []
    for combo in combinations(range(m), 2):
        wi, wj = combo
        pair_label = labels[wi] + labels[wj]
        grid_pair = tuple(sorted([grid_map[wi], grid_map[wj]]))
        grid_label = labels[grid_pair[0]] + labels[grid_pair[1]]

        # Tile geometry
        (x1, y1), (x2, y2) = TILES[wi], TILES[wj]
        skip1, skip2 = x1-y1-1, x2-y2-1
        share_vertex = len({x1,y1} & {x2,y2}) > 0
        same_row = y1 == y2
        skip_sum = skip1 + skip2

        # Which edges does this pair create?
        edges = set()
        sl_count = 0
        for mask in range(1 << m):
            cn1 = canon_map[mask]
            flip_mask = mask ^ (1 << wi) ^ (1 << wj)
            cn2 = canon_map[flip_mask]
            if cn1 == cn2:
                sl_count += 1
            else:
                edges.add((min(cn1,cn2), max(cn1,cn2)))

        new_edges = edges - d1_edges
        n_total = len(edges)
        n_new = len(new_edges)

        pair_data.append({
            'pair': pair_label, 'tiles': combo,
            'grid_pair': grid_label, 'share_vertex': share_vertex,
            'same_row': same_row, 'skip_sum': skip_sum,
            'total_edges': n_total, 'new_edges': n_new,
            'sl': sl_count, 'edges': edges, 'new': new_edges
        })

    # Sort by new edges (most productive first)
    pair_data.sort(key=lambda x: -x['new_edges'])

    print(f"\n  LETTER PAIRS (sorted by new edges beyond d=1):")
    print(f"  {'Pair':>5} {'Grid':>5} {'Share':>5} {'SameR':>5} {'SkipS':>5} {'Total':>6} {'New':>5} {'SL':>5}")
    for pd in pair_data:
        print(f"  {pd['pair']:>5} {pd['grid_pair']:>5} {'Y' if pd['share_vertex'] else 'N':>5} {'Y' if pd['same_row'] else 'N':>5} {pd['skip_sum']:>5} {pd['total_edges']:>6} {pd['new_edges']:>5} {pd['sl']:>5}")

    # Grid-symmetry pairing
    print(f"\n  GRID-SYMMETRY PAIRS (same edge set?):")
    for pd in pair_data:
        if pd['pair'] < pd['grid_pair']:
            partner = next((p for p in pair_data if p['pair'] == pd['grid_pair']), None)
            if partner:
                jaccard = len(pd['edges'] & partner['edges']) / len(pd['edges'] | partner['edges']) if pd['edges'] | partner['edges'] else 0
                print(f"    {pd['pair']} <-> {pd['grid_pair']}: J={jaccard:.3f}, new={pd['new_edges']}/{partner['new_edges']}")

    # Which geometric features predict productivity?
    sharing = [pd['new_edges'] for pd in pair_data if pd['share_vertex']]
    not_sharing = [pd['new_edges'] for pd in pair_data if not pd['share_vertex']]
    print(f"\n  PRODUCTIVITY BY GEOMETRY:")
    print(f"    Share vertex: avg new = {sum(sharing)/len(sharing):.1f}" if sharing else "    No sharing pairs")
    print(f"    No shared vertex: avg new = {sum(not_sharing)/len(not_sharing):.1f}" if not_sharing else "    All share")

    same_row_new = [pd['new_edges'] for pd in pair_data if pd['same_row']]
    diff_row_new = [pd['new_edges'] for pd in pair_data if not pd['same_row']]
    print(f"    Same row: avg new = {sum(same_row_new)/len(same_row_new):.1f}" if same_row_new else "")
    print(f"    Diff row: avg new = {sum(diff_row_new)/len(diff_row_new):.1f}" if diff_row_new else "")

    # ================================================================
    # d=3 LETTER TRIPLES (for completeness at n=5)
    # ================================================================
    if n == 5:
        d2_all_edges = set()
        for pd in pair_data:
            d2_all_edges |= pd['edges']

        remaining_after_d2 = set()
        # Full metagraph edges
        all_edges = set()
        for mask in range(1 << m):
            cn1 = canon_map[mask]
            for mask2 in range(mask+1, 1 << m):
                cn2 = canon_map[mask2]
                if cn1 != cn2:
                    all_edges.add((min(cn1,cn2), max(cn1,cn2)))

        # Actually: just check d=3 triples
        print(f"\n  d=3: {len(list(combinations(range(m), 3)))} letter triples")
        d12_edges = d1_edges | d2_all_edges

        triple_data = []
        for combo in combinations(range(m), 3):
            triple_label = ''.join(labels[i] for i in combo)
            edges = set()
            sl = 0
            for mask in range(1 << m):
                cn1 = canon_map[mask]
                flip_mask = mask
                for wi in combo:
                    flip_mask ^= (1 << wi)
                cn2 = canon_map[flip_mask]
                if cn1 == cn2:
                    sl += 1
                else:
                    edges.add((min(cn1,cn2), max(cn1,cn2)))
            new = edges - d12_edges
            triple_data.append({
                'triple': triple_label, 'total': len(edges),
                'new': len(new), 'sl': sl
            })

        triple_data.sort(key=lambda x: -x['new'])
        print(f"  TOP TRIPLES (most new edges):")
        for td in triple_data[:10]:
            print(f"    {td['triple']}: total={td['total']}, new={td['new']}, SL={td['sl']}")

        # Which triples are needed for the 6 edges that d<=2 misses?
        needed = set()
        for td in triple_data:
            if td['new'] > 0:
                needed.add(td['triple'])
        print(f"  Triples with new edges: {len(needed)}/{len(triple_data)}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
