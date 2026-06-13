#!/usr/bin/env python3
"""
complement_flip_overlap_s20ez.py — Classes where complement = single flip neighbor
kind-pasteur-2026-03-24-S20ez

At n=6: 5 pairs of classes are BOTH complement pairs AND single-flip neighbors.
At n=4,5: 0 such pairs. Investigate these at n=6,7 and characterize them.

For each overlap pair {C, C'} where C' = complement(C) AND C' is a single-flip neighbor:
  - What are H, |Aut|, scores?
  - Which tile(s) connect them? (which wiggly class?)
  - Are they SC or NS?
  - What is their position in the metagraph (degree, spine/rib/sea)?
  - Can we predict which classes have this property?
"""

import sys
from math import factorial
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  COMPLEMENT-FLIP OVERLAP: Special Cases")
print("  kind-pasteur-2026-03-24-S20ez")
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
    m = len(TILES)
    N = n
    VERTS = list(range(n, 0, -1))
    all_perms = list(permutations(range(N)))

    print(f"\n{'#'*70}")
    print(f"  n = {n}, m = {m} tiles")
    print(f"{'#'*70}")

    if n >= 8:
        print(f"  Skipping (too large for brute force)")
        continue

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

    def count_ham_paths(A):
        nn = len(A)
        dp = [[0]*nn for _ in range(1 << nn)]
        for v in range(nn): dp[1 << v][v] = 1
        for mask in range(1, 1 << nn):
            for v in range(nn):
                if not (mask & (1 << v)): continue
                if dp[mask][v] == 0: continue
                for u in range(nn):
                    if mask & (1 << u): continue
                    if A[v][u]: dp[mask | (1 << u)][u] += dp[mask][v]
        return sum(dp[(1 << nn) - 1])

    def score_seq(A):
        return tuple(sorted(sum(A[i]) for i in range(len(A))))

    # Build data
    canon_map = {}
    adj_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bits_to_adj(bits)
        adj_map[mask] = A
        canon_map[mask] = canonicalize(A)

    # Unmerged classes with properties
    class_data = {}
    for cn in set(canon_map.values()):
        for mask, c in canon_map.items():
            if c == cn and cn not in class_data:
                A = adj_map[mask]
                comp_cn = canonicalize(adj_complement(A))
                size = sum(1 for m2, c2 in canon_map.items() if c2 == cn)
                aut = factorial(N) // size
                H = count_ham_paths(A)
                scores = score_seq(A)
                class_data[cn] = {
                    'comp': comp_cn, 'aut': aut, 'H': H,
                    'scores': scores, 'size': size,
                    'sc': comp_cn == cn, 'mask': mask
                }

    # Find wiggly edges (unmerged)
    wiggly_edges = defaultdict(set)  # (cn1, cn2) -> set of tile indices
    for mask in range(1 << m):
        cn1 = canon_map[mask]
        for wi in range(m):
            cn2 = canon_map[mask ^ (1 << wi)]
            if cn1 != cn2:
                edge = (min(cn1, cn2), max(cn1, cn2))
                wiggly_edges[edge].add(wi)

    # Find complement pairs (unmerged, non-SC only)
    comp_pairs = set()
    for cn, data in class_data.items():
        if not data['sc']:
            pair = (min(cn, data['comp']), max(cn, data['comp']))
            comp_pairs.add(pair)

    # Find overlap: complement pairs that are also wiggly neighbors
    overlap = comp_pairs & set(wiggly_edges.keys())

    V = len(class_data)
    E_wiggly = len(wiggly_edges)
    SC = sum(1 for d in class_data.values() if d['sc'])

    print(f"  V = {V}, SC = {SC}, complement pairs = {len(comp_pairs)}")
    print(f"  Wiggly edges = {E_wiggly}")
    print(f"  OVERLAP (complement AND single-flip): {len(overlap)}")

    if overlap:
        print(f"\n  OVERLAP PAIRS IN DETAIL:")
        for pair in sorted(overlap):
            cn1, cn2 = pair
            d1 = class_data[cn1]
            d2 = class_data[cn2]
            tiles_connecting = wiggly_edges[pair]
            tile_names = [f"({TILES[wi][0]},{TILES[wi][1]})" for wi in sorted(tiles_connecting)]
            skips = [TILES[wi][0] - TILES[wi][1] - 1 for wi in tiles_connecting]

            print(f"\n    Pair: {cn1[:20]}... <-> {cn2[:20]}...")
            print(f"      H: {d1['H']} <-> {d2['H']}")
            print(f"      |Aut|: {d1['aut']} <-> {d2['aut']}")
            print(f"      Scores: {d1['scores']} <-> {d2['scores']}")
            print(f"      SC: {d1['sc']} <-> {d2['sc']}")
            print(f"      Class size: {d1['size']} <-> {d2['size']}")
            print(f"      Connecting tiles: {tile_names}")
            print(f"      Skips: {skips}")
            print(f"      # tiles connecting: {len(tiles_connecting)}")

            # Hamming distance between representative tilings
            mask1 = d1['mask']
            # Find a mask in class cn2
            for mask2, c in canon_map.items():
                if c == cn2:
                    bits1 = [(mask1 >> k) & 1 for k in range(m)]
                    bits2 = [(mask2 >> k) & 1 for k in range(m)]
                    hamming = sum(b1 != b2 for b1, b2 in zip(bits1, bits2))
                    # Also compute full adjacency Hamming
                    A1 = adj_map[mask1]
                    A2 = adj_map[mask2]
                    adj_hamming = sum(A1[i][j] != A2[i][j] for i in range(N) for j in range(N) if i != j)
                    print(f"      Tiling Hamming (rep): {hamming}/{m}")
                    print(f"      Adjacency Hamming (rep): {adj_hamming}/{N*(N-1)}")
                    break

        # Statistics on overlap pairs
        print(f"\n  OVERLAP STATISTICS:")
        overlap_H = [(class_data[p[0]]['H'], class_data[p[1]]['H']) for p in overlap]
        print(f"    H values: {overlap_H}")
        overlap_scores = [(class_data[p[0]]['scores'], class_data[p[1]]['scores']) for p in overlap]
        print(f"    Score pairs: {overlap_scores}")
        overlap_tiles = [len(wiggly_edges[p]) for p in overlap]
        print(f"    # tiles connecting per pair: {overlap_tiles}")

        # Are the overlapping classes clustered in any part of the metagraph?
        overlap_nodes = set()
        for p in overlap:
            overlap_nodes.add(p[0])
            overlap_nodes.add(p[1])
        overlap_H_vals = [class_data[cn]['H'] for cn in overlap_nodes]
        all_H_vals = [d['H'] for d in class_data.values()]
        print(f"\n    Overlap nodes H range: {min(overlap_H_vals)}-{max(overlap_H_vals)}")
        print(f"    All nodes H range: {min(all_H_vals)}-{max(all_H_vals)}")
        print(f"    Overlap nodes: {len(overlap_nodes)}/{V} ({len(overlap_nodes)/V*100:.1f}%)")

    else:
        print(f"  No overlap at n={n}.")

    # Growth prediction
    if len(overlap) > 0:
        frac = len(overlap) / len(comp_pairs)
        print(f"\n  Overlap/complement_pairs = {len(overlap)}/{len(comp_pairs)} = {frac:.4f}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
