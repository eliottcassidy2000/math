#!/usr/bin/env python3
"""
creative_traversal_s20fx.py — Creative traversal and compression of tournament space
kind-pasteur-2026-03-24-S20fx

5 creative techniques for traversing/compressing Q_m:

1. WIGGLY SPIRAL: cycle through wiggly classes A,B,C,D,E,F,A,B,C,...
   Does this "spiral" have better H-properties than random Gray codes?

2. DESCENT MINIMIZATION: find the Gray code with fewest H-descents.
   This is the "minimal backtracking" traversal.

3. SCORE-BLOCK ORDERING: group tilings by score sequence, traverse
   within each block, then jump between blocks.

4. ENTROPY COMPRESSION: how many bits are needed per iso class?
   Shannon entropy of pi gives the theoretical minimum.

5. TWO-COORDINATE SYSTEM: use (H, Hamming_weight) as 2D coordinates.
   How much of the metagraph does this capture?
"""

import sys
import numpy as np
from math import factorial, comb, log2
from itertools import permutations
from collections import defaultdict, Counter
import time
import random

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  CREATIVE TRAVERSAL AND COMPRESSION")
print("  kind-pasteur-2026-03-24-S20fx")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [4, 5]:
    t0 = time.time()
    TILES = get_tiles(n)
    m = len(TILES)
    N = n
    VERTS = list(range(n, 0, -1))
    all_perms = list(permutations(range(N)))
    tv = [(VERTS.index(x), VERTS.index(y)) for x, y in TILES]

    def b2a(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(m):
            xi, yi = tv[i]
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

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

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    # Build all data
    H_map = {}
    canon_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        H_map[mask] = count_hp(A)
        canon_map[mask] = canonicalize(A)

    classes = sorted(set(canon_map.values()))
    V = len(classes)
    tilings = np.array([sum(1 for c in canon_map.values() if c == cn) for cn in classes], dtype=float)
    pi = tilings / tilings.sum()

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m}, 2^m = {2**m}, V = {V}")
    print(f"{'#'*60}")

    # ================================================================
    # 1. WIGGLY SPIRAL
    # ================================================================
    # Build a path by cycling through tiles: flip tile 0, then 1, then 2, ...
    # When a tile would revisit, skip to next
    print(f"\n  1. WIGGLY SPIRAL:")

    best_spiral_descents = m * (2**m)  # worst case
    for start in range(1 << m):
        path = [start]
        visited = {start}
        pos = start
        wi = 0
        while len(path) < (1 << m):
            next_pos = pos ^ (1 << (wi % m))
            if next_pos not in visited:
                visited.add(next_pos)
                path.append(next_pos)
                pos = next_pos
            wi += 1
            if wi > m * (2**m):  # safety
                break

        if len(path) == (1 << m):
            H_seq = [H_map[v] for v in path]
            descents = sum(1 for i in range(len(H_seq)-1) if H_seq[i] > H_seq[i+1])
            if descents < best_spiral_descents:
                best_spiral_descents = descents
                best_spiral_seq = H_seq

    if best_spiral_descents < m * (2**m):
        total_steps = (1 << m) - 1
        print(f"    Best spiral descents: {best_spiral_descents}/{total_steps} = {best_spiral_descents/total_steps*100:.1f}%")
        if n == 4:
            print(f"    H sequence: {best_spiral_seq}")
    else:
        print(f"    Spiral didn't complete (some tilings unreachable)")

    # ================================================================
    # 2. DESCENT MINIMIZATION (exact at n=4)
    # ================================================================
    if n == 4:
        print(f"\n  2. MINIMUM DESCENT GRAY CODE:")
        adj = defaultdict(list)
        for mask in range(1 << m):
            for wi in range(m):
                adj[mask].append(mask ^ (1 << wi))

        # Find ALL Hamiltonian paths and pick the one with min descents
        all_paths = []
        def dfs_all(path, visited):
            if len(path) == (1 << m):
                all_paths.append(list(path))
                return
            for nb in adj[path[-1]]:
                if nb not in visited:
                    visited.add(nb)
                    path.append(nb)
                    dfs_all(path, visited)
                    path.pop()
                    visited.remove(nb)

        for start in range(1 << m):
            dfs_all([start], {start})

        min_desc = (1 << m)
        min_path = None
        for path in all_paths:
            H_seq = [H_map[v] for v in path]
            desc = sum(1 for i in range(len(H_seq)-1) if H_seq[i] > H_seq[i+1])
            if desc < min_desc:
                min_desc = desc
                min_path = path

        print(f"    Minimum descents: {min_desc}/{(1<<m)-1}")
        if min_path:
            H_seq = [H_map[v] for v in min_path]
            print(f"    H sequence: {H_seq}")
            ascents = sum(1 for i in range(len(H_seq)-1) if H_seq[i] < H_seq[i+1])
            levels = sum(1 for i in range(len(H_seq)-1) if H_seq[i] == H_seq[i+1])
            print(f"    Ascents: {ascents}, Levels: {levels}, Descents: {min_desc}")

    # ================================================================
    # 3. ENTROPY COMPRESSION
    # ================================================================
    print(f"\n  3. ENTROPY COMPRESSION:")
    entropy = -sum(p * log2(p) for p in pi if p > 0)
    print(f"    Shannon entropy of pi: {entropy:.4f} bits")
    print(f"    log2(V) = {log2(V):.4f} bits (uniform coding)")
    print(f"    Compression ratio: {log2(V)/entropy:.3f}x")
    print(f"    Savings: {(1 - entropy/log2(V))*100:.1f}% fewer bits than uniform")

    # Huffman-like: assign shorter codes to common classes
    sorted_pi = sorted(enumerate(pi), key=lambda x: -x[1])
    print(f"    Top 3 classes: {[(i, f'{p:.3f}') for i, p in sorted_pi[:3]]}")
    print(f"    Bottom 3: {[(i, f'{p:.4f}') for i, p in sorted_pi[-3:]]}")

    # ================================================================
    # 4. TWO-COORDINATE SYSTEM (H, Hamming_weight)
    # ================================================================
    print(f"\n  4. TWO-COORDINATE SYSTEM (H, HammingWeight):")

    # For each tiling: compute (H, hw)
    coord_map = {}
    for mask in range(1 << m):
        hw = bin(mask).count('1')
        H = H_map[mask]
        cn = canon_map[mask]
        coord_map[mask] = (H, hw)

    # How many DISTINCT (H, hw) pairs?
    distinct_pairs = len(set(coord_map.values()))
    # How many iso classes does (H, hw) distinguish?
    # Group tilings by (H, hw) -> set of classes
    pair_to_classes = defaultdict(set)
    for mask, (H, hw) in coord_map.items():
        pair_to_classes[(H, hw)].add(canon_map[mask])

    multi_class_pairs = sum(1 for classes in pair_to_classes.values() if len(classes) > 1)

    print(f"    Distinct (H, hw) pairs: {distinct_pairs}")
    print(f"    Pairs that distinguish classes uniquely: {distinct_pairs - multi_class_pairs}/{distinct_pairs}")
    print(f"    Pairs with multiple classes: {multi_class_pairs}")
    print(f"    Total classes: {V}")
    print(f"    (H, hw) resolves {(distinct_pairs - multi_class_pairs)/V*100:.1f}% of classes uniquely")

    # ================================================================
    # 5. SCORE-BASED ORDERING
    # ================================================================
    print(f"\n  5. SCORE-BASED BLOCKS:")
    score_blocks = defaultdict(list)
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        sc = tuple(sorted(sum(A[i]) for i in range(N)))
        score_blocks[sc].append(mask)

    print(f"    Score sequences: {len(score_blocks)}")
    for sc, masks in sorted(score_blocks.items(), key=lambda x: len(x[1]), reverse=True)[:5]:
        H_vals_in_block = [H_map[mask] for mask in masks]
        H_range = f"{min(H_vals_in_block)}-{max(H_vals_in_block)}"
        classes_in = len(set(canon_map[mask] for mask in masks))
        print(f"    score={sc}: {len(masks)} tilings, H range={H_range}, {classes_in} classes")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
