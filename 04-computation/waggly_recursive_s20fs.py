#!/usr/bin/env python3
"""
waggly_recursive_s20fs.py — Recursive patterns in waggly letter pairs
kind-pasteur-2026-03-24-S20fs

At each n: m = C(n-1,2) tiles decompose as:
  OLD tiles: C(n-2,2) tiles NOT involving vertex n (= the n-1 staircase)
  NEW tiles: n-2 tiles involving vertex n: (n,1), (n,2), ..., (n,n-2)

For d=2 pairs: OLD-OLD, OLD-NEW, NEW-NEW.
OLD-OLD pairs at n are exactly the d=2 pairs from n-1.

Question: does the PRODUCTIVITY of OLD-OLD pairs at n match their
productivity at n-1? Does the recursion preserve the pair structure?

Compute at n=5,6,7 and compare across levels.
"""

import sys
from math import factorial, comb
from itertools import permutations, combinations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  WAGGLY RECURSIVE PATTERNS")
print("  kind-pasteur-2026-03-24-S20fs")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

# Store results across n for comparison
all_results = {}

for n in [5, 6, 7]:
    t0 = time.time()
    TILES = get_tiles(n)
    m = len(TILES)
    N = n
    VERTS = list(range(n, 0, -1))
    labels = [chr(65+i) if i < 26 else chr(65+i-26)+str(i//26) for i in range(m)]
    all_perms = list(permutations(range(N)))

    tv = [(VERTS.index(x), VERTS.index(y)) for x, y in TILES]

    # Classify tiles: OLD (not involving vertex n) vs NEW (involving vertex n)
    # In VERTS labeling, vertex n has VERTS label n, index VERTS.index(n) = 0 (since VERTS = [n, n-1, ..., 1])
    # NEW tiles: those with x = n (TILES[i][0] == n)
    new_tiles = [i for i in range(m) if TILES[i][0] == n]
    old_tiles = [i for i in range(m) if TILES[i][0] != n]

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m} tiles")
    print(f"  OLD tiles (from n-1): {len(old_tiles)} = C({n-2},2) = {comb(n-2,2)}")
    print(f"  NEW tiles (vertex {n}): {len(new_tiles)} = {n-2}")
    for i in new_tiles:
        print(f"    {labels[i]}: ({TILES[i][0]},{TILES[i][1]})")
    print(f"{'#'*60}")

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

    print(f"  Building canon map...", end=" ", flush=True)
    t_c = time.time()
    canon_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        canon_map[mask] = canonicalize(A)
    print(f"done ({time.time()-t_c:.0f}s)")

    classes = sorted(set(canon_map.values()))
    V = len(classes)

    # d=1 edges
    d1_edges = set()
    for mask in range(1 << m):
        cn1 = canon_map[mask]
        for wi in range(m):
            cn2 = canon_map[mask ^ (1 << wi)]
            if cn1 != cn2:
                d1_edges.add((min(cn1,cn2), max(cn1,cn2)))

    # d=2 pairs classified
    oo_pairs = [(i,j) for i,j in combinations(range(m), 2) if i in old_tiles and j in old_tiles]
    on_pairs = [(i,j) for i,j in combinations(range(m), 2) if (i in old_tiles) != (j in old_tiles)]
    nn_pairs = [(i,j) for i,j in combinations(range(m), 2) if i in new_tiles and j in new_tiles]

    print(f"  d=2 pairs: {len(oo_pairs)} OLD-OLD, {len(on_pairs)} OLD-NEW, {len(nn_pairs)} NEW-NEW")

    def compute_pair_stats(pairs, label):
        total_new = 0
        total_edges = 0
        total_sl = 0
        pair_new_counts = []

        for wi, wj in pairs:
            edges = set()
            sl = 0
            for mask in range(1 << m):
                cn1 = canon_map[mask]
                cn2 = canon_map[mask ^ (1 << wi) ^ (1 << wj)]
                if cn1 == cn2: sl += 1
                else: edges.add((min(cn1,cn2), max(cn1,cn2)))
            new = len(edges - d1_edges)
            total_new += new
            total_edges += len(edges)
            total_sl += sl
            pair_new_counts.append(new)

        if pairs:
            avg_new = total_new / len(pairs)
            min_new = min(pair_new_counts)
            max_new = max(pair_new_counts)
        else:
            avg_new = min_new = max_new = 0

        print(f"  {label}: {len(pairs)} pairs, avg_new={avg_new:.1f}, range=[{min_new},{max_new}], total_new={total_new}, avg_SL={total_sl/max(len(pairs),1):.1f}")
        return pair_new_counts

    oo_counts = compute_pair_stats(oo_pairs, "OLD-OLD")
    on_counts = compute_pair_stats(on_pairs, "OLD-NEW")
    nn_counts = compute_pair_stats(nn_pairs, "NEW-NEW")

    # Store for cross-n comparison
    all_results[n] = {
        'oo': oo_counts, 'on': on_counts, 'nn': nn_counts,
        'oo_pairs': oo_pairs, 'on_pairs': on_pairs, 'nn_pairs': nn_pairs,
        'd1_edges': len(d1_edges), 'V': V
    }

    # Top 5 most productive pairs in each category
    for category, pairs, counts in [("OLD-OLD", oo_pairs, oo_counts),
                                      ("OLD-NEW", on_pairs, on_counts),
                                      ("NEW-NEW", nn_pairs, nn_counts)]:
        if not pairs: continue
        sorted_pairs = sorted(zip(counts, pairs), reverse=True)
        print(f"\n  Top 5 {category}:")
        for new, (wi, wj) in sorted_pairs[:5]:
            t1, t2 = TILES[wi], TILES[wj]
            print(f"    {labels[wi]}{labels[wj]} = ({t1[0]},{t1[1]})+({t2[0]},{t2[1]}): new={new}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.0f}s")

# Cross-n comparison
print(f"\n{'='*60}")
print(f"  CROSS-n COMPARISON")
print(f"{'='*60}")

for n in [6, 7]:
    if n in all_results and n-1 in all_results:
        oo_n = all_results[n]['oo']
        all_prev = all_results[n-1]['oo'] + all_results[n-1]['on'] + all_results[n-1]['nn']

        # OLD-OLD at n should correspond to ALL pairs at n-1
        # (since old tiles at n = all tiles at n-1)
        print(f"\n  n={n}: OLD-OLD pairs = {len(oo_n)}, ALL pairs at n-1 = {len(all_prev)}")
        print(f"  Count match? {len(oo_n) == len(all_prev)}")

        if len(oo_n) == len(all_prev):
            # Compare productivity
            avg_oo = sum(oo_n)/len(oo_n) if oo_n else 0
            avg_prev = sum(all_prev)/len(all_prev) if all_prev else 0
            print(f"  Avg new edges: OLD-OLD@n={avg_oo:.1f}, ALL@n-1={avg_prev:.1f}")
            print(f"  Ratio: {avg_oo/avg_prev:.2f}" if avg_prev > 0 else "")

        # How much do NEW tiles contribute?
        on_n = all_results[n]['on']
        nn_n = all_results[n]['nn']
        total_new_from_new = sum(on_n) + sum(nn_n)
        total_new_from_old = sum(oo_n)
        total_all = total_new_from_new + total_new_from_old
        print(f"  New edges from OLD-OLD: {total_new_from_old} ({total_new_from_old/total_all*100:.1f}%)")
        print(f"  New edges from OLD-NEW: {sum(on_n)} ({sum(on_n)/total_all*100:.1f}%)")
        print(f"  New edges from NEW-NEW: {sum(nn_n)} ({sum(nn_n)/total_all*100:.1f}%)")

print("\nDONE.")
print("=" * 80)
