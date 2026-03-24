#!/usr/bin/env python3
"""
waggly_coverage_s20fq.py — Waggly coverage by Hamming distance d
kind-pasteur-2026-03-24-S20fq

WAGGLY lines at distance d: flip exactly d tiles simultaneously.
  d=1: wiggly (single tile flip)
  d=2: double flip
  ...
  d=m: blue/black (complement tiling, flip all tiles)

COMPLETENESS ORDER k*: minimum d such that d=1..d covers all metagraph edges.
THEOREM (opus-S299): k* = diameter of the metagraph.

Compute for n=5,6,7:
  - Coverage at each d (fraction of edges reachable)
  - Cumulative coverage (d<=k covers how many edges?)
  - New edges per distance level
  - The waggly weight matrix W_d for each d
  - Self-loop counts at each d
"""

import sys
from math import factorial, comb
from itertools import permutations, combinations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  WAGGLY COVERAGE BY HAMMING DISTANCE")
print("  kind-pasteur-2026-03-24-S20fq")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [5, 6, 7]:
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

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    # Build full canon map (expensive at n=7 but already proven feasible)
    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m}, 2^m = {2**m}")
    print(f"{'#'*60}")

    print(f"  Building canon map...", end=" ", flush=True)
    t_canon = time.time()
    canon_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        canon_map[mask] = canonicalize(A)
    print(f"done ({time.time()-t_canon:.0f}s)")

    classes = sorted(set(canon_map.values()))
    V = len(classes)

    # All metagraph edges (unmerged, from d=1)
    all_edges_d1 = set()
    for mask in range(1 << m):
        cn1 = canon_map[mask]
        for wi in range(m):
            cn2 = canon_map[mask ^ (1 << wi)]
            if cn1 != cn2:
                all_edges_d1.add((min(cn1,cn2), max(cn1,cn2)))

    # Actually: the FULL set of metagraph edges includes those reachable
    # at ANY distance. But by definition, if two classes are ever connected
    # (at any d), they're connected. The "metagraph" has all such edges.
    # We want to find which d first creates each edge.

    # Coverage by distance
    max_d = min(m, 8)  # limit to d<=8 for tractability
    edges_by_d = {}  # d -> set of new edges at this distance
    cumul_edges = set()
    sl_by_d = {}  # d -> self-loop count

    for d in range(1, max_d + 1):
        t_d = time.time()
        new_edges = set()
        sl_count = 0
        total_pairs = 0

        # For each tiling: generate all C(m,d) distance-d neighbors
        if comb(m, d) > 500 and d > 3:
            # Sample for large d
            n_samples = min(5000, 2**m)
            masks = [__import__('random').randint(0, 2**m - 1) for _ in range(n_samples)]
            n_combos = min(200, comb(m, d))
            for mask in masks:
                cn1 = canon_map[mask]
                sampled = 0
                for combo in combinations(range(m), d):
                    flip_mask = mask
                    for wi in combo:
                        flip_mask ^= (1 << wi)
                    cn2 = canon_map[flip_mask]
                    if cn1 == cn2:
                        sl_count += 1
                    else:
                        edge = (min(cn1, cn2), max(cn1, cn2))
                        if edge not in cumul_edges:
                            new_edges.add(edge)
                    sampled += 1
                    if sampled >= n_combos:
                        break
                total_pairs += sampled
        else:
            # Full enumeration
            for mask in range(1 << m):
                cn1 = canon_map[mask]
                for combo in combinations(range(m), d):
                    flip_mask = mask
                    for wi in combo:
                        flip_mask ^= (1 << wi)
                    cn2 = canon_map[flip_mask]
                    if cn1 == cn2:
                        sl_count += 1
                    else:
                        edge = (min(cn1, cn2), max(cn1, cn2))
                        if edge not in cumul_edges:
                            new_edges.add(edge)
                    total_pairs += 1

        edges_by_d[d] = new_edges
        cumul_edges |= new_edges
        sl_by_d[d] = sl_count

        # Total possible edges
        max_edges = V * (V-1) // 2

        elapsed_d = time.time() - t_d
        print(f"  d={d}: {len(new_edges):5d} new edges, cumul={len(cumul_edges):5d}/{max_edges} ({len(cumul_edges)/max_edges*100:5.1f}%), SL={sl_count}, C(m,d)={comb(m,d):>6d} ({elapsed_d:.1f}s)")

        # Check if complete
        if len(cumul_edges) == max_edges:
            print(f"  *** COMPLETE at d={d}! k* = {d} = diameter ***")
            break

    # Also check d=m (complement = blue/black)
    if m <= max_d:
        pass  # already computed
    else:
        d = m
        t_d = time.time()
        new_edges_m = set()
        sl_m = 0
        for mask in range(1 << m):
            cn1 = canon_map[mask]
            cn2 = canon_map[mask ^ ((1 << m) - 1)]  # flip all
            if cn1 == cn2:
                sl_m += 1
            else:
                edge = (min(cn1, cn2), max(cn1, cn2))
                if edge not in cumul_edges:
                    new_edges_m.add(edge)
        edges_by_d[m] = new_edges_m
        print(f"  d={m} (complement): {len(new_edges_m):5d} new edges, SL={sl_m} ({time.time()-t_d:.1f}s)")

    # Summary
    print(f"\n  SUMMARY:")
    print(f"    V = {V}, max edges = {V*(V-1)//2}")
    print(f"    {'d':>3} {'C(m,d)':>8} {'new':>6} {'cumul':>6} {'%':>6} {'SL':>8}")
    for d in sorted(edges_by_d.keys()):
        c = len(edges_by_d[d])
        cumul = sum(len(edges_by_d[dd]) for dd in edges_by_d if dd <= d)
        pct = cumul / (V*(V-1)//2) * 100
        print(f"    {d:3d} {comb(m,d):8d} {c:6d} {cumul:6d} {pct:5.1f}% {sl_by_d.get(d,0):8d}")

    elapsed = time.time() - t0
    print(f"\n  Total time: {elapsed:.0f}s")

print("\nDONE.")
print("=" * 80)
