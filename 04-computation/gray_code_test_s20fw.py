#!/usr/bin/env python3
"""
gray_code_test_s20fw.py — Does a monotone-H Gray code exist?
kind-pasteur-2026-03-24-S20fw

PREDICTION 5: A Gray code on Q_m where H is weakly non-decreasing
exists at n=4 (m=3, 8 tilings) but NOT at n=5 (m=6, 64 tilings).

A Gray code = Hamiltonian path on Q_m (visit each vertex once,
each step flips one bit). Monotone-H = H values along the path
are weakly non-decreasing.

At n=4: 8 tilings, Q_3 = 3-cube. Gray codes on Q_3 are well-known.
Check if ANY of them gives monotone H.

At n=5: 64 tilings, Q_6. Too many Gray codes to enumerate, but we
can check if the H-value structure ALLOWS monotone traversal by
checking the "descent" structure of H on Q_m.
"""

import sys
from math import factorial
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  MONOTONE-H GRAY CODE TEST")
print("  kind-pasteur-2026-03-24-S20fw")
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

    # Compute H for all tilings
    H_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        H_map[mask] = count_hp(A)

    H_vals = sorted(set(H_map.values()))
    print(f"\n  n = {n}, m = {m}, 2^m = {2**m}")
    print(f"  H values: {H_vals}")
    print(f"  H distribution: {dict(sorted([(h, sum(1 for v in H_map.values() if v == h)) for h in H_vals]))}")

    if n == 4:
        # n=4: m=3, Q_3 has 8 vertices. Find ALL Hamiltonian paths on Q_3.
        # Q_3 adjacency: vertices differ in 1 bit
        adj = defaultdict(list)
        for mask in range(1 << m):
            for wi in range(m):
                adj[mask].append(mask ^ (1 << wi))

        # DFS to find all Hamiltonian paths
        paths = []
        def dfs(path, visited):
            if len(path) == (1 << m):
                paths.append(list(path))
                return
            for nb in adj[path[-1]]:
                if nb not in visited:
                    visited.add(nb)
                    path.append(nb)
                    dfs(path, visited)
                    path.pop()
                    visited.remove(nb)

        for start in range(1 << m):
            dfs([start], {start})

        print(f"  Total Hamiltonian paths on Q_{m}: {len(paths)}")

        # Check monotonicity
        monotone_inc = 0
        monotone_dec = 0
        best_ascents = 0
        best_path = None

        for path in paths:
            H_seq = [H_map[v] for v in path]
            is_inc = all(H_seq[i] <= H_seq[i+1] for i in range(len(H_seq)-1))
            is_dec = all(H_seq[i] >= H_seq[i+1] for i in range(len(H_seq)-1))
            ascents = sum(1 for i in range(len(H_seq)-1) if H_seq[i] < H_seq[i+1])
            if is_inc: monotone_inc += 1
            if is_dec: monotone_dec += 1
            if ascents > best_ascents:
                best_ascents = ascents
                best_path = path

        print(f"  Monotone non-decreasing paths: {monotone_inc}")
        print(f"  Monotone non-increasing paths: {monotone_dec}")
        print(f"  Best path (most ascents): {best_ascents} out of {(1<<m)-1}")
        if best_path:
            H_seq = [H_map[v] for v in best_path]
            print(f"  Best H sequence: {H_seq}")

        # Also check: paths that are monotone on ISO CLASS level (not just H)
        # (since different tournaments can have same H)

    elif n == 5:
        # n=5: 64 tilings on Q_6. Too many paths to enumerate.
        # Instead: check necessary conditions.

        # For monotone H: at each H-level, can we enter from below and exit upward?
        # Count: for each H-value h, how many tilings have ALL neighbors with H >= h?
        # If any h has 0 such "upward-exitable" tilings, monotone is impossible.

        H_levels = sorted(set(H_map.values()))
        print(f"\n  NECESSARY CONDITIONS for monotone Gray code:")

        for h in H_levels:
            tilings_at_h = [mask for mask, H in H_map.items() if H == h]
            # For each tiling at level h: how many neighbors have H > h? H = h? H < h?
            up_exits = 0
            for mask in tilings_at_h:
                has_up = any(H_map[mask ^ (1 << wi)] > h for wi in range(m))
                if has_up: up_exits += 1

            down_entries = 0
            for mask in tilings_at_h:
                has_down = any(H_map[mask ^ (1 << wi)] < h for wi in range(m))
                if has_down: down_entries += 1

            print(f"    H={h:3d}: {len(tilings_at_h):3d} tilings, {up_exits:3d} can go up, {down_entries:3d} can come from below")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
