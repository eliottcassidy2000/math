#!/usr/bin/env python3
"""
tiling_profiles_s20fn.py — How tilings within an iso class see different neighborhoods
kind-pasteur-2026-03-24-S20fn

For each tiling T in class C:
  profile(T) = [class(T flip A), class(T flip B), class(T flip C), ...]
  = one class per wiggly position

Different tilings of C have different profiles. But there's hidden order:

1. Grid-symmetry: T and grid(T) have paired profiles
   (wiggly class X gives the same target as class grid(X))

2. The Aut(C) group acts on the profiles: if |Aut|=1, all profiles
   are distinct; if |Aut|=k, profiles come in orbits of size k.

3. The MULTISET of target classes (ignoring tile labels) might be
   constant across all tilings of C — or might vary.

Study this for important classes: transitive, regular, near-regular.
"""

import sys
from math import factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  TILING NEIGHBORHOOD PROFILES WITHIN ISO CLASSES")
print("  kind-pasteur-2026-03-24-S20fn")
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

    # Grid symmetry map
    tile_idx_map = {(t[0],t[1]): i for i, t in enumerate(TILES)}
    grid_map = []
    for i, (x, y) in enumerate(TILES):
        gx, gy = N+1-y, N+1-x
        grid_map.append(tile_idx_map.get((gx, gy), i))

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

    def count_aut(A):
        return sum(1 for p in all_perms if all(A[p[i]][p[j]] == A[i][j] for i in range(N) for j in range(N)))

    # Build all tilings
    canon_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        canon_map[mask] = canonicalize(A)

    classes = sorted(set(canon_map.values()))
    V = len(classes)
    class_short = {cn: i for i, cn in enumerate(classes)}
    class_tilings = defaultdict(list)
    for mask, cn in canon_map.items():
        class_tilings[cn].append(mask)

    # Class properties
    ci = {}
    for cn in classes:
        mask = class_tilings[cn][0]
        A = b2a([(mask >> k) & 1 for k in range(m)])
        ci[cn] = {'H': count_hp(A), 'aut': count_aut(A), 'tilings': len(class_tilings[cn])}

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m}, V = {V}")
    print(f"{'#'*60}")

    # Study specific important classes
    # Sort by H to find transitive (H=1), near-transitive (H=3), regular (max H)
    by_H = sorted(classes, key=lambda c: ci[c]['H'])

    # Pick a few representative classes
    study_classes = []
    seen_H = set()
    for cn in by_H:
        H = ci[cn]['H']
        if H not in seen_H and ci[cn]['tilings'] >= 1:
            study_classes.append(cn)
            seen_H.add(H)
        if len(study_classes) >= 8:
            break

    for cn in study_classes:
        tilings_list = class_tilings[cn]
        n_tilings = len(tilings_list)
        H = ci[cn]['H']
        aut = ci[cn]['aut']

        print(f"\n  === CLASS: H={H}, |Aut|={aut}, tilings={n_tilings} ===")

        # For each tiling: compute the full profile
        # profile = tuple of class indices for each wiggly position
        profiles = []
        for mask in tilings_list:
            profile = []
            for wi in range(m):
                ncn = canon_map[mask ^ (1 << wi)]
                profile.append(class_short[ncn])
            profiles.append(tuple(profile))

        # How many DISTINCT profiles?
        distinct_profiles = len(set(profiles))
        print(f"    Distinct profiles: {distinct_profiles}/{n_tilings}")
        print(f"    Ratio (should be tilings/|Aut| if Aut acts freely on profiles): {distinct_profiles}, tilings/|Aut| = {n_tilings/aut:.1f}")

        # How many distinct MULTISETS (ignoring tile labels)?
        multisets = [tuple(sorted(p)) for p in profiles]
        distinct_multisets = len(set(multisets))
        print(f"    Distinct multisets (ignoring tile position): {distinct_multisets}/{n_tilings}")

        # Grid symmetry: does grid(profile) give the paired profile?
        grid_consistent = 0
        for mask in tilings_list:
            # Grid of mask
            bits = [(mask >> k) & 1 for k in range(m)]
            grid_bits = [0]*m
            for i in range(m):
                grid_bits[grid_map[i]] = bits[i]
            grid_mask = sum(b << k for k, b in enumerate(grid_bits))

            if grid_mask in canon_map:
                # Check: profile of grid_mask should be grid-permuted profile of mask
                p1 = []
                p2 = []
                for wi in range(m):
                    p1.append(class_short[canon_map[mask ^ (1 << wi)]])
                    p2.append(class_short[canon_map[grid_mask ^ (1 << wi)]])

                # p2[wi] should equal p1[grid_map[wi]] (paired tile gives same class)
                grid_ok = all(p2[wi] == p1[grid_map[wi]] for wi in range(m))
                if grid_ok:
                    grid_consistent += 1

        print(f"    Grid-symmetry consistent: {grid_consistent}/{n_tilings} ({'ALL' if grid_consistent == n_tilings else 'PARTIAL'})")

        # Show a few profiles
        if n_tilings <= 20:
            print(f"    PROFILES (class index per wiggly position {' '.join(labels[:m])}):")
            for i, mask in enumerate(tilings_list[:10]):
                profile = profiles[i]
                self_count = sum(1 for x in profile if x == class_short[cn])
                target_multiset = Counter(profile)
                print(f"      T{i}: [{','.join(str(x) for x in profile)}] self={self_count} distinct={len(set(profile))}")
        else:
            # Summarize profile statistics
            self_counts = [sum(1 for x in p if x == class_short[cn]) for p in profiles]
            distinct_counts = [len(set(p)) for p in profiles]
            print(f"    Self-neighbor counts: {Counter(self_counts)}")
            print(f"    Distinct-neighbor counts: {Counter(distinct_counts)}")

        # The MULTISET distribution: which multisets appear and how often?
        ms_counter = Counter(multisets)
        print(f"    Multiset distribution: {len(ms_counter)} distinct multisets")
        for ms, count in ms_counter.most_common(5):
            self_in_ms = ms.count(class_short[cn])
            print(f"      count={count}: targets={ms}, self={self_in_ms}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
