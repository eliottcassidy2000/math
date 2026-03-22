#!/usr/bin/env python3
"""
morse_deep_s20x.py -- kind-pasteur-2026-03-22-S20x

Deep analysis of the Morse theory of H on tournament space.

Key finding: at n=6, H has TWO critical values for local maxima:
  H=45 (global max, 480 tournaments) and H=37 (local max, 720 tournaments)

Questions:
1. What score class are the H=37 local maxima in?
2. Are they SC or NSC?
3. What is the barrier height between H=37 and H=45 basins?
4. What happens at n=7? (too large for exhaustive, but we can sample)

Author: kind-pasteur-2026-03-22-S20x
"""
import sys
import numpy as np
from math import comb
from collections import defaultdict
from itertools import permutations
import random
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def count_hc(A, n):
    full = (1 << n) - 1
    dp = defaultdict(int)
    dp[(1, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[(full, v)] for v in range(n) if A[v][0])

def is_sc(A, n):
    A_comp = np.zeros_like(A)
    for i in range(n):
        for j in range(n):
            if i != j: A_comp[i][j] = 1 - A[i][j]
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if A[perm[i]][perm[j]] != A_comp[i][j]:
                    match = False
                    break
            if not match: break
        if match: return True
    return False

print("=" * 65)
print("  MORSE THEORY OF H: DEEP ANALYSIS AT n=6")
print("=" * 65)

n = 6
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# Compute H for all tournaments
print(f"\n  Computing H for all {2**m} tournaments at n={n}...")
H_map = {}
for bits in range(2**m):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H_map[bits] = count_hp(A, n)
    if (bits + 1) % 10000 == 0:
        print(f"    {bits+1}/{2**m}...")

# Find local maxima
print(f"\n  Finding local maxima...")
local_maxima = []
for bits in range(2**m):
    h = H_map[bits]
    is_max = True
    for k in range(m):
        if H_map[bits ^ (1 << k)] > h:
            is_max = False
            break
    if is_max:
        local_maxima.append(bits)

print(f"  Found {len(local_maxima)} local maxima")

# Analyze each local maximum
print(f"\n  LOCAL MAXIMA ANALYSIS:")
print(f"  {'H':>4s} {'HC':>4s} {'L':>4s} {'score':>25s} {'c3':>4s} {'SC':>4s} {'count':>6s}")

max_by_type = defaultdict(list)
for bits in local_maxima:
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    s = A.sum(axis=1).astype(int)
    S2 = int(sum(s*s))
    score = tuple(sorted(s))
    H = H_map[bits]
    HC = count_hc(A, n)
    L = H - n * HC
    c3 = comb(n,3) - (S2 - comb(n,2)) // 2

    # SC check only for SC-eligible scores
    can_sc = all(score[i] + score[n-1-i] == n-1 for i in range(n//2))
    sc = is_sc(A, n) if can_sc else False

    max_by_type[(H, HC, L, score, c3, sc)].append(bits)

for key, bits_list in sorted(max_by_type.items(), key=lambda x: -x[0][0]):
    H, HC, L, score, c3, sc = key
    print(f"  {H:>4d} {HC:>4d} {L:>4d} {str(list(score)):>25s} {c3:>4d} {str(sc):>4s} {len(bits_list):>6d}")

# ================================================================
# THE BARRIER
# ================================================================
print(f"\n  BARRIER ANALYSIS: PATH FROM H=37 LOCAL MAX TO H=45 GLOBAL MAX")
print()

# Find one H=37 local max and one H=45 local max
max37 = [b for b in local_maxima if H_map[b] == 37][0] if any(H_map[b] == 37 for b in local_maxima) else None
max45 = [b for b in local_maxima if H_map[b] == 45][0] if any(H_map[b] == 45 for b in local_maxima) else None

if max37 and max45:
    # Hamming distance between them
    xor = max37 ^ max45
    hamming = bin(xor).count('1')
    print(f"  H=37 max: bits={max37:015b}")
    print(f"  H=45 max: bits={max45:015b}")
    print(f"  Hamming distance: {hamming}")
    print(f"  Differing arcs: ", end="")
    for k in range(m):
        if (xor >> k) & 1:
            print(f"{pairs[k]}", end=" ")
    print()

    # BFS from H=37 max to H=45 max: what's the minimum H along any path?
    # This is expensive for full BFS, so let's just check the greedy path
    # (flip arcs that increase H the most, or decrease it the least)
    print(f"\n  Greedy path from H=37 max toward H=45 max:")
    current = max37
    path = [(current, H_map[current])]
    visited = {current}
    for step in range(hamming + 10):  # allow some extra steps
        h_curr = H_map[current]
        # Among neighbors, pick the one closest to max45 (shortest Hamming)
        # that we haven't visited
        best = None
        best_score = -1e9
        for k in range(m):
            nb = current ^ (1 << k)
            if nb in visited: continue
            h_nb = H_map[nb]
            # Score: higher H is better, closer to max45 is better
            hamming_to_goal = bin(nb ^ max45).count('1')
            score_val = h_nb - 2 * hamming_to_goal  # heuristic
            if score_val > best_score:
                best_score = score_val
                best = nb
        if best is None: break
        current = best
        visited.add(current)
        path.append((current, H_map[current]))
        if current == max45: break

    min_H_on_path = min(h for _, h in path)
    print(f"  Path length: {len(path)} steps")
    print(f"  Min H on path: {min_H_on_path}")
    print(f"  Barrier height: {37 - min_H_on_path}")
    print(f"  H along path: {[h for _, h in path]}")

# ================================================================
# THE n=5 COMPARISON: WHY NO LOCAL MAX EXCEPT GLOBAL?
# ================================================================
print(f"\n  WHY n=5 HAS NO NON-GLOBAL LOCAL MAXIMA:")
print()

# At n=5, H=13 is the second-highest value. Do H=13 tournaments always
# have a neighbor with H=15?
n5 = 5
pairs5 = [(i,j) for i in range(n5) for j in range(i+1, n5)]
m5 = len(pairs5)

H_map5 = {}
for bits in range(2**m5):
    A = np.zeros((n5,n5), dtype=np.int8)
    for k, (i,j) in enumerate(pairs5):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H_map5[bits] = count_hp(A, n5)

# Check each non-maximum H value: does every tournament at that H have a neighbor with higher H?
for H_val in sorted(set(H_map5.values())):
    if H_val == 15: continue
    bits_at_H = [b for b in range(2**m5) if H_map5[b] == H_val]
    stuck = 0
    for bits in bits_at_H:
        has_higher = False
        for k in range(m5):
            if H_map5[bits ^ (1 << k)] > H_val:
                has_higher = True
                break
        if not has_higher:
            stuck += 1
    print(f"  H={H_val}: {len(bits_at_H)} tournaments, {stuck} are local maxima")

# At n=6, same analysis for H=37 and nearby
print(f"\n  n=6 STUCKNESS ANALYSIS:")
for H_val in [33, 37, 41, 43, 45]:
    bits_at_H = [b for b in range(2**m) if H_map[b] == H_val]
    stuck = 0
    max_neighbor = []
    for bits in bits_at_H:
        has_higher = False
        max_nb_H = 0
        for k in range(m):
            nb_H = H_map[bits ^ (1 << k)]
            if nb_H > H_val: has_higher = True
            max_nb_H = max(max_nb_H, nb_H)
        if not has_higher: stuck += 1
        max_neighbor.append(max_nb_H)
    avg_max_nb = sum(max_neighbor) / len(max_neighbor)
    print(f"  H={H_val}: {len(bits_at_H)} tournaments, {stuck} local maxima, avg max neighbor H = {avg_max_nb:.1f}")

# ================================================================
# THE MORSE INDEX (number of downward directions)
# ================================================================
print(f"\n  MORSE INDEX OF LOCAL MAXIMA AT n=6:")
print()

for bits in local_maxima[:20]:
    H = H_map[bits]
    down = 0
    same = 0
    for k in range(m):
        nb_H = H_map[bits ^ (1 << k)]
        if nb_H < H: down += 1
        elif nb_H == H: same += 1
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    score = tuple(sorted(A.sum(axis=1).astype(int)))
    print(f"  H={H}, down={down}, same={same}, up={m-down-same}, score={score}")

# ================================================================
# SUMMARY TABLE
# ================================================================
print(f"\n" + "=" * 65)
print(f"  MORSE THEORY SUMMARY TABLE")
print(f"=" * 65)
print()

for n_val, m_val in [(3, 3), (4, 6), (5, 10), (6, 15)]:
    if n_val <= 5:
        H_map_temp = {}
        pairs_temp = [(i,j) for i in range(n_val) for j in range(i+1, n_val)]
        for bits in range(2**m_val):
            A = np.zeros((n_val,n_val), dtype=np.int8)
            for k, (i,j) in enumerate(pairs_temp):
                if (bits >> k) & 1: A[i][j] = 1
                else: A[j][i] = 1
            H_map_temp[bits] = count_hp(A, n_val)

        lmax = 0; lmin = 0; saddle = 0
        lmax_H = set()
        for bits in range(2**m_val):
            h = H_map_temp[bits]
            has_up = any(H_map_temp[bits ^ (k)] > h for k in range(m_val))
            has_down = any(H_map_temp[bits ^ (k)] < h for k in range(m_val))
            if not has_up: lmax += 1; lmax_H.add(h)
            elif not has_down: lmin += 1
            else: saddle += 1

        H_max = max(H_map_temp.values())
        H_min = min(H_map_temp.values())
        n_H = len(set(H_map_temp.values()))
        print(f"  n={n_val}: m={m_val}, |T|={2**m_val}, H in [{H_min},{H_max}], {n_H} H-values")
        print(f"    Local maxima: {lmax} (at H={sorted(lmax_H)})")
        print(f"    Local minima: {lmin}")
        print(f"    Saddle points: {saddle}")
        print(f"    chi = max-saddle+min = {lmax - saddle + lmin}")
        print()
    else:
        lmax = len(local_maxima)
        lmin_count = sum(1 for b in range(2**m) if all(H_map[b ^ (1<<k)] >= H_map[b] for k in range(m)) and all(H_map[b ^ (1<<k)] > H_map[b] for k in range(m) if True))
        # Recount properly
        lmin = 0; saddle = 0; lmax2 = 0
        lmax_H = set()
        for bits in range(2**m):
            h = H_map[bits]
            has_up = any(H_map[bits ^ (1 << k)] > h for k in range(m))
            has_down = any(H_map[bits ^ (1 << k)] < h for k in range(m))
            if not has_up: lmax2 += 1; lmax_H.add(h)
            elif not has_down: lmin += 1
            else: saddle += 1

        H_max = max(H_map.values())
        H_min = min(H_map.values())
        n_H = len(set(H_map.values()))
        print(f"  n={n_val}: m={m_val}, |T|={2**m_val}, H in [{H_min},{H_max}], {n_H} H-values")
        print(f"    Local maxima: {lmax2} (at H={sorted(lmax_H)})")
        print(f"    Local minima: {lmin}")
        print(f"    Saddle points: {saddle}")
        print(f"    chi = max-saddle+min = {lmax2 - saddle + lmin}")
        print()
