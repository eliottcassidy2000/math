#!/usr/bin/env python3
"""
delta_gap_s20x.py -- kind-pasteur-2026-03-22-S20x

INVESTIGATE THE DELTA GAP: at n=5, delta=+-10 is missing from arc flip deltas.
10 = C(5,2) = the dimension of tournament space.

Conjecture: delta=+-C(n,2) is always missing?
Or: delta=+-10 is missing because 10 = 2*(H_max-H_min)/3?
Or: some other structural reason?

Also: explore the basin structure at n=6 via sampling.

Author: kind-pasteur-2026-03-22-S20x
"""
import sys
import numpy as np
from math import comb
from collections import defaultdict
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

print("=" * 65)
print("  THE DELTA GAP: WHY IS +/-10 MISSING?")
print("=" * 65)

# ================================================================
# 1. Delta distribution at n=3, 4, 5 (exhaustive)
# ================================================================
print()
for n in [3, 4, 5]:
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    H_map = {}
    for bits in range(2**m):
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        H_map[bits] = count_hp(A, n)

    delta_set = set()
    for bits in range(2**m):
        for k in range(m):
            bits2 = bits ^ (1 << k)
            delta_set.add(H_map[bits2] - H_map[bits])

    deltas = sorted(delta_set)
    # Missing even values
    min_d, max_d = min(deltas), max(deltas)
    even_range = set(range(min_d, max_d + 1, 2))
    missing = sorted(even_range - set(deltas))

    print(f"  n={n}: m=C({n},2)={m}")
    print(f"    Deltas: {deltas}")
    print(f"    Missing even: {missing}")
    print(f"    m in missing: {m in [abs(x) for x in missing]}")
    print()

# ================================================================
# 2. Delta at n=6 via sampling (exhaustive too expensive: 32768*15)
# Actually 32768*15 = 491520, that's fine for computing deltas
# The bottleneck is HP count: 32768 tournaments already computed above
# ================================================================
print("  n=6: EXHAUSTIVE DELTA ANALYSIS")
print("  (Computing H for all 32768 tournaments...)")

n = 6
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)  # 15

H_map = {}
for bits in range(2**m):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H_map[bits] = count_hp(A, n)
    if (bits + 1) % 5000 == 0:
        print(f"    {bits+1}/{2**m}...")

print("  Computing deltas...")
delta_set = set()
delta_dist = defaultdict(int)
for bits in range(2**m):
    for k in range(m):
        bits2 = bits ^ (1 << k)
        d = H_map[bits2] - H_map[bits]
        delta_set.add(d)
        delta_dist[d] += 1

deltas = sorted(delta_set)
min_d, max_d = min(deltas), max(deltas)
even_range = set(range(min_d, max_d + 1, 2))
missing = sorted(even_range - set(deltas))

print(f"\n  n={n}: m=C({n},2)={m}")
print(f"  Delta range: [{min_d}, {max_d}]")
print(f"  Distinct deltas: {len(deltas)}")
print(f"  Missing even values: {missing}")
print(f"  m={m} in missing: {m in [abs(x) for x in missing]}")
print()

# Full distribution
print(f"  Delta distribution:")
for d in deltas:
    count = delta_dist[d]
    pct = 100 * count / (2**m * m)
    bar = '#' * int(pct * 2)
    print(f"    {d:>+4d}: {count:>8d} ({pct:>5.2f}%) {bar}")

# ================================================================
# 3. THE LEVEL SET STRUCTURE AT n=6
# ================================================================
print(f"\n  n=6: LEVEL SET STRUCTURE")
print()

H_vals = sorted(set(H_map.values()))
print(f"  {'H':>4s} {'count':>6s} {'avg_same':>9s} {'avg_up':>8s} {'avg_down':>8s} {'local_max':>9s} {'local_min':>9s}")

local_max_list = []
local_min_list = []

for H_val in H_vals:
    bits_at_H = [b for b, h in H_map.items() if h == H_val]
    total_same = 0; total_up = 0; total_down = 0
    lmax = 0; lmin = 0
    for bits in bits_at_H:
        same = 0; up = 0; down = 0
        for k in range(m):
            bits2 = bits ^ (1 << k)
            h2 = H_map[bits2]
            if h2 == H_val: same += 1
            elif h2 > H_val: up += 1
            else: down += 1
        total_same += same
        total_up += up
        total_down += down
        if up == 0: lmax += 1; local_max_list.append(bits)
        if down == 0: lmin += 1; local_min_list.append(bits)
    cnt = len(bits_at_H)
    print(f"  {H_val:>4d} {cnt:>6d} {total_same/cnt:>9.2f} {total_up/cnt:>8.2f} {total_down/cnt:>8.2f} {lmax:>9d} {lmin:>9d}")

print(f"\n  Total local maxima: {len(local_max_list)}")
print(f"  Total local minima: {len(local_min_list)}")

# Are all local maxima at the same H?
max_H_vals = set(H_map[b] for b in local_max_list)
min_H_vals = set(H_map[b] for b in local_min_list)
print(f"  Max H values at local maxima: {sorted(max_H_vals)}")
print(f"  Min H values at local minima: {sorted(min_H_vals)}")

# ================================================================
# 4. BASIN ANALYSIS AT n=6 (sampling for speed)
# ================================================================
print(f"\n  BASIN ANALYSIS (steepest ascent from random starts)")
print()

random.seed(42)

def steepest_ascent(bits, H_map, m):
    current = bits
    for _ in range(1000):  # max steps
        h_curr = H_map[current]
        best_h = h_curr
        best_next = current
        for k in range(m):
            bits2 = current ^ (1 << k)
            if H_map[bits2] > best_h:
                best_h = H_map[bits2]
                best_next = bits2
        if best_next == current:
            break
        current = best_next
    return current

# Sample 1000 random starting points
basin_peaks = defaultdict(int)
for _ in range(2000):
    start = random.randint(0, 2**m - 1)
    peak = steepest_ascent(start, H_map, m)
    basin_peaks[H_map[peak]] += 1

print(f"  Gradient ascent destinations (2000 random starts):")
for h_val, count in sorted(basin_peaks.items(), key=lambda x: -x[1]):
    print(f"    H={h_val}: {count} ({100*count/2000:.1f}%)")

# ================================================================
# 5. THE FORBIDDEN DELTA PATTERN
# ================================================================
print(f"\n  THE FORBIDDEN DELTA PATTERN")
print()

for n_val, missing_list in [(3, []), (4, []), (5, [-10, 10])]:
    m_val = comb(n_val, 2)
    H_max = {3: 3, 4: 5, 5: 15}[n_val]
    H_min = 1
    print(f"  n={n_val}: m={m_val}, H_range=[{H_min},{H_max}], missing={missing_list}")

print(f"  n=6: m={m}, H_range=[1,{max(H_vals)}], missing={missing}")
print()

# Check: is the missing value related to H_max - H_min?
# n=5: missing 10, H_max-H_min = 14
# n=6: H_max-H_min = 44
# missing at n=6: let's see what they are

# Check: are the missing values related to the forbidden H values?
print(f"  Forbidden H at n=5: [7]")
print(f"  Missing delta at n=5: [10] (=7+3=forbidden+H_min+2?)")
print(f"  Forbidden H at n=6: {sorted(set(range(1, max(H_vals)+1, 2)) - set(H_vals))}")
print(f"  Missing delta at n=6: {missing}")
print()

# Check primality of missing values
from math import gcd as _gcd

# Are the deltas that DO appear always of the form H_i - H_j
# for ADJACENT tournaments (differing in 1 arc)?
# That's the definition. So the question is: which differences H_i - H_j
# are achievable by adjacent tournaments?

# H values at n=5: 1, 3, 5, 9, 11, 13, 15
# All pairwise differences: 2,4,8,10,12,14 and their negatives
# Of these, 10 and 14 do NOT appear as arc-flip deltas
# (12 appears, -12 appears)

H_n5 = [1, 3, 5, 9, 11, 13, 15]
all_diffs_n5 = sorted(set(abs(a-b) for a in H_n5 for b in H_n5 if a != b))
print(f"  n=5 H values: {H_n5}")
print(f"  All pairwise |H_i - H_j|: {all_diffs_n5}")
print(f"  Missing from arc-flip deltas: {sorted(set(all_diffs_n5) - set(range(0, 13, 2)))}")
print()
print(f"  n=5 achievable deltas (positive): {sorted(d for d in deltas if d > 0)}" if n == 5 else "")

# Correct: let me recompute for n=5 stored data
n5_deltas_pos = [2, 4, 6, 8, 12]
print(f"  n=5 achievable positive deltas: {n5_deltas_pos}")
print(f"  n=5 pairwise positive diffs: {all_diffs_n5}")
print(f"  n=5 non-achievable diffs: {sorted(set(all_diffs_n5) - set(n5_deltas_pos))}")
print(f"  So: diff=10 and diff=14 cannot be achieved by single arc flip")
print(f"  14 is trivially impossible (max delta=12)")
print(f"  10 is the non-trivial gap: TWO adjacent H values never differ by 10")
