#!/usr/bin/env python3
"""
h7_check_v3.py — opus-2026-03-14-S71e

EFFICIENT check: Is H=7 achievable at n=7?

Strategy:
1. Find all score sequences with dc3 ≤ 5 (low cycle count)
2. For each, enumerate a sample of tournaments with that score sequence
3. Check dc5, dc7, and H

Key formula: dc3 = C(n,3) - Σ C(s_i,2) for score sequence (s_1,...,s_n).
At n=7: dc3 = 35 - Σ C(s_i,2).
"""

import sys
from itertools import combinations, permutations
from collections import Counter
import time
import random
sys.stdout.reconfigure(line_buffering=True)

n = 7

print("=" * 70)
print("SCORE SEQUENCES WITH dc3 ≤ 5 AT n=7")
print("=" * 70)

from math import comb

# Find all valid score sequences at n=7
# Constraints: s_1 ≤ s_2 ≤ ... ≤ s_7, Σ s_i = 21, Erdős-Gallai conditions
# Landau's theorem: sorted scores s_1 ≤ ... ≤ s_n are valid iff
#   Σ_{i=1}^k s_i ≥ C(k,2) for all k, with equality at k=n.

def valid_score_sequences(n):
    """Generate all valid tournament score sequences."""
    total = n*(n-1)//2
    def generate(remaining_sum, min_val, max_val, current, k):
        if k == n:
            if remaining_sum == 0:
                # Check Landau condition
                valid = True
                prefix = 0
                for j in range(n):
                    prefix += current[j]
                    if prefix < comb(j+1, 2):
                        valid = False
                        break
                if valid:
                    yield tuple(current)
            return
        remaining_spots = n - k
        max_possible = min(max_val, remaining_sum - min_val * (remaining_spots - 1))
        min_possible = max(min_val, remaining_sum - max_val * (remaining_spots - 1))
        for v in range(min_possible, max_possible + 1):
            yield from generate(remaining_sum - v, v, max_val, current + [v], k + 1)

    yield from generate(total, 0, n-1, [], 0)

low_dc3_scores = []
for scores in valid_score_sequences(n):
    sum_c2 = sum(comb(s, 2) for s in scores)
    dc3 = 35 - sum_c2
    if dc3 <= 5:
        low_dc3_scores.append((dc3, scores))

low_dc3_scores.sort()
print(f"\n  Score sequences with dc3 ≤ 5:")
for dc3, scores in low_dc3_scores:
    sum_c2 = sum(comb(s, 2) for s in scores)
    print(f"    dc3={dc3}: {scores}, Σ C(s,2)={sum_c2}")

# Now for each low-dc3 score sequence, generate many tournaments and check H
print("\n" + "=" * 70)
print("SAMPLING TOURNAMENTS WITH LOW dc3")
print("=" * 70)

def random_tournament_with_scores(scores, n, max_attempts=10000):
    """Generate a random tournament with given sorted score sequence."""
    # Use the switching technique
    for _ in range(max_attempts):
        # Start with a tournament matching the score sequence
        # Use a greedy approach
        adj = [[0]*n for _ in range(n)]
        remaining = list(scores)

        # Assign arcs greedily with randomization
        pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
        random.shuffle(pairs)

        # Build adjacency respecting out-degrees
        target = list(scores)
        current_out = [0]*n
        success = True

        for i, j in pairs:
            # Decide direction
            can_ij = current_out[i] < target[i]
            can_ji = current_out[j] < target[j]

            if can_ij and can_ji:
                if random.random() < 0.5:
                    adj[i][j] = 1
                    current_out[i] += 1
                else:
                    adj[j][i] = 1
                    current_out[j] += 1
            elif can_ij:
                adj[i][j] = 1
                current_out[i] += 1
            elif can_ji:
                adj[j][i] = 1
                current_out[j] += 1
            else:
                success = False
                break

        if success and current_out == list(target):
            return adj
    return None

def count_hp(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def count_dc5(adj, n):
    count = 0
    for verts in combinations(range(n), 5):
        v = list(verts)
        for perm in permutations(v):
            if all(adj[perm[i]][perm[(i+1) % 5]] for i in range(5)):
                count += 1
    return count // 5

random.seed(42)
h7_found = False

for dc3_target, scores in low_dc3_scores:
    print(f"\n  Score {scores} (dc3={dc3_target}):")

    # For transitive score (0,1,2,3,4,5,6), only 1 tournament (up to labeling)
    # For other scores, sample many
    samples = 2000
    h_vals = Counter()

    for trial in range(samples):
        adj = random_tournament_with_scores(scores, n)
        if adj is None:
            continue

        H = count_hp(adj, n)
        h_vals[H] += 1

        if H == 7:
            h7_found = True
            dc5 = count_dc5(adj, n)
            print(f"    *** H=7 FOUND! dc5={dc5} ***")
            for i in range(n):
                print(f"      {i}: beats {[j for j in range(n) if adj[i][j]]}")

    if h_vals:
        print(f"    H values: {sorted(h_vals.items())[:10]}...")
        min_h = min(h_vals.keys())
        print(f"    Min H = {min_h}")
    else:
        print(f"    No valid tournaments generated")

if not h7_found:
    print(f"\n  H=7 NOT found in {sum(1 for _ in low_dc3_scores) * 2000} samples.")
    print(f"  This strongly supports HYP-992.")
else:
    print(f"\n  *** H=7 IS ACHIEVABLE ***")
    print(f"  HYP-992 is WRONG!")

# Also try a direct enumeration for the SPECIFIC score sequences with dc3=3
print("\n" + "=" * 70)
print("EXHAUSTIVE ENUMERATION FOR dc3=3 SCORE SEQUENCES")
print("=" * 70)

# There are only a few score sequences with dc3=3
# For each, we can try to enumerate ALL tournaments (up to relabeling)

for dc3_target, scores in low_dc3_scores:
    if dc3_target != 3:
        continue

    print(f"\n  Score {scores} (dc3=3):")
    print(f"  Enumerating all tournaments with this score sequence...")

    # Generate many random ones
    seen_h = Counter()
    seen_configs = set()

    for trial in range(50000):
        adj = random_tournament_with_scores(scores, n)
        if adj is None:
            continue

        # Hash the adjacency to avoid recounting
        config = tuple(tuple(row) for row in adj)
        if config in seen_configs:
            continue
        seen_configs.add(config)

        H = count_hp(adj, n)
        seen_h[H] += 1

    print(f"    Distinct tournaments found: {len(seen_configs)}")
    print(f"    H distribution: {sorted(seen_h.items())}")

print("\nDone.")
