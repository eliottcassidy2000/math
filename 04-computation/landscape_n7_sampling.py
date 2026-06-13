#!/usr/bin/env python3
"""
landscape_n7_sampling.py — opus-2026-03-13-S67j

TESTING HYP-741: Does the H landscape have no spurious local maxima at n=7 (odd)?

Strategy: We can't enumerate all 2^21 ≈ 2M tournaments at n=7, but we can:
1. Sample 10000 random tournaments
2. Do greedy ascent on EXACT H (using dynamic programming for efficiency)
3. Check if all greedy paths reach the same H value

For exact H at n=7, we use DP for Hamiltonian path counting:
  dp[mask][v] = number of Hamiltonian paths ending at v using vertex set 'mask'
  O(2^n * n^2) per tournament = 2^7 * 49 ≈ 6272 operations (vs 5040 for brute force)

For greedy ascent, we try all m=21 arc flips and pick the best.
"""

import numpy as np
import random
import time

def ham_path_count_dp(A):
    """Count Hamiltonian paths using DP. O(2^n * n^2)."""
    n = A.shape[0]
    # dp[mask][v] = number of paths using vertices in mask, ending at v
    dp = np.zeros((1 << n, n), dtype=np.int64)

    # Base case: single vertex
    for v in range(n):
        dp[1 << v][v] = 1

    # Fill DP table
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            if not (mask & (1 << v)):
                continue
            # Try extending to each unvisited vertex
            for u in range(n):
                if mask & (1 << u):
                    continue  # already visited
                if A[v][u] == 1:  # arc from v to u
                    dp[mask | (1 << u)][u] += dp[mask][v]

    # Sum over all endpoints for full mask
    full_mask = (1 << n) - 1
    return int(np.sum(dp[full_mask]))

def bits_to_tournament(bits, n, edges):
    A = np.zeros((n, n), dtype=np.int8)
    for k, (i, j) in enumerate(edges):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    return A

def score_seq(A):
    return tuple(sorted(A.sum(axis=1)))

# =====================================================================
# PART 1: VERIFY DP AT SMALL n
# =====================================================================
print("=" * 70)
print("VERIFYING DP HAMILTONIAN PATH COUNT")
print("=" * 70)

from itertools import permutations

for n in [3, 4, 5]:
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    # Check a few random tournaments
    correct = 0
    total = min(50, 2**m)
    for bits in random.sample(range(2**m), total):
        A = bits_to_tournament(bits, n, edges)

        h_dp = ham_path_count_dp(A)

        h_brute = 0
        for perm in permutations(range(n)):
            valid = True
            for i in range(n-1):
                if A[perm[i]][perm[i+1]] != 1:
                    valid = False
                    break
            if valid:
                h_brute += 1

        if h_dp == h_brute:
            correct += 1

    print(f"  n={n}: {correct}/{total} correct")

# =====================================================================
# PART 2: n=7 LANDSCAPE SAMPLING
# =====================================================================
print("\n" + "=" * 70)
print("n=7 LANDSCAPE SAMPLING (m=21, 2^21 = 2,097,152 tournaments)")
print("=" * 70)

n = 7
edges = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(edges)
N = 2**m

n_samples = 5000
random.seed(42)

print(f"  Sampling {n_samples} random tournaments...")
t0 = time.time()

# Cache for H values (avoid recomputation)
H_cache = {}

def get_H(bits):
    if bits not in H_cache:
        A = bits_to_tournament(bits, n, edges)
        H_cache[bits] = ham_path_count_dp(A)
    return H_cache[bits]

# Greedy ascent results
final_H_values = []
final_bits = []
steps_list = []
initial_H_values = []

for trial in range(n_samples):
    bits = random.randint(0, N - 1)
    initial_H_values.append(get_H(bits))

    current = bits
    steps = 0
    while True:
        h_current = get_H(current)
        best_k = -1
        best_delta = 0
        for k in range(m):
            neighbor = current ^ (1 << k)
            h_neighbor = get_H(neighbor)
            delta = h_neighbor - h_current
            if delta > best_delta:
                best_delta = delta
                best_k = k

        if best_k == -1:
            break
        current = current ^ (1 << best_k)
        steps += 1

    final_H_values.append(get_H(current))
    final_bits.append(current)
    steps_list.append(steps)

    if trial % 1000 == 0 and trial > 0:
        elapsed = time.time() - t0
        rate = trial / elapsed
        print(f"    {trial}/{n_samples} done ({rate:.0f}/s), cache size: {len(H_cache)}")

elapsed = time.time() - t0
print(f"  Done in {elapsed:.1f}s (cache size: {len(H_cache)})")

# Analyze results
unique_final_H = sorted(set(final_H_values))
print(f"\n  Unique final H values: {unique_final_H}")

max_final_H = max(final_H_values)
min_final_H = min(final_H_values)

count_global_max = sum(1 for h in final_H_values if h == max_final_H)

print(f"  Maximum H reached: {max_final_H}")
print(f"  Minimum final H: {min_final_H}")
print(f"  Fraction reaching max H: {count_global_max}/{n_samples} = {100*count_global_max/n_samples:.1f}%")

if len(unique_final_H) > 1:
    print(f"\n  *** SPURIOUS LOCAL MAXIMA DETECTED at n=7 ***")
    for h_val in unique_final_H:
        count = sum(1 for h in final_H_values if h == h_val)
        print(f"    H={h_val}: {count} ({100*count/n_samples:.1f}%)")

        # Find a representative and check its score sequence
        idx = final_H_values.index(h_val)
        bits = final_bits[idx]
        A = bits_to_tournament(bits, n, edges)
        ss = score_seq(A)
        print(f"      Score sequence: {ss}")
else:
    print(f"\n  *** NO SPURIOUS LOCAL MAXIMA at n=7 ***")
    print(f"  All {n_samples} random starts converge to H={max_final_H}")

# Average steps
print(f"\n  Steps statistics:")
print(f"    Mean: {np.mean(steps_list):.1f}")
print(f"    Max: {max(steps_list)}")
print(f"    Min: {min(steps_list)}")

# Check: is the maximum H the value for regular tournaments?
# Regular tournament at n=7: all scores = 3
# Build one and check
print(f"\n  Checking regular tournament H value:")
# Paley P_7
QR = set()
for k in range(1, 7):
    QR.add((k*k) % 7)
A_paley = np.zeros((7,7), dtype=np.int8)
for i in range(7):
    for j in range(7):
        if i != j and ((j-i)%7) in QR:
            A_paley[i][j] = 1

H_paley = ham_path_count_dp(A_paley)
print(f"    H(P_7) = {H_paley}")
print(f"    Score sequence: {score_seq(A_paley)}")
print(f"    Is this the max? {H_paley == max_final_H}")

# =====================================================================
# PART 3: QUICK CHECK n=8 (even)
# =====================================================================
print("\n" + "=" * 70)
print("n=8 LANDSCAPE QUICK CHECK (m=28)")
print("=" * 70)

n = 8
edges_8 = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(edges_8)
N_8 = 2**m

H_cache_8 = {}

def get_H_8(bits):
    if bits not in H_cache_8:
        A = bits_to_tournament(bits, n, edges_8)
        H_cache_8[bits] = ham_path_count_dp(A)
    return H_cache_8[bits]

n_samples_8 = 500
random.seed(123)

print(f"  Sampling {n_samples_8} random tournaments...")
t0 = time.time()

final_H_8 = []
for trial in range(n_samples_8):
    bits = random.randint(0, N_8 - 1)
    current = bits
    steps = 0
    while steps < 100:  # safety limit
        h_current = get_H_8(current)
        best_k = -1
        best_delta = 0
        for k in range(m):
            neighbor = current ^ (1 << k)
            h_neighbor = get_H_8(neighbor)
            delta = h_neighbor - h_current
            if delta > best_delta:
                best_delta = delta
                best_k = k
        if best_k == -1:
            break
        current = current ^ (1 << best_k)
        steps += 1

    final_H_8.append(get_H_8(current))

    if trial % 100 == 0 and trial > 0:
        elapsed = time.time() - t0
        print(f"    {trial}/{n_samples_8} done ({trial/elapsed:.0f}/s)")

elapsed = time.time() - t0
print(f"  Done in {elapsed:.1f}s")

unique_H_8 = sorted(set(final_H_8))
max_H_8 = max(final_H_8)
print(f"  Unique final H values: {unique_H_8}")
print(f"  Max H: {max_H_8}")
print(f"  Fraction at max: {sum(1 for h in final_H_8 if h == max_H_8)}/{n_samples_8}")

if len(unique_H_8) > 1:
    print(f"  SPURIOUS LOCAL MAXIMA at n=8 (even)")
    for h_val in unique_H_8:
        count = sum(1 for h in final_H_8 if h == h_val)
        print(f"    H={h_val}: {count} ({100*count/n_samples_8:.1f}%)")

print("\n\nDONE — landscape_n7_sampling.py complete")
