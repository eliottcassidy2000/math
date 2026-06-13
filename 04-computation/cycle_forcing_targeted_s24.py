#!/usr/bin/env python3
"""
cycle_forcing_targeted_s24.py — opus-2026-04-05-S24

Targeted verification of cycle-forcing thresholds at large n.

Random sampling misses near-transitive tournaments (very low c₃).
This script CONSTRUCTS tournaments with specific small c₃ values
and checks whether c₅, c₇, c₉ are forced > 0.

Strategy: Start from the transitive tournament and flip a controlled
number of arcs to create tournaments with specific c₃ values.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter
import random
import time
import sys

sys.stdout.reconfigure(encoding='utf-8')

PERMS = {}
for k in [3, 5, 7, 9]:
    PERMS[k] = list(permutations(range(k - 1)))

def transitive_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def score_c3(A, n):
    total = n * (n-1) * (n-2) // 6
    for v in range(n):
        s = sum(A[v])
        total -= s * (s-1) // 2
    return total

def count_kcycles(A, n, k):
    count = 0
    perms = PERMS[k]
    for combo in combinations(range(n), k):
        start = combo[0]
        rest = combo[1:]
        for p in perms:
            prev = start
            valid = True
            for idx in p:
                nxt = rest[idx]
                if not A[prev][nxt]:
                    valid = False
                    break
                prev = nxt
            if valid and A[prev][start]:
                count += 1
    return count

def flip_arc(A, i, j):
    B = [row[:] for row in A]
    B[i][j] = 1 - B[i][j]
    B[j][i] = 1 - B[j][i]
    return B

def construct_with_target_c3(n, target_c3, attempts=10000):
    """Construct tournaments with specific c₃ by random walks from transitive."""
    found = []
    for _ in range(attempts):
        A = transitive_tournament(n)
        # Randomly flip arcs until we hit or pass target_c3
        for _ in range(n * n):  # plenty of flips
            c3 = score_c3(A, n)
            if c3 == target_c3:
                found.append([row[:] for row in A])
                break
            if c3 < target_c3:
                # Need more 3-cycles: flip a random arc
                i = random.randint(0, n-2)
                j = random.randint(i+1, n-1)
                A = flip_arc(A, i, j)
            else:
                # Too many: restart
                break
        if len(found) >= 100:
            break
    return found

print("=" * 70)
print("TARGETED CYCLE-FORCING VERIFICATION AT LARGE n")
print("=" * 70)

# Verify the 3*(2^k - 1) pattern
print("\n--- Verifying c₃ → c_{2k+1} threshold = 3*(2^k - 1) at n = 2k+1 ---")
print()

for k in range(1, 5):
    cycle_len = 2*k + 1
    n = 2*k + 1
    predicted_threshold = 3 * (2**k - 1)

    if n > 11:
        break

    m = n * (n-1) // 2
    print(f"k={k}: c₃ → c_{cycle_len} at n={n}, predicted threshold = {predicted_threshold}")

    # Check at threshold - 1: should be possible to have c_{cycle_len} = 0
    for target_c3 in [predicted_threshold - 1, predicted_threshold]:
        tournaments = construct_with_target_c3(n, target_c3, attempts=50000)
        if not tournaments:
            print(f"  c₃={target_c3}: could NOT construct any tournaments (might not exist)")
            continue

        n_checked = len(tournaments)
        min_cq = float('inf')
        max_cq = 0
        zero_count = 0

        for A in tournaments:
            cq = count_kcycles(A, n, cycle_len)
            min_cq = min(min_cq, cq)
            max_cq = max(max_cq, cq)
            if cq == 0:
                zero_count += 1

        label = "BELOW threshold" if target_c3 < predicted_threshold else "AT threshold"
        print(f"  c₃={target_c3} ({label}): {n_checked} tournaments, "
              f"c_{cycle_len} range [{min_cq}, {max_cq}], "
              f"c_{cycle_len}=0 in {zero_count}/{n_checked}")

    print()

# Also verify for n > 2k+1 (threshold may change)
print("\n--- How does c₃ → c₅ threshold change with n? ---")
print()

for n in range(5, 12):
    m = n * (n-1) // 2
    if n > 9 and m > 40:
        print(f"n={n}: skipping (too large for cycle counting)")
        continue

    # Find the threshold by testing c₃ = 0, 1, 2, 3, 4, 5, ...
    threshold = None
    for target_c3 in range(0, min(20, n*(n-1)*(n-2)//6)):
        tournaments = construct_with_target_c3(n, target_c3, attempts=20000)
        if not tournaments:
            continue

        has_zero_c5 = False
        for A in tournaments:
            c5 = count_kcycles(A, n, 5)
            if c5 == 0:
                has_zero_c5 = True
                break

        if has_zero_c5:
            last_zero = target_c3
        else:
            if threshold is None:
                threshold = target_c3
            # Keep checking to make sure it's truly monotone

    if threshold is not None:
        print(f"  n={n}: c₃ → c₅ threshold = {threshold} (max c₃ with c₅=0: {last_zero})")
    else:
        print(f"  n={n}: threshold not found in range")

print("\n--- How does c₃ → c₇ threshold change with n? ---")
print()

for n in range(7, 11):
    threshold = None
    last_zero = None
    for target_c3 in range(0, 20):
        tournaments = construct_with_target_c3(n, target_c3, attempts=20000)
        if not tournaments:
            continue

        has_zero_c7 = False
        for A in tournaments:
            c7 = count_kcycles(A, n, 7)
            if c7 == 0:
                has_zero_c7 = True
                break

        if has_zero_c7:
            last_zero = target_c3
        else:
            if threshold is None:
                threshold = target_c3

    if threshold is not None:
        print(f"  n={n}: c₃ → c₇ threshold = {threshold} (max c₃ with c₇=0: {last_zero})")
    else:
        print(f"  n={n}: threshold not found in range tested")
