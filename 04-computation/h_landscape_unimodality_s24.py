#!/usr/bin/env python3
"""
h_landscape_unimodality_s24.py — opus-2026-04-05-S24

VERIFY: Is the H-landscape unimodal under single-arc-flip gradient ascent?

Conjecture: For every tournament T on n vertices, greedy H-gradient ascent
(flip the arc that maximally increases H) reaches the GLOBAL H-maximum.
Equivalently: there are no local maxima of H other than the global max.

If true, this means the H-landscape has a single "mountain" — a remarkable
structural property of tournaments.

Method: exhaustive at n=4,5; heavy sampling at n=6.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
import time
import sys

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def adj_to_bits(A, n):
    bits = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] == 1:
                bits |= (1 << idx)
            idx += 1
    return bits

def count_ham_paths(A, n):
    m = 1 << n
    dp = [[0]*n for _ in range(m)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(m):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            if not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u] == 1:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = m - 1
    return sum(dp[full][v] for v in range(n))

def flip_arc(A, i, j):
    B = A.copy()
    B[i][j] = 1 - B[i][j]
    B[j][i] = 1 - B[j][i]
    return B

def h_gradient_ascent(A, n, max_steps=100):
    """Perform greedy H-gradient ascent. Return (final_H, steps, trajectory)."""
    current = A.copy()
    H = count_ham_paths(current, n)
    trajectory = [H]
    for step in range(max_steps):
        best_delta = 0
        best_ij = None
        for i in range(n):
            for j in range(i+1, n):
                B = flip_arc(current, i, j)
                new_H = count_ham_paths(B, n)
                delta = new_H - H
                if delta > best_delta:
                    best_delta = delta
                    best_ij = (i, j)
        if best_ij is None:
            # Local maximum reached
            break
        current = flip_arc(current, best_ij[0], best_ij[1])
        H += best_delta
        trajectory.append(H)
    return H, step + 1, trajectory

def find_all_local_maxima(n):
    """Find all LOCAL maxima of H under single-arc flips. Exhaustive."""
    m = n * (n - 1) // 2
    total = 2 ** m

    local_max_H = Counter()
    global_max_H = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)
        global_max_H = max(global_max_H, H)

        # Check if this is a local max: no single arc flip increases H
        is_local_max = True
        for i in range(n):
            for j in range(i+1, n):
                B = flip_arc(A, i, j)
                new_H = count_ham_paths(B, n)
                if new_H > H:
                    is_local_max = False
                    break
            if not is_local_max:
                break

        if is_local_max:
            local_max_H[H] += 1

    return local_max_H, global_max_H

def h_gradient_ascent_exhaustive(n):
    """Run gradient ascent from EVERY tournament. Check all converge to global max."""
    m = n * (n - 1) // 2
    total = 2 ** m

    terminal_H = Counter()
    max_steps_needed = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        final_H, steps, _ = h_gradient_ascent(A, n)
        terminal_H[final_H] += 1
        max_steps_needed = max(max_steps_needed, steps)

    return terminal_H, max_steps_needed

print("=" * 70)
print("H-LANDSCAPE UNIMODALITY ANALYSIS")
print("=" * 70)

for n in [4, 5]:
    m = n * (n - 1) // 2
    total = 2 ** m
    print(f"\n{'='*50}")
    print(f"n = {n}, m = {m}, total tournaments = {total}")
    print(f"{'='*50}")

    # Find all local maxima (exhaustive)
    t0 = time.time()
    local_max_H, global_max_H = find_all_local_maxima(n)
    t1 = time.time()
    print(f"\n  Local maxima search: {t1-t0:.1f}s")
    print(f"  Global max H = {global_max_H}")
    print(f"  Local maxima by H value: {dict(local_max_H)}")

    non_global = {h: c for h, c in local_max_H.items() if h < global_max_H}
    if non_global:
        print(f"  *** NON-GLOBAL LOCAL MAXIMA FOUND: {non_global} ***")
    else:
        print(f"  ALL local maxima are at H = {global_max_H}. H-LANDSCAPE IS UNIMODAL.")

    # Verify: gradient ascent from every tournament reaches global max
    t0 = time.time()
    terminal_H, max_steps = h_gradient_ascent_exhaustive(n)
    t1 = time.time()
    print(f"\n  Gradient ascent from all {total} tournaments: {t1-t0:.1f}s")
    print(f"  Terminal H values: {dict(terminal_H)}")
    print(f"  Max steps needed: {max_steps}")

    if len(terminal_H) == 1 and list(terminal_H.keys())[0] == global_max_H:
        print(f"  CONFIRMED: Every tournament reaches H = {global_max_H} by gradient ascent.")
    else:
        print(f"  *** GRADIENT ASCENT DOES NOT ALWAYS REACH GLOBAL MAX ***")

# n=6: sample-based analysis
n = 6
m = n * (n - 1) // 2
total = 2 ** m
print(f"\n{'='*50}")
print(f"n = {n}, m = {m}, total tournaments = {total}")
print(f"{'='*50}")

import random

# First find global max by sampling
print("\n  Finding global max by sampling...")
max_H_seen = 0
sample_size = 5000
for _ in range(sample_size):
    bits = random.randint(0, total - 1)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    max_H_seen = max(max_H_seen, H)

# Also find it by gradient ascent from random starts
t0 = time.time()
terminal_H = Counter()
num_trials = 2000
max_steps_seen = 0
for trial in range(num_trials):
    bits = random.randint(0, total - 1)
    A = bits_to_adj(bits, n)
    final_H, steps, trajectory = h_gradient_ascent(A, n)
    terminal_H[final_H] += 1
    max_steps_seen = max(max_steps_seen, steps)
t1 = time.time()

print(f"\n  Gradient ascent from {num_trials} random tournaments: {t1-t0:.1f}s")
print(f"  Terminal H values: {dict(terminal_H)}")
print(f"  Max steps needed: {max_steps_seen}")

global_max_n6 = max(terminal_H.keys())
non_max = {h: c for h, c in terminal_H.items() if h < global_max_n6}
if non_max:
    print(f"\n  *** NON-GLOBAL TERMINALS FOUND: {non_max} ***")
    print(f"  This means local maxima exist at n=6!")
    # Analyze the local maxima
    for h_val, count in sorted(non_max.items()):
        print(f"    H={h_val}: {count} tournaments ({count/num_trials*100:.1f}%)")
else:
    print(f"\n  All {num_trials} trials reached H={global_max_n6}.")

# Now check: are there local maxima by direct search on random tournaments?
print("\n  Direct local max search on 5000 random tournaments...")
local_max_non_global = 0
local_max_H_vals = Counter()
for _ in range(5000):
    bits = random.randint(0, total - 1)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    is_local_max = True
    for i in range(n):
        for j in range(i+1, n):
            B = flip_arc(A, i, j)
            new_H = count_ham_paths(B, n)
            if new_H > H:
                is_local_max = False
                break
        if not is_local_max:
            break

    if is_local_max and H < global_max_n6:
        local_max_non_global += 1
        local_max_H_vals[H] += 1

print(f"  Non-global local maxima found: {local_max_non_global}/5000")
if local_max_H_vals:
    print(f"  H-values of local maxima: {dict(local_max_H_vals)}")
else:
    print(f"  STRONG EVIDENCE: H-landscape is unimodal at n=6 too.")

# Analysis: number of H-maximal tournaments
print(f"\n  H-max tournaments at each n:")
for nn in [3, 4, 5]:
    mm = nn * (nn - 1) // 2
    tt = 2 ** mm
    max_H = 0
    max_count = 0
    for bits in range(tt):
        A = bits_to_adj(bits, nn)
        H = count_ham_paths(A, nn)
        if H > max_H:
            max_H = H
            max_count = 1
        elif H == max_H:
            max_count += 1
    print(f"    n={nn}: max(H)={max_H}, count={max_count}/{tt} ({max_count/tt*100:.1f}%)")
