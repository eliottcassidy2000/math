#!/usr/bin/env python3
"""Sample test of β₁ deletion bound at n=7 (random tournaments)."""

import numpy as np
import random
import sys
sys.path.insert(0, '/home/e/Documents/claude/math/04-computation')

# Import core functions from the analysis script
from beta2_filler_vertex_analysis import (
    path_betti_numbers, delete_vertex, count_3cycles, find_3cycles_through
)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

n = 7
N_SAMPLES = 2000
random.seed(42)

print(f"Sampling {N_SAMPLES} random tournaments at n={n}")
print(f"Testing: β₁(T)=0 ⟹ Σ_v β₁(T\\v) ≤ 3")
print("="*60)

from collections import Counter
sum_dist = Counter()
max_sum = 0
total_b1_zero = 0
total_tested = 0
crit_size_dist = Counter()

for trial in range(N_SAMPLES):
    A = random_tournament(n)
    total_tested += 1

    betti = path_betti_numbers(A, n, max_dim=1)
    if betti[1] != 0:
        continue

    total_b1_zero += 1
    crit = []
    for v in range(n):
        B, _ = delete_vertex(A, n, v)
        b_del = path_betti_numbers(B, n-1, max_dim=1)
        if b_del[1] > 0:
            crit.append(v)

    s = len(crit)
    sum_dist[s] += 1
    crit_size_dist[s] += 1

    if s > max_sum:
        max_sum = s
        print(f"\n  NEW MAX at trial {trial}: Σ = {s}, crit = {crit}")
        print(f"  Score: {sorted([sum(A[i]) for i in range(n)])}, t₃ = {count_3cycles(A, n)}")

        # Check 3-cycle structure
        if len(crit) >= 3:
            from itertools import combinations
            for triple in combinations(crit, 3):
                cyc = find_3cycles_through(A, n, list(triple))
                if cyc:
                    print(f"  Triple {triple} in 3-cycle: {cyc}")

    if s > 3:
        print(f"\n  *** COUNTEREXAMPLE at trial {trial}! Σ = {s} > 3 ***")
        print(f"  crit = {crit}")
        for i in range(n):
            row = [A[i][j] for j in range(n)]
            print(f"  A[{i}] = {row}")

    if (trial+1) % 500 == 0:
        print(f"  ... {trial+1}/{N_SAMPLES} done, {total_b1_zero} with β₁=0, max Σ = {max_sum}")

print(f"\n{'='*60}")
print(f"RESULTS: n={n}, {N_SAMPLES} samples")
print(f"  β₁(T)=0: {total_b1_zero} / {total_tested}")
print(f"  Σ distribution: {dict(sorted(sum_dist.items()))}")
print(f"  MAX Σ = {max_sum}")
print(f"  Claim (Σ ≤ 3): {'HOLDS' if max_sum <= 3 else 'FAILS!'}")
