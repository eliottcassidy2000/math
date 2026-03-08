#!/usr/bin/env python3
"""
n=9 tournament path homology: FAST version.
Only compute up to max_dim=5 (enough to detect β_5).
No need for dim 6,7,8 — those won't appear in new dimensions.
"""
import numpy as np
from collections import Counter
import random
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import path_betti_numbers, count_3cycles, ham_path_count

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

print("=" * 70)
print("n=9 FAST: β up to dim 5 (hunting β_5)")
print("=" * 70)

n = 9
betti_dist = Counter()
NUM_TRIALS = 100

for trial in range(NUM_TRIALS):
    A = random_tournament(n)
    # max_dim=5 is enough to see β_5
    betti = path_betti_numbers(A, n, max_dim=5)
    bt = tuple(betti)
    betti_dist[bt] += 1

    t3 = count_3cycles(A, n)

    # Classify
    nonzero_dims = [p for p in range(1, len(betti)) if betti[p] > 0]
    phase = "P"
    if 5 in nonzero_dims: phase = "S5"
    elif 3 in nonzero_dims: phase = "S3"
    elif 1 in nonzero_dims: phase = "C"

    print(f"  Trial {trial}: t3={t3}, β={betti}, phase={phase}", flush=True)

print(f"\n  Betti distribution (n={n}, {sum(betti_dist.values())} samples):")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")

# Check β_2=0
any_b2 = any(bt[2] > 0 for bt in betti_dist.keys())
print(f"\n  β_2 > 0 ever? {any_b2}")
any_b4 = any(len(bt) > 4 and bt[4] > 0 for bt in betti_dist.keys())
print(f"  β_4 > 0 ever? {any_b4}")
any_b5 = any(len(bt) > 5 and bt[5] > 0 for bt in betti_dist.keys())
print(f"  β_5 > 0 ever? {any_b5}")

# Mutual exclusion check
for bt in betti_dist:
    nonzero = [p for p in range(1, len(bt)) if bt[p] > 0]
    if len(nonzero) > 1:
        print(f"  MULTIPLE nonzero: β={list(bt)}, dims={nonzero}")

print("\nDone.")
