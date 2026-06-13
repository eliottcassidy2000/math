#!/usr/bin/env python3
"""
n=9 tournament path homology: hunting for β_5

Based on the pattern:
- n≤6: only β_1 appears
- n=7: β_3 appears (8.4%)
- n=8: β_3 (16%), β_4 (rare)
- n=9: β_5 predicted?

Also: β_2=0 verification at n=8,9 and deeper geometric analysis.
"""
import numpy as np
from itertools import combinations
from collections import Counter
import random
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    path_betti_numbers, count_3cycles, ham_path_count
)

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

# ===== n=9 exploration =====
print("=" * 70)
print("n=9 TOURNAMENT PATH HOMOLOGY")
print("=" * 70)

n = 9
betti_dist = Counter()
phase_H = {'P': [], 'C': [], 'S3': [], 'S5': [], 'other': []}

for trial in range(30):
    A = random_tournament(n)
    # Only compute up to dim 8 (full for n=9)
    # This will be slow — each computation at n=9 may take minutes
    betti = path_betti_numbers(A, n, max_dim=8)
    bt = tuple(betti)
    betti_dist[bt] += 1

    t3 = count_3cycles(A, n)
    H = ham_path_count(A, n)

    # Classify phase
    if betti[5] > 0:
        phase = 'S5'
    elif betti[3] > 0:
        phase = 'S3'
    elif betti[1] > 0:
        phase = 'C'
    elif all(betti[p] == 0 for p in range(1, len(betti))):
        phase = 'P'
    else:
        phase = 'other'

    phase_H[phase].append(H)
    print(f"  Trial {trial}: t3={t3}, H={H}, β={betti}, phase={phase}", flush=True)

print(f"\n  Betti distribution (n={n}, {sum(betti_dist.values())} samples):")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")

print(f"\n  Phase distribution:")
for phase in ['P', 'C', 'S3', 'S5', 'other']:
    if phase_H[phase]:
        print(f"    {phase}: {len(phase_H[phase])} tournaments, mean H={np.mean(phase_H[phase]):.1f}")

# Check β_2=0
any_b2 = any(bt[2] > 0 for bt in betti_dist.keys())
print(f"\n  β_2 > 0 ever? {any_b2}")

# Check mutual exclusion
for bt in betti_dist:
    nonzero = [p for p in range(1, len(bt)) if bt[p] > 0]
    if len(nonzero) > 1:
        print(f"  MULTIPLE nonzero: β={list(bt)}, dims={nonzero}")

print("\nDone.")
