#!/usr/bin/env python3
"""
beta2_min_surplus_n8_heavy.py - Heavy sampling to pin down min surplus at n=8

From 3000 samples: min surplus = 27.
Need more samples to confirm this is the true minimum.

Also: analyze WHICH tournaments achieve min surplus.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os, random
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def compute_surplus_data(A, n):
    allowed = {}
    for p in range(5):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if not allowed[p]:
            break

    omega_basis = {}
    for p in range(4):
        if p not in allowed or not allowed[p]:
            omega_basis[p] = np.zeros((0, 0))
            continue
        basis = compute_omega_basis(A, n, p, allowed[p],
                                     allowed[p-1] if p >= 1 and p-1 in allowed else [])
        omega_basis[p] = basis

    dim_O2 = omega_basis[2].shape[1] if omega_basis[2].ndim == 2 else 0
    dim_O3 = omega_basis[3].shape[1] if omega_basis[3].ndim == 2 else 0

    if dim_O2 == 0:
        Z2 = 0
    else:
        bd2 = build_full_boundary_matrix(allowed[2], allowed.get(1, []))
        bd2_om = bd2 @ omega_basis[2]
        S_v = np.linalg.svd(bd2_om, compute_uv=False)
        rk2 = int(np.sum(np.abs(S_v) > 1e-8))
        Z2 = dim_O2 - rk2

    surplus = dim_O3 - Z2
    return surplus, dim_O2, dim_O3, Z2


# Also try specific "extreme" tournaments to find lower surplus
def near_transitive(n, num_flips):
    """Start from transitive and flip a few arcs."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    for _ in range(num_flips):
        i = random.randint(0, n-1)
        j = random.randint(0, n-1)
        while i == j:
            j = random.randint(0, n-1)
        if A[i][j] == 1:
            A[i][j] = 0
            A[j][i] = 1
    return A

def tournament_with_score(n, target_scores):
    """Try to build a tournament with given score sequence (greedy)."""
    from itertools import permutations
    A = [[0]*n for _ in range(n)]
    remaining = list(target_scores)
    for i in range(n):
        for j in range(i+1, n):
            if remaining[i] > 0:
                A[i][j] = 1
                remaining[i] -= 1
            else:
                A[j][i] = 1
                remaining[j] -= 1
    return A

n = 8
print("=" * 70)
print(f"MINIMUM SURPLUS SEARCH AT n={n}")
print("=" * 70)

N_SAMPLES = 20000
t0 = time.time()
min_surplus = float('inf')
min_examples = []

for trial in range(N_SAMPLES):
    if trial % 5000 == 0 and trial > 0:
        dt = time.time() - t0
        print(f"  ... {trial}/{N_SAMPLES} ({dt:.0f}s), min={min_surplus}")

    A = random_tournament(n)
    surplus, O2, O3, Z2 = compute_surplus_data(A, n)

    if surplus < min_surplus:
        min_surplus = surplus
        score = tuple(sorted(sum(row) for row in A))
        t3 = 0
        for i in range(n):
            for j in range(n):
                if i==j: continue
                for k in range(n):
                    if k==i or k==j: continue
                    if A[i][j] and A[j][k] and A[k][i]:
                        t3 += 1
        t3 //= 3
        min_examples = [(score, O2, O3, Z2, t3)]
        print(f"  New min surplus={surplus}: O2={O2}, O3={O3}, Z2={Z2}, score={score}, t3={t3}")
    elif surplus == min_surplus:
        score = tuple(sorted(sum(row) for row in A))
        if score not in [x[0] for x in min_examples]:
            t3 = 0
            for i in range(n):
                for j in range(n):
                    if i==j: continue
                    for k in range(n):
                        if k==i or k==j: continue
                        if A[i][j] and A[j][k] and A[k][i]:
                            t3 += 1
            t3 //= 3
            min_examples.append((score, O2, O3, Z2, t3))
            print(f"  Another min surplus={surplus}: O2={O2}, O3={O3}, Z2={Z2}, score={score}, t3={t3}")

dt = time.time() - t0
print(f"\n  n={n}: min surplus = {min_surplus} ({dt:.0f}s)")
print(f"  Examples at minimum:")
for score, O2, O3, Z2, t3 in min_examples:
    print(f"    score={score}, O2={O2}, O3={O3}, Z2={Z2}, t3={t3}")

# Try near-transitive tournaments (few flips from transitive)
print(f"\n{'='*70}")
print("NEAR-TRANSITIVE TOURNAMENTS")
print(f"{'='*70}")

for nf in [1, 2, 3, 4, 5, 10]:
    min_s = float('inf')
    for _ in range(2000):
        A = near_transitive(n, nf)
        surplus, O2, O3, Z2 = compute_surplus_data(A, n)
        if surplus < min_s:
            min_s = surplus
    print(f"  {nf} flips from transitive: min surplus = {min_s}")

# Try tournaments with extreme score sequences
print(f"\n{'='*70}")
print("EXTREME SCORE TOURNAMENTS")
print(f"{'='*70}")

# Near-regular (all degrees close)
for _ in range(2000):
    A = random_tournament(n)
    scores = sorted(sum(row) for row in A)
    if max(scores) - min(scores) <= 1:  # Almost regular
        surplus, O2, O3, Z2 = compute_surplus_data(A, n)
        print(f"  Near-regular {scores}: surplus={surplus}")
        break

# Very unequal scores
for _ in range(2000):
    A = random_tournament(n)
    scores = sorted(sum(row) for row in A)
    if scores[0] == 0 and scores[-1] == n-1:
        surplus, O2, O3, Z2 = compute_surplus_data(A, n)
        if surplus < 30:
            print(f"  Extreme scores {scores}: surplus={surplus}, O2={O2}, O3={O3}, Z2={Z2}")

print("\nDone.")
