#!/usr/bin/env python3
"""
beta3_analysis.py - Analyze tournaments with beta3 > 0

DISCOVERY: beta3 can be nonzero for tournaments!
At n=6: beta3=1 for ~1.35% of tournaments.
Betti vector: [1, 0, 0, 1, 0, 0]

Questions:
1. Which tournaments have beta3 > 0?
2. What is their score structure?
3. Can beta3 > 1?
4. At what n does beta3 first appear?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, path_betti_numbers
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


def count_c3(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    c3 += 1
    return c3


# ============================================================
# n=4,5: Check for beta3 (should be 0 based on n=5 exhaustive)
# ============================================================
for n in [4, 5]:
    print(f"n={n}: checking all {1 << (n*(n-1)//2)} tournaments for beta3 > 0")
    total = 1 << (n*(n-1)//2)
    b3_nonzero = 0
    for bits in range(total):
        A = build_adj(n, bits)
        betti = path_betti_numbers(A, n, max_dim=n-1)
        if len(betti) > 3 and betti[3] > 0:
            b3_nonzero += 1
    print(f"  beta3 > 0: {b3_nonzero}/{total}\n")


# ============================================================
# n=6: EXHAUSTIVE beta3 check
# ============================================================
print("=" * 70)
print("BETA3 ANALYSIS: n=6 (exhaustive)")
print("=" * 70)

n = 6
total = 1 << (n*(n-1)//2)

b3_examples = []
betti_full_dist = Counter()
t0 = time.time()

for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s), beta3>0: {len(b3_examples)}")

    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=n-1)
    bt = tuple(betti)
    betti_full_dist[bt] += 1

    if betti[3] > 0:
        c3 = count_c3(A, n)
        scores = tuple(sorted([sum(row) for row in A]))
        b3_examples.append((bits, scores, c3, list(betti)))

dt = time.time() - t0
print(f"\nDone in {dt:.0f}s")

print(f"\nFull Betti vector distribution:")
for bt, cnt in sorted(betti_full_dist.items()):
    print(f"  {list(bt)}: {cnt}")

print(f"\nbeta3 > 0: {len(b3_examples)}/{total}")
if b3_examples:
    print(f"\nExamples with beta3 > 0:")
    score_dist = Counter()
    for bits, scores, c3, betti in b3_examples:
        score_dist[scores] += 1

    print(f"\n  Score distribution:")
    for sc, cnt in sorted(score_dist.items()):
        print(f"    {sc}: {cnt}")

    print(f"\n  First 10 examples:")
    for bits, scores, c3, betti in b3_examples[:10]:
        print(f"    bits={bits}, scores={scores}, c3={c3}, betti={betti}")

    # Max beta3
    max_b3 = max(b[3] for _, _, _, b in b3_examples)
    print(f"\n  Max beta3: {max_b3}")


# ============================================================
# Analyze structure of beta3=1 tournaments at n=6
# ============================================================
if b3_examples:
    print(f"\n{'='*70}")
    print("STRUCTURE OF BETA3=1 TOURNAMENTS")
    print("=" * 70)

    # Pick a specific example
    bits0, scores0, c3_0, betti0 = b3_examples[0]
    A = build_adj(n, bits0)

    print(f"\nExample: bits={bits0}, scores={scores0}, c3={c3_0}")
    print(f"Adjacency matrix:")
    for i in range(n):
        print(f"  {A[i]}")

    # Compute all Om dimensions
    for p in range(n):
        ap = enumerate_allowed_paths(A, n, p)
        apm1 = enumerate_allowed_paths(A, n, p-1) if p > 0 else []
        omp = compute_omega_basis(A, n, p, ap, apm1)
        dim_om = omp.shape[1] if omp.ndim == 2 else 0
        print(f"  |A_{p}|={len(ap)}, dim(Om_{p})={dim_om}")

    # Check: is the tournament a specific well-known structure?
    # Count cycles of all lengths
    def count_directed_cycles(A, n, k):
        """Count directed k-cycles."""
        from itertools import permutations
        count = 0
        for perm in permutations(range(n), k):
            if all(A[perm[i]][perm[(i+1)%k]] for i in range(k)):
                count += 1
        return count // k

    for k in [3, 4, 5, 6]:
        ck = count_directed_cycles(A, n, k)
        print(f"  c_{k} = {ck}")

    # Is it a regular tournament?
    is_reg = all(sum(row) == (n-1)//2 for row in A)
    print(f"  Regular: {is_reg}")

    # Check if beta3=1 is related to the Paley construction
    # At n=6 there's no Paley (6 is not prime)


# ============================================================
# n=7: sampled beta3 check
# ============================================================
print(f"\n{'='*70}")
print("BETA3 ANALYSIS: n=7 (sampled)")
print("=" * 70)

import random
random.seed(42)
n = 7
n_samples = 2000

b3_nonzero7 = 0
b3_dist7 = Counter()
t0 = time.time()

for trial in range(n_samples):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=4)
    if betti[3] > 0:
        b3_nonzero7 += 1
        b3_dist7[betti[3]] += 1
    if len(betti) > 4 and betti[4] > 0:
        scores = tuple(sorted([sum(row) for row in A]))
        print(f"  beta4>0: bits={bits}, scores={scores}, betti={betti}")

    if trial % 500 == 499:
        dt = time.time() - t0
        print(f"  ... {trial+1}/{n_samples} ({dt:.0f}s)")

dt = time.time() - t0
print(f"\nn=7: {n_samples} samples in {dt:.0f}s")
print(f"beta3 > 0: {b3_nonzero7}/{n_samples}")
print(f"beta3 distribution: {dict(b3_dist7)}")


print("\n\nDone.")
