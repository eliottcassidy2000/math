#!/usr/bin/env python3
"""
STABILITY TEST: Path homology of C_n^S as n varies

Tang-Yau Theorem 4.5: For "no wrap-around" S, β stabilizes for large primes.
Test this computationally with specific S values.

Also test: Betti numbers of complete tournaments at n=8.
"""
import numpy as np
from itertools import combinations
import random
import sys
sys.path.insert(0, '04-computation')
from path_homology_v2 import (
    path_betti_numbers, circulant_digraph, count_3cycles, ham_path_count
)

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

# ===== Test 1: S={1,3} stability =====
print("=" * 70)
print("STABILITY TEST: C_n^{1,3}")
print("=" * 70)

for n in range(5, 24):
    A = circulant_digraph(n, [1, 3])
    betti = path_betti_numbers(A, n, max_dim=min(n-1, 5))
    is_prime = all(n % i != 0 for i in range(2, int(n**0.5)+1)) and n > 1
    p_str = " (prime)" if is_prime else ""
    print(f"  n={n:2d}{p_str}: β={betti}")

# ===== Test 2: S={1,2} stability =====
print("\n" + "=" * 70)
print("STABILITY TEST: C_n^{1,2}")
print("=" * 70)

for n in range(4, 20):
    A = circulant_digraph(n, [1, 2])
    betti = path_betti_numbers(A, n, max_dim=min(n-1, 5))
    is_prime = all(n % i != 0 for i in range(2, int(n**0.5)+1)) and n > 1
    p_str = " (prime)" if is_prime else ""
    print(f"  n={n:2d}{p_str}: β={betti}")

# ===== Test 3: S={1,2,4} stability =====
print("\n" + "=" * 70)
print("STABILITY TEST: C_n^{1,2,4}")
print("=" * 70)

for n in range(5, 18):
    A = circulant_digraph(n, [1, 2, 4])
    betti = path_betti_numbers(A, n, max_dim=min(n-1, 5))
    is_prime = all(n % i != 0 for i in range(2, int(n**0.5)+1)) and n > 1
    p_str = " (prime)" if is_prime else ""
    print(f"  n={n:2d}{p_str}: β={betti}")

# ===== Test 4: n=8 tournaments =====
print("\n\n" + "=" * 70)
print("n=8 TOURNAMENTS: HUNTING HIGHER BETTI NUMBERS")
print("=" * 70)

n = 8
betti_dist = {}
for trial in range(200):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=7)
    bt = tuple(betti)
    betti_dist[bt] = betti_dist.get(bt, 0) + 1

    if trial < 5:
        t3 = count_3cycles(A, n)
        print(f"  Trial {trial}: t3={t3}, β={betti}")

print(f"\n  Betti distribution (n=8, 200 samples):")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")

# Any new topology types?
for bt in betti_dist:
    b = list(bt)
    if any(b[p] > 0 for p in range(4, len(b))):
        print(f"\n  *** HIGHER BETTI: β={b} ***")

# ===== Test 5: Complement duality (quick, small n) =====
print("\n\n" + "=" * 70)
print("COMPLEMENT DUALITY (quick test)")
print("=" * 70)

# n=5 exhaustive
n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)
matches = 0
total = 0
for mask in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (mask >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    # Complement
    A_op = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                A_op[i][j] = 1 - A[i][j]

    b1 = path_betti_numbers(A, n, max_dim=n-1)
    b2 = path_betti_numbers(A_op, n, max_dim=n-1)
    total += 1
    if b1 == b2:
        matches += 1

print(f"  n=5 exhaustive: β(T) = β(T^op) for {matches}/{total} tournaments")

# n=6 sample
n = 6
matches = 0
total = 0
for trial in range(200):
    A = random_tournament(n)
    A_op = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                A_op[i][j] = 1 - A[i][j]
    b1 = path_betti_numbers(A, n, max_dim=5)
    b2 = path_betti_numbers(A_op, n, max_dim=5)
    total += 1
    if b1 == b2:
        matches += 1

print(f"  n=6 (200 samples): β(T) = β(T^op) for {matches}/{total}")

print("\n\nDone.")
