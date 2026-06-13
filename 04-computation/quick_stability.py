#!/usr/bin/env python3
"""Quick stability test for circulant path homology"""
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import path_betti_numbers, circulant_digraph
import random

random.seed(42)

def is_prime(n):
    if n < 2: return False
    return all(n % i != 0 for i in range(2, int(n**0.5)+1))

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

# Stability: C_n^{1,3}
print("C_n^{1,3} stability:")
for n in range(5, 24):
    A = circulant_digraph(n, [1, 3])
    betti = path_betti_numbers(A, n, max_dim=min(n-1, 4))
    p = "P" if is_prime(n) else " "
    print(f"  n={n:2d} {p}: β={betti}", flush=True)

# Stability: C_n^{1,2}
print("\nC_n^{1,2} stability:")
for n in range(4, 18):
    A = circulant_digraph(n, [1, 2])
    betti = path_betti_numbers(A, n, max_dim=min(n-1, 4))
    p = "P" if is_prime(n) else " "
    print(f"  n={n:2d} {p}: β={betti}", flush=True)

# Complement duality: n=5 exhaustive
print("\nComplement duality n=5 (exhaustive):")
n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)
matches = 0
for mask in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (mask >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    A_op = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                A_op[i][j] = 1 - A[i][j]
    b1 = path_betti_numbers(A, n, max_dim=n-1)
    b2 = path_betti_numbers(A_op, n, max_dim=n-1)
    if b1 == b2:
        matches += 1
print(f"  β(T) = β(T^op): {matches}/1024", flush=True)

# n=8 tournament sample (just 30)
print("\nn=8 tournaments (30 samples):")
n = 8
betti_dist = {}
for trial in range(30):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=7)
    bt = tuple(betti)
    betti_dist[bt] = betti_dist.get(bt, 0) + 1
    if trial < 3:
        print(f"  Trial {trial}: β={betti}", flush=True)

print(f"\n  Betti distribution:")
for bt in sorted(betti_dist.keys()):
    print(f"    β={list(bt)}: {betti_dist[bt]}")

print("\nDone.")
