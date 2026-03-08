#!/usr/bin/env python3
"""
n9_max_betti_quick.py — Quickly find ONE H=3357 n=9 tournament and compute its Betti

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, random, os
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full])

def random_regular_tournament(n):
    assert n % 2 == 1
    target = (n - 1) // 2
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    for _ in range(1000):
        scores = [sum(A[i]) for i in range(n)]
        if max(abs(s - target) for s in scores) == 0:
            break
        high_v = max(range(n), key=lambda v: scores[v])
        low_v = min(range(n), key=lambda v: scores[v])
        if high_v == low_v:
            break
        if A[high_v][low_v]:
            A[high_v][low_v] = 0
            A[low_v][high_v] = 1
        else:
            candidates = [w for w in range(n) if w != high_v and w != low_v
                         and A[high_v][w] and A[w][low_v]]
            if candidates:
                w = random.choice(candidates)
                A[high_v][w] = 0
                A[w][high_v] = 1
                A[w][low_v] = 0
                A[low_v][w] = 1
    return A

n = 9
print(f"Searching for H=3357 regular tournament on n={n}...")
t0 = time.time()

found_A = None
for trial in range(100000):
    A = random_regular_tournament(n)
    scores = [sum(A[i]) for i in range(n)]
    if max(scores) - min(scores) > 0:
        continue
    H = H_tournament(A, n)
    if H == 3357:
        found_A = [row[:] for row in A]
        print(f"Found H=3357 at trial {trial} ({time.time()-t0:.1f}s)")
        break

if found_A is None:
    print("ERROR: Could not find H=3357")
    sys.exit(1)

# Print adjacency matrix
print("\nAdjacency matrix:")
for row in found_A:
    print("  " + " ".join(map(str, row)))

# Eigenvalues
S = np.array([[found_A[i][j] - found_A[j][i] for j in range(n)] for i in range(n)], dtype=float)
evals = np.linalg.eigvals(S)
pos_imag = sorted([e.imag for e in evals if e.imag > 0.01])
gap = max(pos_imag) - min(pos_imag)
print(f"\nEigenvalues: [{', '.join(f'{e:.4f}' for e in pos_imag)}]")
print(f"Spectral gap: {gap:.4f}")

# Count 3-cycles
c3 = 0
for i in range(n):
    for j in range(i+1, n):
        for k in range(j+1, n):
            if (found_A[i][j] and found_A[j][k] and found_A[k][i]) or \
               (found_A[i][k] and found_A[k][j] and found_A[j][i]):
                c3 += 1
print(f"c3 = {c3}")

# Compute Betti (max_dim=7 is too slow, use 5)
print(f"\nComputing Betti numbers (max_dim=5)...")
t1 = time.time()
beta = path_betti_numbers(found_A, n, max_dim=5)
beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(6)]
print(f"β = {beta_list} ({time.time()-t1:.1f}s)")

# Check vertex-deletion topology
print(f"\nVertex-deletion topology:")
from collections import Counter
del_bettis = Counter()
for v in range(n):
    An = [[found_A[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
    Hn = H_tournament(An, n-1)
    bv = path_betti_numbers(An, n-1, max_dim=5)
    bv_list = tuple(int(bv[k]) if k < len(bv) else 0 for k in range(6))
    del_bettis[bv_list] += 1
    print(f"  v={v}: H(T-v)={Hn}, β={list(bv_list)}")

print(f"\nDeletion summary:")
for bv, cnt in del_bettis.most_common():
    print(f"  β={list(bv)}: {cnt} deletions")

print(f"\nTotal time: {time.time()-t0:.1f}s")
