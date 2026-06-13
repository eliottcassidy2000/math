#!/usr/bin/env python3
"""
h_lambda_n8_vitali.py — opus-2026-03-13-S71c

At n=8: Check H variation within ACTUAL lambda fibers.
Previous (1,1,2,2) reversals don't preserve lambda at n=8!
Need to verify lambda is actually preserved before checking H.
"""

import sys, time
import numpy as np
from itertools import combinations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def lambda_key(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
    return tuple(L[i][j] for i in range(n) for j in range(i+1, n))

n = 8
tb = n*(n-1)//2
np.random.seed(42)

print(f"n={n}: Finding lambda-preserving reversals and checking H variation...")

lambda_preserving = 0
h_changed = 0
dH_dist = defaultdict(int)

t0 = time.time()
for trial in range(1000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    lk_orig = lambda_key(A, n)
    H_orig = count_ham_paths(A, n)

    for combo in combinations(range(n), 4):
        verts = list(combo)
        sub_scores = [sum(A[v][w] for w in verts if w != v) for v in verts]
        if sorted(sub_scores) != [1, 1, 2, 2]:
            continue

        # Reverse edges within this 4-set
        A_new = A.copy()
        for i_idx in range(4):
            for j_idx in range(i_idx+1, 4):
                u, v = verts[i_idx], verts[j_idx]
                A_new[u][v], A_new[v][u] = A_new[v][u], A_new[u][v]

        lk_new = lambda_key(A_new, n)
        if lk_new == lk_orig:
            lambda_preserving += 1
            H_new = count_ham_paths(A_new, n)
            dH = H_new - H_orig
            if dH != 0:
                h_changed += 1
                dH_dist[abs(dH)] += 1
                if h_changed <= 10:
                    print(f"  trial {trial}, {combo}: H={H_orig}→{H_new}, ΔH={dH}")

    if trial % 250 == 0:
        dt = time.time() - t0
        print(f"  ... trial {trial}: {dt:.1f}s, lam-preserving: {lambda_preserving}, dH≠0: {h_changed}")

dt = time.time() - t0
print(f"\nTotal: lambda-preserving (1,1,2,2): {lambda_preserving}")
print(f"  Of which H changed: {h_changed}")
print(f"  |ΔH| distribution: {dict(sorted(dH_dist.items()))}")

print("\nDone.")
