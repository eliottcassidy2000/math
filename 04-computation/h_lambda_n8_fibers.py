#!/usr/bin/env python3
"""
h_lambda_n8_fibers.py — opus-2026-03-13-S71c

At n=7: only 31/46010 lambda fibers show H ambiguity (ΔH=2 always).
QUESTION: What happens at n=8?
- More ambiguous fibers? (c7 has more variation at n=8: |dc7|≤3)
- Larger H gaps?
- α₂ still lambda-determined?
"""

import sys, time
import numpy as np
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

lam_to_H = defaultdict(set)
lam_to_c3c5 = defaultdict(set)

t0 = time.time()
SAMPLES = 20000
for trial in range(SAMPLES):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    lk = lambda_key(A, n)
    H = count_ham_paths(A, n)
    c3 = int(np.trace(A @ A @ A)) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
    lam_to_H[lk].add(H)
    lam_to_c3c5[lk].add((c3, c5))

    if trial % 5000 == 0 and trial > 0:
        dt = time.time() - t0
        ambig = sum(1 for v in lam_to_H.values() if len(v) > 1)
        print(f"  trial {trial}: {dt:.1f}s, {len(lam_to_H)} groups, {ambig} ambiguous")

dt = time.time() - t0
print(f"\nn=8: {len(lam_to_H)} lambda groups from {SAMPLES} samples, {dt:.1f}s")

ambig_H = 0
gap_dist = defaultdict(int)
for lk, H_vals in lam_to_H.items():
    if len(H_vals) > 1:
        ambig_H += 1
        gap = max(H_vals) - min(H_vals)
        gap_dist[gap] += 1

print(f"H ambiguous: {ambig_H}/{len(lam_to_H)} ({100*ambig_H/len(lam_to_H):.3f}%)")
print(f"Gap distribution: {dict(sorted(gap_dist.items()))}")

# (c3,c5) still lambda-determined?
ambig_c3c5 = sum(1 for v in lam_to_c3c5.values() if len(v) > 1)
print(f"(c3,c5) ambiguous: {ambig_c3c5}/{len(lam_to_c3c5)}")

if ambig_H > 0:
    print("\nAmbiguous fiber details:")
    count = 0
    for lk, H_vals in sorted(lam_to_H.items()):
        if len(H_vals) > 1:
            c3c5 = lam_to_c3c5.get(lk, set())
            print(f"  H ∈ {sorted(H_vals)}, (c3,c5) = {sorted(c3c5)}, ΔH = {max(H_vals)-min(H_vals)}")
            count += 1
            if count >= 20: break

print("\nDone.")
