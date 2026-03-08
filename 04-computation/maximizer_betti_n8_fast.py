#!/usr/bin/env python3
"""
maximizer_betti_n8_fast.py — Fast n=8 maximizer Betti check

Strategy: first find high-H tournaments WITHOUT computing Betti,
then compute Betti only for the very best ones.

From memory: max H at n=8 is 661, score (3,3,3,3,4,4,4,4)

Author: kind-pasteur-2026-03-08-S40
"""
import sys, os, time, random
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

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

def score_seq(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)]))

def skew_adj(A, n):
    S = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            S[i][j] = A[i][j] - A[j][i]
    return S

def spectral_gap(S):
    eigs = np.linalg.eigvals(S)
    pos = sorted([abs(e.imag) for e in eigs if abs(e.imag) > 0.01])
    if len(pos) < 2: return 0.0
    return pos[-1] - pos[0]

# Phase 1: Find high-H tournaments quickly (no Betti)
print("=" * 70)
print("PHASE 1: Find high-H n=8 tournaments (no Betti)")
print("=" * 70)

n = 8
best_tours = []
t0 = time.time()

for trial in range(100000):
    A = random_tournament(n)
    H = H_tournament(A, n)
    if H >= 550:
        sc = score_seq(A, n)
        best_tours.append({'H': H, 'A': [row[:] for row in A], 'score': sc})

    if (trial + 1) % 20000 == 0:
        print(f"  {trial+1}/100000 ({time.time()-t0:.1f}s), found {len(best_tours)} with H>=550")

print(f"\nDone in {time.time()-t0:.1f}s")
print(f"Found {len(best_tours)} with H >= 550")

# Show H distribution
H_dist = Counter(d['H'] for d in best_tours)
for H_val in sorted(H_dist.keys(), reverse=True)[:10]:
    print(f"  H={H_val}: {H_dist[H_val]}")

# Phase 2: Compute Betti for top tournaments only
print("\n" + "=" * 70)
print("PHASE 2: Betti numbers for top tournaments")
print("=" * 70)

# Sort by H descending, take top 30
best_tours.sort(key=lambda d: -d['H'])
top = best_tours[:30]

for i, d in enumerate(top):
    A = d['A']
    H = d['H']
    S = skew_adj(A, n)
    gap = spectral_gap(S)

    t1 = time.time()
    try:
        beta = path_betti_numbers(A, n, max_dim=5)
        beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(6)]
    except Exception as e:
        beta_list = f"FAILED: {e}"

    dt = time.time() - t1
    print(f"  [{i+1}] H={H}, score={d['score']}, gap={gap:.4f}, beta={beta_list} ({dt:.1f}s)")

# Phase 3: Also check a few near-max for comparison
print("\n--- Near-max comparison ---")
near_max = [d for d in best_tours if d['H'] >= 600 and d['H'] < max(dd['H'] for dd in best_tours)]
for d in near_max[:5]:
    A = d['A']
    H = d['H']
    S = skew_adj(A, n)
    gap = spectral_gap(S)
    try:
        beta = path_betti_numbers(A, n, max_dim=5)
        beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(6)]
    except Exception as e:
        beta_list = f"FAILED: {e}"
    print(f"  H={H}, score={d['score']}, gap={gap:.4f}, beta={beta_list}")
