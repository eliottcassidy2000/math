#!/usr/bin/env python3
"""
maximizer_betti_n8.py — Check Betti numbers for H-maximizers at n=8

From memory: max H at n=8 is 661, score (3,3,3,3,4,4,4,4)
Check: do maximizers have beta_5 > 0 (= beta_{n-3})?

Strategy: targeted search for H>=600 tournaments, compute Betti.
n=8 has 2^28 = 268M tournaments, so we sample.

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
    if len(pos) < 2:
        return 0.0
    return pos[-1] - pos[0]

# ===== Search for H-maximizers at n=8 =====
print("=" * 70)
print("n=8 H-MAXIMIZER SEARCH + BETTI NUMBERS")
print("=" * 70)

n = 8
high_H = []
t0 = time.time()

# H_tournament at n=8 needs 2^8 * 8 = 2048 entries per DP, manageable
# But total search space is huge, so we sample
print("Sampling random n=8 tournaments...")

for trial in range(20000):
    A = random_tournament(n)
    H = H_tournament(A, n)

    if H >= 500:  # Near-max
        sc = score_seq(A, n)
        S = skew_adj(A, n)
        gap = spectral_gap(S)
        try:
            beta = path_betti_numbers(A, n, max_dim=5)
            beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(6)]
        except:
            beta_list = [1, 0, 0, 0, 0, 0]

        high_H.append({
            'H': H, 'beta': beta_list, 'score': sc,
            'gap': gap
        })

    if (trial + 1) % 5000 == 0:
        print(f"  {trial+1}/20000 ({time.time()-t0:.1f}s), found {len(high_H)} with H>=500")

elapsed = time.time() - t0
print(f"\nDone in {elapsed:.1f}s")
print(f"Found {len(high_H)} with H >= 500")

# Results by H
print("\nH >= 500 tournaments:")
H_groups = {}
for d in high_H:
    H = d['H']
    if H not in H_groups:
        H_groups[H] = []
    H_groups[H].append(d)

for H_val in sorted(H_groups.keys(), reverse=True)[:15]:
    entries = H_groups[H_val]
    betti_counts = Counter(tuple(d['beta']) for d in entries)
    scores = Counter(d['score'] for d in entries)
    gaps = [d['gap'] for d in entries]
    print(f"\n  H={H_val} ({len(entries)} found):")
    for b, cnt in betti_counts.most_common(5):
        print(f"    beta={list(b)}: {cnt}")
    for s, cnt in scores.most_common(3):
        print(f"    score={s}: {cnt}")
    print(f"    gap: mean={np.mean(gaps):.4f}, min={min(gaps):.4f}, max={max(gaps):.4f}")

# Also try to construct near-regular tournaments more likely to be maximizers
print("\n" + "=" * 70)
print("TARGETED: Near-regular n=8 tournaments")
print("=" * 70)

# Generate tournaments with score close to (3,3,3,3,4,4,4,4) or (3,4,3,4,3,4,3,4)
# Use rejection sampling
near_reg = []
random.seed(1234)

for trial in range(50000):
    A = random_tournament(n)
    sc = score_seq(A, n)
    # Check if near-regular (low variance)
    scores = [sum(A[i]) for i in range(n)]
    var = np.var(scores)
    if var <= 0.5:  # Very regular
        H = H_tournament(A, n)
        try:
            beta = path_betti_numbers(A, n, max_dim=5)
            beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(6)]
        except:
            beta_list = [1, 0, 0, 0, 0, 0]

        S = skew_adj(A, n)
        gap = spectral_gap(S)
        near_reg.append({
            'H': H, 'beta': beta_list, 'score': sc,
            'gap': gap
        })

    if (trial + 1) % 10000 == 0:
        print(f"  {trial+1}/50000, found {len(near_reg)} near-regular")

print(f"\nFound {len(near_reg)} near-regular tournaments")

if near_reg:
    for d in sorted(near_reg, key=lambda x: -x['H'])[:20]:
        print(f"  H={d['H']}, beta={d['beta']}, score={d['score']}, gap={d['gap']:.4f}")

    # Best H found
    max_found = max(d['H'] for d in near_reg)
    max_tours = [d for d in near_reg if d['H'] == max_found]
    print(f"\nBest H found: {max_found} ({len(max_tours)} instances)")
    for d in max_tours[:5]:
        print(f"  beta={d['beta']}, gap={d['gap']:.4f}")
