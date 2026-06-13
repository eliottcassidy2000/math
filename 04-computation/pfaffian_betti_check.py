#!/usr/bin/env python3
"""
pfaffian_betti_check.py - Verify: β₃>0 iff |Pf(S)| = 7 at n=6?

DISCOVERY from spectral_betti_connection.py: ALL 6 sampled β₃>0 at n=6
had |det(S)| = 49 = Pf(S)² = 7².

Is this universal? If so, it's a SPECTRAL CRITERION for path homology!

Note: for tournament on n=6, the skew-adjacency S has det(S) = Pf(S)².
Possible |Pf| values for 6-tournament: each S[i][j] = ±1 (off-diag),
so Pf is sum of 15 terms each ±1, giving Pf ∈ {-15,-13,...,13,15}.
Actually Pf(S) = Σ_{perfect matchings} sgn(M)·∏S[i_k][j_k].
For n=6: 15 perfect matchings of K_6.

Author: opus-2026-03-07-S46e
"""
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import path_betti_numbers
import numpy as np
from collections import Counter
import random

random.seed(123)  # Different seed for independence

def skew_adj(A, n):
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            S[i][j] = A[i][j] - A[j][i]
    return S

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

# === EXHAUSTIVE n=6 ===
print("=" * 60)
print("n=6: EXHAUSTIVE Pfaffian vs β₃")
print("=" * 60)

n = 6
m = n*(n-1)//2  # 15
total = 1 << m  # 32768

pf_by_beta = {'b0': Counter(), 'b1': Counter(), 'b3': Counter()}
count_by_beta = Counter()

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0

    S = skew_adj(A, n)
    det_val = round(np.linalg.det(S))
    pf_sq = abs(det_val)

    if b3 > 0:
        key = 'b3'
    elif b1 > 0:
        key = 'b1'
    else:
        key = 'b0'

    pf_by_beta[key][pf_sq] += 1
    count_by_beta[key] += 1

    if (bits + 1) % 10000 == 0:
        print(f"  ... {bits+1}/{total}")

print(f"\nTotal: {sum(count_by_beta.values())}")
for key in ['b0', 'b1', 'b3']:
    print(f"\n{key}: {count_by_beta[key]} tournaments")
    print(f"  |Pf(S)|² distribution: {dict(sorted(pf_by_beta[key].items()))}")

# Check: is β₃>0 iff |Pf|² = 49?
if 'b3' in pf_by_beta and pf_by_beta['b3']:
    b3_pf_values = set(pf_by_beta['b3'].keys())
    print(f"\nβ₃>0 Pfaffian² values: {b3_pf_values}")
    if b3_pf_values == {49}:
        print("CONFIRMED: β₃>0 iff |Pf(S)|² = 49 (|Pf| = 7)")

        # But does |Pf|² = 49 imply β₃>0?
        pf49_b0 = pf_by_beta['b0'].get(49, 0)
        pf49_b1 = pf_by_beta['b1'].get(49, 0)
        pf49_b3 = pf_by_beta['b3'].get(49, 0)
        print(f"\n|Pf|² = 49 tournaments: β₃={pf49_b3}, β₁={pf49_b1}, β=0={pf49_b0}")
        if pf49_b0 == 0 and pf49_b1 == 0:
            print("→ β₃>0 IFF |Pf| = 7 at n=6!")
        else:
            print("→ |Pf| = 7 is necessary but NOT sufficient for β₃>0")

print("\n" + "=" * 60)
print("PFAFFIAN VALUE ANALYSIS")
print("=" * 60)

# What are all possible Pf values for 6-tournaments?
all_pf_sq = set()
for key in pf_by_beta:
    all_pf_sq |= set(pf_by_beta[key].keys())
print(f"All |Pf|² values: {sorted(all_pf_sq)}")
print(f"All |Pf| values: {sorted(set(int(round(np.sqrt(v))) for v in all_pf_sq))}")
