#!/usr/bin/env python3
"""
beta3_maximizer_n7.py - Does β₃ > 0 correlate with high H at n=7?

At n=6: H-maximizers (H=45) split between β₁=1 (14) and β₃=1 (21).
Question: At n=7, do β₃>0 tournaments have higher H than average?

Also: check the β₁/β₃ split among high-H tournaments at n=7.

Author: opus-2026-03-07-S46e
"""
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import path_betti_numbers
from itertools import permutations
import random

random.seed(42)

def ham_count(A, n):
    count = 0
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[i+1]] for i in range(n-1)):
            count += 1
    return count

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

# === n=7: β vs H correlation ===
print("=" * 60)
print("n=7: β₃ vs H (sampling 300 tournaments)")
print("=" * 60)

n = 7
results = {'b0': [], 'b1': [], 'b3': []}

for trial in range(300):
    A = random_tournament(n)
    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0
    H = ham_count(A, n)

    if b3 > 0:
        results['b3'].append(H)
    elif b1 > 0:
        results['b1'].append(H)
    else:
        results['b0'].append(H)

    if (trial + 1) % 50 == 0:
        print(f"  ... {trial+1}/300 done")

print(f"\nResults:")
for key in ['b0', 'b1', 'b3']:
    Hs = results[key]
    if Hs:
        mean_H = sum(Hs) / len(Hs)
        print(f"  {key}: n={len(Hs)}, mean H={mean_H:.1f}, range=[{min(Hs)}, {max(Hs)}]")
    else:
        print(f"  {key}: n=0")

# Expected H for random tournament at n=7: n!/2^(n-1) = 5040/64 = 78.75
print(f"\nExpected H (random): {5040/64:.2f}")
print(f"Max possible H at n=7: {max(max(results[k]) for k in results if results[k])}")

# === Targeted search for high-H at n=7 ===
print("\n" + "=" * 60)
print("n=7: HIGH-H TOURNAMENTS β classification")
print("=" * 60)

max_H_found = 0
high_H_data = []

for trial in range(2000):
    A = random_tournament(n)
    H = ham_count(A, n)

    if H > max_H_found:
        max_H_found = H

    if H >= 100:  # Focus on high-H
        beta = path_betti_numbers(A, n)
        b1 = int(beta[1]) if len(beta) > 1 else 0
        b3 = int(beta[3]) if len(beta) > 3 else 0
        high_H_data.append((H, b1, b3))

print(f"Max H found: {max_H_found}")
print(f"Tournaments with H >= 100: {len(high_H_data)}")

if high_H_data:
    from collections import Counter
    by_H = {}
    for H, b1, b3 in high_H_data:
        if H not in by_H:
            by_H[H] = Counter()
        label = f"β₁={b1},β₃={b3}"
        by_H[H][label] += 1

    for H in sorted(by_H.keys(), reverse=True):
        print(f"  H={H}: {dict(by_H[H])}")

print("\n" + "=" * 60)
print("INTERPRETATION")
print("=" * 60)
print("""
At n=6: β₃>0 ↔ H-maximizer was a strong signal (60% of H=45 had β₃=1).
At n=7: Does this pattern persist? If β₃>0 correlates with high H,
this suggests path homology detects "cycle-optimal" tournament structure.

Connection to OCF: H = I(Ω(T), 2) = 1 + 2α₁ + 4α₂ + ...
High H means many vertex-disjoint odd cycle collections.
β₃ > 0 means nontrivial 3-dimensional directed topology.
Both may be measuring "richness" of odd-cycle interactions.
""")
