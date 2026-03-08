#!/usr/bin/env python3
"""
beta3_maximizer_check.py - Is β₃ > 0 related to H-maximizers?

DISCOVERY at n=6: ALL 3 sampled β₃>0 tournaments had H=45 (the maximum!).
They all had t₃=8, t₅=12, α₂=1.

HYPOTHESIS: At n=6, β₃ > 0 iff T is an H-maximizer.
If true, this is a deep connection between path homology and H-maximization.

Author: opus-2026-03-07-S46e
"""
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import path_betti_numbers
from itertools import combinations, permutations
from collections import Counter
from math import factorial

def ham_count(A, n):
    count = 0
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[i+1]] for i in range(n-1)):
            count += 1
    return count

def count_t3(A, n):
    return sum(1 for i,j,k in combinations(range(n), 3)
               if (A[i][j] and A[j][k] and A[k][i]) or
                  (A[i][k] and A[k][j] and A[j][i]))

# Exhaustive n=5 check: does β₁ = 1 correlate with H-max?
print("=" * 60)
print("n=5: β₁ vs H")
print("=" * 60)

n = 5
m = n*(n-1)//2
b1_h = {}
for bits in range(1 << m):
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
    H = ham_count(A, n)
    if b1 not in b1_h:
        b1_h[b1] = []
    b1_h[b1].append(H)

for b1_val in sorted(b1_h.keys()):
    Hs = b1_h[b1_val]
    H_dist = Counter(Hs)
    print(f"  β₁={b1_val}: {len(Hs)} tournaments, H range [{min(Hs)}, {max(Hs)}]")
    print(f"    H distribution: {dict(sorted(H_dist.items()))}")

max_H = max(max(Hs) for Hs in b1_h.values())
max_H_count = sum(1 for H_list in b1_h.values() for H in H_list if H == max_H)
max_H_b1 = {b1: sum(1 for H in Hs if H == max_H) for b1, Hs in b1_h.items()}
print(f"\nH-maximizers (H={max_H}): {max_H_count} total")
for b1_val, count in max_H_b1.items():
    if count > 0:
        print(f"  β₁={b1_val}: {count}")

# Now check n=6: is β₃ = 1 concentrated on H-maximizers?
print("\n" + "=" * 60)
print("n=6: β₃ vs H (sampling H=43-45 range)")
print("=" * 60)

# At n=6, max H = 45. Let me check β for high-H tournaments.
import random
random.seed(42)

n = 6
high_H_betas = {}
samples_by_H = Counter()
total_checked = 0

# Targeted search: generate many and check β for high-H ones
for trial in range(2000):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    H = ham_count(A, n)
    total_checked += 1

    if H >= 39:  # Check β only for high-H tournaments
        beta = path_betti_numbers(A, n)
        b3 = int(beta[3]) if len(beta) > 3 else 0
        b1 = int(beta[1]) if len(beta) > 1 else 0
        key = (H, b1, b3)
        high_H_betas[key] = high_H_betas.get(key, 0) + 1
        samples_by_H[H] += 1

print(f"Checked {total_checked} tournaments, {sum(samples_by_H.values())} with H >= 39")
print("\nH vs (β₁, β₃):")
for key in sorted(high_H_betas.keys(), reverse=True):
    H, b1, b3 = key
    count = high_H_betas[key]
    print(f"  H={H}, β₁={b1}, β₃={b3}: {count}")

# Direct check: build the C6 regular tournament and check
print("\n" + "=" * 60)
print("n=6: SPECIFIC HIGH-H TOURNAMENTS")
print("=" * 60)

# Known H=45 tournaments at n=6 are SC with score (2,2,3,3,3,2)?
# Actually score (2,2,2,3,3,3) for even n=6.
# Let me try some specific constructions.

# Circulant tournament C6^{1,2,3}:
A_circ = [[0]*6 for _ in range(6)]
for i in range(6):
    for d in [1, 2, 3]:
        A_circ[i][(i+d)%6] = 1

# Wait, this is a tournament? For n=6, each vertex sends to 3 others (outdeg=3=(n-1)/2).
# But wait, for i→(i+1)%6 and also i→(i+4)%6? No:
# Connection set {1,2,3}: i→(i+1), i→(i+2), i→(i+3).
# Then (i+4)→i is the complement. But is (i+3)→i or i→(i+3)?
# i→(i+3)%6 if 3 ∈ S. Yes. But also (i+3)→i iff (i+3)+d ≡ i mod 6 for some d ∈ {1,2,3}.
# d = 6-3 = 3. So (i+3)→(i+6)%6 = i. Both i→(i+3) AND (i+3)→i?
# NO: in a tournament, this can't happen. The circulant S={1,2,3} at n=6 is NOT a tournament
# because 3 and 6-3=3 are the same mod 6. For d=3: i→(i+3) and (i+3)→i would both be in S.
# This is a bidirected edge, not a tournament!

# For n=6, valid circulant tournament connection sets S must have |S| = 3 = (n-1)/2,
# and {d, n-d} must have exactly one element in S for each d.
# d=1: choose 1 or 5
# d=2: choose 2 or 4
# d=3: choose 3... but n-d = 3 too! So 3 is in S iff 3 is in S — no constraint.
# Actually for d=3 = n/2 at even n=6: the arc i→(i+3) and (i+3)→i are different arcs.
# To be a tournament, exactly one must be present. So 3 ∈ S or 3 ∉ S (since n-3=3=d).
# If 3 ∈ S: i→(i+3) for all i. Then also (i+3)→((i+3)+3)%6 = (i+6)%6 = i. Wait no:
# (i+3)→(i+3+d) for d∈S. If 3∈S then (i+3)→(i+6)%6 = i. So (i+3)→i AND i→(i+3).
# This is bidirectional! So 3 CANNOT be in S at n=6.

# Valid S at n=6: choose from {1 or 5}, {2 or 4}, and 3 not included. |S|=2, not 3!
# But we need |S| = (n-1)/2 = 5/2 which is not integer. n=6 is even, so
# circulant tournaments at even n need connection set of size (n-1)/2 which is non-integer.
# So there are NO regular circulant tournaments at even n.

# Let me just build a self-complementary tournament with H=45 at n=6.
# From the OEIS: a(6) = 45. SC tournaments achieve this.
# The Paley construction doesn't apply at n=6.
# Let me generate and find one.

max_found_H = 0
max_A = None
for trial in range(5000):
    A = [[0]*6 for _ in range(6)]
    for i in range(6):
        for j in range(i+1, 6):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    H = ham_count(A, 6)
    if H > max_found_H:
        max_found_H = H
        max_A = [row[:] for row in A]

print(f"Best H found: {max_found_H}")
if max_A:
    beta = path_betti_numbers(max_A, 6)
    t3 = count_t3(max_A, 6)
    print(f"  β = {[int(b) for b in beta]}")
    print(f"  t₃ = {t3}")
    print(f"  Adjacency:")
    for i in range(6):
        out = [j for j in range(6) if max_A[i][j]]
        print(f"    {i} → {out}")
