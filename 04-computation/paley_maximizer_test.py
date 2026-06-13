#!/usr/bin/env python3
"""
Test the Paley maximizer conjecture at n=7 (P(7) should give max H=189=a(7))
and investigate the gap at n=8.

Also: check whether a(8)=661 is achieved by ANY tournament (not just SC).
We can't enumerate all 2^28 = 268M tournaments at n=8, but we can check:
1. All SC tournaments with different anti-automorphisms
2. Random regular tournaments 
3. Whether 661 = a(8) from OEIS

Instance: opus-2026-03-05-S9
"""
from itertools import combinations, permutations
from collections import Counter
import random

def count_ham_dp(T, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1][v] for v in range(n))

# PART 1: Verify a(7)=189 and P(7) achieves it
print("="*60)
print("PART 1: n=7 — exhaustive check")
print("="*60)

# P(7)
P7 = [[0]*7 for _ in range(7)]
qr7 = {1, 2, 4}
for i in range(7):
    for j in range(7):
        if i != j and (j - i) % 7 in qr7:
            P7[i][j] = 1

h_p7 = count_ham_dp(P7, 7)
print(f"H(P(7)) = {h_p7}")

# Enumerate ALL tournaments on 7 vertices (2^21 = 2M)
n7 = 7
num_arcs_7 = n7*(n7-1)//2  # 21
max_h7 = 0
max_count = 0
paley_is_max = True

arc_pairs = [(i,j) for i in range(n7) for j in range(i+1, n7)]

print(f"Enumerating all {2**num_arcs_7} tournaments on 7 vertices...")
for bits in range(2**num_arcs_7):
    T = [[0]*n7 for _ in range(n7)]
    for k, (i, j) in enumerate(arc_pairs):
        if (bits >> k) & 1:
            T[i][j] = 1
        else:
            T[j][i] = 1
    h = count_ham_dp(T, n7)
    if h > max_h7:
        max_h7 = h
        max_count = 1
    elif h == max_h7:
        max_count += 1

print(f"Max H over ALL n=7 tournaments: {max_h7}")
print(f"Number achieving max: {max_count}")
print(f"P(7) achieves max: {h_p7 == max_h7}")
print(f"OEIS a(7) = 189, matches: {max_h7 == 189}")

# PART 2: n=8 random sampling to see if anything beats 661
print(f"\n{'='*60}")
print("PART 2: n=8 — random sampling")
print("="*60)

n8 = 8
best_random = 0
random.seed(42)
num_samples = 100000

for _ in range(num_samples):
    T = [[0]*n8 for _ in range(n8)]
    for i in range(n8):
        for j in range(i+1, n8):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    h = count_ham_dp(T, n8)
    if h > best_random:
        best_random = h

print(f"Best H among {num_samples} random n=8 tournaments: {best_random}")
print(f"OEIS a(8) = 661")

# PART 3: Check the Paley conjecture values
print(f"\n{'='*60}")
print("PART 3: Paley conjecture check")
print("="*60)

# Build P(3)
P3 = [[0,1,0],[0,0,1],[1,0,0]]
h_p3 = count_ham_dp(P3, 3)
print(f"H(P(3)) = {h_p3}, a(3) = 3, match: {h_p3 == 3}")

# P(7) already computed
print(f"H(P(7)) = {h_p7}, a(7) = 189, match: {h_p7 == 189}")

# Build P(11)
P11 = [[0]*11 for _ in range(11)]
qr11 = {1, 3, 4, 5, 9}  # QRs mod 11
for i in range(11):
    for j in range(11):
        if i != j and (j - i) % 11 in qr11:
            P11[i][j] = 1
h_p11 = count_ham_dp(P11, 11)
print(f"H(P(11)) = {h_p11}, a(11) = 95095, match: {h_p11 == 95095}")

# P(19) — might be slow but let's try
P19 = [[0]*19 for _ in range(19)]
qr19 = {1, 4, 5, 6, 7, 9, 11, 16, 17}  # QRs mod 19
for i in range(19):
    for j in range(19):
        if i != j and (j - i) % 19 in qr19:
            P19[i][j] = 1

# P(19) with n=19 DP is O(2^19 * 19) = ~10M, should be fine
print("Computing H(P(19))...")
h_p19 = count_ham_dp(P19, 19)
print(f"H(P(19)) = {h_p19}")
