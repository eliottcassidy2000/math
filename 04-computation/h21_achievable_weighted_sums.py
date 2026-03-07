"""
h21_achievable_weighted_sums.py — Which values of alpha_1 + 2*alpha_2 are achievable?

H = 1 + 2*(alpha_1 + 2*alpha_2 + 4*alpha_3 + ...)
For H=7: need alpha_1+2*alpha_2+4*alpha_3+... = 3
For H=21: need alpha_1+2*alpha_2+4*alpha_3+... = 10

Key question: which values of this "weighted sum" w = (H-1)/2 are achievable?
If 3 and 10 are NEVER achievable, what's the pattern?

At each n, compute the set of achievable w values.

Author: opus-2026-03-07-S43
"""
import random
from itertools import combinations
from math import comb

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def hamiltonian_paths_held_karp(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

# Exhaustive for n <= 7
print("=== Achievable w = (H-1)/2 values ===")

for n in range(3, 8):
    w_set = set()
    h_set = set()
    nbits = n*(n-1)//2
    total = 1 << nbits

    for T in range(total):
        # Build adjacency matrix from bits
        A = [[0]*n for _ in range(n)]
        pos = 0
        for i in range(n):
            for j in range(i+1, n):
                if (T >> pos) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                pos += 1

        H = hamiltonian_paths_held_karp(A, n)
        h_set.add(H)
        w = (H - 1) // 2
        w_set.add(w)

    max_w = max(w_set)
    missing = [w for w in range(max_w + 1) if w not in w_set]

    print(f"\nn={n}: max_H={(max(h_set))}, max_w={max_w}")
    print(f"  Achievable H: {sorted(h_set)}")
    print(f"  Achievable w: {sorted(w_set)}")
    print(f"  Missing w in [0..{max_w}]: {missing}")

# n=8 sampling
print("\n\n=== n=8: sampling for achievable w values ===")
n = 8
w_set = set()
h_set = set()
trials = 2000000
for _ in range(trials):
    A = random_tournament(n)
    H = hamiltonian_paths_held_karp(A, n)
    h_set.add(H)
    w_set.add((H - 1) // 2)

max_w = max(w_set)
missing_small = [w for w in range(min(100, max_w + 1)) if w not in w_set]
print(f"n=8: max_H={max(h_set)}, max_w={max_w}")
print(f"  Missing w in [0..99]: {missing_small}")
print(f"  Number of achievable w in [0..99]: {sum(1 for w in range(100) if w in w_set)}")

# Check: is w=3 (H=7) missing at all n?
# Check: is w=10 (H=21) missing at all n?
print(f"\n  w=3 (H=7) achieved? {'YES' if 3 in w_set else 'NO'}")
print(f"  w=10 (H=21) achieved? {'YES' if 10 in w_set else 'NO'}")

# Also check w=31 (H=63)
print(f"  w=31 (H=63) achieved? {'YES' if 31 in w_set else 'NO'}")
