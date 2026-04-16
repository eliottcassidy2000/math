#!/usr/bin/env python3
"""
n17_scan.py -- Find the H-maximizing circulant on Z_17 and compute alpha decomposition.

17 is prime but 17 ≡ 1 mod 4, so no Paley tournament.
We need |S| = 8 (pick one from each pair {k, 17-k} for k=1..8).

Strategy:
  1. Quick H scan over SC circulants (S closed under complement mod 17)
  2. Then test top candidates more carefully
  3. Full alpha decomposition on the winner
"""

from itertools import combinations
from math import factorial, gcd
import time

def circulant(n, S):
    adj = [0]*n
    for v in range(n):
        for k in S:
            adj[v] |= 1 << ((v+k) % n)
    return adj

def H_dp(adj):
    n = len(adj); N = 1 << n
    dp = [[0]*n for _ in range(N)]
    for v in range(n): dp[1 << v][v] = 1
    full = N-1
    for mask in range(1, N):
        row = dp[mask]
        for v in range(n):
            val = row[v]
            if not val: continue
            outs = adj[v] & ~mask & full
            while outs:
                ub = outs & -outs; u = ub.bit_length()-1
                dp[mask|ub][u] += val; outs ^= ub
    return sum(dp[full])

# ── Quick scan: SC circulants on Z_17 ─────────────────────────────────────────
# SC means S = {17-k : k in S}, i.e. pick full pairs {k, 17-k}.
# There are 8 complementary pairs; we pick 4 of them.
# C(8,4) = 70 circulants to try.

n = 17
pairs = [(k, n-k) for k in range(1, n//2+1)]  # [(1,16),(2,15),...,(8,9)]

print(f"Scanning SC circulants on Z_{n}  [C(8,4)={len(list(combinations(range(8),4)))} total]")
print(f"{'S':>50}  H")

best_H, best_S = 0, None
results = []

for chosen_pairs in combinations(range(8), 4):
    S = []
    for i in chosen_pairs:
        S.extend(pairs[i])
    S.sort()
    adj = circulant(n, S)
    H = H_dp(adj)
    results.append((H, S))
    if H > best_H:
        best_H, best_S = H, S

results.sort(reverse=True)
for H, S in results[:10]:
    marker = " <-- MAX" if H == best_H else ""
    print(f"  S={S!s:>45}  {H:>12,}{marker}")

print(f"\nBest SC circulant: S={best_S}, H={best_H:,}")

# ── Also try a few non-SC candidates ─────────────────────────────────────────
# n=9 best was {1,2,3,5} (not SC). Check some greedy picks.
print(f"\nSpot-checking non-SC candidates:")

extras = [
    [1,2,3,4,5,6,7,8],   # first half
    [1,3,5,7,9,11,13,15], # all-odd
    [2,4,6,8,10,12,14,16], # all-even
    [1,2,4,8,9,13,15,16], # QR-like pattern
]

for S in extras:
    S = sorted(k for k in S if 0 < k < n)
    if len(S) != 8:
        continue
    adj = circulant(n, S)
    H = H_dp(adj)
    marker = " <-- NEW MAX" if H > best_H else ""
    print(f"  S={S!s:<45}  {H:>12,}{marker}")
    if H > best_H:
        best_H, best_S = H, S

print(f"\nOverall best: S={best_S}, H={best_H:,}")
