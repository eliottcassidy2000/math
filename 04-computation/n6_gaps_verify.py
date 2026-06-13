#!/usr/bin/env python3
"""
Exhaustive n=6 H-value verification.
opus-2026-03-14-S84

n=6 has 2^15 = 32768 tournaments. Enumerate ALL to find exact H distribution.
Previous sampling found gaps at {7, 21, 35, 39}. Let's verify exhaustively.
"""

from itertools import permutations
from collections import Counter
import sys

n = 6
m = n * (n - 1) // 2  # 15
N = 1 << m  # 32768

print(f"Exhaustive n={n} H computation: {N} tournaments, {m} arcs")
print(f"This may take a few minutes...")

arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
all_perms = list(permutations(range(n)))

H_dist = Counter()

for bits in range(N):
    if bits % 5000 == 0:
        print(f"  Progress: {bits}/{N} ({100*bits/N:.1f}%)", file=sys.stderr)

    # Build adjacency
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    # Count Hamiltonian paths
    H = 0
    for p in all_perms:
        valid = True
        for i in range(n-1):
            if adj[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            H += 1

    H_dist[H] += 1

# Results
print(f"\n{'='*60}")
print(f"EXHAUSTIVE n=6 RESULTS")
print(f"{'='*60}")

max_H = max(H_dist.keys())
all_odd = sorted(range(1, max_H + 1, 2))
achievable = sorted(H_dist.keys())
missing = sorted(set(all_odd) - set(achievable))

print(f"\nTotal tournaments: {sum(H_dist.values())} (should be {N})")
print(f"Max H = {max_H}")
print(f"Distinct H values: {len(achievable)}")
print(f"Achievable H: {achievable}")
print(f"Missing odd values: {missing}")

print(f"\nH distribution:")
total_H = 0
total_H2 = 0
for H in sorted(H_dist.keys()):
    c = H_dist[H]
    total_H += H * c
    total_H2 += H * H * c
    print(f"  H={H:3d}: {c:5d} tournaments ({c/N:.6f})")

mean_H = total_H / N
mean_H2 = total_H2 / N
var_H = mean_H2 - mean_H**2

print(f"\nMean(H) = {mean_H:.6f} (expected n!/2^(n-1) = {720/32:.6f})")
print(f"Var(H) = {var_H:.6f}")
print(f"Var/Mean^2 = {var_H/mean_H**2:.10f}")

# Analyze gaps
print(f"\nGap analysis:")
for h in missing:
    # What's special about this value?
    factors = []
    temp = h
    for p in [2, 3, 5, 7, 11, 13]:
        while temp % p == 0:
            factors.append(p)
            temp //= p
    if temp > 1:
        factors.append(temp)
    print(f"  H={h}: factors={factors}, h mod 6 = {h%6}, h mod 8 = {h%8}")

# Check: are all multiples of 7 that are odd and <= max_H missing?
print(f"\nMultiples of 7 analysis:")
for k in range(1, max_H // 7 + 1):
    h = 7 * k
    if h % 2 == 1 and h <= max_H:
        status = "MISSING" if h in missing else f"PRESENT ({H_dist.get(h, 0)} times)"
        print(f"  7*{k} = {h}: {status}")
