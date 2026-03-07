#!/usr/bin/env python3
"""
H=21 gap mechanism at n=7.

At n=7, H = 1 + sw(1) + sw(2) + sw(3) where:
  sw(1) = 2*t3 (directed 3-cycles)
  sw(2) = 2*t5 + 4*p33 (5-cycles + disjoint 3-cycle pairs)
  sw(3) = 2*t7 + 4*p53 + 8*p333 (7-cycles + disjoint (5,3)-pairs + triple 3-cycle triples)

For H=21: sw(1) + sw(2) + sw(3) = 20

Let's enumerate all tournaments at n=7 and find the achievable
(sw1, sw2, sw3) triples near the target sum of 20.

But n=7 has 2^21 = 2M tournaments — too many for full sw analysis
(need cycle counting). Let me sample and use a faster approach.

Actually, let me just focus on the t3 constraint since it determines sw(1)=2*t3,
and check what range of H values each t3 allows.

From exhaustive n=7 data (h21_n7_fast.py ran all 2M tournaments):
H=19 exists, H=21 does not, H=23 exists.

opus-2026-03-07-S39
"""
from collections import defaultdict
import random


def held_karp(n, adj_bits):
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S & (1<<v)): continue
            if dp[S][v] == 0: continue
            for u in range(n):
                if S & (1<<u): continue
                if adj_bits[v] & (1<<u):
                    dp[S|(1<<u)][u] += dp[S][v]
    return sum(dp[(1<<n)-1][v] for v in range(n))


def count_t3(n, adj_bits):
    """Count directed 3-cycles."""
    t3 = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                ab = (adj_bits[a] >> b) & 1
                bc = (adj_bits[b] >> c) & 1
                ca = (adj_bits[c] >> a) & 1
                if ab and bc and ca:
                    t3 += 1
                ac = (adj_bits[a] >> c) & 1
                cb = (adj_bits[c] >> b) & 1
                ba = (adj_bits[b] >> a) & 1
                if ac and cb and ba:
                    t3 += 1
    return t3


# === Exhaustive n=7: H by t3 ===
print("=== n=7 exhaustive: H values by t3 ===")
n = 7
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)

t3_to_H = defaultdict(set)
H_to_t3 = defaultdict(set)
H_count = defaultdict(int)

for bits in range(1 << m):
    adj = [0]*n
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            adj[j] |= (1 << i)
        else:
            adj[i] |= (1 << j)

    H = held_karp(n, adj)
    t3 = count_t3(n, adj)
    t3_to_H[t3].add(H)
    H_to_t3[H].add(t3)
    H_count[H] += 1

    if bits % 500000 == 0 and bits > 0:
        print(f"  ... processed {bits}/{1<<m}")

print(f"\nt3 -> achievable H (values near 21):")
for t3 in sorted(t3_to_H.keys()):
    H_near = sorted([h for h in t3_to_H[t3] if 15 <= h <= 27])
    if H_near:
        print(f"  t3={t3}: H near 21 = {H_near}")

print(f"\nH values near 21:")
for h in range(15, 28):
    if h in H_to_t3:
        print(f"  H={h}: count={H_count[h]}, t3 in {sorted(H_to_t3[h])}")
    else:
        print(f"  H={h}: NOT ACHIEVABLE")

# === What t3 values exist? ===
print(f"\nAll achievable t3 values: {sorted(t3_to_H.keys())}")
print(f"Max t3: {max(t3_to_H.keys())}")
