#!/usr/bin/env python3
"""
paley11_count_only.py — opus-2026-03-13-S71

Just COUNT directed paths from vertex 0 in Paley P_11.
No storage, no matrices — just counting.
This tells us |A_m|/11 for each m.
"""

import time

p = 11
QR = {a*a % p for a in range(1, p)}
print(f"P_{p}: QR = {sorted(QR)}")

# Adjacency as set for fast lookup
adj = {}
for i in range(p):
    adj[i] = set()
    for s in QR:
        adj[i].add((i + s) % p)

# Count paths from vertex 0 by DFS
def count_paths_from_0(max_m):
    counts = [0] * (max_m + 1)

    def dfs(last, depth, visited_mask):
        counts[depth] += 1
        if depth >= max_m:
            return
        for v in adj[last]:
            if not (visited_mask & (1 << v)):
                dfs(v, depth + 1, visited_mask | (1 << v))

    dfs(0, 0, 1 << 0)
    return counts

t0 = time.time()
counts = count_paths_from_0(10)
t1 = time.time()

print(f"\nPath counts from vertex 0 (time={t1-t0:.1f}s):")
for m in range(11):
    total = counts[m] * p
    print(f"  m={m:2d}: from_0={counts[m]:10d}, total={total:12d}")

print(f"\n|A_m|/11 = {counts}")
print(f"|A_m|    = {[c*11 for c in counts]}")

# Compare with P_7
print(f"\nP_7 |A_m|: [7, 21, 63, 147, 273, 315, 189]")
print(f"P_7 |A_m|/7: [1, 3, 9, 21, 39, 45, 27]")

print("\nDONE.")
