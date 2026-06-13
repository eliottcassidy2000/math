#!/usr/bin/env python3
"""
Fast check: does H(T)=21 occur at n=7?
Only compute H via Held-Karp. 2^21 = 2,097,152 tournaments.

opus-2026-03-07-S38
"""
from collections import Counter
import time

def held_karp(n, adj):
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S & (1<<v)): continue
            if dp[S][v] == 0: continue
            for u in range(n):
                if S & (1<<u): continue
                if adj[v] & (1<<u):
                    dp[S|(1<<u)][u] += dp[S][v]
    full = (1<<n)-1
    return sum(dp[full][v] for v in range(n))

n = 7
edges = [(i, j) for i in range(n) for j in range(i+1, n)]
m = len(edges)  # 21
total = 1 << m  # 2,097,152

H_counter = Counter()
t0 = time.time()

for bits in range(total):
    adj = [0] * n
    for k, (i, j) in enumerate(edges):
        if bits & (1 << k):
            adj[j] |= (1 << i)
        else:
            adj[i] |= (1 << j)
    H = held_karp(n, adj)
    H_counter[H] += 1

    if bits % 200000 == 0 and bits > 0:
        elapsed = time.time() - t0
        rate = bits / elapsed
        eta = (total - bits) / rate
        print(f"  {bits}/{total} ({100*bits/total:.1f}%), "
              f"elapsed {elapsed:.0f}s, ETA {eta:.0f}s", flush=True)

elapsed = time.time() - t0
print(f"\nDone in {elapsed:.1f}s")
print(f"Total tournaments: {total}")

# Print all achievable H values
print(f"\nAll achievable H values at n=7:")
for h in sorted(H_counter.keys()):
    print(f"  H={h}: {H_counter[h]}")

print(f"\nH=21 achievable: {21 in H_counter}")

# Print values near 21
print(f"\nValues near 21:")
for h in range(17, 30):
    if h in H_counter:
        print(f"  H={h}: {H_counter[h]}")
    else:
        print(f"  H={h}: GAP")
