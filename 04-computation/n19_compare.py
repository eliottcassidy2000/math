#!/usr/bin/env python3
"""
n19_compare.py -- Compare Paley T_19 vs cyclic interval C_19 (and others).
19 ≡ 3 mod 4, so Paley exists. THM-135 claims C_19 beats T_19.
"""

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

n = 19

# Paley T_19: QR mod 19
paley_S = sorted({i*i % n for i in range(1, n)} - {0})
print(f"Paley S = {paley_S}")

# Cyclic interval C_19: S = {1,...,9}
cyclic_S = list(range(1, (n+1)//2))
print(f"Cyclic S = {cyclic_S}")

# All-odd
allodd_S = [k for k in range(1, n, 2)]
print(f"All-odd S = {allodd_S}")

candidates = [
    ("Paley T_19",    paley_S),
    ("Cyclic C_19",   cyclic_S),
    ("All-odd n=19",  allodd_S),
]

results = []
for name, S in candidates:
    t0 = time.time()
    adj = circulant(n, S)
    H = H_dp(adj)
    dt = time.time() - t0
    results.append((H, name, S, dt))
    print(f"  {name}: H={H:,}  ({dt:.1f}s)")

results.sort(reverse=True)
print(f"\nRanking:")
for H, name, S, dt in results:
    print(f"  {H:>20,}  {name}")

best_H, best_name = results[0][0], results[0][1]
paley_H = next(H for H, name, *_ in results if name == "Paley T_19")
ratio = best_H / paley_H
print(f"\n{best_name} / Paley = {ratio:.4f}×")
