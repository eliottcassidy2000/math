#!/usr/bin/env python3
"""
paley_vs_cyclic.py -- At which n does Paley stop being optimal?

Tests all odd n=3..17 (including non-primes):
  Paley: only for p ≡ 3 mod 4 prime
  Cyclic interval C_n: S = {1,...,(n-1)/2}
"""

from sympy import isprime
import time

def qr(p):
    return sorted({i*i % p for i in range(1, p)} - {0})

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

print(f"{'n':>4} {'mod4':>5} {'Paley':>22} {'Cyclic':>22} {'winner':>10} {'ratio':>8}")
print("-"*80)

for n in range(3, 20, 2):
    t0 = time.time()
    # Cyclic interval
    cyc_S = list(range(1, (n+1)//2))
    H_cyc = H_dp(circulant(n, cyc_S))

    # Paley (if applicable)
    has_paley = isprime(n) and n % 4 == 3
    if has_paley:
        H_pal = H_dp(circulant(n, qr(n)))
        ratio = H_cyc / H_pal
        winner = "Cyclic" if H_cyc > H_pal else ("Paley" if H_pal > H_cyc else "TIE")
    else:
        H_pal = None
        ratio = float('nan')
        winner = "Cyclic(only)"

    pal_str = f"{H_pal:,}" if H_pal else "N/A"
    print(f"{n:>4} {n%4:>5}   {pal_str:>22} {H_cyc:>22,} {winner:>10} {ratio:>8.5f}"
          if H_pal else
          f"{n:>4} {n%4:>5}   {'N/A':>22} {H_cyc:>22,} {winner:>10}")
