#!/usr/bin/env python3
"""
Compute H(P(p)) for Paley primes p = 3 mod 4.
Paley primes: 3, 7, 11, 19, 23, 31, 43, 47, ...

For p=23 the DP is O(2^23 * 23) ~ 193M operations — might take a few minutes.
p=31 would be O(2^31 * 31) ~ 66B — too much.

Instance: opus-2026-03-05-S9
"""

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

def quadratic_residues(p):
    """Return set of quadratic residues mod p (excluding 0)."""
    return {(k*k) % p for k in range(1, p)} - {0}

def build_paley(p):
    """Build Paley tournament on p vertices."""
    qr = quadratic_residues(p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                T[i][j] = 1
    return T

# Paley primes p ≡ 3 (mod 4)
paley_primes = [3, 7, 11, 19, 23]

oeis_known = {3: 3, 7: 189, 11: 95095}

for p in paley_primes:
    print(f"Computing H(P({p}))...", flush=True)
    T = build_paley(p)
    h = count_ham_dp(T, p)
    match_str = ""
    if p in oeis_known:
        match_str = f" [a({p})={oeis_known[p]}, match={h==oeis_known[p]}]"
    print(f"  H(P({p})) = {h}{match_str}")
    
    # Check: n!/2^(n-1) lower bound
    import math
    lb = math.factorial(p) / (2**(p-1))
    ratio = h / lb
    print(f"  n!/2^(n-1) = {lb:.1f}, ratio = {ratio:.4f}")
