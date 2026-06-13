#!/usr/bin/env python3
"""
Problem Set 7 — Problem 3: Optimized implementation.

Key optimizations:
  1. GCD precomputation table
  2. LCD scaling (avoid Fraction entirely)
  3. Bucket accumulation (group by exponent, multiply by 2^t once)
  4. Integer-only arithmetic throughout
"""

from math import factorial, gcd
import time

def enumerate_odd_partitions(n):
    results = []
    parts = list(range(n if n % 2 else n-1, 0, -2))
    def rec(idx, rem, cur):
        if rem == 0: results.append(list(cur)); return
        if idx >= len(parts): return
        k = parts[idx]
        if k > rem: rec(idx+1, rem, cur); return
        rec(idx+1, rem, cur)
        for m in range(1, rem//k + 1):
            rec(idx+1, rem-m*k, cur + [(k,m)])
    rec(0, n, [])
    return results

def T_optimized(n):
    """Optimized computation with bucket accumulation."""
    if n <= 1:
        return 1

    partitions = enumerate_odd_partitions(n)

    # GCD table
    G = [[0]*(n+1) for _ in range(n+1)]
    for i in range(1, n+1):
        for j in range(1, n+1):
            G[i][j] = gcd(i, j)

    # Compute z_lambda and c(lambda) for each partition
    data = []  # (z, c) pairs
    for lam in partitions:
        z = 1
        for l, m in lam:
            z *= l**m * factorial(m)
        c = 0
        for i, (li, mi) in enumerate(lam):
            c += mi * (li // 2)
            c += mi * (mi - 1) // 2 * li
            for j in range(i+1, len(lam)):
                lj, mj = lam[j]
                c += mi * mj * G[li][lj]
        data.append((z, c))

    # LCD scaling: L = lcm of all z values
    from math import lcm
    L = 1
    for z, c in data:
        L = lcm(L, z)

    # Bucket accumulation
    buckets = {}  # exponent -> accumulated integer coefficient
    for z, c in data:
        coeff = L // z
        buckets[c] = buckets.get(c, 0) + coeff

    # Final sum: Σ bucket[t] * 2^t, then divide by L
    total = sum(coeff * (2**t) for t, coeff in buckets.items())
    return total // L

# Verify
KNOWN = {0:1, 1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880,
         9:191536, 10:9733056, 11:903753248, 12:154108311168,
         13:48542114686912, 14:28401423719122304,
         15:31021002160355166848}

print("Optimized computation — verification:")
for n in sorted(KNOWN):
    t0 = time.time()
    result = T_optimized(n)
    dt = time.time() - t0
    ok = "✓" if result == KNOWN[n] else "✗"
    print(f"  T({n:2d}) = {result:>25d}  {ok}  ({dt:.4f}s)")

print("\nBenchmark at larger n:")
for n in [20, 30, 40, 50, 60]:
    t0 = time.time()
    result = T_optimized(n)
    dt = time.time() - t0
    ndigits = len(str(result))
    nparts = len(enumerate_odd_partitions(n))
    print(f"  T({n:3d}): {ndigits:5d} digits, {nparts:7d} partitions, {dt:.3f}s")
    if dt > 30:
        break
