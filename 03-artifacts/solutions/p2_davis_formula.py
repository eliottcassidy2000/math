#!/usr/bin/env python3
"""
Problem Set 7 — Problem 2: Implementing Davis's formula.

T(n) = Σ_{λ ⊢ n, all odd} 2^{c(λ)} / z_λ

where c(λ) = Σ_i m_i floor(l_i/2) + Σ_i C(m_i,2) l_i + Σ_{i<j} m_i m_j gcd(l_i,l_j)
and z_λ = Π l_i^{m_i} m_i!
"""

from math import factorial, gcd
from fractions import Fraction
import time

def enumerate_odd_partitions(n):
    """All partitions of n into odd parts, as [(part, mult), ...]."""
    results = []
    parts = list(range(n if n % 2 else n-1, 0, -2))

    def recurse(idx, remaining, current):
        if remaining == 0:
            results.append(list(current))
            return
        if idx >= len(parts):
            return
        k = parts[idx]
        if k > remaining:
            recurse(idx + 1, remaining, current)
            return
        recurse(idx + 1, remaining, current)
        for m in range(1, remaining // k + 1):
            recurse(idx + 1, remaining - m * k, current + [(k, m)])

    recurse(0, n, [])
    return results

def cycle_exponent(partition):
    """c(λ) = orbit count on unordered pairs."""
    c = 0
    for i, (li, mi) in enumerate(partition):
        c += mi * (li // 2)           # within-cycle orbits
        c += mi * (mi - 1) // 2 * li  # same-type pair orbits
        for j in range(i + 1, len(partition)):
            lj, mj = partition[j]
            c += mi * mj * gcd(li, lj) # cross-type pair orbits
    return c

def centralizer_order(partition):
    """z_λ = Π l_i^{m_i} m_i!"""
    z = 1
    for l, m in partition:
        z *= l**m * factorial(m)
    return z

def T(n):
    """Compute T(n) using Davis's formula with exact rational arithmetic."""
    if n <= 1:
        return 1
    total = Fraction(0)
    for lam in enumerate_odd_partitions(n):
        c = cycle_exponent(lam)
        z = centralizer_order(lam)
        total += Fraction(2**c, z)
    return int(total)

# Verify and benchmark
KNOWN = {0:1, 1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880,
         9:191536, 10:9733056, 11:903753248, 12:154108311168}

print("Davis's formula — verification:")
for n in range(13):
    t0 = time.time()
    result = T(n)
    dt = time.time() - t0
    ok = "✓" if result == KNOWN[n] else "✗"
    print(f"  T({n:2d}) = {result:>18d}  {ok}  ({dt:.4f}s)")

print("\nBenchmark:")
for n in [20, 30, 40, 50]:
    t0 = time.time()
    result = T(n)
    dt = time.time() - t0
    nparts = len(enumerate_odd_partitions(n))
    print(f"  T({n}) = {str(result)[:30]}...  ({nparts} partitions, {dt:.3f}s)")
