#!/usr/bin/env python3
"""
a000568_enum.py — A000568 via direct enumeration of odd partitions.

Instead of DP, enumerate all partitions of n into odd parts and sum
2^{t(λ)} / z_λ directly.

For n=100, there are only 444,793 such partitions.

Author: opus-2026-03-07-S46f
"""

import gmpy2
from gmpy2 import mpz
from time import perf_counter
from math import gcd


def euler_totient(n):
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result


def a000568_enum(n):
    """Compute A000568(n) by enumerating all odd partitions."""
    if n <= 1:
        return 1

    # Precompute gcd table and phi
    phi_cache = {}
    for i in range(1, n + 1):
        phi_cache[i] = euler_totient(i)

    gcd_tab = {}
    for i in range(1, n + 1, 2):
        for j in range(i, n + 1, 2):
            gcd_tab[(i, j)] = gcd(i, j)
            gcd_tab[(j, i)] = gcd_tab[(i, j)]

    # Precompute factorials
    fact = [mpz(1)]
    for i in range(1, n + 1):
        fact.append(fact[-1] * i)

    # Compute LCD = Π_{k odd, k≤n} k^{⌊n/k⌋} · ⌊n/k⌋!
    LCD = mpz(1)
    for k in range(1, n + 1, 2):
        mk = n // k
        LCD *= mpz(k) ** mk * fact[mk]

    # Enumerate partitions of n into odd parts (non-increasing)
    # For each partition, compute t(λ) and z_λ, add LCD/z_λ * 2^t to accumulator
    total_sum = mpz(0)
    count = 0

    odd_parts = list(range(1, n + 1, 2))

    def enumerate_partitions(remaining, max_part_idx, parts_so_far):
        """Recursively enumerate partitions, computing contribution on the fly."""
        nonlocal total_sum, count

        if remaining == 0:
            # Compute t(λ) and z_λ from parts_so_far = [(k, m), ...]
            count += 1

            t_val = 0
            z_val = mpz(1)

            for i, (k_i, m_i) in enumerate(parts_so_far):
                # Self-contribution
                t_val += m_i * (m_i - 1) * k_i // 2 + m_i * (k_i - 1) // 2
                z_val *= mpz(k_i) ** m_i * fact[m_i]

                # Cross-contribution with later parts
                for j in range(i + 1, len(parts_so_far)):
                    k_j, m_j = parts_so_far[j]
                    t_val += m_i * m_j * gcd_tab[(k_i, k_j)]

            # Contribution = LCD / z_λ * 2^t
            contrib = (LCD // z_val) * (mpz(1) << t_val)
            total_sum += contrib
            return

        for pi in range(max_part_idx, -1, -1):
            k = odd_parts[pi]
            if k > remaining:
                continue
            max_m = remaining // k
            for m in range(1, max_m + 1):
                parts_so_far.append((k, m))
                enumerate_partitions(remaining - m * k, pi - 1, parts_so_far)
                parts_so_far.pop()

    enumerate_partitions(n, len(odd_parts) - 1, [])

    assert total_sum % LCD == 0
    return int(total_sum // LCD)


known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168}

if __name__ == "__main__":
    print("=" * 70)
    print("A000568 via direct enumeration of odd partitions")
    print("=" * 70, flush=True)

    print("\nValidation:", flush=True)
    for n_val in range(13):
        t0 = perf_counter()
        result = a000568_enum(n_val)
        dt = perf_counter() - t0
        expected = known.get(n_val, "?")
        match = "OK" if result == expected else f"FAIL (got {result}, expected {expected})"
        print(f"  a({n_val:2d}) = {result:>15}  [{dt:.4f}s]  {match}", flush=True)

    print("\nBenchmark:", flush=True)
    for n_test in [20, 30, 40, 50, 60]:
        t0 = perf_counter()
        result = a000568_enum(n_test)
        dt = perf_counter() - t0
        s = str(result)
        print(f"  a({n_test:3d}) ({len(s):4d} digits)  [{dt:8.3f}s]", flush=True)
