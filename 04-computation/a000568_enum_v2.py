#!/usr/bin/env python3
"""
a000568_enum_v2.py — A000568 via enumeration with incremental t/z computation.

Key optimization: maintain running t_val and z_val as we build the partition,
rather than recomputing from scratch for each complete partition.

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


def a000568_enum_v2(n):
    """Compute A000568(n) with incremental t/z tracking."""
    if n <= 1:
        return 1

    gcd_tab = {}
    for i in range(1, n + 1, 2):
        for j in range(i, n + 1, 2):
            g = gcd(i, j)
            gcd_tab[(i, j)] = g
            gcd_tab[(j, i)] = g

    fact = [mpz(1)]
    for i in range(1, n + 1):
        fact.append(fact[-1] * i)

    # LCD = Π_{k odd, k≤n} k^{⌊n/k⌋} · ⌊n/k⌋!
    LCD = mpz(1)
    for k in range(1, n + 1, 2):
        mk = n // k
        LCD *= mpz(k) ** mk * fact[mk]

    odd_parts = list(range(1, n + 1, 2))
    total_sum = mpz(0)

    # Precompute LCD_over_z[k][m] = LCD / (k^m * m!) ... but that depends on the full partition.
    # Instead, precompute z-factors: z_factor[k][m] = k^m * m!
    z_factor_cache = {}
    for k in range(1, n + 1, 2):
        mk = n // k
        for m in range(1, mk + 1):
            z_factor_cache[(k, m)] = mpz(k) ** m * fact[m]

    def enumerate_partitions(remaining, max_part_idx, t_val, z_val):
        """Enumerate with incremental t_val and z_val.

        t_val: current value of t(λ) for the partial partition
        z_val: current value of z_λ for the partial partition
        """
        nonlocal total_sum

        if remaining == 0:
            # Contribution: LCD / z_val * 2^t_val
            total_sum += (LCD // z_val) * (mpz(1) << t_val)
            return

        for pi in range(max_part_idx, -1, -1):
            k = odd_parts[pi]
            if k > remaining:
                continue
            max_m = remaining // k

            for m in range(1, max_m + 1):
                # Incremental t update:
                # Adding (k, m) contributes:
                #   self: m(m-1)k/2 + m(k-1)/2
                #   cross: m * Σ_{(k', m') already placed} m' * gcd(k, k')
                # But we don't have access to previously placed parts here...
                # We need to track the cross contribution somehow.
                pass

        # Actually, we need to know the "weighted gcd sum" against existing parts.
        # Let's track that as a parameter.

    # Revised approach: track a "gcd profile" vector
    # For each odd d that divides some future k, track Σ_{placed k: d|k} m_k
    # Then the cross contribution of adding (k, m) is m * Σ_{d|k} φ(d) * S_d
    # ...but this is exactly the DP state! We've come full circle.

    # Alternative: track the gcd contribution directly as a number.
    # When adding (k, m), the cross contribution is m * Σ_{(k', m') placed} m' * gcd(k, k')
    # We can maintain a "gcd accumulator" for each possible future k.

    # Let's just track the list of placed parts and compute cross incrementally.

    def enumerate_v2(remaining, max_part_idx, parts_so_far, t_val, z_val):
        nonlocal total_sum

        if remaining == 0:
            total_sum += (LCD // z_val) * (mpz(1) << t_val)
            return

        for pi in range(max_part_idx, -1, -1):
            k = odd_parts[pi]
            if k > remaining:
                continue
            max_m = remaining // k

            for m in range(1, max_m + 1):
                # Self contribution to t
                dt_self = m * (m - 1) * k // 2 + m * (k - 1) // 2

                # Cross contribution to t
                dt_cross = 0
                for k_prev, m_prev in parts_so_far:
                    dt_cross += m * m_prev * gcd_tab[(k, k_prev)]

                new_t = t_val + dt_self + dt_cross
                new_z = z_val * z_factor_cache[(k, m)]

                parts_so_far.append((k, m))
                enumerate_v2(remaining - m * k, pi - 1, parts_so_far, new_t, new_z)
                parts_so_far.pop()

    enumerate_v2(n, len(odd_parts) - 1, [], 0, mpz(1))

    assert total_sum % LCD == 0
    return int(total_sum // LCD)


known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168}

if __name__ == "__main__":
    print("=" * 70)
    print("A000568 via enumeration with incremental t/z (v2)")
    print("=" * 70, flush=True)

    print("\nValidation:", flush=True)
    for n_val in range(13):
        t0 = perf_counter()
        result = a000568_enum_v2(n_val)
        dt = perf_counter() - t0
        expected = known.get(n_val, "?")
        match = "OK" if result == expected else f"FAIL (got {result}, expected {expected})"
        print(f"  a({n_val:2d}) = {result:>15}  [{dt:.4f}s]  {match}", flush=True)

    print("\nBenchmark:", flush=True)
    for n_test in [20, 30, 40, 50, 60, 70, 76, 80, 90, 100]:
        t0 = perf_counter()
        result = a000568_enum_v2(n_test)
        dt = perf_counter() - t0
        s = str(result)
        print(f"  a({n_test:3d}) ({len(s):4d} digits)  [{dt:8.3f}s]", flush=True)
