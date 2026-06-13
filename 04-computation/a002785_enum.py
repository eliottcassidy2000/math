#!/usr/bin/env python3
"""
A002785: Self-complementary oriented graphs (self-converse tournaments).

PARI formula from OEIS (Andrew Howroyd, Sep 2018):
  permcount(v) = s!/m where s=sum(v), m=prod(v_i * freq_running)
  edges(v) = 2*sum_{i<j} gcd(v_i,v_j) + sum(v_i)
  a(n) = (1/n!) * sum_{p partition of floor(n/2), all odd parts}
         permcount(2*p) * 2^edges(p) * (n*2^#p if n odd, 1 if even)

Compressed (part, multiplicity) form:
  edges = sum_r m_r^2 * r + 2 * sum_{r<s} m_r * m_s * gcd(r,s)
  z_{2λ} = prod (2r)^m_r * m_r!
  permcount(2λ) = (2k)! / z_{2λ}

Uses same partition enumeration as A000568 but over partitions of n/2.

Author: opus-2026-03-08-S48
"""

from math import gcd, factorial
from fractions import Fraction
from time import perf_counter


def a002785(n):
    """Compute A002785(n)."""
    if n <= 1:
        return 1, 1

    k = n // 2
    is_odd_n = n % 2 == 1

    odd_parts = [i for i in range(1, k + 1, 2)]

    total_sum = Fraction(0)
    count = 0
    fact_n = factorial(n)
    fact_2k = factorial(2 * k)

    def enumerate_partitions(remaining, max_part_idx, parts_so_far):
        nonlocal total_sum, count
        if remaining == 0:
            count += 1

            # edges = sum_r m_r^2 * r + 2 * sum_{r<s} m_r * m_s * gcd(r,s)
            edge_val = 0
            for i, (r, mr) in enumerate(parts_so_far):
                edge_val += mr * mr * r
                for j in range(i + 1, len(parts_so_far)):
                    s, ms = parts_so_far[j]
                    edge_val += 2 * mr * ms * gcd(r, s)

            # z_{2λ} = prod (2*r)^m_r * m_r!
            z_2lam = 1
            num_parts_total = 0
            for r, mr in parts_so_far:
                z_2lam *= (2 * r) ** mr * factorial(mr)
                num_parts_total += mr

            # permcount(2λ) = (2k)! / z_{2λ}
            extra = n * (1 << num_parts_total) if is_odd_n else 1

            contrib = Fraction(fact_2k * (1 << edge_val) * extra,
                               z_2lam * fact_n)
            total_sum += contrib
            return

        for pi in range(max_part_idx, -1, -1):
            part = odd_parts[pi]
            if part > remaining:
                continue
            max_m = remaining // part
            for m in range(1, max_m + 1):
                parts_so_far.append((part, m))
                enumerate_partitions(remaining - m * part, pi - 1, parts_so_far)
                parts_so_far.pop()

    enumerate_partitions(k, len(odd_parts) - 1, [])

    result = int(total_sum)
    return result, count


# Known OEIS values (from OEIS A002785, verified)
known = {1: 1, 2: 1, 3: 2, 4: 2, 5: 8, 6: 12, 7: 88, 8: 176, 9: 2752,
         10: 8784, 11: 279968, 12: 1492288, 13: 95458560, 14: 872687552}

if __name__ == "__main__":
    print("A002785: Self-complementary oriented graphs")
    print("=" * 60)

    for n in range(1, 101):
        t0 = perf_counter()
        result, nparts = a002785(n)
        dt = perf_counter() - t0
        expected = known.get(n, "?")
        match = "OK" if result == expected else f"FAIL (expected {expected})"
        if n <= 20 or n % 10 == 0:
            print(f"  a({n:2d}) = {result}  [{dt:.4f}s, {nparts} parts]  {match}")

    print("\nDone.")
