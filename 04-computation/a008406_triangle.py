#!/usr/bin/env python3
"""
A008406: Triangle T(n,k) = number of simple graphs on n nodes with k edges.

Uses the Burnside/Polya approach:
  T(n,k) = [x^k] (1/n!) * sum_{lambda ⊢ n} permcount(lambda) * (1+x)^{edges(lambda)}

where edges(lambda) = sum_{i<j} gcd(v_i,v_j) + sum_i floor(v_i/2)
for expanded cycle type (v_1,...,v_p).

In compressed form [(r_a, m_a)]:
  edges = sum_a C(m_a,2)*r_a + sum_a m_a*floor(r_a/2) + sum_{a<b} m_a*m_b*gcd(r_a,r_b)

The polynomial (1+x)^edges has degree C(n,2), and we track it exactly.

Author: opus-2026-03-08-S48
"""

from math import gcd, factorial
from fractions import Fraction
from time import perf_counter
import sys


def binomial(n, k):
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    k = min(k, n - k)
    result = 1
    for i in range(k):
        result = result * (n - i) // (i + 1)
    return result


def compute_edges_simple(compressed):
    """Compute edges(lambda) for simple graph enumeration (A000088).
    edges = sum_{i<j} gcd(v_i, v_j) + sum floor(v_i/2)
    In compressed form: same-part pairs + cross-part pairs + self terms.
    """
    result = 0
    for r, m in compressed:
        # Self term: m * floor(r/2)
        result += m * (r // 2)
        # Same-part pairs: C(m,2) * gcd(r,r) = C(m,2) * r
        result += m * (m - 1) // 2 * r
    # Cross-part pairs
    d = len(compressed)
    for i in range(d):
        ri, mi = compressed[i]
        for j in range(i + 1, d):
            rj, mj = compressed[j]
            result += mi * mj * gcd(ri, rj)
    return result


def poly_mul_scalar_shift(poly, edges):
    """Multiply polynomial by (1+x)^edges.

    Uses the fact that (1+x)^edges has binomial coefficients.
    The result has degree = old_degree + edges.
    """
    max_deg = len(poly) - 1 + edges
    result = [0] * (max_deg + 1)

    for d in range(len(poly)):
        if poly[d] == 0:
            continue
        # Multiply poly[d] * x^d by (1+x)^edges = sum C(edges, j) x^j
        for j in range(edges + 1):
            result[d + j] += poly[d] * binomial(edges, j)

    return result


def a008406_row(n):
    """Compute row n of A008406: T(n,k) for k=0..C(n,2)."""
    max_edges = n * (n - 1) // 2
    fact_n = factorial(n)

    # Polynomial accumulator: poly[k] = sum of (n!/z_lambda) * C(edges, k) / n!
    # We work with Fraction to stay exact
    poly = [Fraction(0)] * (max_edges + 1)

    count = 0

    def enumerate_partitions(remaining, max_part, parts_so_far):
        nonlocal count
        if remaining == 0:
            count += 1
            edges = compute_edges_simple(parts_so_far)
            z = 1
            for r, m in parts_so_far:
                z *= r ** m * factorial(m)
            coeff = Fraction(1, z)  # permcount/n! = 1/z

            # Add coeff * (1+x)^edges to poly
            for k in range(edges + 1):
                poly[k] += coeff * binomial(edges, k)
            return

        for p in range(min(remaining, max_part), 0, -1):
            max_m = remaining // p
            for m in range(1, max_m + 1):
                parts_so_far.append((p, m))
                enumerate_partitions(remaining - m * p, p - 1, parts_so_far)
                parts_so_far.pop()

    enumerate_partitions(n, n, [])

    # Convert to integers
    row = []
    for k in range(max_edges + 1):
        val = poly[k]
        assert val.denominator == 1, f"Non-integer at k={k}: {val}"
        row.append(int(val))

    return row, count


if __name__ == "__main__":
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 12

    # Known A008406 rows (from OEIS)
    known = {
        1: [1],
        2: [1, 1],
        3: [1, 1, 1, 1],
        4: [1, 1, 2, 3, 2, 1, 1],
        5: [1, 1, 2, 4, 6, 6, 6, 4, 2, 1, 1],
    }

    print("A008406: Triangle T(n,k) = graphs on n nodes with k edges")
    print("=" * 60)

    for n in range(1, max_n + 1):
        t0 = perf_counter()
        row, nparts = a008406_row(n)
        dt = perf_counter() - t0

        expected = known.get(n)
        if expected:
            match = "OK" if row == expected else f"FAIL (expected {expected})"
        else:
            match = ""

        print(f"  n={n}: {row}  [{dt:.3f}s, {nparts} parts]  {match}")

        if dt > 60:
            print("  Stopping: too slow")
            break

    # Generate b-file format
    print("\n=== B-file format ===")
    idx = 0
    for n in range(1, min(max_n + 1, 30)):
        row, _ = a008406_row(n) if n > max_n else (a008406_row(n)[0] if n <= max_n else None, 0)
        # Already computed above, but redo for formatting
        pass

    print("\nDone.")
