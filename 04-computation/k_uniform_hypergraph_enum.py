#!/usr/bin/env python3
"""
General k-uniform hypergraph enumeration on n unlabeled nodes.

Formula:
  a(n,k) = (1/n!) * sum_{lambda ⊢ n} permcount(lambda) * 2^{c_k(lambda)}

where c_k(lambda) = number of orbits of <sigma> on k-element subsets
of {1,...,n}, for sigma with cycle type lambda.

c_k is computed via Burnside on the cyclic group <sigma>:
  c_k(lambda) = (1/L) * sum_{d|L} f_k(d, lambda)
where L = lcm(lambda) and f_k(d, lambda) = number of k-subsets fixed by sigma^d.

A k-subset is fixed by sigma^d iff it's a union of complete cycles of sigma^d.
For partition lambda in compressed form (r_a, m_a):
  sigma^d has m_a * gcd(d, r_a) cycles of length r_a / gcd(d, r_a) from part r_a.
  GF: P_d(x) = prod_a (1 + x^{r_a/gcd(d,r_a)})^{m_a * gcd(d,r_a)}
  f_k(d) = [x^k] P_d(x)

Covers:
  k=2: A000088 (simple graphs)
  k=3: A000665 (3-uniform hypergraphs)
  k=4: A051240 (4-uniform hypergraphs)
  k=5: column k=5 of A309858

Author: opus-2026-03-08-S48
"""

from math import gcd, factorial
from fractions import Fraction
from functools import reduce
from time import perf_counter
import sys


def lcm(a, b):
    return a * b // gcd(a, b)


def divisors(n):
    """Return all divisors of n."""
    divs = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            divs.append(i)
            if i != n // i:
                divs.append(n // i)
        i += 1
    return sorted(divs)


def euler_phi(n):
    """Euler's totient function."""
    result = n
    p = 2
    m = n
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            result -= result // p
        p += 1
    if m > 1:
        result -= result // m
    return result


def poly_mul_trunc(p1, p2, max_deg):
    """Multiply two polynomials (as lists), truncate to degree max_deg."""
    result = [0] * (max_deg + 1)
    for i in range(min(len(p1), max_deg + 1)):
        if p1[i] == 0:
            continue
        for j in range(min(len(p2), max_deg + 1 - i)):
            if p2[j] == 0:
                continue
            result[i + j] += p1[i] * p2[j]
    return result


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


def compute_ck(compressed_partition, k):
    """Compute c_k(lambda) = number of orbits on k-subsets.

    compressed_partition: list of (part_value, multiplicity) pairs.
    k: the uniformity parameter.

    Uses Burnside: c_k = (1/L) * sum_{t=0}^{L-1} f_k(sigma^t)
    where f_k(sigma^t) = [x^k] prod_a (1 + x^{r_a/gcd(t,r_a)})^{m_a*gcd(t,r_a)}
    """
    if not compressed_partition:
        return 0

    # Compute L = lcm of all parts
    L = 1
    for r, m in compressed_partition:
        L = lcm(L, r)

    # Optimization: group t values by their gcd pattern
    # Two t values give the same f_k if they have the same (gcd(t, r_a)) for all a.
    # Iterate over all t from 0 to L-1 (L is manageable for n < 80).
    total = 0
    for t in range(L):
        poly = [0] * (k + 1)
        poly[0] = 1

        for r, m in compressed_partition:
            g = gcd(t, r)
            length = r // g  # cycle length in sigma^t
            cnt = m * g      # number of such cycles

            # (1 + x^length)^cnt truncated to degree k
            factor = [0] * (k + 1)
            for j in range(cnt + 1):
                deg = j * length
                if deg > k:
                    break
                factor[deg] = binomial(cnt, j)

            poly = poly_mul_trunc(poly, factor, k)

        f_t = poly[k] if k < len(poly) else 0
        total += f_t

    assert total % L == 0, f"c_k not integer: {total}/{L}"
    return total // L


def a_k_uniform(n, k):
    """Compute the number of k-uniform hypergraphs on n unlabeled nodes."""
    if n < k:
        return 1, 1  # only the empty hypergraph

    total = Fraction(0)
    fact_n = factorial(n)
    count = 0

    def enumerate_partitions(remaining, max_part, parts_so_far):
        nonlocal total, count
        if remaining == 0:
            count += 1
            ck = compute_ck(parts_so_far, k)
            # permcount in compressed form: n! / prod(r^m * m!)
            z = 1
            for r, m in parts_so_far:
                z *= r ** m * factorial(m)
            pc = fact_n // z
            total += Fraction(pc * (1 << ck), fact_n)
            return

        for p in range(min(remaining, max_part), 0, -1):
            max_m = remaining // p
            for m in range(1, max_m + 1):
                parts_so_far.append((p, m))
                enumerate_partitions(remaining - m * p, p - 1, parts_so_far)
                parts_so_far.pop()

    enumerate_partitions(n, n, [])

    result = int(total)
    assert total == result, f"Non-integer result: {total}"
    return result, count


# Known values from OEIS
known_k2 = {0: 1, 1: 1, 2: 2, 3: 4, 4: 11, 5: 34, 6: 156, 7: 1044}  # A000088
known_k3 = {0: 1, 1: 1, 2: 1, 3: 2, 4: 5, 5: 34, 6: 2136, 7: 7013320, 8: 1788782616656}  # A000665
known_k4 = {0: 1, 1: 1, 2: 1, 3: 1, 4: 2, 5: 6, 6: 156, 7: 7013320,
            8: 29281354514767168}  # A051240


if __name__ == "__main__":
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    max_n = int(sys.argv[2]) if len(sys.argv) > 2 else 15

    known = {2: known_k2, 3: known_k3, 4: known_k4}.get(k, {})
    seq_id = {2: "A000088", 3: "A000665", 4: "A051240"}.get(k, f"k={k}")

    print(f"{seq_id}: {k}-uniform hypergraphs on n unlabeled nodes")
    print("=" * 60)

    for n in range(0, max_n + 1):
        t0 = perf_counter()
        result, nparts = a_k_uniform(n, k)
        dt = perf_counter() - t0

        expected = known.get(n, "?")
        if expected != "?":
            match = "OK" if result == expected else f"FAIL (expected {expected})"
        else:
            match = ""

        print(f"  a({n:3d}) = {result}  [{dt:.4f}s, {nparts} parts]  {match}")

        if dt > 300:
            print(f"  Stopping: took {dt:.1f}s")
            break

    print("\nDone.")
