#!/usr/bin/env python3
"""
A000665: Number of 3-uniform hypergraphs on n unlabeled nodes.

Formula (Andrew Howroyd, OEIS):
  a(n) = (1/n!) * sum_{p partition of n} permcount(p) * 2^edges(p)

where edges(p) counts the number of orbits of 3-element subsets
fixed by a permutation with cycle type p:

  edges(p) = sum_i ceil((p_i-1)*(p_i-2)/6)
           + sum_{i>j} [gcd(c,d)*(c+d-2+((c-d)/gcd(c,d))%2)/2
                        + sum_{k<j} c*d*p_k/lcm(c,d,p_k)]

where the sums are over all parts of p (with repetition).

OEIS has 29 terms (n=0..28). We aim to extend significantly.

Author: opus-2026-03-08-S48
"""

from math import gcd, factorial, ceil
from fractions import Fraction
from time import perf_counter
import sys


def lcm(a, b):
    return a * b // gcd(a, b)


def lcm3(a, b, c):
    return lcm(lcm(a, b), c)


def permcount(v):
    """Permutation count: s!/m where s=sum(v), m=prod(v_i * running_freq)."""
    m = 1
    s = 0
    k = 0
    for i, t in enumerate(v):
        if i > 0 and t == v[i - 1]:
            k += 1
        else:
            k = 1
        m *= t * k
        s += t
    return factorial(s) // m


def edges_3uniform(p):
    """Count fixed 3-element subsets for a permutation with cycle type p.

    p is the full partition vector (parts listed with repetition, sorted).
    """
    n = len(p)
    # Single-cycle contributions: 3-subsets within one cycle of length c
    result = sum(ceil((p[i] - 1) * (p[i] - 2) / 6) for i in range(n))

    # Pair and triple contributions
    for i in range(1, n):
        for j in range(i):
            c, d = p[i], p[j]
            g = gcd(c, d)
            # Pair contribution: 2-subsets from cycle i, 1 from cycle j (and vice versa)
            pair_term = g * (c + d - 2 + ((c - d) // g) % 2) // 2
            result += pair_term

            # Triple contribution: 1 from each of three distinct cycles
            for k in range(j):
                result += c * d * p[k] // lcm3(c, d, p[k])

    return result


def a000665(n):
    """Compute A000665(n) via Burnside/Polya enumeration."""
    if n <= 2:
        return 1, 1

    total = Fraction(0)
    count = 0
    fact_n = factorial(n)

    def enumerate_partitions(remaining, max_part, parts):
        nonlocal total, count
        if remaining == 0:
            count += 1
            pc = permcount(parts)
            e = edges_3uniform(parts)
            total += Fraction(pc * (1 << e), fact_n)
            return

        for p in range(min(remaining, max_part), 0, -1):
            parts.append(p)
            enumerate_partitions(remaining - p, p, parts)
            parts.pop()

    enumerate_partitions(n, n, [])

    result = int(total)
    assert total == result, f"Non-integer result: {total}"
    return result, count


# Known OEIS values (n=0..28)
known = {
    0: 1, 1: 1, 2: 1, 3: 2, 4: 5, 5: 34, 6: 2136, 7: 7013320,
    8: 1788782616656
}

if __name__ == "__main__":
    print("A000665: 3-uniform hypergraphs on n unlabeled nodes")
    print("=" * 60)

    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 30

    for n in range(0, max_n + 1):
        t0 = perf_counter()
        result, nparts = a000665(n)
        dt = perf_counter() - t0

        expected = known.get(n, "?")
        if expected != "?":
            match = "OK" if result == expected else f"FAIL (expected {expected})"
        else:
            match = ""

        print(f"  a({n:3d}) = {result}  [{dt:.4f}s, {nparts} parts]  {match}")

        if dt > 300:  # 5 minute limit per term
            print(f"  Stopping: n={n} took {dt:.1f}s")
            break

    print("\nDone.")
