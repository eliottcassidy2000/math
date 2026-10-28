#!/usr/bin/env python3
"""Audit the exact abundance classification in the seven-proper-divisor layer."""

import argparse
import json
from math import isqrt


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def factor(n):
    factors = []
    p = 2
    while p * p <= n:
        a = 0
        while n % p == 0:
            n //= p
            a += 1
        if a:
            factors.append((p, a))
        p += 1 if p == 2 else 2
    if n > 1:
        factors.append((n, 1))
    return factors


def proper_divisors(n):
    ds = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            ds.update((d, n // d))
    ds.discard(n)
    return sorted(ds)


def prime_sieve(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit + 1:p] = b"\x00" * ((limit - start) // p + 1)
    return [p for p in range(2, limit + 1) if sieve[p]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=100000)
    args = parser.parse_args()
    limit = args.limit
    require(1000 <= limit <= 2000000, "limit must lie in [1000,2000000]")

    sums, counts = [0] * (limit + 1), [0] * (limit + 1)
    for d in range(1, limit // 2 + 1):
        for n in range(d * 2, limit + 1, d):
            sums[n] += d
            counts[n] += 1
    primes = prime_sieve(limit)
    exceptions = [24, 40, 54, 56, 70, 88, 104]
    family = [6 * q for q in primes if 5 <= q <= limit // 6]
    require(not set(exceptions) & set(family), "claimed disjoint sets overlap")
    expected = sorted(exceptions + family)
    actual = [n for n in range(1, limit + 1) if counts[n] == 7 and sums[n] > n]
    perfect = [n for n in range(1, limit + 1) if counts[n] == 7 and sums[n] == n]
    require(actual == expected, "seven-term abundance classification mismatch")
    require(not perfect, "unexpected perfect seven-term source")

    square_edge_sources = [n for n in actual if isqrt(n + sums[n]) ** 2 == n + sums[n]]
    prime_set = set(primes)
    quadratic_sources = []
    for k in range(1, isqrt((limit + 6) // 72) + 1):
        q = 12 * k * k - 1
        if q in prime_set:
            n = 6 * q
            require(n <= limit and sums[n] == 72 * k * k + 6, f"quadratic target {k}")
            require(n + sums[n] == (12 * k) ** 2, f"quadratic square total {k}")
            require(sums[2 * q] == 12 * k * k + 2, f"semiprime target {k}")
            require(n + sums[n] == 4 * (2 * q + sums[2 * q]), f"square transfer {k}")
            quadratic_sources.append(n)
    require(square_edge_sources == sorted([70] + quadratic_sources), "square-edge classification")
    require(min(square_edge_sources) == 66, "minimum square-edge source")
    require(min(sums[n] for n in square_edge_sources) == 74, "minimum square-edge target")

    shape_counts = {"p^7": 0, "p^3*q": 0, "p*q*r": 0}
    for n in range(1, limit + 1):
        if counts[n] == 7:
            shape = tuple(sorted(a for _, a in factor(n)))
            require(shape in ((7,), (1, 3), (1, 1, 1)), f"shape at {n}")
            shape_counts[{(7,): "p^7", (1, 3): "p^3*q", (1, 1, 1): "p*q*r"}[shape]] += 1
        if n <= min(2000, limit):
            ds = proper_divisors(n)
            require(sum(ds) == sums[n] and len(ds) == counts[n], f"pair oracle at {n}")

    # Check every cutoff, not only a selected final count.
    prime_prefix = [0] * (limit + 1)
    for p in primes:
        prime_prefix[p] = 1
    for n in range(1, limit + 1):
        prime_prefix[n] += prime_prefix[n - 1]
    abundance_prefix = 0
    for cutoff in range(limit + 1):
        if cutoff > 0 and counts[cutoff] == 7 and sums[cutoff] > cutoff:
            abundance_prefix += 1
        formula = max(prime_prefix[cutoff // 6] - 2, 0) + sum(n <= cutoff for n in exceptions)
        require(abundance_prefix == formula, f"counting formula at {cutoff}")
        if cutoff >= 104:
            require(abundance_prefix == prime_prefix[cutoff // 6] + 5, f"simplified count {cutoff}")

    named = {}
    for n in (6, 12, 18, 24, 30, 40, 42, 54, 56, 66, 70, 88, 104, 105, 110, 128, 136, 210):
        ds = proper_divisors(n)
        require(sum(ds) == sums[n] and len(ds) == counts[n], f"named oracle at {n}")
        named[str(n)] = {"factorization": factor(n), "proper_divisors": ds,
                         "s": sums[n], "terms": counts[n], "excess": sums[n] - n}
    require(named["70"]["terms"] == 7 and named["70"]["s"] == 74, "corrected 70 control")
    require(named["210"]["terms"] == 15 and named["210"]["s"] == 366, "rejected 210 control")
    require(named["54"]["excess"] == 12 and 54 not in family, "54 is not six times a prime")
    require(named["12"]["terms"] == named["18"]["terms"] == 5, "coprimality term guard")
    require(named["12"]["excess"] == 4 and named["18"]["excess"] == 3, "coprimality excess guard")

    adjunction_tests = 0
    for m in (6, 28, 496):
        dm = proper_divisors(m)
        require(sum(dm) == m, f"perfect seed {m}")
        for q in primes:
            if q > 100:
                break
            if m % q == 0:
                continue
            dn = proper_divisors(m * q)
            require(len(dn) == 2 * len(dm) + 1, f"adjunction count {m},{q}")
            require(sum(dn) - m * q == 2 * m, f"adjunction excess {m},{q}")
            adjunction_tests += 1

    print(json.dumps({
        "status": "passed", "source_universe": [1, limit],
        "paired_divisor_universe": [1, min(2000, limit)],
        "all_integer_counting_cutoffs": [0, limit], "seven_term_shape_counts": shape_counts,
        "abundant_count": len(actual), "abundant_exceptions": exceptions,
        "six_times_prime_family_count": len(family), "perfect_sources": perfect,
        "square_edge_sources": square_edge_sources,
        "smallest_square_edge_source": {"source": 66, "target": 78},
        "smallest_square_edge_target": {"source": 70, "target": 74},
        "prime_adjunction_tests": adjunction_tests, "named_controls": named,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
