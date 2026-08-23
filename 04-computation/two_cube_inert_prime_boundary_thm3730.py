#!/usr/bin/env python3
"""Finite exact hostile probe for the inert-prime unique two-cube family.

The all-scale proof is in RESEARCH-two-cube-boundary.md. This probe checks the
valuation, injectivity, and singleton-fibre mechanism on a bounded window.
"""

from __future__ import annotations

from hashlib import sha256
from math import isqrt


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def primes_up_to(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for value in range(2, isqrt(limit) + 1):
        if sieve[value]:
            start = value * value
            sieve[start:limit + 1:value] = b"\x00" * (
                (limit - start) // value + 1
            )
    return tuple(index for index, flag in enumerate(sieve) if flag)


def divisors(value):
    low = []
    high = []
    for divisor in range(1, isqrt(value) + 1):
        if value % divisor == 0:
            low.append(divisor)
            if divisor * divisor != value:
                high.append(value // divisor)
    return tuple(low + high[::-1])


def positive_distinct_representations(value):
    answer = []
    for divisor in divisors(value):
        numerator = divisor**3 - value
        denominator = 3 * divisor
        if numerator % denominator:
            continue
        product = numerator // denominator
        discriminant = divisor * divisor - 4 * product
        if discriminant <= 0:
            continue
        gap = isqrt(discriminant)
        if gap * gap != discriminant or (divisor - gap) % 2:
            continue
        left = (divisor - gap) // 2
        right = (divisor + gap) // 2
        if 0 < left < right and left**3 + right**3 == value:
            answer.append((left, right, divisor))
    return tuple(answer)


def stable_digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


def main():
    prime_limit = 5000
    inert_primes = tuple(
        prime for prime in primes_up_to(prime_limit)
        if prime % 3 == 2 and prime >= 5
    )
    census = []
    constructed = {}
    checked_unique = 0

    for prime in inert_primes:
        row_values = []
        for left in range(1, (prime - 1) // 2 + 1):
            right = prime - left
            cofactor = left * left - left * right + right * right
            value = prime * cofactor
            require(
                cofactor == 3 * left * left - 3 * prime * left + prime**2,
                ("cofactor identity", prime, left),
            )
            require(cofactor % prime != 0, ("q mod p", prime, left))
            require(value % prime == 0 and (value // prime) % prime != 0,
                    ("p-adic valuation", prime, left, value))
            require(value == left**3 + right**3,
                    ("cube identity", prime, left))
            row_values.append(value)

            if prime <= 101:
                representations = positive_distinct_representations(value)
                require(
                    representations == ((left, right, prime),),
                    ("singleton divisor fibre", prime, left,
                     representations),
                )
                checked_unique += 1

            previous = constructed.setdefault(value, (prime, left, right))
            require(previous == (prime, left, right),
                    ("cross-prime collision", value, previous,
                     (prime, left, right)))

        require(len(row_values) == (prime - 1) // 2,
                ("row cardinality", prime, len(row_values)))
        require(all(a > b for a, b in zip(row_values, row_values[1:])),
                ("strict row monotonicity", prime))
        census.append((prime, len(row_values), row_values[0], row_values[-1]))

    expected_count = sum((prime - 1) // 2 for prime in inert_primes)
    require(len(constructed) == expected_count,
            ("global cardinality", len(constructed), expected_count))
    semantic = (
        prime_limit,
        len(inert_primes),
        checked_unique,
        tuple(census[:12]),
        tuple(census[-12:]),
        len(constructed),
    )
    print("two-cube inert-prime boundary hostile probe")
    print("prime_limit=" + repr(prime_limit))
    print("inert_prime_count=" + repr(len(inert_primes)))
    print("direct_singleton_fibre_checks=" + repr(checked_unique))
    print("first_census=" + repr(tuple(census[:12])))
    print("last_census=" + repr(tuple(census[-12:])))
    print("constructed_value_count=" + repr(len(constructed)))
    print("semantic_sha256=" + stable_digest(semantic))
    print("scope=finite_exact_probe_only;all_scale_proof_is_in_scratch_report")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
