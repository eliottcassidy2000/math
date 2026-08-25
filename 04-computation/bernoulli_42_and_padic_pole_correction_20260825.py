#!/usr/bin/env python3
"""Exact finite controls for the bernoulli-hurwitz correction.

The von Staudt--Clausen denominator formula is evaluated for every even index
2<=n<=1000.  The p-adic interpolation value zeta_7(-5) is also reduced exactly.
No quaternion ramification claim is inferred from these denominator data.
"""

from fractions import Fraction
from math import isqrt


LIMIT = 1000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primes_through(n: int) -> list[int]:
    sieve = [True] * (n + 1)
    sieve[0:2] = [False, False]
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            for multiple in range(p * p, n + 1, p):
                sieve[multiple] = False
    return [p for p, is_prime in enumerate(sieve) if is_prime]


def von_staudt_denominator(index: int, primes: list[int]) -> int:
    value = 1
    for p in primes:
        if p - 1 > index:
            break
        if index % (p - 1) == 0:
            value *= p
    return value


def valuation(numerator: int, prime: int) -> int:
    value = 0
    while numerator % prime == 0 and numerator:
        numerator //= prime
        value += 1
    return value


def main() -> None:
    primes = primes_through(LIMIT + 1)
    denominators = {
        index: von_staudt_denominator(index, primes)
        for index in range(2, LIMIT + 1, 2)
    }
    hits = [index for index, denominator in denominators.items() if denominator == 42]

    require(denominators[2] == 6, "B_2 denominator control failed")
    require(denominators[4] == 30, "B_4 denominator control failed")
    require(hits[:2] == [6, 114], "the first two denominator-42 indices changed")
    require(len(hits) == 15, "unexpected denominator-42 census count")

    zeta_7_minus_5 = Fraction(7**5 - 1, 252)
    v7 = valuation(zeta_7_minus_5.numerator, 7) - valuation(zeta_7_minus_5.denominator, 7)
    require(zeta_7_minus_5 == Fraction(2801, 42), "zeta_7(-5) interpolation value failed")
    require(v7 == -1, "zeta_7(-5) valuation failed")

    print("BERNOULLI-42 AND P-ADIC POLE CORRECTION")
    print(f"universe: even Bernoulli indices 2<=n<={LIMIT}")
    print(f"denom(B_n)=42 indices = {hits}")
    print(f"count                    = {len(hits)}")
    print("first two                = [6, 114]")
    print(f"zeta_7(-5)               = {zeta_7_minus_5}")
    print(f"v_7(zeta_7(-5))          = {v7}")
    print("negative valuation is a finite denominator, not a pole")
    print("FINITE-EXACT CORRECTION: PASS")


if __name__ == "__main__":
    main()
