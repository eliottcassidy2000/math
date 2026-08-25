#!/usr/bin/env python3
"""Exact hostile for THM-4036's pointwise polylog-loss frontier."""

from math import comb


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(label)


def first_odd_primes(k: int) -> tuple[int, ...]:
    primes = []
    candidate = 3
    while len(primes) < k:
        if all(candidate % divisor for divisor in range(2, int(candidate**0.5) + 1)):
            primes.append(candidate)
        candidate += 2
    return tuple(primes)


def divisors_of_squarefree(primes: tuple[int, ...]) -> tuple[int, ...]:
    divisors = [1]
    for prime in primes:
        divisors += [prime * divisor for divisor in divisors]
    return tuple(divisors)


def triangular(index: int) -> int:
    return comb(index, 2)


def main() -> None:
    for k in range(1, 11):
        primes = first_odd_primes(k)
        product = 1
        for prime in primes:
            product *= prime
        quartic_index = 24 * product
        difference = comb(quartic_index, 4)
        pairs = set()
        for divisor in divisors_of_squarefree(primes):
            cofactor = 2 * difference // divisor
            require(divisor * cofactor == 2 * difference, (k, divisor, "factor"))
            require(divisor % 2 == 1 and cofactor % 2 == 0, (k, divisor, "parity"))
            require(cofactor - divisor >= 3, (k, divisor, "positive indices"))
            upper = (divisor + cofactor + 1) // 2
            lower = (cofactor - divisor + 1) // 2
            require(upper >= 2 and lower >= 2, (k, divisor, "domain"))
            require(
                triangular(upper) - triangular(lower) == difference,
                (k, divisor, "triangular difference"),
            )
            pairs.add((upper, lower))
        require(len(pairs) == 2**k, (k, "pair count"))
        require(
            difference == comb(quartic_index, 4) - comb(3, 4),
            (k, "structured tail difference"),
        )
        print(
            f"k={k};P={product};t={quartic_index};h={difference};"
            f"constructed_R2={len(pairs)};"
            f"max_target={max(triangular(upper) for upper, _lower in pairs)}"
        )
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
