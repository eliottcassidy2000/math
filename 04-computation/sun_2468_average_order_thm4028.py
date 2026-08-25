#!/usr/bin/env python3
"""Exact finite controls and constants for THM-4028."""

from __future__ import annotations

import math
from fractions import Fraction


DEGREES = (2, 4, 6, 8)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def binomial_values(r: int, first_top: int, limit: int) -> list[int]:
    values = []
    top = first_top
    while math.comb(top, r) <= limit:
        values.append(math.comb(top, r))
        top += 1
    return values


def triangular_count(limit: int) -> int:
    """Number of w>=2 for which C(w,2)<=limit."""
    if limit < 1:
        return 0
    root = math.isqrt(8 * limit + 1)
    top = (root + 1) // 2
    while math.comb(top + 1, 2) <= limit:
        top += 1
    while math.comb(top, 2) > limit:
        top -= 1
    return top - 1


def cumulative_count(limit: int) -> int:
    c4 = binomial_values(4, 3, limit - 1)
    c6 = binomial_values(6, 5, limit - 1)
    c8 = binomial_values(8, 7, limit - 1)
    total = 0
    for t8 in c8:
        for t6 in c6:
            if t8 + t6 >= limit:
                break
            base = limit - t8 - t6
            for t4 in c4:
                if t4 >= base:
                    break
                total += triangular_count(base - t4)
    return total


def mixed_2469_hostile(target: int) -> tuple[int, int]:
    """Count minimal-domain 2-4-6-9 representations, allowing T=0."""
    c4 = binomial_values(4, 3, target)
    c6 = binomial_values(6, 5, target)
    c9 = binomial_values(9, 8, target)
    triples = 0
    representations = 0
    for t9 in c9:
        for t6 in c6:
            if t9 + t6 > target:
                break
            for t4 in c4:
                if t9 + t6 + t4 > target:
                    break
                triples += 1
                remainder = target - t9 - t6 - t4
                root = math.isqrt(8 * remainder + 1)
                if root * root == 8 * remainder + 1:
                    representations += 1
    return triples, representations


def constants() -> tuple[Fraction, float, float]:
    reciprocal_sum = sum(Fraction(1, r) for r in DEGREES)
    numerator = 1.0
    for r in DEGREES:
        numerator *= math.gamma(1 + 1 / r) * math.factorial(r) ** (1 / r)
    volume = numerator / math.gamma(1 + float(reciprocal_sum))
    shell = float(reciprocal_sum) * volume
    return reciprocal_sum, volume, shell


def main() -> None:
    reciprocal_sum, volume, shell = constants()
    require(reciprocal_sum == Fraction(25, 24), "reciprocal-sum regression")
    print(f"degrees={DEGREES}")
    print(f"reciprocal_sum={reciprocal_sum} mean_exponent={reciprocal_sum - 1}")
    print(f"cumulative_constant={volume:.15f}")
    print(f"formal_shell_constant={shell:.15f}")
    for limit in (10**4, 10**5, 10**6, 10**7, 10**8, 10**9):
        count = cumulative_count(limit)
        ratio = count / limit ** float(reciprocal_sum)
        mean_scaled = (count / limit) / limit ** float(reciprocal_sum - 1)
        require(abs(ratio - mean_scaled) < 1e-12, "mean normalization drift")
        print(
            f"X={limit} cumulative={count} ratio_to_X^(25/24)={ratio:.12f} "
            f"ratio_to_constant={ratio / volume:.12f}"
        )
    hostile_n = 1_061_619
    hostile_triples, hostile_count = mixed_2469_hostile(hostile_n)
    require(hostile_triples == 28_037, "2-4-6-9 universe regression")
    require(hostile_count == 0, "2-4-6-9 hostile unexpectedly represented")
    print(
        f"hostile_2_4_6_9 N={hostile_n} reciprocal_sum={Fraction(37, 36)} "
        f"triples={hostile_triples} representations={hostile_count}"
    )
    print("finite_controls=PASS exact_integer_counts")


if __name__ == "__main__":
    main()
