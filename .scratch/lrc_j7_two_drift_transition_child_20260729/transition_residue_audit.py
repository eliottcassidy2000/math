#!/usr/bin/env python3
"""Exact finite referee for the normalized two-drift handoff residue.

This only audits the arithmetic lemma and its hostile non-near-dilate ray.
It does not search for or assert an LRC cover.
"""

from __future__ import annotations

import math
from fractions import Fraction as F


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def positive_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue if residue else modulus


def overlap_numerator(a: int, b: int, n: int, m: int) -> int:
    return b * (14 * n + 1) - a * (14 * m - 1)


def predicted_numerator(a: int, b: int) -> int:
    g = math.gcd(a, b)
    return g * positive_residue(a // g + b // g, 14)


def bezout_witness_for_minimum(a: int, b: int) -> tuple[int, int, int]:
    """Construct tooth addresses attaining the predicted positive residue."""
    g = math.gcd(a, b)
    alpha = a // g
    beta = b // g
    target = predicted_numerator(a, b)
    quotient = (target - a - b) // (14 * g)
    require(
        target - a - b == 14 * g * quotient,
        "predicted residue is outside the endpoint lattice",
    )
    if alpha == 1:
        n = 0
    else:
        n = (quotient * pow(beta, -1, alpha)) % alpha
    m = (beta * n - quotient) // alpha
    require(beta * n - alpha * m == quotient, "Bezout witness failed")
    return target, n, m


def main() -> None:
    pairs = 0
    seam_pairs = 0
    for a in range(1, 180):
        for b in range(a + 1, 240):
            expected = predicted_numerator(a, b)
            g = math.gcd(a, b)
            actual, n, m = bezout_witness_for_minimum(a, b)
            require(
                overlap_numerator(a, b, n, m) == actual == expected,
                f"{a,b}: residue witness failed",
            )
            require(
                all(
                    (candidate - a - b) % (14 * g) != 0
                    for candidate in range(1, expected)
                ),
                f"{a,b}: smaller positive lattice numerator exists",
            )
            alpha = a // g
            beta = b // g
            if (alpha + beta) % 14 == 0:
                seam_pairs += 1
                require(
                    F(actual, 14 * a * b) == F(1, math.lcm(a, b)),
                    f"{a,b}: seam-compatible strict quantum failed",
                )
            pairs += 1

    for n in range(1, 10_001):
        a = 1
        b = 14 * n + 1
        numerator = overlap_numerator(a, b, 0, n)
        require(numerator == 2, f"N={n}: hostile ray numerator changed")
        require(
            F(numerator, 14 * a * b) == F(1, 7 * b),
            f"N={n}: hostile ray width changed",
        )

    print("LRC14 normalized two-drift transition residue scratch audit")
    print(f"finite_pairs={pairs};seam_compatible_pairs={seam_pairs}")
    print(
        "law=min_positive_numerator="
        "gcd(a,b)*least_positive_residue(a/g+b/g mod 14)"
    )
    print(
        "hostile_ray=(a,b)=(1,14N+1);addresses=(0,N);"
        "numerator=2;width=1/(7(14N+1));ratio_unbounded"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
