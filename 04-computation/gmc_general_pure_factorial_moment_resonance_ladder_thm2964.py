#!/usr/bin/env python3
"""Exact regression for THM-2964's pure factorial-moment wall ladder.

For order k and first-gap width M, put

    S_k(n) = sum_j (-1)^j binom(k,j)
             (k n + 1)_{j M} / (n + 1)_M^j.

The denominator-cleared pure coefficient is

    P_{k,M}(n) = D_{k,M}(n) S_k(n),

where D has exponent k-1 at -1,...,-(M-1) and exponent k-2 at -M.
This companion constructs P without importing any earlier GMC script and
checks the all-order theorem on a finite hostile bank.
"""

from __future__ import annotations

from hashlib import sha256
from math import comb

from flint import fmpz_poly


X = fmpz_poly([0, 1])


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def rising(slope: int, intercept: int, length: int) -> fmpz_poly:
    answer = fmpz_poly(1)
    for offset in range(length):
        answer *= slope * X + intercept + offset
    return answer


def pure_coefficient(order: int, width: int) -> fmpz_poly:
    """Return D_{k,M} S_k by one exact polynomial division."""

    base = rising(1, 1, width)
    numerator = fmpz_poly(0)
    for selected in range(order + 1):
        numerator += (
            (-1) ** selected
            * comb(order, selected)
            * rising(order, 1, selected * width)
            * base ** (order - selected)
        )

    # base^order / D_{k,M} = base * (X+M).
    clearing_quotient = base * (X + width)
    require(
        numerator % clearing_quotient == 0,
        f"denominator clearing failed: k={order}, M={width}",
    )
    return numerator // clearing_quotient


def valuation_at(polynomial: fmpz_poly, root: int) -> int:
    factor = X + root
    value = polynomial
    answer = 0
    while value(-root) == 0:
        value //= factor
        answer += 1
    return answer


def predicted(order: int, width: int, root: int) -> bool:
    d = width - root
    return (
        d >= 1
        and width == (order + 1) * d + 1
        and root == order * d + 1
        and ((order - 1) * d) % 2 == 0
    )


def main() -> None:
    records: list[str] = []
    total_walls = 0
    total_roots = 0

    for order in range(2, 8):
        roots: list[tuple[int, int]] = []
        for width in range(3, 36):
            coefficient = pure_coefficient(order, width)
            for root in range(1, width):
                total_walls += 1
                valuation = valuation_at(coefficient, root)
                expectation = predicted(order, width, root)
                require(
                    (valuation > 0) == expectation,
                    "resonance classification failed: "
                    f"k={order},M={width},r={root},v={valuation}",
                )
                require(
                    valuation in (0, 1),
                    "a predicted pure wall was not simple: "
                    f"k={order},M={width},r={root},v={valuation}",
                )
                if valuation:
                    roots.append((width, root))
                    total_roots += 1
        record = f"k={order};roots={roots}"
        records.append(record)
        print(record)

    digest = sha256("\n".join(records).encode()).hexdigest()
    print(f"orders=2..7;widths=3..35;walls={total_walls};roots={total_roots}")
    print(f"record_digest={digest}")
    print("general_resonance_regression=PASS")


if __name__ == "__main__":
    main()
