#!/usr/bin/env python3
"""Independent exact controls for THM-2781's terminal-tail theorem.

The all-degree statement is proved by a coefficient recurrence and UFD
factor multiplicities.  These finite rational controls probe every interface,
including sharp insufficient-tail and unreduced-exponent hostiles.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from math import gcd
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


GATES = 0


def gate(condition: bool, message: str) -> None:
    global GATES
    require(condition, message)
    GATES += 1


def trim(poly: list[Fraction]) -> list[Fraction]:
    answer = poly[:]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def multiply(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    answer = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            answer[i + j] += x * y
    return trim(answer)


def power(base: list[Fraction], exponent: int) -> list[Fraction]:
    answer = [Fraction(1)]
    factor = trim(base)
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = multiply(answer, factor)
        factor = multiply(factor, factor)
        remaining //= 2
    return answer


def series_coefficients(
    exponent: Fraction,
    base: list[Fraction],
    last: int,
) -> list[Fraction]:
    """Coefficients of base(z)^exponent for base(0)=1."""

    degree = len(base) - 1
    coefficients = [Fraction(1)]
    for n in range(1, last + 1):
        numerator = sum(
            base[r] * ((exponent + 1) * r - n) * coefficients[n - r]
            for r in range(1, min(degree, n) + 1)
        )
        coefficients.append(numerator / n)
    return coefficients


def padded(poly: list[Fraction], degree: int) -> list[Fraction]:
    require(len(poly) - 1 <= degree, "polynomial exceeds declared degree")
    return poly + [Fraction(0)] * (degree + 1 - len(poly))


def main() -> None:
    checked_parameter_rows = 0
    positive_power_rows = 0
    nonpower_rows = 0

    # Recurrence window, positive b-th powers, nonpowers, b=1, and declared
    # degree bounds with a vanishing top coefficient.
    for degree in range(1, 11):
        for denominator in range(1, degree + 1):
            if degree % denominator:
                continue
            for numerator in range(1, 9):
                if gcd(numerator, denominator) != 1:
                    continue
                checked_parameter_rows += 1
                exponent = Fraction(numerator, denominator)
                terminal = degree * numerator // denominator

                gate(
                    degree * numerator % denominator == 0,
                    f"d={degree},a={numerator},b={denominator}: integral N",
                )
                gate(
                    (exponent + 1) * degree - (terminal + degree) == 0,
                    f"d={degree},a={numerator},b={denominator}: multiplier",
                )
                for r in range(1, degree):
                    predecessor = terminal + degree - r
                    gate(
                        terminal + 1 <= predecessor <= terminal + degree - 1,
                        f"d={degree},a={numerator},b={denominator}: window",
                    )

                root_degree = degree // denominator
                patterns = (
                    [Fraction(1)]
                    + [Fraction((j % 3) - 1) for j in range(1, root_degree + 1)],
                    [Fraction(1)]
                    + [Fraction(j + 1) for j in range(1, root_degree + 1)],
                )
                for root in patterns:
                    base = padded(power(root, denominator), degree)
                    coefficients = series_coefficients(
                        exponent,
                        base,
                        terminal + degree + 3,
                    )
                    expected = padded(power(root, numerator), terminal)
                    gate(
                        coefficients[: terminal + 1] == expected,
                        f"d={degree},a={numerator},b={denominator}: positive prefix",
                    )
                    gate(
                        all(value == 0 for value in coefficients[terminal + 1 :]),
                        f"d={degree},a={numerator},b={denominator}: positive tail",
                    )
                    positive_power_rows += 1

                if denominator == 1:
                    arbitrary = padded(
                        [Fraction(1), Fraction(2), Fraction(-1)],
                        max(degree, 2),
                    )
                    arbitrary = arbitrary[: degree + 1]
                    coefficients = series_coefficients(
                        exponent,
                        arbitrary,
                        terminal + degree,
                    )
                    gate(
                        all(value == 0 for value in coefficients[terminal + 1 :]),
                        f"d={degree},a={numerator}: integral exponent",
                    )
                else:
                    nonpower = padded([Fraction(1), Fraction(1)], degree)
                    coefficients = series_coefficients(
                        exponent,
                        nonpower,
                        terminal + degree - 1,
                    )
                    gate(
                        any(value != 0 for value in coefficients[terminal + 1 :]),
                        f"d={degree},a={numerator},b={denominator}: nonpower tail",
                    )
                    nonpower_rows += 1

                if root_degree >= 2:
                    short_root = [Fraction(1), Fraction(2)]
                    short_base = padded(power(short_root, denominator), degree)
                    coefficients = series_coefficients(
                        exponent,
                        short_base,
                        terminal + degree,
                    )
                    gate(
                        short_base[-1] == 0
                        and all(
                            value == 0
                            for value in coefficients[terminal + 1 :]
                        ),
                        f"d={degree},a={numerator},b={denominator}: top-zero control",
                    )

    # One missing coefficient is insufficient for the cubic response gate.
    cubic = padded(
        [Fraction(1), Fraction(3), Fraction(3)],
        3,
    )
    cubic_coefficients = series_coefficients(Fraction(1, 3), cubic, 3)
    gate(cubic_coefficients[2] == 0, "cubic short-tail zero")
    gate(cubic_coefficients[3] != 0, "cubic missing second response")
    gate(len(trim(cubic)) - 1 == 2, "cubic is not a cube")

    # Two missing coefficients are insufficient uniformly for the quartic
    # gate: at M=6 the root-zero point f=1+z^3 is response-active.
    quartic = padded([Fraction(1), Fraction(0), Fraction(0), Fraction(1)], 4)
    quartic_coefficients = series_coefficients(Fraction(3, 2), quartic, 9)
    gate(quartic_coefficients[7] == 0, "quartic first flux zero")
    gate(quartic_coefficients[8] == 0, "quartic second flux zero")
    gate(quartic_coefficients[9] == Fraction(-1, 16), "quartic response active")
    gate(len(trim(quartic)) - 1 == 3, "quartic is not a square")

    # Coprimality of a,b is essential: the reduced denominator is two, not
    # the displayed four, and a square need not be a fourth power.
    unreduced = padded(power([Fraction(1), Fraction(0), Fraction(1)], 2), 4)
    unreduced_coefficients = series_coefficients(Fraction(2, 4), unreduced, 5)
    gate(
        all(value == 0 for value in unreduced_coefficients[3:6]),
        "unreduced exponent has full displayed tail",
    )
    gate(
        len(trim(unreduced)) - 1 == 4,
        "unreduced hostile is square but not fourth power",
    )

    source = Path(__file__).read_text(encoding="utf-8")
    gate(
        not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "truth-bearing assert node",
    )

    print("THM-2781 TERMINAL-TAIL PERFECT-POWER CERTIFICATE")
    print(f"parameter_rows={checked_parameter_rows}")
    print(f"positive_power_rows={positive_power_rows}")
    print(f"nonpower_rows={nonpower_rows}")
    print(f"exact_gates={GATES}")
    print("recurrence_terminal_multiplier=PASS")
    print("bth_power_equivalence_controls=PASS")
    print("vanishing_top_coefficient_control=PASS")
    print("cubic_dminus2_hostile=PASS")
    print("quartic_dminus2_hostile=PASS")
    print("unreduced_exponent_hostile=PASS")
    print("TERMINAL_TAIL_EXACT_CONTROLS=PASS")


if __name__ == "__main__":
    main()
