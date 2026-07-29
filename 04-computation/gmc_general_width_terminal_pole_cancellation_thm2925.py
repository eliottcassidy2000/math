#!/usr/bin/env python3
"""Exact companion for THM-2925.

The proof of the general denominator law is symbolic and valid at every
width.  This companion supplies hostile exact controls: direct Laurent
principal parts at the terminal pole, exhaustive coefficient divisibility
on every four-slot shape through width eight and orders two through five,
and numerical specialization checks against the normalized moment
constructor.
"""

from __future__ import annotations

import importlib.util
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import factorial
from pathlib import Path

from flint import fmpz_poly


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


SOURCE = Path(__file__).with_name(
    "gmc_diameter_four_nonconsecutive_macaulay_newton_thm2921.py"
)
SOURCE_BYTES = SOURCE.read_bytes().replace(b"\r\n", b"\n")
require(
    sha256(SOURCE_BYTES).hexdigest()
    == "42e9b5ceddd677d1f2601a9d5d668c9437281596b65999ddcb8549d4e0b9bf64",
    "THM-2921 constructor dependency hash changed",
)
SPEC = importlib.util.spec_from_file_location("thm2921_exact", SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load THM-2921")
t = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(t)


def general_denominator_poly(width: int, order: int) -> fmpz_poly:
    answer = fmpz_poly(1)
    for shift in range(1, width):
        answer *= (t.POLY_X + shift) ** (order - 1)
    answer *= (t.POLY_X + width) ** (order - 2)
    return answer


def general_denominator(depth: int, width: int, order: int) -> int:
    answer = 1
    for shift in range(1, width):
        answer *= (depth + shift) ** (order - 1)
    answer *= (depth + width) ** (order - 2)
    return answer


def signed_coefficient_fraction(
    width: int,
    directions: tuple[int, ...],
) -> tuple[fmpz_poly, fmpz_poly]:
    answer = (fmpz_poly(0), fmpz_poly(1))
    for mask in range(1 << len(directions)):
        selected = tuple(
            width if mask & (1 << position) else direction
            for position, direction in enumerate(directions)
        )
        numerator, denominator = t.normalized_tensor_poly(selected)
        sign = -1 if mask.bit_count() % 2 else 1
        answer = t.add_fraction(
            answer,
            (sign * numerator, denominator),
        )
    return answer


def scaled_coefficient(
    width: int,
    monomial: tuple[int, int, int],
    interior_offsets: tuple[int, int, int],
) -> fmpz_poly:
    order = sum(monomial)
    directions = tuple(
        interior_offsets[index]
        for index, count in enumerate(monomial)
        for _ in range(count)
    )
    numerator, denominator = signed_coefficient_fraction(
        width,
        directions,
    )
    numerator, denominator = t.reduced_fraction(
        t.multinomial(monomial) * numerator,
        denominator,
    )
    common = general_denominator_poly(width, order)
    require(
        common % denominator == 0,
        f"denominator escape: M={width}, alpha={monomial}, b={interior_offsets}",
    )
    answer = (common // denominator) * numerator
    require(
        answer.degree() <= (order - 1) * width - 1,
        f"degree escape: M={width}, alpha={monomial}, b={interior_offsets}",
    )
    return answer


def terminal_principal_part(
    offsets: tuple[int, ...],
    width: int,
) -> tuple[int, Fraction]:
    """Return pole order and leading Laurent coefficient at n=-M."""
    numerator, denominator = t.normalized_tensor_poly(offsets)
    factor = t.POLY_X + width
    numerator_order = 0
    denominator_order = 0
    while numerator(-width) == 0:
        numerator //= factor
        numerator_order += 1
    while denominator(-width) == 0:
        denominator //= factor
        denominator_order += 1
    return (
        denominator_order - numerator_order,
        Fraction(int(numerator(-width)), int(denominator(-width))),
    )


def polynomial_digest(records: list[str]) -> str:
    return sha256(("|".join(records) + "\n").encode()).hexdigest()


def main() -> None:
    terminal_records: list[str] = []
    terminal_controls = 0
    for width in range(1, 13):
        for order in range(2, 8):
            constant = Fraction(
                factorial(order * width - 1),
                factorial(width - 1) ** order,
            )
            pole, leading = terminal_principal_part(
                (width,) * order,
                width,
            )
            all_top_contribution = (-1) ** order * leading
            require(
                pole == order - 1
                and all_top_contribution == -order * constant,
                f"all-top principal part changed: M={width}, m={order}",
            )
            total = all_top_contribution
            for lower_offset in range(width):
                pole, leading = terminal_principal_part(
                    (lower_offset,) + (width,) * (order - 1),
                    width,
                )
                contribution = (-1) ** (order - 1) * leading
                require(
                    pole == order - 1 and contribution == constant,
                    "one-nontop principal part changed: "
                    f"M={width}, m={order}, b={lower_offset}",
                )
                terminal_records.append(
                    f"{width}:{order}:{lower_offset}:{pole}:"
                    f"{contribution.numerator}/{contribution.denominator}"
                )
                terminal_controls += 1

            # The signed tensor has one such contribution for each of its
            # order positions.  The leading value is independent of the
            # retained lower offset.
            total += order * constant
            require(
                total == 0,
                f"terminal residue did not cancel: M={width}, m={order}",
            )

    coefficient_records: list[str] = []
    coefficient_checks = 0
    numeric_checks = 0
    support_checks = 0
    for width in range(3, 9):
        supports = tuple(combinations(range(1, width), 2))
        for first, second in supports:
            support_checks += 1
            interior_offsets = (0, first, second)
            forms: dict[
                int,
                dict[tuple[int, int, int], fmpz_poly],
            ] = {}
            for order in range(2, 6):
                form = {}
                for monomial in t.MONOMIALS[order]:
                    polynomial = scaled_coefficient(
                        width,
                        monomial,
                        interior_offsets,
                    )
                    form[monomial] = polynomial
                    coefficient_records.append(
                        f"{width}:{first}:{second}:{order}:{monomial}:"
                        + ",".join(map(str, polynomial.coeffs()))
                    )
                    coefficient_checks += 1
                forms[order] = form

            # One support per width gets a full specialization audit at
            # three hostile depths.
            if (first, second) == supports[-1]:
                offsets = (0, first, second, width)
                for depth in (0, 1, 3):
                    for order in range(2, 6):
                        exact = t.moment_form(depth, order, offsets)
                        scale = general_denominator(depth, width, order)
                        for monomial, polynomial in forms[order].items():
                            expected = exact[monomial] * scale
                            require(
                                expected.denominator == 1
                                and int(polynomial(depth))
                                == expected.numerator,
                                "specialization mismatch: "
                                f"M={width}, n={depth}, m={order}",
                            )
                            numeric_checks += 1

    require(terminal_controls == 468, "terminal control count changed")
    require(support_checks == 56, "support control count changed")
    require(coefficient_checks == 2912, "coefficient control count changed")
    require(numeric_checks == 936, "numeric control count changed")

    print("THM-2925 GENERAL-WIDTH TERMINAL-POLE CANCELLATION")
    print(
        "terminal_grid=M:1..12,m:2..7,b:0..M-1;"
        f"one_nontop_controls={terminal_controls};PASS"
    )
    print(
        "terminal_identity=all_top:-m*C;"
        "each_one_nontop:+C;sum_over_m_positions:0"
    )
    print(
        "divisibility_grid=M:3..8,m:2..5;"
        f"supports={support_checks};coefficients={coefficient_checks};PASS"
    )
    print(f"numeric_specialization_controls={numeric_checks};PASS")
    print(
        "terminal_digest_sha256="
        + polynomial_digest(terminal_records)
    )
    print(
        "scaled_coefficient_digest_sha256="
        + polynomial_digest(coefficient_records)
    )
    print(
        "proved_denominator="
        "prod_(r=1)^(M-1)(n+r)^(m-1)*(n+M)^(m-2)"
    )
    print("proved_scaled_degree_bound=(m-1)*M-1")
    print("fixed_20Q_10C_6F_minor_degree_bound=58*M-36")
    print(
        "scope=denominator and degree law only;"
        "no uniform minor nonvanishing or arbitrary-width SFC4"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
