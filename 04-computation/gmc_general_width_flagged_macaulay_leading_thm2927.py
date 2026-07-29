#!/usr/bin/env python3
"""Exact companion for THM-2927.

This audits the universal pure-power identity for the fixed 36-row Macaulay
chart, its flag-boundary ranks, and the Cauchy--Binet/strict-total-positivity
specialization to every four-slot width through twelve.
"""

from __future__ import annotations

import importlib.util
from hashlib import sha256
from itertools import combinations
from math import factorial
from pathlib import Path

import sympy as sp
from flint import fmpz_mat


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


def multinomial(exponents: tuple[int, int, int]) -> int:
    answer = factorial(sum(exponents))
    for exponent in exponents:
        answer //= factorial(exponent)
    return answer


def pure_form(
    vector: tuple[object, object, object],
    degree: int,
) -> dict[tuple[int, int, int], object]:
    return {
        exponents: multinomial(exponents)
        * vector[0] ** exponents[0]
        * vector[1] ** exponents[1]
        * vector[2] ** exponents[2]
        for exponents in t.MONOMIALS[degree]
    }


def selected_matrix(
    u: tuple[object, object, object],
    v: tuple[object, object, object],
    w: tuple[object, object, object],
) -> list[list[object]]:
    rows: list[list[object]] = []
    for degree, form in (
        (2, pure_form(u, 2)),
        (3, pure_form(v, 3)),
        (4, pure_form(w, 4)),
    ):
        for multiplier in t.MONOMIALS[t.TARGET_DEGREE - degree]:
            row: list[object] = [0] * len(t.TARGET_MONOMIALS)
            for monomial, coefficient in form.items():
                target = tuple(
                    multiplier[index] + monomial[index]
                    for index in range(3)
                )
                row[t.TARGET_INDEX[target]] = coefficient
            rows.append(row)
    return [rows[index] for index in t.SELECTED_ROWS]


def det2(u: tuple[int, int, int], v: tuple[int, int, int]) -> int:
    return u[0] * v[1] - u[1] * v[0]


def det3(
    u: tuple[int, int, int],
    v: tuple[int, int, int],
    w: tuple[int, int, int],
) -> int:
    return int(fmpz_mat([list(u), list(v), list(w)]).det())


def flagged_formula(
    u: tuple[int, int, int],
    v: tuple[int, int, int],
    w: tuple[int, int, int],
) -> int:
    return (
        3
        * u[0] ** 14
        * v[0] ** 4
        * det2(u, v) ** 2
        * det3(u, v, w) ** 24
    )


def selected_determinant(
    u: tuple[int, int, int],
    v: tuple[int, int, int],
    w: tuple[int, int, int],
) -> int:
    return int(fmpz_mat(selected_matrix(u, v, w)).det())


def response_row(base: int, offset: int, width: int) -> int:
    return base**width - base**offset


def response_vector(
    base: int,
    first: int,
    second: int,
    width: int,
) -> tuple[int, int, int]:
    return tuple(
        response_row(base, offset, width)
        for offset in (0, first, second)
    )


def generalized_vandermonde(
    bases: tuple[int, ...],
    exponents: tuple[int, ...],
) -> int:
    return int(
        fmpz_mat(
            [
                [(base - 1) * base**exponent for exponent in exponents]
                for base in bases
            ]
        ).det()
    )


def cauchy_binet_2(first: int, width: int) -> tuple[int, int]:
    answer = 0
    terms = 0
    for low in range(first):
        for high in range(first, width):
            summand = generalized_vandermonde((2, 3), (low, high))
            require(summand > 0, "TP2 summand lost positivity")
            answer += summand
            terms += 1
    return answer, terms


def cauchy_binet_3(
    first: int,
    second: int,
    width: int,
) -> tuple[int, int]:
    answer = 0
    terms = 0
    for low in range(first):
        for middle in range(first, second):
            for high in range(second, width):
                summand = generalized_vandermonde(
                    (2, 3, 4),
                    (low, middle, high),
                )
                require(summand > 0, "TP3 summand lost positivity")
                answer += summand
                terms += 1
    return answer, terms


def main() -> None:
    print("THM-2927 GENERAL-WIDTH FLAGGED MACAULAY LEADING COEFFICIENT")
    print(
        "selected_rows="
        + ",".join(map(str, t.SELECTED_ROWS))
        + ";allocation=(20Q,10C,6F)"
    )
    print(
        "identity=3*u0^14*v0^4*det2(u,v)^2*det3(u,v,w)^24"
    )

    b, d = sp.symbols("b d")
    symbolic_matrix = sp.Matrix(
        selected_matrix((1, b, 0), (d, 1, 0), (0, 0, 1))
    )
    symbolic_determinant = sp.factor(
        symbolic_matrix.det(method="domain-ge")
    )
    require(
        sp.expand(
            symbolic_determinant - 3 * d**4 * (b * d - 1) ** 26
        )
        == 0,
        "two-parameter flag slice changed",
    )
    print("symbolic_flag_slice=3*d^4*(b*d-1)^26;PASS")

    generic_controls = 0
    generic_values = []
    for seed in range(1, 41):
        u = (seed + 1, 2 * seed - 3, 3 - seed)
        v = (2 * seed + 1, seed + 2, 5 - 3 * seed)
        w = (3 * seed - 1, 7 - seed, seed + 4)
        determinant = selected_determinant(u, v, w)
        require(
            determinant == flagged_formula(u, v, w),
            f"generic flagged identity failed at seed {seed}",
        )
        generic_values.append(determinant)
        generic_controls += 1
    print(
        f"generic_integer_identity_controls={generic_controls};"
        f"digest={sha256(','.join(map(str, generic_values)).encode()).hexdigest()}"
    )

    boundary_data = (
        ("u0", (0, 1, 0), (1, 0, 0), (0, 0, 1), 29),
        ("v0", (1, 0, 0), (0, 1, 0), (0, 0, 1), 34),
        ("det2", (1, 0, 0), (1, 0, 1), (0, 1, 0), 35),
        ("generic", (1, 2, 3), (2, 1, 4), (3, 5, 1), 36),
    )
    boundary_summary = []
    for name, u, v, w, expected_rank in boundary_data:
        rank = fmpz_mat(selected_matrix(u, v, w)).rank()
        require(rank == expected_rank, f"{name} boundary rank changed")
        boundary_summary.append(f"{name}:{rank}")
    print("flag_boundary_ranks=" + ",".join(boundary_summary))

    support_count = 0
    cb2_terms = 0
    cb3_terms = 0
    leading_values = []
    for width in range(3, 13):
        for first, second in combinations(range(1, width), 2):
            r2 = response_vector(2, first, second, width)
            r3 = response_vector(3, first, second, width)
            r4 = response_vector(4, first, second, width)

            minor2, terms2 = cauchy_binet_2(first, width)
            minor3, terms3 = cauchy_binet_3(first, second, width)
            require(minor2 == det2(r2, r3) > 0, "CB2 identity failed")
            require(minor3 == det3(r2, r3, r4) > 0, "CB3 identity failed")

            determinant = selected_determinant(r2, r3, r4)
            formula = flagged_formula(r2, r3, r4)
            require(
                determinant == formula > 0,
                f"width/support leading identity failed: "
                f"{width},{first},{second}",
            )
            support_count += 1
            cb2_terms += terms2
            cb3_terms += terms3
            leading_values.append(formula)

    require(support_count == 220, "width/support census changed")
    print(
        f"widths=3..12;supports={support_count};"
        f"positive_CB2_terms={cb2_terms};positive_CB3_terms={cb3_terms}"
    )
    print(
        "leading_value_digest="
        + sha256(",".join(map(str, leading_values)).encode()).hexdigest()
    )
    print(
        "consequence=degree 58M-36 is exact and the fixed chart is "
        "eventually positive for every 0<a<b<M"
    )
    print(
        "scope=leading coefficient only;finite-depth nonvanishing and "
        "all-width SFC(4) remain open"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
