#!/usr/bin/env python3
"""Exact referee for THM-2194; the symbolic proof supplies the quantifiers.

The first path derives the five Faber boundary/flux rows over Q[c,L,O].
The second path is an independent exact-Fraction Laurent solver on a bounded
hostile bank.  The bounded bank is a regression test, not the proof.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction

import sympy as sp


def symbolic_faber_bank() -> None:
    x, c, lam, omega, mu = sp.symbols("X c Lambda Omega mu")
    k = 1 + x + c * x**2 / 2
    t_hat = sp.expand(k**2 + lam * x**3 + omega * x**4)
    t_coeff = [None] + [t_hat.coeff(x, i) for i in range(1, 5)]

    expected = {
        1: (
            sp.Rational(1, 2),
            (2 * c - 1) / 8,
            lam / 4,
        ),
        2: (
            c / 2,
            lam / 2,
            -(lam - 2 * omega) / 4,
        ),
        3: (
            (12 * lam + 6 * c - 1) / 16,
            3 * (-16 * lam + 32 * omega + 4 * c**2 - 4 * c + 1) / 128,
            -3 * lam * (2 * c - 1) / 32,
        ),
        5: (
            (
                80 * lam * c
                - 40 * lam
                + 160 * omega
                + 60 * c**2
                - 20 * c
                + 3
            )
            / 256,
            5
            * (
                32 * lam**2
                - 32 * lam * c
                + 16 * lam
                + 64 * omega * c
                - 32 * omega
                + 8 * c**3
                - 12 * c**2
                + 6 * c
                - 1
            )
            / 1024,
            -5
            * lam
            * (16 * lam - 32 * omega + 4 * c**2 - 4 * c + 1)
            / 512,
        ),
        6: (
            (3 * lam**2 + 6 * omega * c + c**3) / 8,
            -3 * lam * (lam - 2 * omega) / 8,
            3 * ((lam - 2 * omega) ** 2 + lam**2 * (1 - 2 * c)) / 32,
        ),
    }

    for degree, wanted in expected.items():
        alpha = sp.Rational(degree, 4)
        coeff = [sp.Integer(1)]
        for index in range(1, degree + 3):
            value = sum(
                t_coeff[i]
                * ((alpha + 1) * i - index)
                * coeff[index - i]
                for i in range(1, min(4, index) + 1)
            ) / index
            coeff.append(sp.expand(value))
        found = (
            coeff[degree],
            coeff[degree + 1],
            coeff[degree + 2] + coeff[degree + 1] / 2,
        )
        assert all(sp.expand(a - b) == 0 for a, b in zip(found, wanted))

    b6 = expected[6][0]
    cusp = sp.expand(8 * b6.subs(omega, mu - c**2 / 4))
    assert sp.expand(cusp - (3 * lam**2 - c**3 / 2 + 6 * mu * c)) == 0
    print("symbolic_faber_bank=PASS degrees=1,2,3,5,6 observables=15")
    print("cubic_cusp_identity=PASS")


Poly = dict[int, Fraction]


def poly_add(*polys: Poly) -> Poly:
    out: defaultdict[int, Fraction] = defaultdict(Fraction)
    for poly in polys:
        for exponent, coefficient in poly.items():
            out[exponent] += coefficient
    return {e: a for e, a in out.items() if a}


def poly_scale(poly: Poly, scalar: Fraction, shift: int = 0) -> Poly:
    return {
        exponent + shift: coefficient * scalar
        for exponent, coefficient in poly.items()
        if coefficient * scalar
    }


def poly_mul(left: Poly, right: Poly) -> Poly:
    out: defaultdict[int, Fraction] = defaultdict(Fraction)
    for e, a in left.items():
        for f, b in right.items():
            out[e + f] += a * b
    return {e: a for e, a in out.items() if a}


ONE: Poly = {0: Fraction(1)}


def faber_coefficients(degree: int, t_coeff: list[Poly | None]) -> list[Poly]:
    alpha = Fraction(degree, 4)
    coeff = [ONE]
    for index in range(1, degree + 3):
        value: Poly = {}
        for i in range(1, min(4, index) + 1):
            term = poly_mul(t_coeff[i] or {}, coeff[index - i])
            value = poly_add(
                value,
                poly_scale(term, (alpha + 1) * i - index),
            )
        coeff.append(poly_scale(value, Fraction(1, index)))
    return coeff


def laurent_system_is_consistent(
    degree: int,
    pole_depth: int,
    t_coeff: list[Poly | None],
) -> bool:
    degrees = [j for j in range(1, degree + 1) if j % 4]
    variable_count = len(degrees) - 1
    bank = {j: faber_coefficients(j, t_coeff) for j in degrees}
    rows: defaultdict[tuple[str, int], list[Fraction]] = defaultdict(
        lambda: [Fraction(0)] * (variable_count + 1)
    )

    for index, current in enumerate(degrees):
        column = index if current < degree else variable_count
        observables = (
            ("B", current, bank[current][current]),
            ("F", current + 1, bank[current][current + 1]),
            (
                "G",
                current + 2,
                poly_add(
                    bank[current][current + 2],
                    poly_scale(bank[current][current + 1], Fraction(1, 2)),
                ),
            ),
        )
        for label, power, polynomial in observables:
            for exponent, coefficient in polynomial.items():
                laurent_exponent = exponent - power * pole_depth
                if laurent_exponent < 0:
                    rows[(label, laurent_exponent)][column] += (
                        coefficient * Fraction(1, 2**power)
                    )

    matrix = [row for row in rows.values() if any(row)]
    rank = 0
    for column in range(variable_count):
        pivot = next(
            (i for i in range(rank, len(matrix)) if matrix[i][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [entry / pivot_value for entry in matrix[rank]]
        for i, row in enumerate(matrix):
            if i == rank or not row[column]:
                continue
            multiple = row[column]
            matrix[i] = [
                a - multiple * b for a, b in zip(row, matrix[rank])
            ]
        rank += 1

    return not any(
        not any(row[:-1]) and row[-1]
        for row in matrix
    )


def exact_hostile_bank() -> None:
    coefficients = (
        Fraction(-2),
        Fraction(-1),
        Fraction(-1, 2),
        Fraction(0),
        Fraction(1, 2),
        Fraction(1),
        Fraction(2),
    )
    checked = 0
    survivors = 0

    for pole_depth in range(1, 7):
        for c_order in range(1, 2 * pole_depth):
            for cubic_order in range(1, 4 * pole_depth + 1):
                for cubic_coefficient in coefficients:
                    t_coeff: list[Poly | None] = [
                        None,
                        {0: Fraction(2)},
                        poly_add(ONE, {c_order: Fraction(1)}),
                        (
                            {cubic_order: cubic_coefficient}
                            if cubic_coefficient
                            else {}
                        ),
                        {},
                    ]
                    checked += 1
                    survivors += laurent_system_is_consistent(
                        6, pole_depth, t_coeff
                    )

    assert checked == 4508
    assert survivors == 0

    hostile = [
        None,
        {0: Fraction(2)},
        poly_add(ONE, {3: Fraction(1)}),
        {3: Fraction(1)},
        {},
    ]
    assert not laurent_system_is_consistent(6, 2, hostile)

    square_control = [
        None,
        {0: Fraction(2)},
        ONE,
        {},
        {},
    ]
    assert laurent_system_is_consistent(6, 1, square_control)

    print(f"finite_hostile_bank=PASS cases={checked} survivors={survivors}")
    print("hostile_split_degree6_control=PASS")
    print("exact_square_positive_control=PASS")
    print("scope=FINITE_REFEREE_ONLY symbolic_valuation_proof_supplies_quantifiers")


if __name__ == "__main__":
    symbolic_faber_bank()
    exact_hostile_bank()
