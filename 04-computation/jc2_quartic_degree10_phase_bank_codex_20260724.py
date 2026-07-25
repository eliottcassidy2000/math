#!/usr/bin/env python3
"""Exact referee for THM-2202; the valuation proof supplies the quantifiers.

The symbolic path derives the degree-ten boundary/two-flux bank over
Q[c,Lambda,Omega].  The independent exact-Fraction Laurent solver checks a
bounded hostile universe.  That finite universe is a regression test only.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction

import sympy as sp


def require(condition: bool, message: str) -> None:
    """Keep every referee check active under both normal Python and -O."""
    if not condition:
        raise RuntimeError(message)


def symbolic_faber_bank() -> None:
    x, c, lam, omega, mu = sp.symbols("X c Lambda Omega mu")
    k = 1 + x + c * x**2 / 2
    t_hat = sp.expand(k**2 + lam * x**3 + omega * x**4)
    t_coeff = [None] + [t_hat.coeff(x, i) for i in range(1, 5)]

    expected_top = (
        (
            -10 * lam**3
            + 30 * lam**2 * omega
            + 30 * omega**2 * c
            + 10 * omega * c**3
            + c**5
        )
        / 32,
        -5
        * lam
        * (
            lam**2 * c
            - 2 * lam**2
            + 6 * lam * omega
            - 6 * omega**2
        )
        / 32,
        5
        * (
            -lam**4
            + 6 * lam**3 * c
            - 4 * lam**3
            - 12 * lam**2 * omega * c
            + 12 * lam**2 * omega
            - 12 * lam * omega**2
            + 8 * omega**3
        )
        / 128,
    )

    bank: dict[int, tuple[sp.Expr, sp.Expr, sp.Expr]] = {}
    for degree in (1, 2, 3, 5, 6, 7, 9, 10):
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
        bank[degree] = (
            coeff[degree],
            coeff[degree + 1],
            coeff[degree + 2] + coeff[degree + 1] / 2,
        )

    require(
        all(
            sp.expand(a - b) == 0
            for a, b in zip(bank[10], expected_top)
        ),
        "degree-ten top Faber row mismatch",
    )

    expected_even = {
        2: (
            c / 2,
            lam / 2,
            -(lam - 2 * omega) / 4,
        ),
        6: (
            (3 * lam**2 + 6 * omega * c + c**3) / 8,
            -3 * lam * (lam - 2 * omega) / 8,
            3
            * ((lam - 2 * omega) ** 2 + lam**2 * (1 - 2 * c))
            / 32,
        ),
    }
    for degree, wanted in expected_even.items():
        require(
            all(
                sp.expand(a - b) == 0
                for a, b in zip(bank[degree], wanted)
            ),
            f"degree-{degree} even Faber row mismatch",
        )

    odd_residues = {
        1: (
            sp.Rational(1, 2),
            sp.Rational(-1, 8),
            sp.Rational(1, 4),
        ),
        3: (
            sp.Rational(-1, 16),
            sp.Rational(3, 128),
            sp.Rational(3, 32),
        ),
        5: (
            sp.Rational(3, 256),
            sp.Rational(-5, 1024),
            sp.Rational(-5, 512),
        ),
        7: (
            sp.Rational(-5, 2048),
            sp.Rational(35, 32768),
            sp.Rational(7, 4096),
        ),
        9: (
            sp.Rational(35, 65536),
            sp.Rational(-63, 262144),
            sp.Rational(-45, 131072),
        ),
    }
    zero_substitution = {c: 0, lam: 0, omega: 0}
    for degree, wanted in odd_residues.items():
        b, f, g = bank[degree]
        found = (
            b.subs(zero_substitution),
            f.subs(zero_substitution),
            sp.diff(g, lam).subs(zero_substitution),
        )
        require(found == wanted, f"degree-{degree} odd residue mismatch")
        require(
            sp.expand(g.subs(lam, 0)) == 0,
            f"degree-{degree} second flux is not Lambda-divisible",
        )

    b10_mu = sp.expand(32 * bank[10][0].subs(omega, mu - c**2 / 4))
    wanted_b10_mu = (
        -10 * lam**3
        - sp.Rational(15, 2) * lam**2 * c**2
        + 30 * lam**2 * mu
        + sp.Rational(3, 8) * c**5
        - 5 * c**3 * mu
        + 30 * c * mu**2
    )
    require(
        sp.expand(b10_mu - wanted_b10_mu) == 0,
        "degree-ten mu-substitution identity mismatch",
    )

    b6_mu = sp.expand(8 * bank[6][0].subs(omega, mu - c**2 / 4))
    require(
        sp.expand(
            b6_mu - (3 * lam**2 - c**3 / 2 + 6 * mu * c)
        )
        == 0,
        "degree-six mu-substitution identity mismatch",
    )

    print("symbolic_faber_bank=PASS degrees=1,2,3,5,6,7,9,10")
    print("degree10_top_row=PASS observables=3")
    print("odd_residue_bank=PASS rows=5")
    print("mu_substitution_identities=PASS degrees=6,10")


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


def faber_coefficients(
    degree: int,
    t_coeff: list[Poly | None],
) -> list[Poly]:
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
                    poly_scale(
                        bank[current][current + 1],
                        Fraction(1, 2),
                    ),
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
                        10,
                        pole_depth,
                        t_coeff,
                    )

    require(checked == 4508, "hostile-bank cardinality mismatch")
    require(survivors == 0, "hostile bank contains a survivor")

    split_lambda_zero = [
        None,
        {0: Fraction(2)},
        poly_add(ONE, {3: Fraction(1)}),
        {3: Fraction(1)},
        {},
    ]
    require(
        not laurent_system_is_consistent(10, 2, split_lambda_zero),
        "split Lambda-zero negative control survived",
    )

    square_control = [
        None,
        {0: Fraction(2)},
        ONE,
        {},
        {},
    ]
    require(
        laurent_system_is_consistent(10, 1, square_control),
        "exact-square positive control failed",
    )

    print(f"finite_hostile_bank=PASS cases={checked} survivors={survivors}")
    print("split_Lambda_zero_degree10_control=PASS")
    print("exact_square_positive_control=PASS")
    print("scope=FINITE_REFEREE_ONLY symbolic_valuation_proof_supplies_quantifiers")


if __name__ == "__main__":
    symbolic_faber_bank()
    exact_hostile_bank()
