#!/usr/bin/env python3
"""Exact low-degree companion for the THM-3494 discriminant atlas.

The proof of THM-3494 is all-degree and uses the trace-form discriminant
identity together with THM-3438's S_n monodromy.  This companion supplies a
separate exact resultant route in the first three grades n=3,4,5.  Over
Q(P,Q,C) it constructs the canonical weighted seeds, their inverse equations,
the x- and y-coordinate eliminants, and checks

    Disc(E_coordinate) / Disc(T) = (rational index)^2.

It also computes the z-coordinate cubic in the n=3 row, checks that the branch
carrier and every index core are coprime, and freezes a flat-coordinate hostile
whose eliminant has zero discriminant.  All gates use ``require`` so optimized
Python cannot disable them.
"""

from __future__ import annotations

from hashlib import sha256
from math import isqrt

import sympy as sp


EXPECTED_SEMANTIC_SHA256 = "5d327ad97829f806a16fbe9116329c32c3941fe1def2de6234155a9f8e31e075"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def polynomial_digest(expression: sp.Expr, generators: tuple[sp.Symbol, ...]) -> str:
    polynomial = sp.Poly(sp.cancel(expression), *generators, domain=sp.QQ)
    payload = tuple(
        (monomial, int(coefficient.p), int(coefficient.q))
        for monomial, coefficient in polynomial.terms()
    )
    return sha256(repr(payload).encode("ascii")).hexdigest()


def polynomial_stats(expression: sp.Expr, generators: tuple[sp.Symbol, ...]) -> tuple[object, ...]:
    polynomial = sp.Poly(sp.cancel(expression), *generators, domain=sp.QQ)
    return (
        len(polynomial.terms()),
        polynomial.total_degree(),
        tuple(polynomial.degree(generator) for generator in generators),
        polynomial_digest(expression, generators),
    )


def primitive_integer_part(
    expression: sp.Expr, generators: tuple[sp.Symbol, ...]
) -> tuple[sp.Expr, sp.Rational]:
    polynomial = sp.Poly(expression, *generators, domain=sp.QQ)
    _, integral = polynomial.clear_denoms(convert=True)
    _, primitive = integral.primitive()
    if primitive.LC() < 0:
        primitive = -primitive
    primitive_expression = primitive.as_expr()
    scalar = sp.cancel(expression / primitive_expression)
    require(not scalar.free_symbols, ("nonconstant primitive scalar", scalar))
    return primitive_expression, sp.Rational(scalar)


def rational_square_root(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    require(value >= 0, ("negative rational square", value))
    numerator_root = isqrt(int(value.p))
    denominator_root = isqrt(int(value.q))
    require(
        numerator_root * numerator_root == value.p
        and denominator_root * denominator_root == value.q,
        ("nonsquare rational", value),
    )
    return sp.Rational(numerator_root, denominator_root)


def perfect_square_root(
    expression: sp.Expr, generators: tuple[sp.Symbol, ...]
) -> sp.Expr:
    expression = sp.cancel(expression)
    numerator, denominator = expression.as_numer_denom()
    numerator_coefficient, numerator_factors = sp.factor_list(numerator, *generators)
    denominator_coefficient, denominator_factors = sp.factor_list(denominator, *generators)
    coefficient_root = rational_square_root(
        sp.Rational(numerator_coefficient) / sp.Rational(denominator_coefficient)
    )
    root = coefficient_root
    for factor, exponent in numerator_factors:
        require(exponent % 2 == 0, ("odd numerator exponent", factor, exponent))
        root *= factor ** (exponent // 2)
    for factor, exponent in denominator_factors:
        require(exponent % 2 == 0, ("odd denominator exponent", factor, exponent))
        root /= factor ** (exponent // 2)
    root = sp.factor(root)
    require(sp.cancel(root**2 - expression) == 0, ("square reconstruction", expression, root))
    return root


w, P, Q, C, X, Y, Z = sp.symbols("w P Q C X Y Z")


def canonical_seed(d: int) -> sp.Expr:
    require(d >= 2, d)
    delta = sp.Rational(6, d * (d + 1))
    return sp.expand(
        2 * w
        - 3 * w**2
        + w ** (d - 1)
        - w**d
        - delta * w
        + delta * w**2
    )


def coordinate_row(
    label: str,
    inverse_equation: sp.Expr,
    branch_discriminant: sp.Expr,
    branch_primitive: sp.Expr,
    relation: sp.Expr,
    coordinate: sp.Symbol,
    degree: int,
    c_root_exponent: int,
) -> tuple[object, ...]:
    eliminant = sp.cancel(sp.resultant(inverse_equation, relation, w))
    eliminant_polynomial = sp.Poly(eliminant, coordinate, P, Q, C, domain=sp.QQ)
    require(eliminant_polynomial.degree(coordinate) == degree, (label, "degree", eliminant))
    coordinate_discriminant = sp.factor(sp.discriminant(eliminant, coordinate))
    require(coordinate_discriminant != 0, (label, "zero discriminant"))
    ratio = sp.cancel(coordinate_discriminant / branch_discriminant)
    index_root = perfect_square_root(ratio, (P, Q, C))
    index_core = sp.cancel(index_root / C**c_root_exponent)
    require(not index_core.has(C), (label, "C survives", index_core))
    index_core_polynomial = sp.Poly(index_core, P, Q, domain=sp.QQ)
    if index_core_polynomial.LC() < 0:
        index_core = -index_core
        index_root = -index_root
        index_core_polynomial = -index_core_polynomial
    require(sp.cancel(index_root**2 - ratio) == 0, (label, "signed square root"))
    branch_polynomial = sp.Poly(branch_primitive, P, Q, domain=sp.QQ)
    require(
        sp.gcd(branch_polynomial, index_core_polynomial).total_degree() == 0,
        (label, "branch/index common factor"),
    )
    return (
        label,
        len(eliminant_polynomial.terms()),
        polynomial_digest(eliminant, (coordinate, P, Q, C)),
        c_root_exponent,
        polynomial_stats(index_core, (P, Q)),
        sp.sstr(index_core) if degree == 3 else "hash-pinned-in-row",
    )


def degree_row(degree: int) -> tuple[object, ...]:
    d = degree - 1
    seed = canonical_seed(d)
    integral = sp.integrate(seed, w)
    inverse_equation = sp.expand(integral - P * w + Q)
    gamma = P - seed

    require(sp.expand(sp.diff(inverse_equation, w) - (seed - P)) == 0, degree)
    require(sp.Poly(inverse_equation, w).degree() == degree, degree)
    require(seed.subs(w, 0) == 0 and seed.subs(w, 1) == -1, (degree, seed))
    require(integral.subs(w, 1) == 0, (degree, integral.subs(w, 1)))

    branch_discriminant = sp.factor(sp.discriminant(inverse_equation, w))
    resultant_discriminant = sp.resultant(inverse_equation, sp.diff(inverse_equation, w), w)
    leading_coefficient = sp.Poly(inverse_equation, w).LC()
    require(
        sp.cancel(
            branch_discriminant
            - (-1) ** (degree * (degree - 1) // 2)
            * resultant_discriminant
            / leading_coefficient
        )
        == 0,
        (degree, "resultant/discriminant sign"),
    )

    branch_primitive, branch_scalar = primitive_integer_part(branch_discriminant, (P, Q))
    factor_coefficient, factors = sp.factor_list(branch_primitive, P, Q)
    require(
        sp.Rational(factor_coefficient) in (1, -1)
        and len(factors) == 1
        and factors[0][1] == 1,
        (degree, "reducible branch", factor_coefficient, factors),
    )

    c_index_exponent = degree * (degree - 1) // 2
    x_row = coordinate_row(
        "x",
        inverse_equation,
        branch_discriminant,
        branch_primitive,
        X * gamma - C,
        X,
        degree,
        c_index_exponent,
    )
    y_row = coordinate_row(
        "y",
        inverse_equation,
        branch_discriminant,
        branch_primitive,
        C * Y - w + gamma,
        Y,
        degree,
        c_index_exponent,
    )

    flat_eliminant = sp.factor(sp.resultant(inverse_equation, X - P, w))
    require(sp.cancel(flat_eliminant - leading_coefficient ** 0 * (X - P) ** degree) == 0, (degree, flat_eliminant))
    require(sp.discriminant(flat_eliminant, X) == 0, (degree, "flat hostile fires"))

    return (
        degree,
        sp.sstr(seed),
        sp.sstr(inverse_equation),
        branch_scalar,
        polynomial_stats(branch_primitive, (P, Q)),
        sp.sstr(branch_primitive),
        x_row,
        y_row,
        "flat_coordinate_resultant=(X-P)^n_has_zero_discriminant",
    )


def cubic_z_row() -> tuple[object, ...]:
    degree = 3
    d = 2
    seed = canonical_seed(d)
    inverse_equation = sp.expand(sp.integrate(seed, w) - P * w + Q)
    gamma = P - seed
    kappa = sp.diff(seed, w).subs(w, 1)
    a = sp.factor(-(1 + kappa) / (2 + kappa))
    require(a == sp.Rational(-3, 2), a)
    z_numerator = sp.expand(gamma * (gamma * (gamma - 1 + a) - a * w))
    branch_discriminant = sp.factor(sp.discriminant(inverse_equation, w))
    branch_primitive, _ = primitive_integer_part(branch_discriminant, (P, Q))
    return coordinate_row(
        "z",
        inverse_equation,
        branch_discriminant,
        branch_primitive,
        C**2 * Z - z_numerator,
        Z,
        degree,
        degree * (degree - 1),
    )


def main() -> None:
    rows = tuple(degree_row(degree) for degree in (3, 4, 5))
    z_row = cubic_z_row()

    branch_hashes = tuple(row[4][3] for row in rows)
    coordinate_hashes = tuple(
        (row[6][4][3], row[7][4][3])
        for row in rows
    )
    semantic_surface = (
        "universe=Q(P,Q,C);canonical_THM3438_seeds_d=2,3,4",
        "T=integral(p)-Pw+Q;T_prime=p-P=-gamma",
        "Ex=Res_w(T,X*gamma-C);Ey=Res_w(T,CY-w+gamma)",
        "Disc(Ex)/Disc(T)=Ix^2;Disc(Ey)/Disc(T)=Iy^2",
        "branch_and_index_cores_coprime;flat_view_has_zero_discriminant",
        rows,
        z_row,
    )
    semantic_sha256 = sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    require(
        semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
        (semantic_sha256, EXPECTED_SEMANTIC_SHA256),
    )

    print("THM-3494 weighted-lift primitive-coordinate discriminant atlas")
    print("status=FINITE-EXACT_LOW_DEGREE_COMPANION;all_degree_step_is_proof_text")
    print("universe=Q(P,Q,C);degrees=3,4,5;coordinates=x,y;degree3_extra_coordinate=z")
    for row in rows:
        degree, seed, inverse_equation, scalar, branch_stats, branch, x_row, y_row, hostile = row
        print(f"n={degree};p={seed};T={inverse_equation}")
        print(f"n={degree};branch_scalar={scalar};branch_stats={branch_stats};branch={branch}")
        print(f"n={degree};x_row={x_row}")
        print(f"n={degree};y_row={y_row}")
        print(f"n={degree};hostile={hostile}")
    print(f"n=3;z_row={z_row}")
    print(f"branch_hashes={branch_hashes}")
    print(f"coordinate_index_hashes={coordinate_hashes}")
    print("mechanism=common_reduced_branch_carrier_times_coordinate_projection_index_square")
    print("scope=no_arbitrary_map_classification;no_composition_norm;no_LRC_or_tournament_transfer")
    print(f"semantic_sha256={semantic_sha256}")


if __name__ == "__main__":
    main()
