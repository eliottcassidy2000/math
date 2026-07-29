#!/usr/bin/env python3
"""Exact referee for the Gamma positive-cone y=0 exit portal."""

from __future__ import annotations

from hashlib import sha256
from math import comb
from pathlib import Path
import sys

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

BASE = COMP / "gmc_gamma_transverse_null_holotopy_thm2865.py"
require(
    sha256(BASE.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    == "ff32c787fb4fef58b77d06d84f3b032197744515ec03d6dcff8c1be6e70afc52",
    "THM-2865 exact tensor builder changed",
)

import gmc_gamma_transverse_null_holotopy_thm2865 as gamma


ALPHA, X, Y = gamma.ALPHA, gamma.X, gamma.Y

P9_COEFFICIENTS = (
    8,
    3364,
    145397,
    2870386,
    32586944,
    229937674,
    1025612391,
    2810119176,
    4307331852,
    2819687976,
)

P23_COEFFICIENTS = (
    405,
    67862,
    4458573,
    152480084,
    2640768818,
    739093820,
    -1227760066990,
    -36903041464160,
    -665697843962479,
    -8704796296181594,
    -88073949221238223,
    -711703149831960260,
    -4673300030645893432,
    -25164153182568166184,
    -111485361916416386032,
    -405718990297087577216,
    -1204903266566658216448,
    -2884528384378398958208,
    -5458050914397962088448,
    -7916087911348550367232,
    -8378757135131489026048,
    -5946662811236683546624,
    -2388466752014210826240,
    -339470562155929010176,
)


def sign_variations(coefficients: tuple[int, ...]) -> int:
    signs = [1 if coefficient > 0 else -1 for coefficient in coefficients if coefficient]
    return sum(left != right for left, right in zip(signs, signs[1:]))


def bernstein_coefficients(
    expression: sp.Expr,
    left: sp.Rational,
    right: sp.Rational,
) -> tuple[sp.Expr, ...]:
    parameter = sp.symbols("parameter")
    polynomial = sp.Poly(
        sp.expand(
            expression.subs(
                ALPHA, left + (right - left) * parameter
            )
        ),
        parameter,
    )
    degree = polynomial.degree()
    power = tuple(polynomial.nth(index) for index in range(degree + 1))
    return tuple(
        sum(
            power[index]
            * sp.Rational(comb(bernstein_index, index), comb(degree, index))
            for index in range(bernstein_index + 1)
        )
        for bernstein_index in range(degree + 1)
    )


def main() -> None:
    lower = {1: sp.Integer(1), 3: X}
    upper = {2: sp.Integer(1), 3: Y}

    g11 = gamma.multilinear((lower, lower))
    g12 = gamma.multilinear((lower, upper))
    g22 = gamma.multilinear((upper, upper))
    t111 = gamma.multilinear((lower, lower, lower))
    t112 = gamma.multilinear((lower, lower, upper))
    t122 = gamma.multilinear((lower, upper, upper))
    t222 = gamma.multilinear((upper, upper, upper))

    invariant_one = sp.factor(
        3 * t112 * g11 * g22
        - t222 * g11**2
        - 2 * t111 * g12 * g22
    )
    invariant_two = sp.factor(
        3 * t122 * g11 * g22
        - 2 * t222 * g12 * g11
        - t111 * g22**2
    )
    numerator_one, denominator_one = sp.together(
        invariant_one
    ).as_numer_denom()
    numerator_two, denominator_two = sp.together(
        invariant_two
    ).as_numer_denom()
    common_denominator = (
        ALPHA**4
        * (ALPHA + 1) ** 4
        * (ALPHA + 2) ** 4
        * (ALPHA + 3) ** 3
    )
    require(
        sp.expand(denominator_one - common_denominator) == 0
        and sp.expand(denominator_two - common_denominator) == 0,
        "positive invariant denominator changed",
    )
    boundary_one = sp.expand(numerator_one.subs(Y, 0))
    boundary_two = sp.expand(numerator_two.subs(Y, 0))

    p9 = sp.Poly.from_list(P9_COEFFICIENTS, gens=ALPHA).as_expr()
    p23 = sp.Poly.from_list(P23_COEFFICIENTS, gens=ALPHA).as_expr()
    resultant = sp.resultant(boundary_one, boundary_two, X)
    expected_resultant = (
        -sp.Integer(7077888)
        * (ALPHA + 2) ** 6
        * (ALPHA + 3) ** 22
        * (ALPHA + 4) ** 4
        * (ALPHA + 5) ** 2
        * (ALPHA + 8) ** 2
        * p9
        * p23
    )
    require(
        sp.expand(resultant - expected_resultant) == 0,
        "boundary resultant factorization changed",
    )
    require(
        all(coefficient > 0 for coefficient in P9_COEFFICIENTS),
        "degree-nine factor gained a positive root candidate",
    )
    require(
        sign_variations(P23_COEFFICIENTS) == 1,
        "degree-23 Descartes pattern changed",
    )

    alpha_left = sp.Rational(232446, 10000)
    alpha_right = sp.Rational(232447, 10000)
    require(
        p23.subs(ALPHA, alpha_left) < 0
        < p23.subs(ALPHA, alpha_right)
        and sp.Poly(p23, ALPHA).count_roots(
            alpha_left, alpha_right
        )
        == 1
        and sp.gcd(p23, sp.diff(p23, ALPHA)) == 1,
        "positive Gamma-shape isolation changed",
    )

    first = sp.Poly(boundary_one, X).primitive()[1].as_expr()
    second = sp.Poly(boundary_two, X).primitive()[1].as_expr()
    subresultants = sp.subresultants(first, second, X)
    require(
        tuple(sp.degree(polynomial, X) for polynomial in subresultants)
        == (4, 3, 2, 1, 0),
        "boundary subresultant degree profile changed",
    )
    content, primitive_linear = sp.Poly(
        subresultants[-2], X
    ).primitive()
    require(
        sp.expand(
            content
            - 12
            * (ALPHA + 2) ** 3
            * (ALPHA + 3) ** 2
            * (ALPHA + 5) ** 2
            * (ALPHA + 8)
        )
        == 0,
        "linear subresultant content changed",
    )
    denominator_x = sp.Poly(primitive_linear, X).LC()
    numerator_x = -sp.Poly(primitive_linear, X).nth(0)
    require(
        sp.degree(denominator_x, ALPHA) == 26
        and sp.degree(numerator_x, ALPHA) == 26
        and all(
            coefficient > 0
            for coefficient in sp.Poly(
                denominator_x, ALPHA
            ).all_coeffs()
        )
        and all(
            coefficient > 0
            for coefficient in sp.Poly(
                numerator_x, ALPHA
            ).all_coeffs()
        )
        and sp.gcd(denominator_x, p23) == 1,
        "positive specialized x formula changed",
    )

    inverse_denominator_x = sp.invert(
        denominator_x, p23, domain=sp.QQ
    )
    x_mod_p23 = sp.rem(
        numerator_x * inverse_denominator_x, p23, domain=sp.QQ
    )

    def evaluate_x_mod_p23(polynomial: sp.Expr) -> sp.Expr:
        value = sp.Integer(0)
        for coefficient in sp.Poly(polynomial, X).all_coeffs():
            value = sp.rem(
                value * x_mod_p23 + coefficient,
                p23,
                domain=sp.QQ,
            )
        return value

    require(
        evaluate_x_mod_p23(first) == 0
        and evaluate_x_mod_p23(second) == 0,
        "linear subresultant did not restore the common root",
    )

    # A narrow exact bracket for the positive x-coordinate.
    x_left = sp.Rational(71257, 100000)
    x_right = sp.Rational(71259, 100000)
    require(
        all(
            coefficient > 0
            for coefficient in bernstein_coefficients(
                numerator_x - x_left * denominator_x,
                alpha_left,
                alpha_right,
            )
        )
        and all(
            coefficient > 0
            for coefficient in bernstein_coefficients(
                x_right * denominator_x - numerator_x,
                alpha_left,
                alpha_right,
            )
        ),
        "positive x-coordinate bracket changed",
    )

    # The Gram form is strictly positive at every alpha>0 because the two
    # real polynomials d1+x*d3 and d2 are linearly independent.  The exact
    # algebra above is therefore a genuine quadratic/cubic null portal.
    require(
        sp.Poly(first, X).degree() == 4
        and sp.Poly(second, X).degree() == 3
        and all(
            coefficient < 0
            for coefficient in sp.Poly(
                sp.Poly(first, X).LC(), ALPHA
            ).all_coeffs()
        )
        and all(
            coefficient < 0
            for coefficient in sp.Poly(
                sp.Poly(second, X).LC(), ALPHA
            ).all_coeffs()
        ),
        "specialization changed the boundary polynomial degrees",
    )

    print("GAMMA POSITIVE-CONE Y-ZERO EXIT PORTAL -- exact referee")
    print("boundary=U=d1+x*d3; V=d2")
    print(
        "resultant_factors=negative_constant*(alpha+2)^6"
        "*(alpha+3)^22*(alpha+4)^4*(alpha+5)^2*(alpha+8)^2"
        "*P9*P23"
    )
    print("P9_coefficients_positive=10/10; P23_sign_variations=1")
    print(f"alpha_unique_positive_bracket=({alpha_left},{alpha_right})")
    print("subresultant_degrees=(4,3,2,1,0); linear_gcd=PASS")
    print("x_formula_degrees=(26,26); positive_coefficients=(27,27)")
    print(f"x_positive_bracket=({x_left},{x_right})")
    print("common_root_remainders_mod_P23=(0,0)")
    print(
        "consequence=exact unique positive Gamma shape and x-coordinate"
        " on the y=0 quadratic/cubic-null face"
    )
    print(
        "scope=boundary portal only;"
        " connection to the THM-2865 local branch not proved"
    )


if __name__ == "__main__":
    main()
