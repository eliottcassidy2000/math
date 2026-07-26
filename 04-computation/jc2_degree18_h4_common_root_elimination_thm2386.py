#!/usr/bin/env python3
"""Exact companion for THM-2386's H4 common-root elimination.

The local argument in the theorem reduces every common root which is not
already on a closed multiple-P wall to the Q-multiple incidence

    P(1) = Q(1) = Q'(1) = 0.

This script derives that one-parameter family exactly and certifies that its
degree-nine residual cannot have the form H3*S3^2.  It supplies two exact
certificates:

1. a complete discriminant-factor/subresultant gcd-degree atlas; and
2. an independent coefficient-comparison Groebner basis equal to [1].

All arithmetic is over the rationals.  ``require`` is used deliberately so
that ``python`` and ``python -O`` execute identical validity checks.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive_integer_polynomial(
    expression: sp.Expr,
    *variables: sp.Symbol,
) -> sp.Poly:
    """Clear denominators and return the positive-leading primitive part."""

    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    _, integral = polynomial.clear_denoms(convert=True)
    _, primitive = integral.primitive()
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive


def subresultant_of_degree(
    sequence: list[sp.Expr],
    variable: sp.Symbol,
    degree: int,
) -> sp.Expr:
    matches = [item for item in sequence if sp.degree(item, variable) == degree]
    require(
        len(matches) == 1,
        f"subresultant sequence has {len(matches)} entries of degree {degree}",
    )
    return matches[0]


def vanishes_mod_parameter_factor(
    expression: sp.Expr,
    variable: sp.Symbol,
    parameter: sp.Symbol,
    factor: sp.Poly,
) -> bool:
    """Test whether every y-coefficient vanishes in Q[B]/(factor)."""

    polynomial = sp.Poly(expression, variable, domain=sp.QQ[parameter])
    for coefficient in polynomial.all_coeffs():
        coefficient_poly = sp.Poly(coefficient, parameter, domain=sp.QQ)
        if not sp.rem(coefficient_poly, factor).is_zero:
            return False
    return True


def normalized_factor_list(
    polynomial: sp.Poly,
) -> list[tuple[sp.Poly, int]]:
    """Return positive-leading primitive irreducible factors."""

    _, factors = sp.factor_list(polynomial)
    normalized: list[tuple[sp.Poly, int]] = []
    for factor, exponent in factors:
        primitive = primitive_integer_polynomial(
            factor.as_expr(), *polynomial.gens
        )
        normalized.append((primitive, exponent))
    return normalized


def main() -> None:
    y, bvar, cvar, dvar, wvar = sp.symbols("y B C D W")

    quartic_p = (
        245 * y**4
        + 1890 * bvar * y**2
        - 24300 * bvar**2
        + 122472 * dvar
    )
    sextic_q = (
        539 * y**6
        + 11340 * bvar * y**4
        + 183708 * cvar * y**3
        + (72900 * bvar**2 - 367416 * dvar) * y**2
        + (2361960 * bvar * cvar + 2480058 * wvar) * y
    )

    # The common root is moved to y=1 by the weighted action.  The three
    # incidence equations are linear in C,D,W and have a constant nonzero
    # determinant, so this is the complete Q-multiple common-root family.
    incidence_equations = [
        sp.expand(quartic_p.subs(y, 1)),
        sp.expand(sextic_q.subs(y, 1)),
        sp.expand(sp.diff(sextic_q, y).subs(y, 1)),
    ]
    incidence_matrix, incidence_rhs = sp.linear_eq_to_matrix(
        incidence_equations,
        [cvar, dvar, wvar],
    )
    incidence_determinant = sp.factor(incidence_matrix.det())
    require(
        incidence_determinant != 0,
        "Q-multiple incidence system stopped being uniquely solvable",
    )
    incidence_solution = sp.linsolve(
        (incidence_matrix, incidence_rhs),
        [cvar, dvar, wvar],
    )
    require(
        len(incidence_solution.args) == 1,
        "Q-multiple incidence system stopped having one solution",
    )
    solved_c, solved_d, solved_w = next(iter(incidence_solution))
    expected_c = -35 * (81 * bvar + 7) / sp.Integer(26244)
    expected_d = (
        5 * (4860 * bvar**2 - 378 * bvar - 49) / sp.Integer(122472)
    )
    expected_w = (
        4050 * bvar**2 + 395 * bvar + 7
    ) / sp.Integer(39366)
    require(
        all(
            sp.cancel(actual - expected) == 0
            for actual, expected in [
                (solved_c, expected_c),
                (solved_d, expected_d),
                (solved_w, expected_w),
            ]
        ),
        "Q-multiple incidence parametrization changed",
    )

    specialization = {
        cvar: expected_c,
        dvar: expected_d,
        wvar: expected_w,
    }
    specialized_p = sp.factor(quartic_p.subs(specialization))
    specialized_q = sp.factor(sextic_q.subs(specialization))
    expected_specialized_p = (
        35 * (y - 1) * (y + 1) * (54 * bvar + 7 * y**2 + 7)
    )
    expected_specialized_q = (
        7
        * y
        * (y - 1) ** 2
        * (
            1620 * bvar * y
            + 405 * bvar
            + 77 * y**3
            + 154 * y**2
            + 231 * y
            + 63
        )
    )
    require(
        sp.expand(specialized_p - expected_specialized_p) == 0
        and sp.expand(specialized_q - expected_specialized_q) == 0,
        "specialized P,Q factorizations changed",
    )

    p3 = sp.cancel(specialized_p / (y - 1))
    q4 = sp.cancel(specialized_q / (y - 1) ** 2)
    require(
        sp.expand(p3.subs(y, 1) - 140 * (27 * bvar + 7)) == 0
        and sp.expand(q4.subs(y, 1) - 525 * (27 * bvar + 7)) == 0,
        "simple-P incidence evaluations changed",
    )

    residual_r9 = sp.expand(4 * p3**3 + 49 * (y - 1) * q4**2)
    residual_poly = sp.Poly(residual_r9, y, domain=sp.QQ[bvar])
    residual_content, residual_primitive = residual_poly.primitive()
    require(
        residual_content == 343
        and residual_poly.degree() == 9
        and residual_poly.LC() == 73060029
        and residual_primitive.LC() == 213003,
        "degree-nine residual normalization changed",
    )

    factor_27 = sp.Poly(27 * bvar + 7, bvar)
    factor_54 = sp.Poly(54 * bvar + 7, bvar)
    factor_1215 = sp.Poly(1215 * bvar + 91, bvar)
    factor_quadratic = sp.Poly(
        2105352 * bvar**2 + 518049 * bvar + 31997,
        bvar,
    )
    factor_degree_13 = sp.Poly(
        2581517833820205830323568640000000 * bvar**13
        + 27433606764111716665224806400000000 * bvar**12
        + 22675744145697996379882215628800000 * bvar**11
        + 7555894329856417179696508849920000 * bvar**10
        + 1039840519524574652314595234244000 * bvar**9
        - 40352748361733406218104457728875 * bvar**8
        - 35168809165204401406565197022160 * bvar**7
        - 4785349181326208520124913840124 * bvar**6
        - 143144590616746631993379033504 * bvar**5
        + 35604040498012586373759791190 * bvar**4
        + 5212992075701472327189319920 * bvar**3
        + 330046225412419196564414580 * bvar**2
        + 10674830277647551598569056 * bvar
        + 143991306409917519520541,
        bvar,
    )
    expected_discriminant_factors = [
        (factor_27, 6),
        (factor_54, 3),
        (factor_1215, 3),
        (factor_quadratic, 3),
        (factor_degree_13, 1),
    ]

    discriminant = sp.Poly(
        sp.discriminant(residual_primitive.as_expr(), y),
        bvar,
    )
    actual_discriminant_factors = normalized_factor_list(discriminant)
    actual_factor_map = {
        factor.as_expr(): exponent
        for factor, exponent in actual_discriminant_factors
    }
    expected_factor_map = {
        factor.as_expr(): exponent
        for factor, exponent in expected_discriminant_factors
    }
    require(
        actual_factor_map == expected_factor_map,
        "degree-nine residual discriminant factorization changed: "
        f"{[(factor.as_expr(), exponent) for factor, exponent in actual_discriminant_factors]}",
    )
    require(
        factor_quadratic.is_irreducible and factor_degree_13.is_irreducible,
        "nonlinear discriminant factor stopped being irreducible over Q",
    )

    # Constant leading coefficients make subresultant specialization valid at
    # every parameter.  Sres_0 vanishes on each discriminant component;
    # Sres_1 vanishes only on 27B+7; Sres_2 never vanishes there.
    subresultants = sp.subresultants(
        residual_primitive.as_expr(),
        sp.diff(residual_primitive.as_expr(), y),
        y,
    )
    subresultant_profile = [sp.degree(item, y) for item in subresultants]
    require(
        subresultant_profile == list(range(9, -1, -1)),
        f"degree-nine subresultant profile changed: {subresultant_profile}",
    )
    sres_0 = subresultant_of_degree(subresultants, y, 0)
    sres_1 = subresultant_of_degree(subresultants, y, 1)
    sres_2 = subresultant_of_degree(subresultants, y, 2)
    gcd_degree_atlas: list[int] = []
    vanishing_atlas: list[tuple[bool, bool, bool]] = []
    for index, (factor, _) in enumerate(expected_discriminant_factors):
        vanishing = tuple(
            vanishes_mod_parameter_factor(item, y, bvar, factor)
            for item in [sres_0, sres_1, sres_2]
        )
        expected_vanishing = (
            (True, True, False)
            if index == 0
            else (True, False, False)
        )
        require(
            vanishing == expected_vanishing,
            f"subresultant atlas changed on factor {factor.as_expr()}",
        )
        vanishing_atlas.append(vanishing)
        gcd_degree_atlas.append(2 if index == 0 else 1)

    linear_gcds = []
    for factor, _ in expected_discriminant_factors[:3]:
        root = sp.solve(factor.as_expr(), bvar)[0]
        specialized_residual = sp.Poly(
            residual_primitive.as_expr().subs(bvar, root),
            y,
            domain=sp.QQ,
        )
        specialized_gcd = sp.gcd(
            specialized_residual,
            specialized_residual.diff(),
        ).monic()
        linear_gcds.append(specialized_gcd)
    require(
        [item.as_expr() for item in linear_gcds]
        == [
            sp.expand((y - 1) ** 2),
            y,
            y + 1,
        ],
        "rational-point gcd atlas changed: "
        f"{[item.as_expr() for item in linear_gcds]}",
    )
    require(
        max(gcd_degree_atlas) == 2,
        "a degree-three square divisor appeared in the subresultant atlas",
    )

    # Independent coefficient comparison.  Any factorization R9=H3*S3^2
    # can be scaled so that S3 is monic and lc(H3)=lc(R9).
    a0, a1, a2, h0, h1, h2 = sp.symbols("a0 a1 a2 h0 h1 h2")
    square_root = y**3 + a2 * y**2 + a1 * y + a0
    squarefree_carrier = 73060029 * y**3 + h2 * y**2 + h1 * y + h0
    coefficient_difference = sp.Poly(
        sp.expand(squarefree_carrier * square_root**2 - residual_r9),
        y,
    )
    equation_8 = coefficient_difference.coeff_monomial(y**8)
    solved_h2 = sp.solve(equation_8, h2)[0]
    equation_7 = sp.expand(
        coefficient_difference.coeff_monomial(y**7).subs(h2, solved_h2)
    )
    solved_h1 = sp.solve(equation_7, h1)[0]
    equation_6 = sp.expand(
        coefficient_difference.coeff_monomial(y**6).subs(
            {h2: solved_h2, h1: solved_h1}
        )
    )
    solved_h0 = sp.solve(equation_6, h0)[0]
    expected_h2 = -73060029 * (2 * a2 - 3)
    expected_h1 = -453789 * (
        -4320 * bvar
        + 322 * a1
        - 483 * a2**2
        + 966 * a2
        - 966
    )
    expected_h0 = -16807 * (
        233280 * bvar * a2
        - 287550 * bvar
        + 8694 * a0
        - 26082 * a1 * a2
        + 26082 * a1
        + 17388 * a2**3
        - 39123 * a2**2
        + 52164 * a2
        - 38080
    )
    require(
        sp.expand(solved_h2 - expected_h2) == 0
        and sp.expand(solved_h1 - expected_h1) == 0
        and sp.expand(solved_h0 - expected_h0) == 0,
        "top coefficient triangular solve changed",
    )
    coefficient_substitution = {
        h2: solved_h2,
        h1: solved_h1,
        h0: solved_h0,
    }
    remaining_equations = [
        primitive_integer_polynomial(
            sp.expand(
                coefficient_difference.coeff_monomial(y**power).subs(
                    coefficient_substitution
                )
            ),
            bvar,
            a2,
            a1,
            a0,
        )
        for power in range(5, -1, -1)
    ]
    equation_signatures = [
        (len(equation.terms()), equation.total_degree())
        for equation in remaining_equations
    ]
    require(
        equation_signatures
        == [(16, 4), (20, 5), (23, 5), (25, 5), (17, 5), (13, 5)],
        f"coefficient equation signatures changed: {equation_signatures}",
    )

    print("THM-2386 degree-18 H4 common-root elimination exact companion")
    print(f"incidence determinant: {incidence_determinant}")
    print("incidence family: C,D,W are the stated rational functions of B")
    print("specialized orders: ord_1(P)=1, ord_1(Q)=2 when 27B+7!=0")
    print(
        "R9 normalization: degree=9, content=343, "
        "leading coefficient=73060029"
    )
    print("discriminant factors (degree, exponent): 1^6,1^3,1^3,2^3,13^1")
    print("subresultant gcd-degree atlas: 2,1,1,1,1")
    print("rational gcds: (y-1)^2; y; y+1")
    print(
        "coefficient equation signatures (terms,total-degree): "
        "16/4,20/5,23/5,25/5,17/5,13/5"
    )
    print("computing independent exact coefficient-comparison basis...", flush=True)

    coefficient_basis = sp.groebner(
        [equation.as_expr() for equation in remaining_equations],
        a0,
        a1,
        a2,
        bvar,
        order="grevlex",
        domain=sp.QQ,
    )
    require(
        len(coefficient_basis.polys) == 1
        and coefficient_basis.polys[0].is_one,
        "coefficient-comparison ideal stopped being the unit ideal",
    )
    print("coefficient-comparison Groebner basis: [1]")
    print("maximum specialized gcd degree: 2 < 3")
    print(
        "conclusion: the Q-multiple residual has no H3*S3^2 factorization"
    )


if __name__ == "__main__":
    main()
