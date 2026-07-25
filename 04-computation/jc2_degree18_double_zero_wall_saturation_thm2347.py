#!/usr/bin/env python3
"""Exact referee for THM-2347's double-zero-wall saturation.

On 20*B*C+21*W=0, the structured sextic Q gains a second zero at y=0.
In the B!=0, C!=0 chart, normalize B=1 and set

    X=C^2,  J=126*D-25.

A low square class forces gcd(F,F') to have degree at least four.  Three
coefficients from the degree-three and degree-two subresultants suffice:
after stripping their exact J- and C-factors, their grevlex ideal contains
J^2.  Hence the low-square-class locus lies on the common-root wall J=0.

Every executable check remains active under optimized Python.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def exact_quotient(
    expression: sp.Expr,
    divisor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> sp.Expr:
    """Return an exact polynomial quotient, failing on a nonzero remainder."""

    quotient, remainder = sp.div(
        sp.Poly(expression, *variables),
        sp.Poly(divisor, *variables),
    )
    require(remainder.is_zero, f"expected exact division by {divisor}")
    return quotient.as_expr()


def primitive_shifted_even_polynomial(
    expression: sp.Expr,
    cvar: sp.Symbol,
    dvar: sp.Symbol,
    xvar: sp.Symbol,
    jvar: sp.Symbol,
) -> sp.Expr:
    """Replace C^(2k) by X^k and D by (J+25)/126, then primitive-normalize."""

    polynomial = sp.Poly(expression, cvar, dvar)
    require(
        all(c_exponent % 2 == 0 for (c_exponent, _), _ in polynomial.terms()),
        "the C-parity reduction to X=C^2 failed",
    )
    in_x_d = sum(
        coefficient * xvar ** (c_exponent // 2) * dvar**d_exponent
        for (c_exponent, d_exponent), coefficient in polynomial.terms()
    )
    shifted = sp.together(in_x_d.subs(dvar, (jvar + 25) / 126))
    numerator, denominator = shifted.as_numer_denom()
    require(
        denominator.is_Integer and denominator != 0,
        "the shifted polynomial acquired a nonconstant denominator",
    )
    content, primitive = sp.Poly(numerator, xvar, jvar).primitive()
    require(content != 0, "primitive normalization met the zero polynomial")
    return primitive.as_expr()


def main() -> None:
    y = sp.symbols("y")
    bvar, cvar, dvar, wvar = sp.symbols("B C D W")
    xvar, jvar = sp.symbols("X J")

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
    mordell = sp.expand(4 * quartic_p**3 + 49 * sextic_q**2)

    # Normalize B=1 and impose 20*C+21*W=0.  The last coefficient of Q
    # vanishes, so Q has the promised double zero at the origin.
    normalized_q = sp.expand(
        sextic_q.subs(
            {
                bvar: 1,
                wvar: -sp.Rational(20, 21) * cvar,
            }
        )
    )
    expected_q = y**2 * (
        539 * y**4
        + 11340 * y**2
        + 183708 * cvar * y
        + 72900
        - 367416 * dvar
    )
    require(
        sp.expand(normalized_q - expected_q) == 0,
        "the double-zero factorization of Q changed",
    )

    normalized_p = sp.expand(quartic_p.subs(bvar, 1))
    normalized_f = sp.expand(4 * normalized_p**3 + 49 * normalized_q**2)
    normalized_f_poly = sp.Poly(normalized_f, y)
    require(
        normalized_f_poly.degree() == 12
        and normalized_f_poly.LC() == 73060029,
        "the normalized Mordell polynomial changed degree or leading coefficient",
    )

    subresultants = sp.subresultants(
        normalized_f,
        sp.diff(normalized_f, y),
        y,
    )
    require(
        [sp.degree(item, y) for item in subresultants]
        == list(range(12, -1, -1)),
        "the complete subresultant degree profile changed",
    )
    degree_three = sp.Poly(subresultants[-4], y)
    degree_two = sp.Poly(subresultants[-3], y)
    require(
        degree_three.degree() == 3 and degree_two.degree() == 2,
        "the selected terminal subresultants changed degree",
    )

    common_wall = 126 * dvar - 25
    selected = [
        (degree_three.LC(), common_wall**3, 0),
        (degree_three.TC(), common_wall**5, 1),
        (degree_two.LC(), common_wall**6, 0),
    ]

    shifted_generators: list[sp.Expr] = []
    exact_orders: list[tuple[int, int]] = []
    for expression, wall_power, c_power in selected:
        after_wall = exact_quotient(
            expression,
            wall_power,
            (cvar, dvar),
        )
        _, extra_wall_remainder = sp.div(
            sp.Poly(after_wall, cvar, dvar),
            sp.Poly(common_wall, cvar, dvar),
        )
        require(
            not extra_wall_remainder.is_zero,
            "a selected coefficient has a higher common-wall order than recorded",
        )
        after_c = exact_quotient(
            after_wall,
            cvar**c_power,
            (cvar, dvar),
        )
        if c_power:
            _, extra_c_remainder = sp.div(
                sp.Poly(after_c, cvar, dvar),
                sp.Poly(cvar, cvar, dvar),
            )
            require(
                not extra_c_remainder.is_zero,
                "the odd selected coefficient has a higher C-order than recorded",
            )
        shifted_generators.append(
            primitive_shifted_even_polynomial(
                after_c,
                cvar,
                dvar,
                xvar,
                jvar,
            )
        )
        exact_orders.append((sp.degree(wall_power, dvar), c_power))

    require(
        [
            (
                sp.total_degree(generator),
                sp.degree(generator, xvar),
                sp.degree(generator, jvar),
                len(sp.Poly(generator, xvar, jvar).terms()),
            )
            for generator in shifted_generators
        ]
        == [
            (15, 9, 15, 89),
            (13, 8, 13, 70),
            (16, 10, 16, 104),
        ],
        "the shifted generator signatures changed",
    )

    basis = sp.groebner(
        shifted_generators,
        xvar,
        jvar,
        order="grevlex",
        method="f5b",
    )
    large_factor = (
        1233383522515969671 * jvar * xvar**2
        + 27991525644715500 * jvar * xvar
        + 158573241187500 * jvar
        - 659000910244935988096920 * xvar**5
        - 54559736202362496221340 * xvar**4
        - 1758935316808435500000 * xvar**3
        - 27710180945527500000 * xvar**2
        - 213905002500000000 * xvar
        - 648671875000000
    )
    expected_basis = [
        -(30618 * xvar + 361) ** 2 * large_factor,
        jvar
        * (30618 * xvar + 361) ** 3
        * (199644669 * xvar**2 + 7654500 * xvar + 62500),
        jvar**2,
    ]
    actual_basis = [item.as_expr() for item in basis.polys]
    require(
        len(actual_basis) == 3
        and all(
            sp.expand(actual - expected) == 0
            for actual, expected in zip(actual_basis, expected_basis)
        ),
        "the exact three-element Groebner basis changed",
    )
    _, square_remainder = basis.reduce(jvar**2)
    require(
        square_remainder == 0,
        "J^2 stopped belonging to the stripped coefficient ideal",
    )

    # Positive and hostile controls for the actual polynomial gcd.
    def gcd_degree(cvalue: sp.Expr, dvalue: sp.Expr) -> int:
        specialization = sp.Poly(
            normalized_f.subs({cvar: cvalue, dvar: dvalue}),
            y,
            domain=sp.QQ,
        )
        return sp.gcd(specialization, specialization.diff()).degree()

    deep_degree = gcd_degree(1, sp.Rational(25, 126))
    hostile_degrees = (gcd_degree(1, 0), gcd_degree(1, 1))
    require(
        deep_degree == 5 and hostile_degrees == (0, 0),
        "the deep/off-deep gcd controls changed",
    )

    print("THM-2347 exact double-zero-wall saturation referee")
    print("normalized wall: 20*C+21*W=0 and Q=y^2*R4")
    print("Mordell degree/leading coefficient: 12 / 73060029")
    print("subresultant degrees: 12,11,10,9,8,7,6,5,4,3,2,1,0")
    print(f"selected exact (J-order,C-order): {exact_orders}")
    print(
        "shifted generator signatures:"
        " [(15,9,15,89),(13,8,13,70),(16,10,16,104)]"
    )
    print("grevlex variables: (X,J), where X=C^2 and J=126*D-25")
    print("Groebner basis degrees: 7,6,2")
    print(
        "middle basis factors:"
        " J*(30618*X+361)^3"
        "*(199644669*X^2+7654500*X+62500)"
    )
    print("terminal basis element: J^2")
    print(f"deep control gcd degree: {deep_degree}")
    print(f"off-deep hostile gcd degrees: {hostile_degrees}")
    print("VERDICT: gcd degree >=4 forces 126*D-25=0 on B=1, C!=0")


if __name__ == "__main__":
    main()
