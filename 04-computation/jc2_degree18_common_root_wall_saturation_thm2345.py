#!/usr/bin/env python3
"""Exact referee for THM-2345's common-root-wall saturation.

On 126D=25B^2, factor the structured Mordell polynomial as y^2 times a
degree-ten polynomial.  In the B!=0 chart, exact subresultants prove that
gcd(G10,G10') can have degree at least three only on 20C+21W=0.  This is
the deep wall already closed by THM-2338--THM-2342.

Every executable check remains active under optimized Python.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive_after_division(
    polynomial: sp.Expr,
    divisor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> sp.Expr:
    """Divide exactly and discard only the nonzero integer content."""

    quotient, remainder = sp.div(polynomial, divisor, *variables)
    require(remainder == 0, f"expected exact division by {divisor}")
    content, primitive = sp.Poly(quotient, *variables).primitive()
    require(content != 0, "primitive-part extraction met the zero polynomial")
    return primitive.as_expr()


def main() -> None:
    y = sp.symbols("y")
    bvar, cvar, dvar, wvar = sp.symbols("B C D W")

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

    common_root_wall = {dvar: sp.Rational(25, 126) * bvar**2}
    wall_mordell = sp.factor(mordell.subs(common_root_wall))
    residual_ten = sp.Poly(
        sp.cancel(wall_mordell / (1323 * y**2)),
        y,
        domain=sp.QQ[bvar, cvar, wvar],
    ).as_expr()
    require(
        sp.expand(wall_mordell - 1323 * y**2 * residual_ten) == 0
        and sp.degree(residual_ten, y) == 10
        and sp.LC(sp.Poly(residual_ten, y)) == 55223,
        "common-root-wall factorization changed",
    )

    # On B!=0, weighted scaling normalizes B=1.  The normalized equation
    # 20C+21W=0 is equivalent to the invariant equation 20BC+21W=0.
    normalized_ten = sp.expand(residual_ten.subs(bvar, 1))
    deep_linear = 20 * cvar + 21 * wvar

    residual_six = (
        7889 * y**6
        + 211680 * y**4
        + 1047816 * cvar * y**3
        + 1814400 * y**2
        + 22044960 * cvar * y
        + 2916000
        + 178564176 * cvar**2
    )
    require(
        sp.expand(
            normalized_ten.subs(wvar, -sp.Rational(20, 21) * cvar)
            - 7 * y**4 * residual_six
        )
        == 0,
        "deep-wall y^4 factorization changed",
    )

    # For degrees 10 and 9, the final subresultants have degrees 2, 1, 0.
    # A gcd of degree at least three forces every coefficient of the
    # quadratic and linear subresultants to vanish.
    subresultants = sp.subresultants(
        normalized_ten,
        sp.diff(normalized_ten, y),
        y,
    )
    require(
        [sp.degree(item, y) for item in subresultants]
        == list(range(10, -1, -1)),
        "subresultant degree profile changed",
    )
    quadratic_coefficients = sp.Poly(subresultants[-3], y).all_coeffs()
    linear_coefficients = sp.Poly(subresultants[-2], y).all_coeffs()
    require(
        len(quadratic_coefficients) == 3 and len(linear_coefficients) == 2,
        "terminal subresultant shapes changed",
    )

    # Strip the exact powers of the already-known deep factor.  Away from
    # deep_linear=0 this does not change the common zero set.
    residual_coefficients = [
        primitive_after_division(
            linear_coefficients[0],
            deep_linear**3,
            (cvar, wvar),
        ),
        primitive_after_division(
            linear_coefficients[1],
            deep_linear**4,
            (cvar, wvar),
        ),
        primitive_after_division(
            quadratic_coefficients[0],
            deep_linear,
            (cvar, wvar),
        ),
        primitive_after_division(
            quadratic_coefficients[1],
            deep_linear**2,
            (cvar, wvar),
        ),
        primitive_after_division(
            quadratic_coefficients[2],
            deep_linear**2,
            (cvar, wvar),
        ),
    ]

    basis = sp.groebner(
        residual_coefficients,
        wvar,
        cvar,
        order="grevlex",
        method="f5b",
    )
    basis_polynomial = (
        (15309 * cvar**2 + 250)
        * (30618 * cvar**2 + 361) ** 3
        * (
            199644669 * cvar**4
            + 7654500 * cvar**2
            + 62500
        )
    )
    primitive_basis = [
        sp.Poly(item.as_expr(), wvar, cvar).primitive()[1].as_expr()
        for item in basis.polys
    ]
    require(
        len(primitive_basis) == 2
        and sp.expand(primitive_basis[0] - basis_polynomial) == 0
        and sp.expand(primitive_basis[1] - deep_linear) == 0,
        "saturated coefficient Groebner basis changed",
    )
    _, deep_remainder = basis.reduce(deep_linear)
    require(
        deep_remainder == 0,
        "deep linear form stopped belonging to the residual ideal",
    )

    # Positive and hostile controls.  A generic deep point has precisely the
    # inherited y^4 collision (gcd degree 3), whereas two off-deep points are
    # squarefree.
    def gcd_degree(cvalue: sp.Expr, wvalue: sp.Expr) -> int:
        specialization = sp.Poly(
            normalized_ten.subs({cvar: cvalue, wvar: wvalue}),
            y,
            domain=sp.QQ,
        )
        return sp.gcd(specialization, specialization.diff()).degree()

    deep_gcd_degree = gcd_degree(1, -sp.Rational(20, 21))
    off_deep_gcd_degrees = (gcd_degree(1, 0), gcd_degree(0, 1))
    require(
        deep_gcd_degree == 3 and off_deep_gcd_degrees == (0, 0),
        "deep/off-deep gcd controls changed",
    )

    print("THM-2345 exact common-root-wall saturation referee")
    print("wall factorization: F=1323*y^2*G10, deg_y(G10)=10")
    print("deep factorization: G10|_(20C+21W=0)=7*y^4*G6")
    print("subresultant degrees: 10,9,8,7,6,5,4,3,2,1,0")
    print("stripped terminal coefficient count: 5")
    print(
        "Groebner basis:"
        " [(15309*C^2+250)*(30618*C^2+361)^3"
        "*(199644669*C^4+7654500*C^2+62500), 20*C+21*W]"
    )
    print("ideal membership: 20*C+21*W has zero remainder")
    print(f"deep control gcd degree: {deep_gcd_degree}")
    print(f"off-deep control gcd degrees: {off_deep_gcd_degrees}")
    print("VERDICT: gcd degree >=3 forces 20*C+21*W=0 on B=1")


if __name__ == "__main__":
    main()
