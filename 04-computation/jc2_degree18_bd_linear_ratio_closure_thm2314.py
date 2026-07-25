#!/usr/bin/env python3
"""Exact referee for the two rational B--D points in THM-2311.

On the B--D weighted line, put B=1 and C=W=0.  This script verifies the
complete branch-discriminant factorizations at

    D = 4075/85176,   D = 25/126,

proves the residual factors squarefree, checks the geometry of the exceptional
fibres, proves absolute irreducibility by an exact polynomial-root Groebner
test, and freezes the ramification counts used in THM-2314.

The Riemann--Hurwitz and rational-map arguments remain mathematical steps in
the theorem text; the script certifies their exact algebraic inputs.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def spectral_curve(
    u: sp.Symbol,
    y: sp.Symbol,
    dvar: sp.Expr,
) -> sp.Expr:
    """THM-2297's G_0 at (B,C,D,W)=(1,0,dvar,0)."""
    return sp.expand(
        -26040609 * u**3
        + (49601160 + 1607445 * y**2) * u**2
        + (
            -20995200
            - 2857680 * y**2
            - 52907904 * dvar
            - 138915 * y**4
        )
        * u
        + 777600 * y**2
        + 33592320 * dvar
        + 78120 * y**4
        + 1959552 * dvar * y**2
        + 1127 * y**6
    )


def even_root_groebner(
    curve: sp.Expr,
    u: sp.Symbol,
    y: sp.Symbol,
) -> sp.GroebnerBasis:
    """Test for a polynomial root u=a*y^2+b over the algebraic closure."""
    a, b, x = sp.symbols("a b x")
    in_x = sp.expand(curve.subs({u: a * x + b, y**2: x}))
    # Simultaneous substitution of y**2 is safest through Poly's even terms.
    even_curve = sp.Poly(curve, y)
    require(
        all(exponent[0] % 2 == 0 for exponent, _ in even_curve.terms()),
        "B--D curve unexpectedly has an odd y term",
    )
    in_x = sp.expand(
        sum(
            coefficient * x ** (exponent[0] // 2)
            for exponent, coefficient in sp.Poly(
                curve.subs(u, a * y**2 + b),
                y,
            ).terms()
        )
    )
    coefficients = [
        sp.Poly(in_x, x).coeff_monomial(x**index) for index in range(4)
    ]
    return sp.groebner(coefficients, a, b, order="lex", domain=sp.QQ)


def main() -> None:
    u, y, x, v, r = sp.symbols("u y x v r")

    t_ramified = sp.Rational(4075, 85176)
    t_triple = sp.Rational(25, 126)
    g_ramified = spectral_curve(u, y, t_ramified)
    g_triple = spectral_curve(u, y, t_triple)

    delta_ramified = sp.factor(sp.discriminant(g_ramified, u))
    delta_triple = sp.factor(sp.discriminant(g_triple, u))

    q2 = 91 * y**2 + 1215
    h8 = (
        1577224103 * y**8
        + 203464170 * y**6
        - 147517112925 * y**4
        + 1389005982000 * y**2
        - 3037628182500
    )
    h6 = (
        7889 * y**6
        + 211680 * y**4
        + 1814400 * y**2
        + 2916000
    )
    expected_delta_ramified = (
        -sp.Rational(56684737689882624, 4826809) * q2**2 * h8
    )
    expected_delta_triple = -19442865027629740032 * y**6 * h6
    require(
        sp.expand(delta_ramified - expected_delta_ramified) == 0,
        "4075/85176 branch-discriminant factorization mismatch",
    )
    require(
        sp.expand(delta_triple - expected_delta_triple) == 0,
        "25/126 branch-discriminant factorization mismatch",
    )

    require(
        sp.gcd(sp.Poly(h8, y), sp.Poly(sp.diff(h8, y), y)).degree() == 0,
        "h8 is not squarefree",
    )
    require(
        sp.gcd(sp.Poly(h8, y), sp.Poly(q2, y)).degree() == 0,
        "h8 meets the triple-ramification factor",
    )
    require(
        sp.gcd(sp.Poly(h6, y), sp.Poly(sp.diff(h6, y), y)).degree() == 0,
        "h6 is not squarefree",
    )
    require(h6.subs(y, 0) != 0, "h6 meets the ordinary triple point")

    # Every member has the same three distinct, unramified weighted branches
    # at infinity: put u=v*y^2 and r=1/y.
    infinity_polynomial = (
        1127 - 138915 * v + 1607445 * v**2 - 26040609 * v**3
    )
    infinity_discriminant = sp.discriminant(infinity_polynomial, v)
    require(
        infinity_discriminant
        == -153384762202971019112448,
        "weighted infinity cubic is not separable",
    )
    for curve in (g_ramified, g_triple):
        infinity_chart = sp.expand(
            r**6 * curve.subs({u: v / r**2, y: 1 / r})
        )
        require(
            sp.expand(infinity_chart.subs(r, 0) - infinity_polynomial) == 0,
            "weighted infinity chart mismatch",
        )

    # At q2=0 the cubic has a triple u-root, but the total curve is smooth:
    # H_x != 0 for x=y^2.  Thus y-y0 has order three on the unique local
    # branch and each of the two q2-roots contributes e-1=2.
    h_ramified = sp.expand(g_ramified.subs(y**2, x))
    x0 = -sp.Rational(1215, 91)
    u0 = sp.Rational(295, 819)
    require(
        sp.factor(h_ramified.subs(x, x0))
        == -sp.Rational(729, 15379) * (819 * u - 295) ** 3,
        "triple ramification fibre mismatch",
    )
    hx_at_ramified = sp.factor(
        sp.diff(h_ramified, x).subs({x: x0, u: u0})
    )
    require(
        hx_at_ramified == -sp.Rational(16329600, 169),
        "triple ramification fibre is not smooth",
    )

    # At D=25/126 the entire y=0 fibre coalesces.  Shift U=u-40/63.
    # Its degree-three tangent cone has three distinct nonvertical lines, so
    # this is an ordinary triple point (delta=3), with three normalization
    # branches all unramified for projection to y.
    capital_u = sp.symbols("U")
    shifted_triple = sp.expand(g_triple.subs(u, capital_u + sp.Rational(40, 63)))
    expected_shifted = 7 * (
        -3720087 * capital_u**3
        + 229635 * capital_u**2 * y**2
        - 19845 * capital_u * y**4
        - 116640 * capital_u * y**2
        + 161 * y**6
        - 1440 * y**4
    )
    require(
        sp.expand(shifted_triple - expected_shifted) == 0,
        "ordinary triple-point chart mismatch",
    )
    tangent_cone = -5103 * capital_u * (5103 * capital_u**2 + 160 * y**2)
    actual_tangent = sum(
        coefficient * capital_u ** exponent[0] * y ** exponent[1]
        for exponent, coefficient in sp.Poly(
            shifted_triple,
            capital_u,
            y,
        ).terms()
        if sum(exponent) == 3
    )
    require(
        sp.expand(actual_tangent - tangent_cone) == 0,
        "ordinary triple-point tangent cone mismatch",
    )
    require(
        sp.discriminant(5103 * capital_u**2 + 160, capital_u) != 0,
        "ordinary triple-point tangent directions coalesce",
    )

    # Absolute irreducibility.  A polynomial u-root has y-degree <=2.  Since
    # the curve is even in y, any root produces an even root a*y^2+b (either
    # directly, or as the third root beside r(y),r(-y)).  The exact coefficient
    # ideals are the unit ideal at both points.
    groebner_ramified = even_root_groebner(g_ramified, u, y)
    groebner_triple = even_root_groebner(g_triple, u, y)
    require(
        len(groebner_ramified.polys) == 1
        and groebner_ramified.polys[0].as_expr() == 1,
        "4075/85176 curve has a polynomial u-root",
    )
    require(
        len(groebner_triple.polys) == 1
        and groebner_triple.polys[0].as_expr() == 1,
        "25/126 curve has a polynomial u-root",
    )

    # Hostile control: away from the THM-2311 bank, D=1 has no repeated
    # branch value.  This detects an accidentally generic factorization.
    g_control = spectral_curve(u, y, sp.Integer(1))
    delta_control = sp.Poly(sp.discriminant(g_control, u), y)
    require(
        sp.gcd(delta_control, delta_control.diff()).degree() == 0,
        "D=1 hostile control unexpectedly has repeated branch",
    )

    print("normalization=B=1,C=W=0")
    print("ratio_1=D=4075/85176")
    print(f"ratio_1_Delta={expected_delta_ramified}")
    print("ratio_1_h8_squarefree=PASS")
    print("ratio_1_special_fibres=2 smooth total-ramification points, e=3")
    print("ratio_1_simple_branch_points=8")
    print("ratio_1_total_ramification=12")
    print("ratio_1_normalization_genus=4")
    print("ratio_2=D=25/126")
    print(f"ratio_2_Delta={expected_delta_triple}")
    print("ratio_2_h6_squarefree=PASS")
    print(f"ratio_2_shifted_curve={expected_shifted}")
    print(f"ratio_2_tangent_cone={tangent_cone}")
    print("ratio_2_origin=ordinary triple point, delta=3, 3 unramified branches")
    print("ratio_2_simple_branch_points=6")
    print("ratio_2_total_ramification=6")
    print("ratio_2_normalization_genus=1")
    print(f"infinity_cubic={infinity_polynomial}")
    print(f"infinity_discriminant={infinity_discriminant}")
    print("infinity=3 distinct unramified points")
    print("absolute_irreducibility=PASS Groebner bases [1],[1]")
    print("hostile_control_D=1 squarefree=PASS")
    print("riemann_hurwitz_and_deck_contradiction=MATHEMATICAL_PROOF_REQUIRED")
    print("status=THM2314_BD_LINEAR_RATIO_EXACT_REFEREE")


if __name__ == "__main__":
    main()
