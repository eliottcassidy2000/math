#!/usr/bin/env python3
"""Exact referee for THM-2316's full C--D ratio-bank closure.

THM-2311 leaves the rational ratio

    t = D^3 / C^4 = -22143375 / 6397664.

The rational representative C=D=1/t treats the linear point over Q.  The
other three ratios form one irreducible cubic orbit; the representative
(C,D)=(alpha^2,alpha^3) treats them uniformly over Q(alpha).  This script
verifies both branch factorizations, all local exceptional fibres, separable
infinity, and absolute irreducibility by exact polynomial-root Groebner
certificates.

Riemann--Hurwitz and the inherited deck contradiction remain mathematical
steps in the theorem text.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def spectral_curve(
    u: sp.Symbol,
    y: sp.Symbol,
    cvar: sp.Expr,
    dvar: sp.Expr,
) -> sp.Expr:
    """THM-2297's G_0 at (B,C,D,W)=(0,cvar,dvar,0)."""
    return sp.expand(
        -26040609 * u**3
        + 1607445 * u**2 * y**2
        - 52907904 * dvar * u
        - 138915 * u * y**4
        + 1959552 * dvar * y**2
        - 435456 * cvar * y**3
        + 1127 * y**6
    )


def polynomial_root_groebner(
    curve: sp.Expr,
    u: sp.Symbol,
    y: sp.Symbol,
    extension: sp.Expr | None = None,
) -> sp.GroebnerBasis:
    """Test for every possible polynomial root u=a*y^2+b*y+c."""
    a, b, c = sp.symbols("a b c")
    substituted = sp.Poly(sp.expand(curve.subs(u, a * y**2 + b * y + c)), y)
    coefficients = [
        substituted.coeff_monomial(y**degree) for degree in range(7)
    ]
    if extension is None:
        return sp.groebner(
            coefficients,
            a,
            b,
            c,
            order="grevlex",
            domain=sp.QQ,
        )
    return sp.groebner(
        coefficients,
        a,
        b,
        c,
        order="grevlex",
        extension=extension,
    )


def field_reduce(expression: sp.Expr, field: sp.AlgebraicField) -> sp.Expr:
    """Return the canonical power-basis representative in an algebraic field."""
    return sp.factor(field.to_sympy(field.from_sympy(sp.cancel(expression))))


def main() -> None:
    u, y, v, z, t = sp.symbols("u y v z t")

    ratio = -sp.Rational(22143375, 6397664)
    representative = 1 / ratio
    require(
        representative**3 / representative**4 == ratio,
        "C=D=1/t does not represent the weighted ratio",
    )

    curve_rational = spectral_curve(u, y, representative, representative)
    scaled_curve = sp.expand(sp.Rational(91125, 49) * curve_rational)
    expected_scaled_curve = (
        -48427561125 * u**3
        + 2989355625 * u**2 * y**2
        - 258339375 * u * y**4
        + 28427563008 * u
        + 2095875 * y**6
        + 233971712 * y**3
        - 1052872704 * y**2
    )
    require(
        scaled_curve == expected_scaled_curve,
        "integer-scaled spectral curve changed",
    )

    exceptional = 15 * y - 52
    h10 = (
        2864799140625 * y**10
        + 19862607375000 * y**9
        + 103285558350000 * y**8
        + 367474086720000 * y**7
        + 527020367328000 * y**6
        - 762214033152000 * y**5
        - 6205500108288000 * y**4
        - 55515871936512000 * y**3
        - 144341267034931200 * y**2
        - 333588706036285440 * y
        - 578220423796228096
    )
    expected_discriminant = (
        -sp.Rational(464767295827968, 1953125) * exceptional**2 * h10
    )
    actual_discriminant = sp.factor(sp.discriminant(curve_rational, u))
    require(
        sp.expand(actual_discriminant - expected_discriminant) == 0,
        "C--D branch-discriminant factorization mismatch",
    )
    h10_poly = sp.Poly(h10, y)
    require(
        sp.gcd(h10_poly, h10_poly.diff()).degree() == 0,
        "degree-ten branch factor is not squarefree",
    )
    require(
        sp.gcd(h10_poly, sp.Poly(exceptional, y)).degree() == 0,
        "simple branch factor meets exceptional fibre",
    )

    # The squared linear discriminant factor is one smooth e=3 fibre.
    y0 = sp.Rational(52, 15)
    u0 = sp.Rational(2704, 10935)
    require(
        sp.factor(curve_rational.subs(y, y0))
        == -sp.Rational(49, 2460375) * (10935 * u - 2704) ** 3,
        "exceptional fibre is not a triple root",
    )
    require(
        sp.factor(sp.diff(curve_rational, y).subs({u: u0, y: y0}))
        == -sp.Rational(2384639688704, 2278125),
        "exceptional triple fibre is not smooth",
    )

    # The three weighted branches at infinity are distinct and unramified.
    infinity_polynomial = (
        1127 - 138915 * v + 1607445 * v**2 - 26040609 * v**3
    )
    infinity_discriminant = sp.discriminant(infinity_polynomial, v)
    require(
        infinity_discriminant == -153384762202971019112448,
        "weighted infinity cubic is not separable",
    )
    infinity_chart = sp.expand(
        z**6 * curve_rational.subs({u: v / z**2, y: 1 / z})
    )
    require(
        sp.expand(infinity_chart.subs(z, 0) - infinity_polynomial) == 0,
        "weighted infinity chart changed",
    )

    # A reducible cubic with constant leading u-coefficient has a polynomial
    # root.  Its y-degree is at most two by the leading-degree comparison.
    # The coefficient ideal for u=a*y^2+b*y+c is the unit ideal.
    groebner = polynomial_root_groebner(curve_rational, u, y)
    require(
        len(groebner.polys) == 1 and groebner.polys[0].as_expr() == 1,
        "spectral cubic has a polynomial u-root",
    )

    # Away from the THM-2311 ratio bank, C=D=1 has a squarefree branch
    # discriminant.  This detects an accidentally generic specialization.
    control = spectral_curve(u, y, sp.Integer(1), sp.Integer(1))
    control_discriminant = sp.Poly(sp.discriminant(control, u), y)
    require(
        sp.gcd(control_discriminant, control_discriminant.diff()).degree() == 0,
        "C=D=1 hostile control unexpectedly has a repeated branch value",
    )

    rational_simple_ramification = h10_poly.degree()
    rational_total_ramification = rational_simple_ramification + 2
    rational_genus = (rational_total_ramification - 4) // 2
    require(rational_simple_ramification == 10, "simple branch count changed")
    require(rational_total_ramification == 12, "total ramification changed")
    require(rational_genus == 4, "rational normalization genus changed")

    # The remaining C--D ratios form one irreducible cubic orbit.
    ratio_cubic = (
        16544432128 * t**3
        + 54880100352 * t**2
        - 8964338040 * t
        + 387420489
    )
    ratio_mod_5 = sp.Poly(ratio_cubic, t, modulus=5).monic()
    require(
        ratio_mod_5 == sp.Poly(t**3 - t**2 - 2, t, modulus=5),
        "cubic ratio polynomial has the wrong mod-5 reduction",
    )
    require(
        all(ratio_mod_5.eval(residue) != 0 for residue in range(5)),
        "cubic ratio polynomial gained a root modulo five",
    )
    require(
        sp.discriminant(ratio_cubic, t)
        == -24314073653663036249383152766404532371456,
        "cubic ratio discriminant changed",
    )

    alpha = sp.CRootOf(ratio_cubic, 0)
    field = sp.QQ.algebraic_field(alpha)
    curve_cubic = spectral_curve(u, y, alpha**2, alpha**3)
    require(
        field_reduce((alpha**3) ** 3 / (alpha**2) ** 4 - alpha, field) == 0,
        "cubic-field representative has the wrong weighted ratio",
    )

    node_y = (
        2130734520 * alpha
        - 175906971
        + 557240320 * alpha**2
    ) / 71030268
    node_u = (
        11851365288 * alpha
        - 1192022163
        + 5827039232 * alpha**2
    ) / 12375051136
    simple_u = -(
        251919849384 * alpha
        - 20715038739
        + 34173669376 * alpha**2
    ) / 111375460224

    delta_cubic = sp.Poly(
        sp.discriminant(curve_cubic, u),
        y,
        extension=alpha,
    )
    node_square = sp.Poly((y - node_y) ** 2, y, extension=alpha)
    cubic_quotient, cubic_remainder = delta_cubic.div(node_square)
    require(
        cubic_remainder.is_zero
        and field_reduce(cubic_quotient.LC(), field)
        == -153384762202971019112448,
        "cubic-field branch-discriminant factorization mismatch",
    )
    cubic_h10 = cubic_quotient.monic()
    require(cubic_h10.degree() == 10, "cubic-field residual degree changed")
    require(
        sp.gcd(cubic_h10, cubic_h10.diff()).degree() == 0,
        "cubic-field residual h10 is not squarefree",
    )
    require(
        sp.gcd(
            cubic_h10,
            sp.Poly(y - node_y, y, extension=alpha),
        ).degree()
        == 0,
        "cubic-field residual h10 meets the node value",
    )

    # At the repeated value the fibre has one simple and one double u-root.
    # The double point is an ordinary node with two nonvertical tangents.
    node_fibre_check = sp.Poly(
        sp.expand(
            curve_cubic.subs(y, node_y)
            + 26040609 * (u - node_u) ** 2 * (u - simple_u)
        ),
        u,
        extension=alpha,
    )
    require(node_fibre_check.is_zero, "cubic-field node fibre changed")
    require(
        field_reduce(curve_cubic.subs({u: node_u, y: node_y}), field) == 0,
        "cubic-field node left the curve",
    )
    require(
        field_reduce(
            sp.diff(curve_cubic, u).subs({u: node_u, y: node_y}),
            field,
        )
        == 0,
        "cubic-field node has nonzero u-derivative",
    )
    require(
        field_reduce(
            sp.diff(curve_cubic, y).subs({u: node_u, y: node_y}),
            field,
        )
        == 0,
        "cubic-field node has nonzero y-derivative",
    )
    tangent_u2 = field_reduce(
        sp.diff(curve_cubic, u, 2).subs({u: node_u, y: node_y}) / 2,
        field,
    )
    tangent_uy = field_reduce(
        sp.diff(curve_cubic, u, y).subs({u: node_u, y: node_y}),
        field,
    )
    tangent_y2 = field_reduce(
        sp.diff(curve_cubic, y, 2).subs({u: node_u, y: node_y}) / 2,
        field,
    )
    expected_tangent_u2 = -59049 * (
        179291068488 * alpha
        - 15721619103
        + 43308511232 * alpha**2
    ) / 126276032
    expected_tangent_discriminant = 31381059609 * (
        23902870780036903800843
        - 559844552503271371280616 * alpha
        + 3545334969099076026159104 * alpha**2
    ) / 2466893777661329408
    tangent_discriminant = field_reduce(
        tangent_uy**2 - 4 * tangent_u2 * tangent_y2,
        field,
    )
    require(
        field_reduce(tangent_u2 - expected_tangent_u2, field) == 0
        and tangent_u2 != 0,
        "cubic-field node gained a vertical tangent",
    )
    require(
        field_reduce(
            tangent_discriminant - expected_tangent_discriminant,
            field,
        )
        == 0
        and tangent_discriminant != 0,
        "cubic-field node tangent cone is not ordinary",
    )

    cubic_infinity_chart = sp.expand(
        z**6 * curve_cubic.subs({u: v / z**2, y: 1 / z})
    )
    require(
        sp.expand(
            cubic_infinity_chart.subs(z, 0) - infinity_polynomial
        )
        == 0,
        "cubic-field infinity chart changed",
    )
    cubic_groebner = polynomial_root_groebner(
        curve_cubic,
        u,
        y,
        extension=alpha,
    )
    require(
        len(cubic_groebner.polys) == 1
        and cubic_groebner.polys[0].as_expr() == 1,
        "cubic-field spectral curve has a polynomial u-root",
    )
    cubic_total_ramification = cubic_h10.degree()
    cubic_genus = (cubic_total_ramification - 4) // 2
    require(cubic_total_ramification == 10, "cubic ramification changed")
    require(cubic_genus == 3, "cubic-field normalization genus changed")

    print("weighted_ratio=D^3/C^4=-22143375/6397664")
    print("rational_representative=C=D=-6397664/22143375")
    print(f"integer_scaled_curve={expected_scaled_curve}")
    print(f"branch_discriminant={expected_discriminant}")
    print("degree_10_factor=squarefree_and_coprime_to_15y-52")
    print("simple_branch_points=10")
    print("exceptional_fibre=y=52/15,u=2704/10935,smooth_total_ramification_e3")
    print(f"infinity_cubic={infinity_polynomial}")
    print(f"infinity_discriminant={infinity_discriminant}")
    print("infinity=3_distinct_unramified_points")
    print("absolute_irreducibility=PASS_Groebner_basis_[1]")
    print("rational_total_ramification=12")
    print("rational_normalization_genus=4")
    print("cubic_ratio_field=minpoly_degree_3_irreducible_mod_5")
    print("cubic_orbit_size=3")
    print("cubic_Delta=universal_lead*(y-Y(alpha))^2*h10")
    print("cubic_h10=degree_10_squarefree_coprime")
    print("cubic_exceptional_fibre=ordinary_node_two_unramified_branches")
    print("cubic_absolute_irreducibility=PASS_Groebner_basis_[1]")
    print("cubic_total_ramification=10")
    print("cubic_normalization_genus=3")
    print("closed_CD_ratios=4")
    print("hostile_control_C=D=1_squarefree=PASS")
    print("riemann_hurwitz_and_deck_contradiction=MATHEMATICAL_PROOF_REQUIRED")
    print("status=THM2316_FULL_CD_BANK_EXACT_REFEREE")


if __name__ == "__main__":
    main()
