#!/usr/bin/env python3
"""Exact referee for THM-2320's full D--W ratio-bank closure.

The rational linear ratio is treated over Q with D=W=1/t.  The other three
ratios form one irreducible cubic orbit; the representative
(D,W)=(alpha^3,alpha^4) treats all conjugates uniformly.  Every check is
exact and remains active under optimized Python.
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
    wvar: sp.Expr,
) -> sp.Expr:
    """THM-2297's G_0 at (B,C,D,W)=(0,0,dvar,wvar)."""
    return sp.expand(
        -5878656 * wvar * y
        - 26040609 * u**3
        + 1607445 * u**2 * y**2
        - 52907904 * dvar * u
        - 138915 * u * y**4
        + 1959552 * dvar * y**2
        + 1127 * y**6
    )


def field_reduce(expression: sp.Expr, field: sp.AlgebraicField) -> sp.Expr:
    """Return the canonical power-basis representative in an algebraic field."""
    return sp.factor(field.to_sympy(field.from_sympy(sp.cancel(expression))))


def polynomial_root_groebner(
    curve: sp.Expr,
    u: sp.Symbol,
    y: sp.Symbol,
) -> sp.GroebnerBasis:
    """Test every possible polynomial root u=a*y^2+b*y+c over Q."""
    a, b, c = sp.symbols("a b c")
    substituted = sp.Poly(
        sp.expand(curve.subs(u, a * y**2 + b * y + c)),
        y,
    )
    coefficients = [
        substituted.coeff_monomial(y**degree) for degree in range(7)
    ]
    return sp.groebner(
        coefficients,
        a,
        b,
        c,
        order="grevlex",
        domain=sp.QQ,
    )


def main() -> None:
    u, y, v, z, t = sp.symbols("u y v z t")

    rational_ratio = -sp.Rational(935886848, 430565625)
    rational_representative = 1 / rational_ratio
    require(
        rational_representative**4 / rational_representative**5
        == rational_ratio,
        "D=W=1/t has the wrong weighted ratio",
    )
    rational_curve = spectral_curve(
        u,
        y,
        rational_representative,
        rational_representative,
    )
    scaled_curve = sp.expand(sp.Rational(7311616, 49) * rational_curve)
    expected_scaled_curve = (
        -3885692518656 * u**3
        + 239857562880 * u**2 * y**2
        - 20728431360 * u * y**4
        + 3632067084375 * u
        + 168167168 * y**6
        - 134521003125 * y**2
        + 403563009375 * y
    )
    require(scaled_curve == expected_scaled_curve, "scaled rational curve changed")

    exceptional = 104 * y - 405
    h10 = (
        851140463828047757312 * y**10
        + 6629074766353064263680 * y**9
        + 38722720389995062886400 * y**8
        + 201060278948051288064000 * y**7
        + 609926274714380083200000 * y**6
        + 999177077047752345600000 * y**5
        - 1467500022145297296000000 * y**4
        - 26582121110232749880000000 * y**3
        - 59732448939525919050000000 * y**2
        - 282296044395397652343750000 * y
        - 549662971058346390380859375
    )
    expected_rational_discriminant = (
        -sp.Rational(
            1628150074335205281,
            97719251621562548224,
        )
        * exceptional**2
        * h10
    )
    require(
        sp.expand(
            sp.discriminant(rational_curve, u)
            - expected_rational_discriminant
        )
        == 0,
        "rational branch-discriminant factorization changed",
    )
    h10_poly = sp.Poly(h10, y)
    require(
        sp.gcd(h10_poly, h10_poly.diff()).degree() == 0,
        "rational h10 is not squarefree",
    )
    require(
        sp.gcd(h10_poly, sp.Poly(exceptional, y)).degree() == 0,
        "rational h10 meets the exceptional fibre",
    )

    rational_y0 = sp.Rational(405, 104)
    rational_u0 = sp.Rational(3375, 10816)
    require(
        sp.factor(rational_curve.subs(y, rational_y0))
        == -sp.Rational(26040609, 1265319018496)
        * (10816 * u - 3375) ** 3,
        "rational exceptional fibre is not a triple root",
    )
    require(
        sp.factor(
            sp.diff(rational_curve, y).subs(
                {u: rational_u0, y: rational_y0}
            )
        )
        == -sp.Rational(692110561078125, 95051008),
        "rational exceptional fibre is not smooth",
    )
    rational_groebner = polynomial_root_groebner(rational_curve, u, y)
    require(
        len(rational_groebner.polys) == 1
        and rational_groebner.polys[0].as_expr() == 1,
        "rational spectral curve has a polynomial u-root",
    )

    infinity_polynomial = (
        1127 - 138915 * v + 1607445 * v**2 - 26040609 * v**3
    )
    infinity_discriminant = sp.discriminant(infinity_polynomial, v)
    require(
        infinity_discriminant == -153384762202971019112448,
        "weighted infinity cubic is not separable",
    )
    rational_infinity = sp.expand(
        z**6 * rational_curve.subs({u: v / z**2, y: 1 / z})
    )
    require(
        sp.expand(rational_infinity.subs(z, 0) - infinity_polynomial) == 0,
        "rational infinity chart changed",
    )

    rational_total_ramification = 10 + 2
    rational_genus = (rational_total_ramification - 4) // 2
    require(rational_genus == 4, "rational normalization genus changed")

    # The other three ratios are one irreducible cubic orbit.
    ratio_cubic = (
        56162900390625 * t**3
        - 1448500838400000 * t**2
        + 17932072576352256 * t
        + 36028797018963968
    )
    ratio_mod_13 = sp.Poly(ratio_cubic, t, modulus=13).monic()
    require(
        ratio_mod_13
        == sp.Poly(t**3 - 4 * t**2 + 5 * t - 1, t, modulus=13),
        "cubic ratio polynomial has the wrong mod-13 reduction",
    )
    require(
        all(ratio_mod_13.eval(residue) != 0 for residue in range(13)),
        "cubic ratio polynomial gained a root modulo thirteen",
    )
    require(
        sp.discriminant(ratio_cubic, t)
        == -1239334538118889697076867375203545322038089020866560000000000000,
        "cubic ratio discriminant changed",
    )

    alpha = sp.CRootOf(ratio_cubic, 0)
    field = sp.QQ.algebraic_field(alpha)
    cubic_curve = spectral_curve(u, y, alpha**3, alpha**4)
    require(
        field_reduce((alpha**4) ** 4 / (alpha**3) ** 5 - alpha, field) == 0,
        "cubic representative has the wrong weighted ratio",
    )

    cubic_discriminant = sp.Poly(
        sp.discriminant(cubic_curve, u),
        y,
        extension=alpha,
    )
    repeated_linear = sp.gcd(
        cubic_discriminant,
        cubic_discriminant.diff(),
    ).monic()
    require(repeated_linear.degree() == 1, "cubic repeated factor is not linear")
    node_y = field_reduce(-repeated_linear.coeff_monomial(1), field)
    expected_node_y = (
        3790040625 * alpha**2
        - 25880297472 * alpha
        - 343597383680
    ) / 42471522304
    require(
        field_reduce(node_y - expected_node_y, field) == 0,
        "cubic node y-coordinate changed",
    )
    quotient, remainder = cubic_discriminant.div(
        sp.Poly((y - node_y) ** 2, y, extension=alpha)
    )
    require(
        remainder.is_zero
        and field_reduce(quotient.LC(), field)
        == -153384762202971019112448,
        "cubic branch-discriminant factorization changed",
    )
    cubic_h10 = quotient.monic()
    require(cubic_h10.degree() == 10, "cubic residual degree changed")
    require(
        sp.gcd(cubic_h10, cubic_h10.diff()).degree() == 0,
        "cubic h10 is not squarefree",
    )
    require(
        sp.gcd(
            cubic_h10,
            sp.Poly(y - node_y, y, extension=alpha),
        ).degree()
        == 0,
        "cubic h10 meets the node value",
    )

    node_fibre = sp.Poly(cubic_curve.subs(y, node_y), u, extension=alpha)
    repeated_u = sp.gcd(node_fibre, node_fibre.diff()).monic()
    require(repeated_u.degree() == 1, "cubic exceptional root is not double")
    node_u = field_reduce(-repeated_u.coeff_monomial(1), field)
    expected_node_u = (
        -7058559375 * alpha**2
        + 104983822336 * alpha
        + 206158430208
    ) / 25517520000
    require(
        field_reduce(node_u - expected_node_u, field) == 0,
        "cubic node u-coordinate changed",
    )
    simple_u = field_reduce(sp.Rational(5, 81) * node_y**2 - 2 * node_u, field)
    require(
        sp.Poly(
            sp.expand(
                cubic_curve.subs(y, node_y)
                + 26040609 * (u - node_u) ** 2 * (u - simple_u)
            ),
            u,
            extension=alpha,
        ).is_zero,
        "cubic exceptional fibre factorization changed",
    )
    require(
        field_reduce(cubic_curve.subs({u: node_u, y: node_y}), field) == 0,
        "node left the cubic curve",
    )
    require(
        field_reduce(
            sp.diff(cubic_curve, u).subs({u: node_u, y: node_y}),
            field,
        )
        == 0,
        "node has nonzero u-derivative",
    )
    require(
        field_reduce(
            sp.diff(cubic_curve, y).subs({u: node_u, y: node_y}),
            field,
        )
        == 0,
        "node has nonzero y-derivative",
    )
    tangent_u2 = field_reduce(
        sp.diff(cubic_curve, u, 2).subs({u: node_u, y: node_y}) / 2,
        field,
    )
    tangent_uy = field_reduce(
        sp.diff(cubic_curve, u, y).subs({u: node_u, y: node_y}),
        field,
    )
    tangent_y2 = field_reduce(
        sp.diff(cubic_curve, y, 2).subs({u: node_u, y: node_y}) / 2,
        field,
    )
    tangent_discriminant = field_reduce(
        tangent_uy**2 - 4 * tangent_u2 * tangent_y2,
        field,
    )
    expected_tangent_u2 = (
        6561
        * (
            611312821875 * alpha**2
            - 11216021225472 * alpha
            - 19310172962816
        )
        / 202520000
    )
    expected_tangent_discriminant = (
        -sp.Rational(
            945539748965690376192,
            57572832660675048828125,
        )
        * (
            136966554945917353125 * alpha**2
            - 4257922067489564393472 * alpha
            - 7854173444688698146816
        )
    )
    require(
        field_reduce(tangent_u2 - expected_tangent_u2, field) == 0,
        "cubic tangent u^2 coefficient changed",
    )
    require(
        field_reduce(
            tangent_discriminant - expected_tangent_discriminant,
            field,
        )
        == 0,
        "cubic tangent discriminant changed",
    )
    require(tangent_u2 != 0, "cubic node gained a vertical tangent")
    require(tangent_discriminant != 0, "cubic tangent cone is not ordinary")

    # A global coefficient ideal proves irreducibility for all conjugates at
    # once, before adjoining a chosen root of the ratio polynomial.
    a, b, c = sp.symbols("a b c")
    universal_curve = spectral_curve(u, y, t**3, t**4)
    substituted = sp.Poly(
        sp.expand(universal_curve.subs(u, a * y**2 + b * y + c)),
        y,
    )
    root_equations = [
        substituted.coeff_monomial(y**degree) for degree in range(7)
    ]
    global_groebner = sp.groebner(
        root_equations + [ratio_cubic],
        a,
        b,
        c,
        t,
        order="grevlex",
        domain=sp.QQ,
    )
    require(
        len(global_groebner.polys) == 1
        and global_groebner.polys[0].as_expr() == 1,
        "cubic orbit has a polynomial u-root",
    )

    cubic_infinity = sp.expand(
        z**6 * cubic_curve.subs({u: v / z**2, y: 1 / z})
    )
    require(
        sp.expand(cubic_infinity.subs(z, 0) - infinity_polynomial) == 0,
        "cubic infinity chart changed",
    )
    cubic_total_ramification = cubic_h10.degree()
    cubic_genus = (cubic_total_ramification - 4) // 2
    require(cubic_total_ramification == 10, "cubic ramification changed")
    require(cubic_genus == 3, "cubic normalization genus changed")

    # Generic off-bank control.
    control_curve = spectral_curve(u, y, sp.Integer(1), sp.Integer(1))
    control_discriminant = sp.Poly(sp.discriminant(control_curve, u), y)
    require(
        sp.gcd(control_discriminant, control_discriminant.diff()).degree() == 0,
        "D=W=1 hostile control unexpectedly has repeated branch",
    )

    print("rational_ratio=W^4/D^5=-935886848/430565625")
    print("rational_representative=D=W=-430565625/935886848")
    print(f"rational_integer_scaled_curve={expected_scaled_curve}")
    print(f"rational_branch_discriminant={expected_rational_discriminant}")
    print("rational_h10=degree_10_squarefree_coprime")
    print("rational_exceptional_fibre=smooth_total_ramification_e3")
    print("rational_absolute_irreducibility=PASS_Groebner_basis_[1]")
    print("rational_total_ramification=12")
    print("rational_normalization_genus=4")
    print("cubic_ratio_field=degree_3_irreducible_mod_13")
    print("cubic_orbit_size=3")
    print("cubic_Delta=universal_lead*(y-Y(alpha))^2*h10")
    print(
        "cubic_node_Y="
        "(3790040625*alpha^2-25880297472*alpha-343597383680)"
        "/42471522304"
    )
    print(
        "cubic_node_R="
        "(-7058559375*alpha^2+104983822336*alpha+206158430208)"
        "/25517520000"
    )
    print("cubic_h10=degree_10_squarefree_coprime")
    print("cubic_exceptional_fibre=ordinary_node_two_unramified_branches")
    print("cubic_tangent_u2=explicit_nonzero_power_basis_element")
    print("cubic_tangent_discriminant=explicit_nonzero_power_basis_element")
    print("cubic_absolute_irreducibility=PASS_global_Groebner_basis_[1]")
    print("cubic_total_ramification=10")
    print("cubic_normalization_genus=3")
    print(f"infinity_cubic={infinity_polynomial}")
    print(f"infinity_discriminant={infinity_discriminant}")
    print("infinity=3_distinct_unramified_points")
    print("closed_DW_ratios=4")
    print("hostile_control_D=W=1_squarefree=PASS")
    print("riemann_hurwitz_and_deck_contradiction=MATHEMATICAL_PROOF_REQUIRED")
    print("status=THM2320_FULL_DW_BANK_EXACT_REFEREE")


if __name__ == "__main__":
    main()
