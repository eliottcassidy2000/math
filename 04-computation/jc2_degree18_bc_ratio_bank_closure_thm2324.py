#!/usr/bin/env python3
"""Exact referee for THM-2324's full B--C ratio-bank closure.

The two rational linear-factor points are checked separately.  The
quadratic and quintic factors are treated uniformly over their ratio
fields, so one representative proves the statement for every conjugate.
Every check remains active under optimized Python.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def spectral_curve(
    u: sp.Symbol,
    y: sp.Symbol,
    bvar: sp.Expr,
    cvar: sp.Expr,
) -> sp.Expr:
    """THM-2297's G_0 at (B,C,D,W)=(bvar,cvar,0,0)."""
    return sp.expand(
        -26040609 * u**3
        + 49601160 * bvar * u**2
        + 1607445 * u**2 * y**2
        - 20995200 * bvar**2 * u
        - 2857680 * bvar * u * y**2
        - 138915 * u * y**4
        + 777600 * bvar**2 * y**2
        - 5598720 * bvar * cvar * y
        + 78120 * bvar * y**4
        - 435456 * cvar * y**3
        + 1127 * y**6
    )


def field_reduce(expression: sp.Expr, field: sp.AlgebraicField) -> sp.Expr:
    """Return the canonical power-basis representative in an algebraic field."""
    return sp.factor(field.to_sympy(field.from_sympy(sp.cancel(expression))))


def primitive_curve(
    curve: sp.Expr,
    u: sp.Symbol,
    y: sp.Symbol,
) -> sp.Expr:
    polynomial = sp.Poly(curve, u, y, domain=sp.QQ)
    _, integral = polynomial.clear_denoms()
    _, primitive = integral.primitive()
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive.as_expr()


def polynomial_root_groebner(
    curve: sp.Expr,
    u: sp.Symbol,
    y: sp.Symbol,
    extra_equations: list[sp.Expr] | None = None,
    extra_generators: list[sp.Symbol] | None = None,
) -> sp.GroebnerBasis:
    """Rule out every possible root u=a*y^2+b*y+c."""
    a, b, c = sp.symbols("a b c")
    substituted = sp.Poly(
        sp.expand(curve.subs(u, a * y**2 + b * y + c)),
        y,
    )
    equations = [
        substituted.coeff_monomial(y**degree) for degree in range(7)
    ]
    equations.extend(extra_equations or [])
    generators: list[sp.Symbol] = [a, b, c]
    generators.extend(extra_generators or [])
    return sp.groebner(
        equations,
        *generators,
        order="grevlex",
        domain=sp.QQ,
    )


def require_unit_groebner(basis: sp.GroebnerBasis, message: str) -> None:
    require(
        len(basis.polys) == 1 and basis.polys[0].as_expr() == 1,
        message,
    )


def polynomial_power_mod(
    base: sp.Poly,
    exponent: int,
    modulus: sp.Poly,
) -> sp.Poly:
    """Binary polynomial powering in the finite field of modulus.domain."""
    result = sp.Poly(1, *modulus.gens, domain=modulus.domain)
    power = base.rem(modulus)
    while exponent:
        if exponent & 1:
            result = (result * power).rem(modulus)
        power = (power * power).rem(modulus)
        exponent //= 2
    return result


def main() -> None:
    u, y, v, z, t = sp.symbols("u y v z t")

    infinity_polynomial = (
        1127 - 138915 * v + 1607445 * v**2 - 26040609 * v**3
    )
    infinity_discriminant = sp.discriminant(infinity_polynomial, v)
    require(
        infinity_discriminant == -153384762202971019112448,
        "weighted infinity cubic is not separable",
    )

    rational_data = (
        {
            "name": "linear_2000_15309",
            "ratio": -sp.Rational(2000, 15309),
            "primitive": (
                3321506250 * u**3
                - 205031250 * u**2 * y**2
                + 48427561125 * u**2
                + 17718750 * u * y**4
                - 2790065250 * u * y**2
                + 156905298045 * u
                - 143750 * y**6
                + 76271625 * y**4
                - 425152800 * y**3
                - 5811307335 * y**2
                + 41841412812 * y
            ),
            "exceptional": 10 * y - 81,
            "residual": (
                (200 * y**2 - 19683)
                * (
                    18400000000 * y**8
                    + 298080000000 * y**7
                    + 1653372000000 * y**6
                    - 11479125600000 * y**5
                    - 441659357460000 * y**4
                    - 3389154437772000 * y**3
                    + 2287679245496100 * y**2
                    + 111181211331110460 * y
                    + 450283905890997363
                )
            ),
            "disc_scalar": -sp.Rational(
                1628150074335205281,
                3906250000,
            ),
            "y0": sp.Rational(81, 10),
            "u0": -sp.Rational(81, 100),
            "fibre": (
                -sp.Rational(26040609, 1000000)
                * (100 * u + 81) ** 2
                * (100 * u + 891)
            ),
            "tangent_u2": -sp.Rational(2109289329, 10),
            "tangent_disc": -sp.Rational(
                106778435362398485784,
                15625,
            ),
        },
        {
            "name": "linear_125_1134",
            "ratio": -sp.Rational(125, 1134),
            "primitive": (
                1660753125 * u**3
                - 102515625 * u**2 * y**2
                + 28697814000 * u**2
                + 8859375 * u * y**4
                - 1653372000 * u * y**2
                + 110199605760 * u
                - 71875 * y**6
                + 45198000 * y**4
                - 251942400 * y**3
                - 4081466880 * y**2
                + 29386561536 * y
            ),
            "exceptional": 5 * y - 54,
            "residual": (
                (5 * y + 54)
                * (
                    44921875 * y**9
                    + 485156250 * y**8
                    - 455625000 * y**7
                    - 59049000000 * y**6
                    - 1142598150000 * y**5
                    - 2008846980000 * y**4
                    + 77484097800000 * y**3
                    + 234311911747200 * y**2
                    - 1446039226782720 * y
                    - 15617223649253376
                )
            ),
            "disc_scalar": -sp.Rational(
                6668902704477000830976,
                244140625,
            ),
            "y0": sp.Rational(54, 5),
            "u0": -sp.Rational(36, 25),
            "fibre": (
                -sp.Rational(26040609, 3125)
                * (5 * u + 36)
                * (25 * u + 36) ** 2
            ),
            "tangent_u2": -sp.Rational(3749847696, 25),
            "tangent_disc": sp.Rational(
                1199902527419435384832,
                15625,
            ),
        },
    )

    rational_genera: list[int] = []
    for data in rational_data:
        ratio = data["ratio"]
        representative = 1 / ratio
        require(
            representative**2 / representative**3 == ratio,
            f"{data['name']} has the wrong weighted ratio",
        )
        curve = spectral_curve(u, y, representative, representative)
        require(
            primitive_curve(curve, u, y) == data["primitive"],
            f"{data['name']} primitive curve changed",
        )
        discriminant = sp.discriminant(curve, u)
        expected_discriminant = (
            data["disc_scalar"]
            * data["exceptional"] ** 2
            * data["residual"]
        )
        require(
            sp.expand(discriminant - expected_discriminant) == 0,
            f"{data['name']} discriminant changed",
        )
        residual = sp.Poly(data["residual"], y)
        exceptional = sp.Poly(data["exceptional"], y)
        require(
            residual.degree() == 10
            and sp.gcd(residual, residual.diff()).degree() == 0
            and sp.gcd(residual, exceptional).degree() == 0,
            f"{data['name']} residual is not squarefree and coprime",
        )
        require(
            sp.expand(curve.subs(y, data["y0"]) - data["fibre"]) == 0,
            f"{data['name']} exceptional fibre changed",
        )
        point = {u: data["u0"], y: data["y0"]}
        require(
            sp.factor(curve.subs(point)) == 0
            and sp.factor(sp.diff(curve, u).subs(point)) == 0
            and sp.factor(sp.diff(curve, y).subs(point)) == 0,
            f"{data['name']} repeated point is not singular",
        )
        tangent_u2 = sp.factor(sp.diff(curve, u, 2).subs(point) / 2)
        tangent_uy = sp.factor(sp.diff(curve, u, y).subs(point))
        tangent_y2 = sp.factor(sp.diff(curve, y, 2).subs(point) / 2)
        tangent_disc = sp.factor(
            tangent_uy**2 - 4 * tangent_u2 * tangent_y2
        )
        require(
            tangent_u2 == data["tangent_u2"]
            and tangent_disc == data["tangent_disc"],
            f"{data['name']} tangent cone changed",
        )
        require_unit_groebner(
            polynomial_root_groebner(curve, u, y),
            f"{data['name']} gained a polynomial root",
        )
        infinity_chart = sp.expand(
            z**6 * curve.subs({u: v / z**2, y: 1 / z})
        )
        require(
            sp.expand(infinity_chart.subs(z, 0) - infinity_polynomial) == 0,
            f"{data['name']} infinity chart changed",
        )
        total_ramification = residual.degree()
        genus = (total_ramification - 4) // 2
        require(
            total_ramification == 10 and genus == 3,
            f"{data['name']} genus changed",
        )
        rational_genera.append(genus)

    quadratic_ratio = (
        3321125 - 161754894 * t - 2812385772 * t**2
    )
    quadratic_mod_19 = sp.Poly(quadratic_ratio, t, modulus=19).monic()
    require(
        quadratic_mod_19
        == sp.Poly(t**2 - 7 * t + 5, t, modulus=19)
        and all(quadratic_mod_19.eval(residue) != 0 for residue in range(19)),
        "quadratic ratio polynomial is not irreducible modulo nineteen",
    )
    require(
        sp.discriminant(quadratic_ratio, t) == 63525784521085236,
        "quadratic ratio discriminant changed",
    )
    alpha2 = sp.CRootOf(quadratic_ratio, 0)
    field2 = sp.QQ.algebraic_field(alpha2)
    curve2 = spectral_curve(u, y, alpha2, alpha2**2)
    require(
        field_reduce(
            (alpha2**2) ** 2 / alpha2**3 - alpha2,
            field2,
        )
        == 0,
        "quadratic representative has the wrong ratio",
    )
    discriminant2 = sp.Poly(
        sp.discriminant(curve2, u),
        y,
        extension=alpha2,
    )
    repeated2 = sp.gcd(discriminant2, discriminant2.diff()).monic()
    require(repeated2.degree() == 1, "quadratic repeated factor changed")
    node_y2 = field_reduce(-repeated2.coeff_monomial(y**0), field2)
    expected_y2 = (69984 * alpha2 - 4075) / 8919
    require(
        field_reduce(node_y2 - expected_y2, field2) == 0,
        "quadratic exceptional y-coordinate changed",
    )
    quotient2, remainder2 = discriminant2.div(
        sp.Poly((y - node_y2) ** 2, y, extension=alpha2)
    )
    require(
        remainder2.is_zero
        and field_reduce(quotient2.LC(), field2)
        == -153384762202971019112448,
        "quadratic discriminant factorization changed",
    )
    residual2 = quotient2.monic()
    require(
        residual2.degree() == 10
        and sp.gcd(residual2, residual2.diff()).degree() == 0
        and sp.gcd(
            residual2,
            sp.Poly(y - node_y2, y, extension=alpha2),
        ).degree()
        == 0,
        "quadratic residual is not squarefree and coprime",
    )
    exceptional_fibre2 = sp.Poly(
        curve2.subs(y, node_y2),
        u,
        extension=alpha2,
    )
    node_u2 = field_reduce(
        5 * (237718152 * alpha2 + 3321125) / 2867360391,
        field2,
    )
    require(
        sp.Poly(
            sp.expand(
                curve2.subs(y, node_y2)
                + 26040609 * (u - node_u2) ** 3
            ),
            u,
            extension=alpha2,
        ).is_zero,
        "quadratic exceptional fibre is not a triple root",
    )
    derivative_y2 = field_reduce(
        sp.diff(curve2, y).subs({u: node_u2, y: node_y2}),
        field2,
    )
    expected_derivative_y2 = field_reduce(
        sp.Rational(128, 60214568211)
        * (16490073897927 * alpha2 - 235328275250),
        field2,
    )
    require(
        derivative_y2 == expected_derivative_y2
        and derivative_y2 != 0,
        "quadratic triple fibre is not smooth",
    )
    universal_curve = spectral_curve(u, y, t, t**2)
    require_unit_groebner(
        polynomial_root_groebner(
            universal_curve,
            u,
            y,
            [quadratic_ratio],
            [t],
        ),
        "quadratic orbit gained a polynomial root",
    )
    require(
        exceptional_fibre2.degree() == 3,
        "quadratic exceptional fibre degree changed",
    )
    quadratic_total_ramification = residual2.degree() + 2
    quadratic_genus = (quadratic_total_ramification - 4) // 2
    require(
        quadratic_total_ramification == 12 and quadratic_genus == 4,
        "quadratic orbit genus changed",
    )

    quintic_ratio = (
        410644531250000
        - 18114791748046875 * t
        - 545436228093750000 * t**2
        - 4951165276923468750 * t**3
        - 18946967714644599000 * t**4
        - 26529827304546537363 * t**5
    )
    quintic_mod_13 = sp.Poly(quintic_ratio, t, modulus=13).monic()
    require(
        quintic_mod_13
        == sp.Poly(
            t**5 - 3 * t**4 + 2 * t**3 - 3 * t**2 - t - 4,
            t,
            modulus=13,
        ),
        "quintic ratio has the wrong mod-thirteen reduction",
    )
    finite_t = sp.Poly(t, t, modulus=13)
    frobenius_1 = polynomial_power_mod(
        finite_t,
        13,
        quintic_mod_13,
    )
    frobenius_5 = polynomial_power_mod(
        finite_t,
        13**5,
        quintic_mod_13,
    )
    require(
        sp.gcd(
            quintic_mod_13,
            frobenius_1 - finite_t,
        ).degree()
        == 0
        and (frobenius_5 - finite_t).rem(quintic_mod_13).is_zero,
        "quintic ratio polynomial is not irreducible modulo thirteen",
    )
    require(
        sp.discriminant(quintic_ratio, t)
        == -199019205840378862392485682250728932226285629024947372230920300533889671177041558900949524968382320366799831390380859375000000000000,
        "quintic ratio discriminant changed",
    )
    alpha5 = sp.CRootOf(quintic_ratio, 0)
    field5 = sp.QQ.algebraic_field(alpha5)
    curve5 = spectral_curve(u, y, alpha5, alpha5**2)
    require(
        field_reduce(
            (alpha5**2) ** 2 / alpha5**3 - alpha5,
            field5,
        )
        == 0,
        "quintic representative has the wrong ratio",
    )
    discriminant5 = sp.Poly(
        sp.discriminant(curve5, u),
        y,
        extension=alpha5,
    )
    repeated5 = sp.gcd(discriminant5, discriminant5.diff()).monic()
    require(repeated5.degree() == 1, "quintic repeated factor changed")
    node_y5 = field_reduce(-repeated5.coeff_monomial(y**0), field5)
    quotient5, remainder5 = discriminant5.div(
        sp.Poly((y - node_y5) ** 2, y, extension=alpha5)
    )
    require(
        remainder5.is_zero
        and field_reduce(quotient5.LC(), field5)
        == -153384762202971019112448,
        "quintic discriminant factorization changed",
    )
    residual5 = quotient5.monic()
    require(
        residual5.degree() == 10
        and sp.gcd(residual5, residual5.diff()).degree() == 0
        and sp.gcd(
            residual5,
            sp.Poly(y - node_y5, y, extension=alpha5),
        ).degree()
        == 0,
        "quintic residual is not squarefree and coprime",
    )
    exceptional_fibre5 = sp.Poly(
        curve5.subs(y, node_y5),
        u,
        extension=alpha5,
    )
    repeated_u5 = sp.gcd(
        exceptional_fibre5,
        exceptional_fibre5.diff(),
    ).monic()
    require(
        repeated_u5.degree() == 1,
        "quintic exceptional root is not double",
    )
    node_u5 = field_reduce(
        -repeated_u5.coeff_monomial(u**0),
        field5,
    )
    simple_u5 = field_reduce(
        sp.Rational(40, 21) * alpha5
        + sp.Rational(5, 81) * node_y5**2
        - 2 * node_u5,
        field5,
    )
    require(
        sp.Poly(
            sp.expand(
                curve5.subs(y, node_y5)
                + 26040609
                * (u - node_u5) ** 2
                * (u - simple_u5)
            ),
            u,
            extension=alpha5,
        ).is_zero,
        "quintic exceptional fibre factorization changed",
    )
    require(
        field_reduce(curve5.subs({u: node_u5, y: node_y5}), field5)
        == 0
        and field_reduce(
            sp.diff(curve5, u).subs({u: node_u5, y: node_y5}),
            field5,
        )
        == 0
        and field_reduce(
            sp.diff(curve5, y).subs({u: node_u5, y: node_y5}),
            field5,
        )
        == 0,
        "quintic repeated point is not singular",
    )
    tangent_u2_5 = field_reduce(
        sp.diff(curve5, u, 2).subs({u: node_u5, y: node_y5}) / 2,
        field5,
    )
    tangent_uy_5 = field_reduce(
        sp.diff(curve5, u, y).subs({u: node_u5, y: node_y5}),
        field5,
    )
    tangent_y2_5 = field_reduce(
        sp.diff(curve5, y, 2).subs({u: node_u5, y: node_y5}) / 2,
        field5,
    )
    tangent_discriminant5 = field_reduce(
        tangent_uy_5**2 - 4 * tangent_u2_5 * tangent_y2_5,
        field5,
    )
    require(
        tangent_u2_5 != 0 and tangent_discriminant5 != 0,
        "quintic tangent cone is not ordinary and nonvertical",
    )
    require_unit_groebner(
        polynomial_root_groebner(
            universal_curve,
            u,
            y,
            [quintic_ratio],
            [t],
        ),
        "quintic orbit gained a polynomial root",
    )
    quintic_total_ramification = residual5.degree()
    quintic_genus = (quintic_total_ramification - 4) // 2
    require(
        quintic_total_ramification == 10 and quintic_genus == 3,
        "quintic orbit genus changed",
    )

    control_curve = spectral_curve(u, y, sp.Integer(1), sp.Integer(1))
    control_discriminant = sp.Poly(sp.discriminant(control_curve, u), y)
    require(
        control_discriminant.degree() == 12
        and sp.gcd(
            control_discriminant,
            control_discriminant.diff(),
        ).degree()
        == 0,
        "B=C=1 hostile control unexpectedly has repeated branch",
    )

    print("rational_linear_orbits=2")
    print("rational_weighted_representatives=B=C=1/t")
    print("rational_absolute_irreducibility=PASS_two_Groebner_bases_[1]")
    print("rational_Delta=(linear_exceptional)^2*h10")
    print("rational_h10=degree_10_squarefree_coprime")
    print("rational_exceptional_fibres=ordinary_nodes_nonvertical")
    print(f"rational_normalization_genera={tuple(rational_genera)}")
    print("quadratic_ratio_field=degree_2_irreducible_mod_19")
    print("quadratic_orbit_size=2")
    print("quadratic_absolute_irreducibility=PASS_global_Groebner_basis_[1]")
    print("quadratic_Delta=universal_lead*(y-Y(alpha))^2*h10")
    print("quadratic_h10=degree_10_squarefree_coprime")
    print("quadratic_exceptional_fibre=smooth_total_ramification_e3")
    print(f"quadratic_normalization_genus={quadratic_genus}")
    print("quintic_ratio_field=degree_5_irreducible_mod_13")
    print("quintic_orbit_size=5")
    print("quintic_absolute_irreducibility=PASS_global_Groebner_basis_[1]")
    print("quintic_Delta=universal_lead*(y-Y(alpha))^2*h10")
    print("quintic_h10=degree_10_squarefree_coprime")
    print("quintic_exceptional_fibre=ordinary_node_two_unramified_branches")
    print(f"quintic_normalization_genus={quintic_genus}")
    print(f"infinity_cubic={infinity_polynomial}")
    print(f"infinity_discriminant={infinity_discriminant}")
    print("infinity=3_distinct_unramified_points")
    print("closed_BC_ratios=9")
    print("hostile_control_B=C=1_squarefree=PASS")
    print("riemann_hurwitz_and_deck_contradiction=MATHEMATICAL_PROOF_REQUIRED")
    print("status=THM2324_FULL_BC_BANK_EXACT_REFEREE")


if __name__ == "__main__":
    main()
