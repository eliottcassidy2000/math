#!/usr/bin/env python3
"""Exact referee for THM-2316's linear C--D ratio closure.

THM-2311 leaves the rational ratio

    t = D^3 / C^4 = -22143375 / 6397664.

The rational representative C=D=1/t turns the normalized spectral equation
into an exact cubic over Q.  This script verifies its complete branch
discriminant, the ten simple branch values, the remaining smooth totally
ramified fibre, the separable infinity fibre, and absolute irreducibility by
an exact polynomial-root Groebner certificate.

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
) -> sp.GroebnerBasis:
    """Test for every possible polynomial root u=a*y^2+b*y+c."""
    a, b, c = sp.symbols("a b c")
    substituted = sp.Poly(
        sp.expand(curve.subs(u, a * y**2 + b * y + c)),
        y,
    )
    coefficients = [
        substituted.coeff_monomial(y**degree) for degree in range(7)
    ]
    return sp.groebner(coefficients, a, b, c, order="grevlex", domain=sp.QQ)


def main() -> None:
    u, y, v, z = sp.symbols("u y v z")

    ratio = -sp.Rational(22143375, 6397664)
    representative = 1 / ratio
    require(
        representative**3 / representative**4 == ratio,
        "C=D=1/t does not represent the weighted ratio",
    )

    curve = spectral_curve(u, y, representative, representative)
    scaled_curve = sp.expand(sp.Rational(91125, 49) * curve)
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
    actual_discriminant = sp.factor(sp.discriminant(curve, u))
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
        sp.factor(curve.subs(y, y0))
        == -sp.Rational(49, 2460375) * (10935 * u - 2704) ** 3,
        "exceptional fibre is not a triple root",
    )
    require(
        sp.factor(sp.diff(curve, y).subs({u: u0, y: y0}))
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
        z**6 * curve.subs({u: v / z**2, y: 1 / z})
    )
    require(
        sp.expand(infinity_chart.subs(z, 0) - infinity_polynomial) == 0,
        "weighted infinity chart changed",
    )

    # A reducible cubic with constant leading u-coefficient has a polynomial
    # root.  Its y-degree is at most two by the leading-degree comparison.
    # The coefficient ideal for u=a*y^2+b*y+c is the unit ideal.
    groebner = polynomial_root_groebner(curve, u, y)
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

    simple_ramification = h10_poly.degree()
    total_ramification = simple_ramification + 2
    genus = (total_ramification - 4) // 2
    require(simple_ramification == 10, "simple branch count changed")
    require(total_ramification == 12, "total ramification changed")
    require(genus == 4, "normalization genus changed")

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
    print("total_ramification=12")
    print("normalization_genus=4")
    print("hostile_control_C=D=1_squarefree=PASS")
    print("riemann_hurwitz_and_deck_contradiction=MATHEMATICAL_PROOF_REQUIRED")
    print("status=THM2316_CD_LINEAR_RATIO_EXACT_REFEREE")


if __name__ == "__main__":
    main()
