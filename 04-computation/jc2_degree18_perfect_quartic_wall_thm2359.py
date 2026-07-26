#!/usr/bin/env python3
"""Exact referee for THM-2359's perfect-quartic wall closure.

On 504*D=115*B^2, normalize B=1.  Then

    P=5*(27+7*y^2)^2

and the Mordell polynomial factors over Q(sqrt(-5)) into two sextics.
A low square class forces gcd(F,F') to have degree at least four.

For C!=0, put X=C^2 and Z=W/C, then center the unique modularly visible
orbit by

    U=8748*X+175,  V=735*Z+682.

Three terminal-subresultant coefficients generate an exact grevlex ideal
containing U^5 and V^5.  The resulting unique orbit has a sixth-order
root, but its residual sextic is squarefree, so its squarefree part has
degree six and cannot be one of THM-2332's degrees 0, 2, 4.

Every executable check remains active under optimized Python.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive(expression: sp.Expr, *variables: sp.Symbol) -> sp.Expr:
    content, polynomial = sp.Poly(expression, *variables).primitive()
    require(content != 0, "primitive normalization met the zero polynomial")
    return polynomial.as_expr()


def reduce_quadratic(
    expression: sp.Expr,
    radical: sp.Symbol,
    radical_square: sp.Expr,
) -> sp.Expr:
    """Reduce a rational expression modulo radical^2-radical_square."""

    numerator, denominator = sp.together(expression).as_numer_denom()
    require(
        not sp.sympify(denominator).has(radical),
        "quadratic reduction met a radical in the denominator",
    )
    remainder = sp.rem(
        sp.Poly(numerator, radical),
        sp.Poly(radical**2 - radical_square, radical),
    ).as_expr()
    return sp.cancel(remainder / denominator)


def ratio_chart_generator(
    expression: sp.Expr,
    expected_c_power: int,
    cvar: sp.Symbol,
    wvar: sp.Symbol,
    xvar: sp.Symbol,
    zvar: sp.Symbol,
    uvar: sp.Symbol,
    vvar: sp.Symbol,
) -> sp.Expr:
    """Pass from (C,W) to the shifted ratio chart (U,V)."""

    in_c_z = sp.Poly(sp.expand(expression.subs(wvar, zvar * cvar)), cvar, zvar)
    c_power = min(c_exponent for (c_exponent, _), _ in in_c_z.terms())
    require(
        c_power == expected_c_power,
        "a selected coefficient changed its exact C-order",
    )
    after_c = sp.cancel(in_c_z.as_expr() / cvar**c_power)
    after_c_poly = sp.Poly(after_c, cvar, zvar)
    require(
        all(c_exponent % 2 == 0 for (c_exponent, _), _ in after_c_poly.terms()),
        "the (C,W)->(-C,-W) parity reduction failed",
    )
    in_x_z = sum(
        coefficient * xvar ** (c_exponent // 2) * zvar**z_exponent
        for (c_exponent, z_exponent), coefficient in after_c_poly.terms()
    )
    shifted = sp.together(
        in_x_z.subs(
            {
                xvar: (uvar - 175) / 8748,
                zvar: (vvar - 682) / 735,
            }
        )
    )
    numerator, denominator = shifted.as_numer_denom()
    require(
        denominator.is_Integer and denominator != 0,
        "the shifted ratio generator acquired a nonconstant denominator",
    )
    return primitive(numerator, uvar, vvar)


def main() -> None:
    y = sp.symbols("y")
    cvar, wvar = sp.symbols("C W")
    xvar, zvar, uvar, vvar = sp.symbols("X Z U V")

    dvalue = sp.Rational(115, 504)
    quartic_p = (
        245 * y**4
        + 1890 * y**2
        - 24300
        + 122472 * dvalue
    )
    sextic_q = (
        539 * y**6
        + 11340 * y**4
        + 183708 * cvar * y**3
        + (72900 - 367416 * dvalue) * y**2
        + (2361960 * cvar + 2480058 * wvar) * y
    )
    mordell = sp.expand(4 * quartic_p**3 + 49 * sextic_q**2)
    line_l = 27 + 7 * y**2

    require(
        sp.expand(quartic_p - 5 * line_l**2) == 0,
        "the perfect-quartic wall factorization changed",
    )
    mordell_poly = sp.Poly(mordell, y)
    require(
        mordell_poly.degree() == 12 and mordell_poly.LC() == 73060029,
        "the Mordell polynomial changed degree or leading coefficient",
    )

    # The split over Q(sqrt(-5)) is a structural sidecar, not the
    # square-class conclusion by itself.
    eta = sp.symbols("eta")
    split_plus = 7 * sextic_q + 10 * eta * line_l**3
    split_minus = 7 * sextic_q - 10 * eta * line_l**3
    require(
        reduce_quadratic(
            sp.expand(split_plus * split_minus - mordell),
            eta,
            -5,
        )
        == 0,
        "the Q(sqrt(-5)) sextic factorization changed",
    )

    q_mod_l = sp.rem(sextic_q, line_l, y)
    expected_q_mod_l = (
        sp.Rational(39366, 7)
        * (y * (294 * cvar + 441 * wvar) + 32)
    )
    common_resultant = sp.factor(sp.resultant(line_l, sextic_q, y))
    common_conic = (
        333396 * cvar**2
        + 1000188 * cvar * wvar
        + 750141 * wvar**2
        + 1024
    )
    require(
        sp.expand(q_mod_l - expected_q_mod_l) == 0
        and sp.expand(
            common_resultant - 3720786376356 * common_conic
        )
        == 0,
        "the L--Q common-root sidecar changed",
    )

    subresultants = sp.subresultants(mordell, sp.diff(mordell, y), y)
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

    # The C=0 edge cannot have gcd degree at least four: two coefficients
    # that would both have to vanish are coprime polynomials in W.
    edge_first = sp.Poly(
        primitive(degree_three.all_coeffs()[0].subs(cvar, 0), wvar),
        wvar,
    )
    edge_second = sp.Poly(
        primitive(degree_three.all_coeffs()[1].subs(cvar, 0), wvar),
        wvar,
    )
    edge_gcd = sp.gcd(edge_first, edge_second)
    require(
        (
            edge_first.degree(),
            len(edge_first.terms()),
            edge_second.degree(),
            len(edge_second.terms()),
            edge_gcd.degree(),
        )
        == (14, 8, 13, 7, 0),
        "the C=0 subresultant edge certificate changed",
    )

    # On C!=0, these three vanishing coefficients suffice.  The constant
    # coefficient of Sres_3 has exact C-order one and may be divided by C.
    shifted_generators = [
        ratio_chart_generator(
            degree_three.LC(),
            0,
            cvar,
            wvar,
            xvar,
            zvar,
            uvar,
            vvar,
        ),
        ratio_chart_generator(
            degree_three.TC(),
            1,
            cvar,
            wvar,
            xvar,
            zvar,
            uvar,
            vvar,
        ),
        ratio_chart_generator(
            degree_two.LC(),
            0,
            cvar,
            wvar,
            xvar,
            zvar,
            uvar,
            vvar,
        ),
    ]
    shifted_signatures = [
        (
            sp.total_degree(generator),
            len(sp.Poly(generator, uvar, vvar).terms()),
            sp.degree(generator, uvar),
            sp.degree(generator, vvar),
        )
        for generator in shifted_generators
    ]
    require(
        shifted_signatures
        == [
            (21, 151, 11, 14),
            (18, 109, 9, 12),
            (27, 219, 12, 18),
        ],
        "the shifted terminal-generator signatures changed",
    )

    basis = sp.groebner(
        shifted_generators,
        uvar,
        vvar,
        order="grevlex",
        method="f5b",
    )
    expected_basis = [
        -212914386392781422592 * uvar**2
        - 776250367057015603200 * uvar * vvar
        + 1873417064878541875 * vvar**4
        + 19823498739475200000 * vvar**3
        - 707519865807175680000 * vvar**2,
        845652816485908992 * uvar**3
        + 25755182917786828800 * uvar**2
        + 93899104387764480000 * uvar * vvar
        + 2724675804767984375 * vvar**3
        + 85585121186764500000 * vvar**2,
        15660237342331648 * uvar**2 * vvar
        - 8718368619626496 * uvar**2
        - 31785718925721600 * uvar * vvar
        - 51227634911866875 * vvar**3
        - 28971358395840000 * vvar**2,
        -8477454513299521536 * uvar**2
        + 3083109226771543200 * uvar * vvar**2
        - 30907386246404505600 * uvar * vvar
        + 6409548769616444375 * vvar**3
        - 28170794755837440000 * vvar**2,
    ]
    actual_basis = [item.as_expr() for item in basis.polys]
    require(
        len(actual_basis) == 4
        and all(
            sp.expand(actual - expected) == 0
            for actual, expected in zip(actual_basis, expected_basis)
        ),
        "the exact four-element Groebner basis changed",
    )
    _, u_remainder = basis.reduce(uvar**5)
    _, v_remainder = basis.reduce(vvar**5)
    require(
        u_remainder == 0 and v_remainder == 0,
        "the U^5/V^5 ideal-membership certificate changed",
    )

    # Audit the unique orbit in Q(sqrt(-21)).
    radical = sp.symbols("r")
    c_exceptional = sp.Rational(5, 162) * radical
    w_exceptional = -sp.Rational(341, 11907) * radical
    require(
        reduce_quadratic(
            8748 * c_exceptional**2 + 175,
            radical,
            -21,
        )
        == 0
        and sp.expand(735 * w_exceptional + 682 * c_exceptional) == 0,
        "the exceptional orbit stopped satisfying U=V=0",
    )

    residual_six = (
        y**6
        - sp.Rational(18, 7) * radical * y**5
        - sp.Rational(8721, 161) * y**4
        + sp.Rational(1356, 49) * radical * y**3
        + sp.Rational(203472, 1127) * y**2
        - sp.Rational(243000, 7889) * radical * y
        - sp.Rational(364500, 7889)
    )
    exceptional_mordell = mordell.subs(
        {
            cvar: c_exceptional,
            wvar: w_exceptional,
        }
    )
    require(
        reduce_quadratic(
            exceptional_mordell
            - 73060029
            * (y + sp.Rational(3, 7) * radical) ** 6
            * residual_six,
            radical,
            -21,
        )
        == 0,
        "the exceptional sixth-order factorization changed",
    )

    residual_discriminant = reduce_quadratic(
        sp.discriminant(residual_six, y),
        radical,
        -21,
    )
    expected_discriminant = sp.Rational(
        1242517308060047668936704000000000000,
        8551083049093133387426609,
    )
    root_value = reduce_quadratic(
        residual_six.subs(y, -sp.Rational(3, 7) * radical),
        radical,
        -21,
    )
    require(
        residual_discriminant == expected_discriminant
        and root_value == -sp.Rational(1093500, 343),
        "the residual-sextic squarefree/coprime certificate changed",
    )

    def gcd_degree(cvalue: sp.Expr, wvalue: sp.Expr) -> int:
        specialization = sp.Poly(
            mordell.subs({cvar: cvalue, wvar: wvalue}),
            y,
            domain=sp.QQ,
        )
        return sp.gcd(specialization, specialization.diff()).degree()

    hostile_degrees = (gcd_degree(1, 0), gcd_degree(0, 1))
    require(
        hostile_degrees == (0, 0),
        "the off-orbit hostile gcd controls changed",
    )

    print("THM-2359 exact perfect-quartic-wall referee")
    print("wall factorization: P=5*(27+7*y^2)^2")
    print("quadratic split: F=(7Q+10*sqrt(-5)*L^3)*(7Q-10*sqrt(-5)*L^3)")
    print(
        "common-root resultant:"
        " 3720786376356"
        "*(333396*C^2+1000188*C*W+750141*W^2+1024)"
    )
    print("subresultant degrees: 12,11,10,9,8,7,6,5,4,3,2,1,0")
    print("C=0 edge signatures/gcd degree: (14,8),(13,7) / 0")
    print("ratio chart: X=C^2, Z=W/C, U=8748*X+175, V=735*Z+682")
    print(
        "shifted generator signatures:"
        " [(21,151,11,14),(18,109,9,12),(27,219,12,18)]"
    )
    print("Groebner basis degrees: 4,3,3,3")
    print("nilpotence certificate: U^5 and V^5 have zero remainder")
    print("unique orbit: C^2=-175/8748, W/C=-682/735")
    print("exceptional factorization: F=73060029*(y+3*r/7)^6*R6, r^2=-21")
    print(f"residual sextic discriminant: {residual_discriminant}")
    print(f"residual value at repeated root: {root_value}")
    print("exceptional squarefree degree: 6")
    print(f"off-orbit hostile gcd degrees: {hostile_degrees}")
    print("VERDICT: the perfect-quartic wall has no square class of degree <=4")


if __name__ == "__main__":
    main()
