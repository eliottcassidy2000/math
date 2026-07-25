#!/usr/bin/env python3
"""Exact referee for THM-2338's deep common-root wall.

The script restricts the structured Mordell polynomial from THM-2332 to
126D=25B^2 and 21W=-20BC, verifies its y^6-times-sextic factorization,
factors the sextic discriminant, and checks the exact gcd at each of the
three multiple-root ratios.  Every executable check remains active under
optimized Python.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


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

    deep_wall = {
        dvar: sp.Rational(25, 126) * bvar**2,
        wvar: -sp.Rational(20, 21) * bvar * cvar,
    }
    wall_p = sp.factor(quartic_p.subs(deep_wall))
    wall_q = sp.factor(sextic_q.subs(deep_wall))
    require(
        wall_p == 35 * y**2 * (54 * bvar + 7 * y**2),
        "wall factorization of P changed",
    )
    require(
        wall_q
        == 7
        * y**3
        * (1620 * bvar * y + 26244 * cvar + 77 * y**3),
        "wall factorization of Q changed",
    )

    residual = (
        7889 * y**6
        + 211680 * bvar * y**4
        + 1047816 * cvar * y**3
        + 1814400 * bvar**2 * y**2
        + 22044960 * bvar * cvar * y
        + 2916000 * bvar**3
        + 178564176 * cvar**2
    )
    wall_mordell = sp.factor(mordell.subs(deep_wall))
    require(
        sp.expand(wall_mordell - 9261 * y**6 * residual) == 0,
        "wall Mordell polynomial is not 9261*y^6 times the residual sextic",
    )

    first_factor = 361 * bvar**3 + 30618 * cvar**2
    second_factor = (
        62500 * bvar**6
        + 7654500 * bvar**3 * cvar**2
        + 199644669 * cvar**4
    )
    residual_discriminant = sp.factor(sp.discriminant(residual, y))
    discriminant_scalar = -(2**29) * (3**30) * (5**9) * (7**15) * 23
    require(
        sp.expand(
            residual_discriminant
            - discriminant_scalar * first_factor**3 * second_factor
        )
        == 0,
        "residual sextic discriminant factorization changed",
    )

    # The B=0 projective edge has C!=0, where both discriminant factors
    # are nonzero.  Thus it has square-class degree six, not two or four.
    require(
        first_factor.subs({bvar: 0, cvar: 1}) == 30618
        and second_factor.subs({bvar: 0, cvar: 1}) == 199644669,
        "B=0 edge control changed",
    )

    # On B!=0, weighted scaling normalizes B=1.  Put t=C^2.  The first
    # discriminant factor has one root and the second has two.
    tvar = sp.symbols("t")
    first_ratio_polynomial = 361 + 30618 * tvar
    second_ratio_polynomial = (
        62500 + 7654500 * tvar + 199644669 * tvar**2
    )
    t_first = -sp.Rational(361, 30618)
    t_plus = (
        -sp.Rational(250, 13041)
        + sp.Rational(500, 117369) * sp.sqrt(3)
    )
    t_minus = (
        -sp.Rational(250, 13041)
        - sp.Rational(500, 117369) * sp.sqrt(3)
    )
    ratio_roots = (t_first, t_plus, t_minus)
    require(
        sp.factor(first_ratio_polynomial.subs(tvar, t_first)) == 0,
        "the first exact ratio stopped vanishing",
    )
    require(
        all(
            sp.simplify(second_ratio_polynomial.subs(tvar, item)) == 0
            for item in (t_plus, t_minus)
        ),
        "the quadratic exact ratios stopped vanishing",
    )
    require(
        sp.factor(second_ratio_polynomial.subs(tvar, t_first))
        == sp.Rational(383, 108),
        "the linear and quadratic ratio factors now intersect",
    )
    require(
        sp.discriminant(second_ratio_polynomial, tvar)
        == 8680203000000,
        "the quadratic ratio discriminant changed",
    )
    require(
        all(
            sp.simplify(ratio_roots[left] - ratio_roots[right]) != 0
            for left in range(3)
            for right in range(left)
        )
        and all(item != 0 for item in ratio_roots),
        "the three nonzero ratio values ceased to be distinct",
    )

    # Reconstruct the unique repeated root over each exact number field.
    # Writing C^2=t and y=kC, the linear-factor branch has k=-486/19.
    # On the quadratic branch, k=-(3168963t+101250)/2500.
    repeated_root_data = (
        (t_first, -sp.Rational(486, 19)),
        (
            t_plus,
            -sp.Rational(81, 5) - sp.Rational(27, 5) * sp.sqrt(3),
        ),
        (
            t_minus,
            -sp.Rational(81, 5) + sp.Rational(27, 5) * sp.sqrt(3),
        ),
    )
    gcd_degrees: list[int] = []
    for ratio, root_scale in repeated_root_data:
        cvalue = sp.sqrt(ratio)
        specialized = sp.Poly(
            residual.subs({bvar: 1, cvar: cvalue}),
            y,
            extension=True,
        )
        common = sp.gcd(specialized, specialized.diff()).monic()
        expected_common = y - root_scale * cvalue
        require(
            sp.degree(common.as_expr(), y) == 1
            and sp.simplify(common.as_expr() - expected_common) == 0,
            f"multiple-root gcd changed at ratio {ratio}",
        )
        gcd_degrees.append(sp.degree(common.as_expr(), y))

    # The quadratic-branch root scale can also be recovered without
    # radicals from its affine-linear relation with t.
    for ratio, root_scale in repeated_root_data[1:]:
        require(
            sp.simplify(
                2500 * root_scale + 3168963 * ratio + 101250
            )
            == 0
            and sp.simplify(
                25 * root_scale**2 + 810 * root_scale + 4374
            )
            == 0,
            "quadratic-branch repeated-root response changed",
        )

    # Generic controls on both projective charts.
    generic_b_chart = sp.Poly(
        residual.subs({bvar: 1, cvar: 0}),
        y,
        domain=sp.QQ,
    )
    generic_c_edge = sp.Poly(
        residual.subs({bvar: 0, cvar: 1}),
        y,
        domain=sp.QQ,
    )
    require(
        sp.gcd(generic_b_chart, generic_b_chart.diff()).degree() == 0
        and sp.gcd(generic_c_edge, generic_c_edge.diff()).degree() == 0,
        "generic squarefree hostile controls changed",
    )

    print("THM-2338 exact deep-common-root-wall referee")
    print("wall factorization: F=9261*y^6*G6")
    print(
        "Disc(G6)=-2^29*3^30*5^9*7^15*23"
        "*(361B^3+30618C^2)^3"
        "*(62500B^6+7654500B^3C^2+199644669C^4)"
    )
    print(f"ratio t0={t_first}")
    print(f"ratio t+={t_plus}")
    print(f"ratio t-={t_minus}")
    print(f"multiple-root gcd degrees={gcd_degrees}")
    print("generic B-chart and B=0 edge: residual sextic squarefree")
    print("VERDICT: exactly three H4 orbits and no H2 orbit on the wall")


if __name__ == "__main__":
    main()
