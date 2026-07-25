#!/usr/bin/env python3
"""Exact referee for THM-2342's first-flux cover obstruction.

The script verifies the homogeneous infinity direction cubic of the
normalized deep-wall plane cubic, reconstructs the recovered first-flux
function Z, and proves that its leading cubic is nonzero on all three
infinity directions.  The line-through-node and double-cover
Riemann--Hurwitz deductions are mathematical steps in the theorem.
Every executable check remains active under optimized Python.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def homogeneous_piece(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    degree: int,
) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), *variables)
    answer = sp.Integer(0)
    for powers, coefficient in polynomial.terms():
        if sum(powers) == degree:
            monomial = coefficient
            for variable, power in zip(variables, powers):
                monomial *= variable**power
            answer += monomial
    return sp.expand(answer)


def main() -> None:
    y, z, h = sp.symbols("y z h")
    cvar = sp.symbols("C")

    wall_p = 35 * y**2 * (54 + 7 * y**2)
    wall_q = 7 * y**3 * (1620 * y + 26244 * cvar + 77 * y**3)
    p_scale = sp.Rational(16, 964467)
    q_scale = sp.Rational(64, 703096443)
    plane_cubic = sp.expand(
        z**3
        + p_scale * (wall_p / y**2) * z
        + q_scale * (wall_q / y**3)
    )
    require(
        sp.Poly(plane_cubic, y, z).total_degree() == 3,
        "the normalized deep-wall equation stopped being a plane cubic",
    )

    leading_plane_cubic = homogeneous_piece(plane_cubic, (y, z), 3)
    direction_cubic = sp.factor(leading_plane_cubic.subs({y: 1, z: h}))
    expected_direction_cubic = (
        h**3
        + sp.Rational(80, 19683) * h
        + sp.Rational(704, 14348907)
    )
    require(
        sp.expand(direction_cubic - expected_direction_cubic) == 0,
        "the infinity direction cubic changed",
    )
    direction_discriminant = sp.discriminant(direction_cubic, h)
    require(
        direction_discriminant == -sp.Rational(94208, 282429536481),
        "the infinity directions stopped being three distinct points",
    )
    require(
        homogeneous_piece(plane_cubic, (y, z), 3).subs({y: 0, z: 1})
        == 1,
        "a vertical infinity direction appeared",
    )

    # Recover original u from the depressed coordinate v=yz.  With
    # B=1, a2=49601160+1607445y^2 and a3=-26040609.
    a2 = 49601160 + 1607445 * y**2
    a3 = -26040609
    original_u = sp.factor(y * z - a2 / (3 * a3))
    require(
        sp.cancel(
            original_u
            - (35 * y**2 + 1701 * y * z + 1080) / 1701
        )
        == 0,
        "the inverse depressed-coordinate translation changed",
    )

    n2 = (
        45927 * original_u**2
        - 58320 * original_u
        - 5670 * original_u * y**2
        + 2160 * y**2
        + 93312 * sp.Rational(25, 126)
        - 15552 * cvar * y
        + 35 * y**4
    )
    first_flux_numerator = (
        139968 * cvar
        + 560 * y**3
        + 34020 * y**2 * z
        - 413343 * y * z**2
        + 12960 * y
    )
    require(
        sp.expand(n2 + sp.Rational(1, 9) * y * first_flux_numerator)
        == 0,
        "the deep-wall N2 factorization changed",
    )
    recovered_z = sp.Rational(2, 45927) * first_flux_numerator
    require(
        sp.expand(recovered_z + 2 * n2 / (5103 * y)) == 0,
        "the recovered first-flux identity Z=-2N2/(5103y) changed",
    )

    leading_flux = homogeneous_piece(recovered_z, (y, z), 3)
    direction_flux = sp.factor(
        (sp.Rational(45927, 2) * leading_flux).subs({y: 1, z: h})
    )
    expected_direction_flux = 560 + 34020 * h - 413343 * h**2
    require(
        sp.expand(direction_flux - expected_direction_flux) == 0,
        "the leading first-flux direction polynomial changed",
    )
    direction_resultant = sp.resultant(
        direction_cubic,
        direction_flux,
        h,
    )
    require(
        direction_resultant == 1024192512
        and sp.factorint(direction_resultant)
        == {2: 12, 3: 6, 7: 3},
        "the first flux acquired a cancellation at infinity",
    )

    # A connected double cover of P1 branched at an even number r of
    # points has genus (r-2)/2.  Three certified odd poles force r>=4.
    certified_odd_poles = 3
    minimum_branch_points = (
        certified_odd_poles
        if certified_odd_poles % 2 == 0
        else certified_odd_poles + 1
    )
    minimum_double_cover_genus = (minimum_branch_points - 2) // 2
    require(
        minimum_branch_points == 4
        and minimum_double_cover_genus == 1,
        "the double-cover genus lower bound changed",
    )

    print("THM-2342 exact deep-wall first-flux-cover referee")
    print(
        "infinity direction cubic="
        "h^3+(80/19683)h+704/14348907"
    )
    print(f"direction discriminant={direction_discriminant}")
    print(
        "Z=(2/45927)*(139968C+560y^3+34020y^2z"
        "-413343yz^2+12960y)"
    )
    print("leading flux direction=560+34020h-413343h^2")
    print(f"direction resultant={direction_resultant}")
    print("certified odd infinity poles=3")
    print(f"minimum double-cover branch points={minimum_branch_points}")
    print(f"minimum double-cover genus={minimum_double_cover_genus}")
    print("VERDICT: neither rational deep-wall curve lifts through Z=T^2")


if __name__ == "__main__":
    main()
