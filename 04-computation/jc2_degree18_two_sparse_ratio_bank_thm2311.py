#!/usr/bin/env python3
"""Exact referee for the degree-18 two-sparse weighted-ratio bank.

Starting from THM-2297's normalized spectral cubic, this script computes the
branch discriminant and the repeated-branch resultant on all six coordinate
planes.  It checks the five nonclosed planes against frozen primitive
weighted-ratio polynomials, verifies squarefreeness and pairwise coprimality,
and counts the resulting 31 algebraic ratio candidates.
"""

from __future__ import annotations

from itertools import combinations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive(expression: sp.Expr, *generators: sp.Symbol) -> sp.Expr:
    polynomial = sp.Poly(expression, *generators, domain=sp.QQ)
    _, integral = polynomial.clear_denoms()
    _, integral = integral.primitive()
    if integral.LC() < 0:
        integral = -integral
    return integral.as_expr()


def lift_ratio(
    polynomial: sp.Expr,
    ratio: sp.Symbol,
    denominator_coordinate: sp.Symbol,
    denominator_power: int,
    numerator_coordinate: sp.Symbol,
    numerator_power: int,
) -> sp.Expr:
    """Homogenize p(t) for t=num^numerator_power/den^denominator_power."""
    data = sp.Poly(polynomial, ratio)
    degree = data.degree()
    return sp.expand(
        sum(
            coefficient
            * denominator_coordinate ** (denominator_power * (degree - exponent))
            * numerator_coordinate ** (numerator_power * exponent)
            for (exponent,), coefficient in data.terms()
        )
    )


def main() -> None:
    u, y, bvar, cvar, dvar, wvar = sp.symbols("u y B C D W")
    ratio = sp.symbols("t")

    spectral = (
        -5878656 * wvar * y
        - 26040609 * u**3
        + 49601160 * bvar * u**2
        + 1607445 * u**2 * y**2
        - 20995200 * bvar**2 * u
        - 2857680 * bvar * u * y**2
        - 52907904 * dvar * u
        - 138915 * u * y**4
        + 777600 * bvar**2 * y**2
        + 33592320 * bvar * dvar
        - 5598720 * bvar * cvar * y
        + 78120 * bvar * y**4
        + 1959552 * dvar * y**2
        - 435456 * cvar * y**3
        + 1127 * y**6
    )
    branch = sp.factor(sp.discriminant(spectral, u))

    ratio_polynomials = {
        "BC": [
            2000 + 15309 * ratio,
            125 + 1134 * ratio,
            3321125 - 161754894 * ratio - 2812385772 * ratio**2,
            410644531250000
            - 18114791748046875 * ratio
            - 545436228093750000 * ratio**2
            - 4951165276923468750 * ratio**3
            - 18946967714644599000 * ratio**4
            - 26529827304546537363 * ratio**5,
        ],
        "BD": [
            4075 - 85176 * ratio,
            25 - 126 * ratio,
            22656250
            - 772734375 * ratio
            + 7600635000 * ratio**2
            - 30805790400 * ratio**3
            + 46376717184 * ratio**4,
        ],
        "BW": [
            5313800000
            + 4508659468656 * ratio
            - 136738899331083 * ratio**2,
            5511577600000000000000000000
            + 4983290602536960000000000000000 * ratio
            - 6564822237254419568640000000000 * ratio**2
            - 3094052863483309848285092659200000 * ratio**3
            - 81862566455344350924421142159812608 * ratio**4
            - 744088924275617882256518828471658624 * ratio**5
            - 2973811237322720333634598763466407943 * ratio**6,
        ],
        "CD": [
            22143375 + 6397664 * ratio,
            387420489
            - 8964338040 * ratio
            + 54880100352 * ratio**2
            + 16544432128 * ratio**3,
        ],
        "DW": [
            935886848 + 430565625 * ratio,
            36028797018963968
            + 17932072576352256 * ratio
            - 1448500838400000 * ratio**2
            + 56162900390625 * ratio**3,
        ],
    }
    multiplicities = {
        "BC": [1, 2, 3, 1],
        "BD": [6, 19, 2],
        "BW": [3, 1],
        "CD": [3, 1],
        "DW": [3, 1],
    }
    plane_data = {
        "BC": (
            {dvar: 0, wvar: 0},
            bvar,
            3,
            cvar,
            2,
            {primitive(bvar, bvar, cvar): 24},
        ),
        "BD": (
            {cvar: 0, wvar: 0},
            bvar,
            2,
            dvar,
            1,
            {},
        ),
        "BW": (
            {cvar: 0, dvar: 0},
            bvar,
            5,
            wvar,
            2,
            {primitive(bvar, bvar, wvar): 6},
        ),
        "CD": (
            {bvar: 0, wvar: 0},
            cvar,
            4,
            dvar,
            3,
            {primitive(dvar, cvar, dvar): 15},
        ),
        "DW": (
            {bvar: 0, cvar: 0},
            dvar,
            5,
            wvar,
            4,
            {primitive(dvar, dvar, wvar): 3},
        ),
    }

    candidate_counts: dict[str, int] = {}
    for name, (
        substitutions,
        denominator_coordinate,
        denominator_power,
        numerator_coordinate,
        numerator_power,
        expected_axes,
    ) in plane_data.items():
        specialized = sp.factor(branch.subs(substitutions))
        resultant = sp.factor(
            sp.resultant(specialized, sp.diff(specialized, y), y)
        )
        _constant, factors = sp.factor_list(resultant)
        actual = {
            primitive(factor, denominator_coordinate, numerator_coordinate): exponent
            for factor, exponent in factors
        }

        expected = dict(expected_axes)
        for polynomial, exponent in zip(
            ratio_polynomials[name],
            multiplicities[name],
            strict=True,
        ):
            lifted = lift_ratio(
                polynomial,
                ratio,
                denominator_coordinate,
                denominator_power,
                numerator_coordinate,
                numerator_power,
            )
            expected[
                primitive(
                    lifted,
                    denominator_coordinate,
                    numerator_coordinate,
                )
            ] = exponent
        require(actual == expected, f"{name} resultant factorization changed")

        polynomials = [sp.Poly(item, ratio) for item in ratio_polynomials[name]]
        for polynomial in polynomials:
            require(polynomial.eval(0) != 0, f"{name} gained a zero ratio")
            require(
                sp.gcd(polynomial, polynomial.diff()).degree() == 0,
                f"{name} ratio polynomial is not squarefree",
            )
        for first, second in combinations(polynomials, 2):
            require(
                sp.gcd(first, second).degree() == 0,
                f"{name} ratio factors overlap",
            )
        candidate_counts[name] = sum(item.degree() for item in polynomials)

    cw_branch = sp.factor(branch.subs({bvar: 0, dvar: 0}))
    require(
        sp.resultant(cw_branch, sp.diff(cw_branch, y), y) == 0,
        "C-W raw resultant should vanish identically",
    )
    require(
        candidate_counts == {"BC": 9, "BD": 6, "BW": 8, "CD": 4, "DW": 4},
        "two-sparse candidate count changed",
    )
    require(sum(candidate_counts.values()) == 31, "total candidate bank changed")

    print("normalized_spectral_cubic=THM2297")
    print("closed_coordinate_plane=CW")
    print("ratio_definitions=BC:C^2/B^3,BD:D/B^2,BW:W^2/B^5,CD:D^3/C^4,DW:W^4/D^5")
    for name in ("BC", "BD", "BW", "CD", "DW"):
        print(f"{name}_candidate_count={candidate_counts[name]}")
        for index, (polynomial, exponent) in enumerate(
            zip(ratio_polynomials[name], multiplicities[name], strict=True),
            start=1,
        ):
            print(f"{name}_factor_{index}_multiplicity={exponent}")
            print(f"{name}_ratio_polynomial_{index}={sp.expand(polynomial)}")
    print("two_sparse_candidate_ratios=31")
    print("ratio_polynomials=squarefree_pairwise_coprime_nonzero_roots")
    print("status=THM2311_DEGREE18_TWO_SPARSE_EXACT_REFEREE")


if __name__ == "__main__":
    main()
