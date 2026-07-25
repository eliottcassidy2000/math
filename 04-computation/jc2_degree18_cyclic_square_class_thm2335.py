#!/usr/bin/env python3
"""Exact referee for THM-2335's cyclic square-class elimination.

For the structured Mordell polynomial F=4P^3+49Q^2 from THM-2332, the
script recursively reconstructs the unique possible degree-six square
root, verifies the six residual coefficient equations, and checks every
exact elimination used to prove that F is a square only at the affine
apex B=C=D=W=0.  Every executable check remains active under optimized
Python.
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
    mordell_poly = sp.Poly(mordell, y)
    require(
        mordell_poly.degree() == 12
        and mordell_poly.LC() == 73060029
        and 73060029 == 69 * 1029**2,
        "the leading square-root datum changed",
    )

    # If mordell=S^2 over C, choose r^2=69 and the sign of S so that
    # S=r*T with LC(T)=1029.  Coefficients y^11,...,y^6 determine T
    # uniquely and successively.
    expected_coefficients = {
        6: sp.Integer(1029),
        5: sp.Integer(0),
        4: sp.Rational(317520, 23) * bvar,
        3: sp.Rational(1571724, 23) * cvar,
        2: -sp.Rational(20412, 529)
        * (1825 * bvar**2 - 12558 * dvar),
        1: sp.Rational(91854, 529)
        * (8060 * bvar * cvar + 5313 * wvar),
        0: -sp.Rational(36741600, 12167)
        * (65 * bvar**3 - 69 * bvar * dvar - 3105 * cvar**2),
    }
    square_root = expected_coefficients[6] * y**6
    for power in range(11, 5, -1):
        index = power - 6
        trial = sp.symbols(f"trial_{index}")
        coefficient = sp.Poly(
            sp.expand(mordell - 69 * (square_root + trial * y**index) ** 2),
            y,
        ).coeff_monomial(y**power)
        solutions = sp.solve(coefficient, trial)
        require(
            len(solutions) == 1
            and sp.cancel(solutions[0] - expected_coefficients[index]) == 0,
            f"square-root recursion changed at y^{power}",
        )
        square_root = sp.expand(
            square_root + expected_coefficients[index] * y**index
        )

    residual = sp.Poly(sp.expand(mordell - 69 * square_root**2), y)
    require(
        all(residual.coeff_monomial(y**power) == 0 for power in range(6, 13)),
        "the reconstructed square root did not cancel degrees 6 through 12",
    )

    e5 = (
        1560 * bvar**2 * cvar
        + 805 * bvar * wvar
        - 9016 * cvar * dvar
    )
    e4 = (
        617375 * bvar**4
        - 6609510 * bvar**2 * dvar
        + 7172550 * bvar * cvar**2
        + 22495725 * cvar * wvar
        + 21595896 * dvar**2
    )
    e3 = (
        233080 * bvar**3 * cvar
        + 192395 * bvar**2 * wvar
        - 1302168 * bvar * cvar * dvar
        - 637560 * cvar**3
        - 1088682 * dvar * wvar
    )
    e2 = (
        1199915000 * bvar**5
        - 12521098800 * bvar**3 * dvar
        + 23960664000 * bvar**2 * cvar**2
        + 41362139700 * bvar * cvar * wvar
        + 34086322368 * bvar * dvar**2
        - 65507551200 * cvar**2 * dvar
        + 25352682075 * wvar**2
    )
    e1 = (
        8060 * bvar * cvar + 5313 * wvar
    ) * (
        65 * bvar**3
        - 69 * bvar * dvar
        - 3105 * cvar**2
    )
    e0 = (
        -105225921875 * bvar**6
        + 1530475458750 * bvar**4 * dvar
        + 445024125000 * bvar**3 * cvar**2
        - 7669002612600 * bvar**2 * dvar**2
        - 472410225000 * bvar * cvar**2 * dvar
        - 10629230062500 * cvar**4
        + 12875106064968 * dvar**3
    )
    residual_equations = (e5, e4, e3, e2, e1, e0)
    residual_contents = (
        sp.Rational(656223346800, 529),
        sp.Rational(19443654720, 12167),
        sp.Rational(1687431463200, 12167),
        sp.Rational(2678462640, 279841),
        sp.Rational(20249177558400, 279841),
        sp.Rational(3673320192, 6436343),
    )
    for power, equation, content in zip(
        range(5, -1, -1),
        residual_equations,
        residual_contents,
    ):
        require(
            sp.expand(residual.coeff_monomial(y**power) - content * equation)
            == 0,
            f"primitive residual equation E{power} changed",
        )

    # B=0: E5 gives CD=0.  Each branch then collapses successively to
    # C=D=W=0.
    require(
        e5.subs(bvar, 0) == -9016 * cvar * dvar,
        "B=0 first split changed",
    )
    require(
        e4.subs({bvar: 0, cvar: 0}) == 21595896 * dvar**2,
        "B=C=0 collapse changed",
    )
    require(
        e2.subs({bvar: 0, cvar: 0, dvar: 0})
        == 25352682075 * wvar**2,
        "B=C=D=0 collapse changed",
    )
    require(
        e4.subs({bvar: 0, dvar: 0}) == 22495725 * cvar * wvar
        and e3.subs({bvar: 0, dvar: 0, wvar: 0})
        == -637560 * cvar**3,
        "B=D=0 collapse changed",
    )

    # B!=0 is normalized to B=1 by the weighted C* action.  E5 then
    # determines W.
    normalized_w = sp.Rational(8, 805) * cvar * (1127 * dvar - 195)
    require(
        sp.solve(e5.subs(bvar, 1), wvar) == [normalized_w],
        "normalized E5 no longer determines W as stated",
    )
    normalized = {bvar: 1, wvar: normalized_w}
    e4n = sp.factor(e4.subs(normalized))
    e3n = sp.factor(e3.subs(normalized))
    e2n = sp.factor(e2.subs(normalized))
    require(
        e4n
        == (
            251952120 * cvar**2 * dvar
            - 36421650 * cvar**2
            + 21595896 * dvar**2
            - 6609510 * dvar
            + 617375
        ),
        "normalized E4 changed",
    )
    g3 = (
        398475 * cvar**2
        + 7620774 * dvar**2
        - 1851500 * dvar
        + 87350
    )
    require(
        e3n == -sp.Rational(8, 5) * cvar * g3,
        "normalized E3 changed",
    )

    # If C=0, E5 also gives W=0.  E4 and E2 are two quadratics in D
    # with nonzero resultant.
    bd_e4 = sp.Poly(e4n.subs(cvar, 0), dvar)
    bd_e2 = sp.Poly(e2n.subs(cvar, 0), dvar)
    bd_resultant = sp.resultant(bd_e4.as_expr(), bd_e2.as_expr(), dvar)
    expected_bd_resultant = 14658253353277217122876416000000
    require(
        bd_resultant == expected_bd_resultant
        and sp.gcd(bd_e4, bd_e2) == 1,
        "the normalized C=W=0 resultant changed",
    )

    # If C!=0, E1 splits into two branches after E5 is imposed.
    first_e1_factor = sp.factor(
        (8060 * bvar * cvar + 5313 * wvar).subs(normalized)
    )
    second_e1_factor = sp.factor(
        (
            65 * bvar**3
            - 69 * bvar * dvar
            - 3105 * cvar**2
        ).subs(bvar, 1)
    )
    require(
        first_e1_factor
        == sp.Rational(4, 5) * cvar * (74382 * dvar - 2795)
        and second_e1_factor
        == 65 - 69 * dvar - 3105 * cvar**2,
        "the normalized E1 branch split changed",
    )

    # First-factor branch: D is fixed.  E4 and E3 then demand two
    # different values of C^2.
    branch_a_d = sp.Rational(2795, 74382)
    c_square = sp.symbols("C_square")
    branch_a_c4 = sp.solve(
        e4n.subs({dvar: branch_a_d, cvar**2: c_square}),
        c_square,
    )
    branch_a_c3 = sp.solve(
        g3.subs({dvar: branch_a_d, cvar**2: c_square}),
        c_square,
    )
    require(
        branch_a_c4 == [sp.Rational(14418035, 972766179)]
        and branch_a_c3 == [-sp.Rational(28025, 391314)],
        "first E1 branch values changed",
    )
    branch_a_gap = sp.factor(branch_a_c4[0] - branch_a_c3[0])
    require(
        branch_a_gap == sp.Rational(20348640145, 235409415318),
        "first E1 branch contradiction vanished",
    )

    # Second-factor branch: substitute C^2=(65-69D)/3105 in E4 and
    # E3/C.  The resulting quadratics have nonzero resultant.
    branch_z_c_square = (65 - 69 * dvar) / 3105
    branch_z_e4 = sp.factor(e4n.subs(cvar**2, branch_z_c_square))
    branch_z_g3 = sp.factor(g3.subs(cvar**2, branch_z_c_square))
    require(
        sp.expand(
            branch_z_e4
            - 5 * (3199392 * dvar**2 - 105156 * dvar - 29015)
        )
        == 0
        and sp.expand(
            branch_z_g3
            - sp.Rational(1, 3)
            * (22862322 * dvar**2 - 5581065 * dvar + 287075)
        )
        == 0,
        "second E1 branch quadratics changed",
    )
    branch_z_f = 3199392 * dvar**2 - 105156 * dvar - 29015
    branch_z_g = 22862322 * dvar**2 - 5581065 * dvar + 287075
    branch_z_resultant = sp.resultant(branch_z_f, branch_z_g, dvar)
    expected_branch_z_resultant = -466513778248576599586500
    require(
        branch_z_resultant == expected_branch_z_resultant
        and sp.gcd(branch_z_f, branch_z_g) == 1,
        "the second E1 branch resultant changed",
    )

    # Positive and hostile controls.  The apex really is a square,
    # while the B-axis already fails E4.
    apex = {bvar: 0, cvar: 0, dvar: 0, wvar: 0}
    require(
        sp.expand(mordell.subs(apex) - 69 * (1029 * y**6) ** 2) == 0,
        "affine-apex positive control stopped being a square",
    )
    require(
        e4.subs({bvar: 1, cvar: 0, dvar: 0, wvar: 0}) == 617375,
        "B-axis hostile control stopped rejecting",
    )

    print("THM-2335 exact cyclic square-class referee")
    print("square-root recursion: degrees 11..6 uniquely reconstructed")
    print("residual equations: E5,E4,E3,E2,E1,E0 verified exactly")
    print(
        "B=0 branch: only (C,D,W)=(0,0,0); "
        f"B!=0,C=0 resultant={bd_resultant}"
    )
    print(
        "B!=0,C!=0 first branch gap="
        f"{branch_a_gap.p}/{branch_a_gap.q}"
    )
    print(
        "B!=0,C!=0 second branch resultant="
        f"{branch_z_resultant}"
    )
    print("positive control: affine apex is a square")
    print("hostile control: normalized B-axis has E4=617375")
    print("VERDICT: 4P^3+49Q^2 is a square iff B=C=D=W=0")


if __name__ == "__main__":
    main()
