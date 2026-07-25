#!/usr/bin/env python3
"""Exact symbolic referee for the nonsplit terminal quartic degree-14 closure.

THM-2245 reduces every surviving degree-fourteen trajectory to

    H(y) = 425 (y-e)^2 ((y+e)^2 + D)

together with a completed-square equation for ``T^2`` and a retained Keller
one-form.  This script checks the two normalization branches used by
THM-2247:

* for ``D != 0``, conic parameterization turns the remaining double cover
  into ``T^2 = 1372*S6(u)/(235298*u^3)``; its exact sextic cannot be a
  scalar square;
* for ``D = 0``, both reducible-conic components give an exact cubic.  The
  only nongeneric cusp has a seventh-power Keller primitive and an
  unavoidable pole in the polynomial approximate root.

The branch-count, Riemann--Hurwitz, rational-primitive, and pole arguments
are mathematical parts of the theorem rather than delegated assertions.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def reduce_c(expression: sp.Expr, c: sp.Symbol) -> sp.Expr:
    """Reduce a polynomial/rational expression using c^2=425."""
    numerator, denominator = sp.fraction(sp.cancel(expression))
    require(
        not denominator.has(c),
        "reduce_c received a denominator involving c",
    )
    remainder = sp.rem(
        sp.Poly(numerator, c),
        sp.Poly(c**2 - 425, c),
    ).as_expr()
    return sp.factor(remainder / denominator)


def main() -> None:
    alpha, beta, gamma, psi0 = sp.symbols(
        "alpha beta gamma Psi0"
    )
    y, e, D, u, c = sp.symbols("y e D u c")

    spectral_h = (
        425 * y**4
        + (840 * beta - 300 * alpha**2) * y**2
        + (
            -1600 * alpha**3
            + 6720 * alpha * beta
            - 15680 * gamma
        )
        * y
        + 13720 * psi0
        + 1200 * alpha**4
        - 6720 * alpha**2 * beta
        + 7840 * alpha * gamma
        + 7056 * beta**2
    )
    singular_h = 425 * (y - e) ** 2 * ((y + e) ** 2 + D)

    beta_singular = sp.factor(
        5 * (17 * D + 12 * alpha**2 - 34 * e**2) / 168
    )
    gamma_singular = sp.factor(
        5
        * (
            68 * D * alpha
            + 17 * D * e
            + 16 * alpha**3
            - 136 * alpha * e**2
        )
        / 1568
    )
    psi_singular = sp.factor(
        -5
        * (
            289 * D**2
            + 136 * D * alpha**2
            + 68 * D * alpha * e
            - 1224 * D * e**2
            + 16 * alpha**4
            - 272 * alpha**2 * e**2
            + 1088 * e**4
        )
        / 10976
    )
    singular_substitution = {
        beta: beta_singular,
        gamma: gamma_singular,
        psi0: psi_singular,
    }
    require(
        sp.factor(
            spectral_h.subs(singular_substitution) - singular_h
        )
        == 0,
        "singular quartic coefficient matching failed",
    )

    spectral_b = (
        219520 * y**3
        - 439040 * alpha**3
        + 1843968 * alpha * beta
        - 4302592 * gamma
    )
    b_singular = sp.factor(
        spectral_b.subs(singular_substitution)
    )
    expected_b = 13720 * (16 * y**3 - 17 * D * e)
    require(
        sp.factor(b_singular - expected_b) == 0,
        "singular spectral B simplification failed",
    )

    # Smooth conic D != 0.  With
    #
    #   u=(y+e)+v,
    #   y=(u-D/u)/2-e,
    #   v=(u+D/u)/2,
    #
    # and W=c(y-e)v, the completed square recovers T^2.
    conic_y = (u - D / u) / 2 - e
    conic_v = (u + D / u) / 2
    cover_numerator = sp.together(
        10976 * c * conic_y * (conic_y - e) * conic_v
        - expected_b.subs(y, conic_y)
    )
    sextic = (
        (c - 20) * u**6
        - 6 * e * (c - 20) * u**5
        + (
            -D * c
            + 60 * D
            + 8 * c * e**2
            - 240 * e**2
        )
        * u**4
        + 10 * e * (-7 * D + 16 * e**2) * u**3
        + D
        * (
            -D * c
            - 60 * D
            + 8 * c * e**2
            + 240 * e**2
        )
        * u**2
        + 6 * D**2 * e * (c + 20) * u
        + D**3 * (c + 20)
    )
    require(
        sp.factor(cover_numerator * u**3 - 1372 * sextic) == 0,
        "conic parameterization did not produce the frozen sextic",
    )
    require(
        reduce_c((c - 20) * (c + 20) - 25, c) == 0,
        "leading/trailing coefficient norm mismatch",
    )

    # If the only branch points are u=0 and u=infinity, every root of the
    # sextic has even multiplicity, so it is a scalar cubic square.  Match
    # the top four coefficients; the u^2 and u residuals are already
    # inconsistent.
    square_b = sp.factor(
        (
            8 * D * c
            + 155 * D
            - 16 * c * e**2
            - 325 * e**2
        )
        / 10
    )
    square_k = sp.factor(
        e
        * (
            10 * D * c
            + 185 * D
            - 16 * c * e**2
            - 335 * e**2
        )
        / 10
    )
    candidate_cubic = (
        u**3 - 3 * e * u**2 + square_b * u + square_k
    )
    square_error = sp.Poly(
        sp.expand(
            reduce_c(
                (c - 20) * candidate_cubic**2 - sextic,
                c,
            )
        ),
        u,
    )
    require(
        all(
            square_error.coeff_monomial(u**degree) == 0
            for degree in (6, 5, 4, 3)
        ),
        "top square-coefficient matching failed",
    )
    residual_u2 = (
        69 * D**2 * c
        + 1420 * D**2
        - 250 * D * c * e**2
        - 7360 * D * e**2
        + 293 * c * e**4
        + 5500 * e**4
    ) / 4
    residual_u1 = e * (
        71 * D**2 * c
        + 1130 * D**2
        - 282 * D * c * e**2
        - 5730 * D * e**2
        + 259 * c * e**4
        + 5380 * e**4
    ) / 2
    require(
        reduce_c(
            square_error.coeff_monomial(u**2) - residual_u2,
            c,
        )
        == 0
        and reduce_c(
            square_error.coeff_monomial(u) - residual_u1,
            c,
        )
        == 0,
        "low square residual mismatch",
    )

    z = sp.symbols("z")
    residual_z2 = (
        (69 * c + 1420) * z**2
        + (-250 * c - 7360) * z
        + 293 * c
        + 5500
    )
    residual_z1 = (
        (71 * c + 1130) * z**2
        + (-282 * c - 5730) * z
        + 259 * c
        + 5380
    )
    residual_resultant = reduce_c(
        sp.resultant(residual_z2, residual_z1, z),
        c,
    )
    expected_resultant = -8670000 * (110 * c + 2281)
    require(
        sp.factor(residual_resultant - expected_resultant) == 0,
        "square-residual resultant mismatch",
    )
    norm_e_zero = 1420**2 - 425 * 69**2
    norm_e_nonzero = 2281**2 - 425 * 110**2
    require(
        norm_e_zero == -7025 and norm_e_nonzero == 60461,
        "nonzero norm controls changed",
    )

    # Reducible conic D=0.  Check both v=sigma(y+e) components.  The
    # generic e != 0 cover is a squarefree cubic.  At e=0 it is the cusp
    # parameterized by y=rho^2, T=lambda*rho^3.
    T, rho, lam, rhop = sp.symbols("T rho lambda rho_prime")
    s_coordinate = (y + 2 * alpha) / 7
    dT_first_flux = (
        343 * T**2
        - 640 * alpha**3
        - 480 * alpha**2 * y
        + 2688 * alpha * beta
        + 1344 * beta * y
        - 6272 * gamma
        + 40 * y**3
    ) / (1960 * y)
    third_flux = -T * (
        -343 * T**2
        + 640 * alpha**3
        - 320 * alpha**2 * y
        - 2688 * alpha * beta
        + 896 * beta * y
        + 6272 * gamma
        + 160 * y**3
    ) / (8960 * y)

    cubic_records: list[tuple[int, sp.Expr]] = []
    cusp_records: list[tuple[int, sp.Expr, sp.Expr, sp.Expr]] = []
    for sigma in (-1, 1):
        r = sigma * c
        cubic = y * ((r - 20) * y**2 - r * e**2)
        component_cover = sp.factor(
            (
                10976
                * c
                * y
                * (y - e)
                * sigma
                * (y + e)
                - expected_b.subs(D, 0)
            )
            / 235298
        )
        require(
            reduce_c(
                component_cover - sp.Rational(16, 343) * cubic,
                c,
            )
            == 0,
            f"sigma={sigma} reducible component cover mismatch",
        )
        cubic_discriminant = reduce_c(
            sp.discriminant(cubic, y),
            c,
        )
        expected_discriminant = reduce_c(
            4 * (r - 20) * r**3 * e**6,
            c,
        )
        require(
            sp.factor(
                cubic_discriminant - expected_discriminant
            )
            == 0,
            f"sigma={sigma} cubic discriminant mismatch",
        )
        cubic_records.append((sigma, cubic_discriminant))

        cusp_t2 = 16 * (r - 20) * y**3 / 343
        cusp_lambda2 = 16 * (r - 20) / 343
        cusp_substitution = {
            D: 0,
            e: 0,
            beta: sp.Rational(5, 14) * alpha**2,
            gamma: sp.Rational(5, 98) * alpha**3,
            T**2: cusp_t2,
        }
        cusp_f = reduce_c(
            third_flux.subs(cusp_substitution),
            c,
        )
        expected_cusp_f = (r - 30) * T * y**2 / 560
        require(
            reduce_c(cusp_f - expected_cusp_f, c) == 0,
            f"sigma={sigma} cusp third flux mismatch",
        )
        cusp_dT = reduce_c(
            dT_first_flux.subs(cusp_substitution),
            c,
        )
        expected_cusp_dT = (2 * r - 35) * y**2 / 245
        require(
            reduce_c(cusp_dT - expected_cusp_dT, c) == 0,
            f"sigma={sigma} cusp dT mismatch",
        )

        # Formal derivative check for the divided Keller one-form.
        cusp_f_rho = expected_cusp_f.subs(
            {T: lam * rho**3, y: rho**2}
        )
        cusp_t_rho = lam * rho**3
        divided_one_form = sp.factor(
            2 * sp.diff(cusp_f_rho, rho) * rhop
            + cusp_f_rho
            * sp.diff(cusp_t_rho, rho)
            * rhop
            / cusp_t_rho
        )
        expected_one_form = (
            17
            * (r - 30)
            * lam
            * sp.diff(rho**7, rho)
            * rhop
            / (560 * 7)
        )
        require(
            reduce_c(
                divided_one_form - expected_one_form,
                c,
            )
            == 0,
            f"sigma={sigma} seventh-power one-form mismatch",
        )

        # Constant coefficient of the polynomial approximate root
        # a(L+s)^2+d.  Terms involving the polynomial E are regular at
        # the forced pole; this is the exact coefficient of rho in
        # a*s^2+d.
        cusp_s = s_coordinate.subs(y, rho**2)
        cusp_a = 1 / (lam * rho**3)
        cusp_d = (
            (2 * r - 35) * rho / (245 * lam)
        )
        approximate_root_constant = sp.expand(
            cusp_a * cusp_s**2 + cusp_d
        )
        polar_coefficient = sp.factor(
            approximate_root_constant.coeff(rho, 1)
        )
        expected_polar = 2 * (r - 15) / (245 * lam)
        require(
            reduce_c(polar_coefficient - expected_polar, c)
            == 0,
            f"sigma={sigma} terminal pole coefficient mismatch",
        )
        cusp_records.append(
            (
                sigma,
                cusp_lambda2,
                expected_cusp_f,
                expected_polar,
            )
        )

    nonzero_norms = {
        "c_minus_20": 20**2 - 425,
        "c_minus_15": 15**2 - 425,
        "c_minus_30": 30**2 - 425,
    }
    require(
        nonzero_norms
        == {
            "c_minus_20": -25,
            "c_minus_15": -200,
            "c_minus_30": 475,
        },
        "hostile component norm controls changed",
    )

    # Explicitly replay both choices c=+/-sqrt(425) and both reducible
    # components sigma=+/-1.  Their product r=sigma*c encounters both
    # algebraic signs twice, and every load-bearing coefficient is nonzero.
    hostile_sign_records: list[tuple[int, int, bool]] = []
    for c_sign in (-1, 1):
        for sigma in (-1, 1):
            r_exact = sigma * c_sign * 5 * sp.sqrt(17)
            all_nonzero = all(
                sp.simplify(r_exact - value) != 0
                for value in (15, 20, 30)
            )
            require(
                all_nonzero,
                "hostile sign replay hit a forbidden coefficient",
            )
            hostile_sign_records.append(
                (c_sign, sigma, all_nonzero)
            )

    print("singular_coefficient_matching=PASS")
    print(f"beta={beta_singular}")
    print(f"gamma={gamma_singular}")
    print(f"Psi0={psi_singular}")
    print(f"spectral_B={expected_b}")
    print(f"smooth_conic_sextic={sextic}")
    print(
        "smooth_conic_cover="
        "T^2=1372*S6(u)/(235298*u^3)"
    )
    print(f"square_residual_u2={sp.factor(residual_u2)}")
    print(f"square_residual_u1={sp.factor(residual_u1)}")
    print(f"square_residual_resultant={expected_resultant}")
    print(
        "square_obstruction_norms="
        f"e_zero:{norm_e_zero},e_nonzero:{norm_e_nonzero}"
    )
    for sigma, discriminant in cubic_records:
        print(
            f"reducible_component_sigma={sigma} "
            f"cubic_discriminant={sp.factor(discriminant)}"
        )
    for sigma, lambda2, cusp_f, polar in cusp_records:
        print(
            f"cusp_sigma={sigma} lambda^2={lambda2} "
            f"F={cusp_f} polar_coefficient={polar}"
        )
    print(f"hostile_nonzero_norms={nonzero_norms}")
    print(f"hostile_sign_replays={hostile_sign_records}")
    print("branch_count_and_riemann_hurwitz=MATHEMATICAL_PROOF_REQUIRED")
    print("rational_primitive_and_pole_endgame=MATHEMATICAL_PROOF_REQUIRED")
    print("status=THM2247_NONSPLIT_TERMINAL_DEGREE14_EXACT_REFEREE")


if __name__ == "__main__":
    main()
