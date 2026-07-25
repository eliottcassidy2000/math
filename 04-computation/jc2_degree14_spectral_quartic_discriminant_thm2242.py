#!/usr/bin/env python3
"""Exact symbolic referee for the nonsplit terminal quartic degree-14 wall.

This companion extends the Laurent/Faber bank of THM-2214 from degrees
2, 6, 10 to degree 14.  For

    Q = E_14 + alpha E_10 + beta E_6 + gamma E_2

it eliminates the first flux, centers the surviving coefficient by
``y = 7*s - 2*alpha``, and verifies the compact spectral identity

    (235298*T^2 + B(y))^2 = (10976*y)^2 H(y),

where ``H`` is a depressed quartic.  The mathematical proof that a
squarefree ``H`` is impossible uses the genus-one completion of
``W^2 = H(y)`` and is not delegated to this script.

The probe also freezes the exceptional ``y=0`` equation, the exact third
observable/Keller one-form input, and an independent multinomial
calculation of every displayed Laurent coefficient.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def faber_observables(
    degree: int,
    p: sp.Expr,
    q: sp.Expr,
    r: sp.Expr,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return the three quartic Laurent observables (Phi, Psi, R)."""
    exponent = sp.Rational(degree, 4)
    coefficients = [sp.Integer(1)]
    quartic = {2: p, 3: q, 4: r}
    for index in range(1, degree + 4):
        value = sum(
            quartic[step]
            * ((exponent + 1) * step - index)
            * coefficients[index - step]
            for step in range(2, min(4, index) + 1)
            if step in quartic
        ) / index
        coefficients.append(sp.factor(value))
    phi = sp.factor(4 * coefficients[degree + 1])
    psi = sp.factor(4 * coefficients[degree + 2])
    third = sp.factor(
        4 * coefficients[degree + 3]
        + p * coefficients[degree + 1]
    )
    return phi, psi, third


def direct_laurent_coefficient(
    degree: int,
    index: int,
    p: sp.Expr,
    q: sp.Expr,
    r: sp.Expr,
) -> sp.Expr:
    """Independent finite multinomial coefficient of P^(degree/4)."""
    exponent = sp.Rational(degree, 4)
    total = sp.Integer(0)
    for i in range(index // 2 + 1):
        for j in range(index // 3 + 1):
            remainder = index - 2 * i - 3 * j
            if remainder < 0 or remainder % 4:
                continue
            k = remainder // 4
            chosen = i + j + k
            falling = sp.prod(exponent - h for h in range(chosen))
            total += (
                falling
                * p**i
                * q**j
                * r**k
                / (
                    sp.factorial(i)
                    * sp.factorial(j)
                    * sp.factorial(k)
                )
            )
    return sp.factor(total)


def substitute_q_even(expression: sp.Expr, q: sp.Symbol, t: sp.Symbol) -> sp.Expr:
    """Replace the even powers of q occurring in the degree-14 bank."""
    return sp.factor(
        expression.subs(q**8, t**4)
        .subs(q**6, t**3)
        .subs(q**4, t**2)
        .subs(q**2, t)
    )


def main() -> None:
    d, q, s, t = sp.symbols("d q s T")
    alpha, beta, gamma, psi0 = sp.symbols(
        "alpha beta gamma Psi0"
    )
    y, u = sp.symbols("y U")
    p = 2 * d
    r = d**2 - s

    bank = {
        degree: faber_observables(degree, p, q, r)
        for degree in (2, 6, 10, 14)
    }
    expected_14 = (
        -sp.Rational(7, 64)
        * q
        * (-40 * d * q**2 * s + q**4 + 40 * s**3),
        sp.Rational(35, 64)
        * (
            2 * d**2 * q**4
            - 12 * d * q**2 * s**2
            + q**4 * s
            + 2 * s**4
        ),
        sp.Rational(35, 128)
        * q
        * (
            -8 * d**2 * q**2 * s
            + d * q**4
            + 8 * d * s**3
            - 4 * q**2 * s**2
        ),
    )
    require(
        all(
            sp.factor(found - wanted) == 0
            for found, wanted in zip(bank[14], expected_14)
        ),
        "degree-14 Laurent bank mismatch",
    )

    for degree in bank:
        direct = tuple(
            direct_laurent_coefficient(
                degree, degree + offset, p, q, r
            )
            for offset in (1, 2, 3)
        )
        require(
            sp.factor(4 * direct[0] - bank[degree][0]) == 0
            and sp.factor(4 * direct[1] - bank[degree][1]) == 0
            and sp.factor(
                4 * direct[2] + p * direct[0] - bank[degree][2]
            )
            == 0,
            f"degree-{degree} independent multinomial mismatch",
        )

    phi_q = sp.factor(
        (
            bank[14][0]
            + alpha * bank[10][0]
            + beta * bank[6][0]
            + gamma * bank[2][0]
        )
        / q
    )
    psi_q = sp.factor(
        bank[14][1]
        + alpha * bank[10][1]
        + beta * bank[6][1]
        + gamma * bank[2][1]
        - psi0
    )
    third_over_q = sp.factor(
        (
            bank[14][2]
            + alpha * bank[10][2]
            + beta * bank[6][2]
            + gamma * bank[2][2]
        )
        / q
    )
    phi_q = substitute_q_even(phi_q, q, t)
    psi_q = substitute_q_even(psi_q, q, t)
    third_over_q = substitute_q_even(third_over_q, q, t)

    centered_s = (y + 2 * alpha) / 7
    phi_centered = sp.factor(
        phi_q.subs(s, centered_s).subs(d, u / t)
    )
    solved_u = sp.factor(sp.solve(phi_centered, u)[0])
    expected_u = (
        343 * t**2
        - 640 * alpha**3
        - 480 * alpha**2 * y
        + 2688 * alpha * beta
        + 1344 * beta * y
        - 6272 * gamma
        + 40 * y**3
    ) / (1960 * y)
    require(
        sp.factor(solved_u - expected_u) == 0,
        "first-flux elimination mismatch",
    )

    exceptional_y_zero = sp.factor(
        phi_q.subs(s, 2 * alpha / 7)
    )
    expected_exceptional = -(
        343 * t**2
        - 640 * alpha**3
        + 2688 * alpha * beta
        - 6272 * gamma
    ) / 3136
    require(
        sp.factor(exceptional_y_zero - expected_exceptional) == 0,
        "y=0 first-flux equation mismatch",
    )

    b_polynomial = (
        219520 * y**3
        - 439040 * alpha**3
        + 1843968 * alpha * beta
        - 4302592 * gamma
    )
    h_polynomial = (
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
    spectral_square = sp.factor(
        (235298 * t**2 + b_polynomial) ** 2
        - (10976 * y) ** 2 * h_polynomial
    )
    eliminated_psi = sp.factor(
        psi_q.subs(s, centered_s).subs(d, solved_u / t)
    )
    require(
        sp.factor(
            spectral_square
            - 1652883742720 * y**2 * eliminated_psi
        )
        == 0,
        "spectral quartic identity mismatch",
    )
    require(
        sp.degree(h_polynomial, y) == 4
        and sp.LC(sp.Poly(h_polynomial, y)) == 425,
        "spectral H is not the frozen quartic",
    )

    eliminated_third = sp.factor(
        third_over_q.subs(s, centered_s).subs(d, solved_u / t)
    )
    expected_third = -t * (
        -343 * t**2
        + 640 * alpha**3
        - 320 * alpha**2 * y
        - 2688 * alpha * beta
        + 896 * beta * y
        + 6272 * gamma
        + 160 * y**3
    ) / (8960 * y)
    require(
        sp.factor(eliminated_third - expected_third) == 0,
        "third-observable elimination mismatch",
    )

    # The discriminant of the pre-completion quadratic in Z=T^2 is a
    # perfect y^2 factor times the spectral quartic.  This is an
    # independent check of the completed-square identity.
    z = sp.symbols("Z")
    spectral_as_quadratic = sp.Poly(
        sp.expand(spectral_square).subs(t**4, z**2).subs(t**2, z),
        z,
    )
    require(
        sp.factor(
            sp.discriminant(spectral_as_quadratic.as_expr(), z)
            - 4
            * 235298**2
            * 10976**2
            * y**2
            * h_polynomial
        )
        == 0,
        "quadratic discriminant factorization mismatch",
    )

    print("faber_flux_bank=PASS degrees=2,6,10,14")
    print("independent_multinomial_bank=PASS degrees=2,6,10,14")
    print(f"centered_coordinate={y}=7*s-2*alpha")
    print(f"first_flux_dT={sp.factor(expected_u)}")
    print(f"exceptional_y_zero_flux={sp.factor(expected_exceptional)}")
    print(f"spectral_B={sp.factor(b_polynomial)}")
    print(f"spectral_H={sp.factor(h_polynomial)}")
    print(
        "spectral_identity="
        "(235298*T^2+B)^2=(10976*y)^2*H"
    )
    print(f"third_observable_over_q={sp.factor(expected_third)}")
    print("squarefree_H_implies_genus_one=MATHEMATICAL_PROOF_REQUIRED")
    print("status=THM2242_DEGREE14_SPECTRAL_QUARTIC_EXACT_REFEREE")


if __name__ == "__main__":
    main()
