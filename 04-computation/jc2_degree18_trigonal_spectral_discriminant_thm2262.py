#!/usr/bin/env python3
"""Exact symbolic referee for the nonsplit terminal quartic degree-18 wall.

For

    P = w^4 + 2*d*w^2 + q*w + (d^2-s),       q^2=T,
    Q = E_18 + alpha*E_14 + beta*E_10
             + gamma*E_6 + delta*E_2,

the first two constant Faber fluxes are reduced in the centered variables

    u=d*T,       Z=T^2,       y=9*s-2*alpha.

Off ``y=0`` the first flux eliminates Z and the second gives an explicit
cubic G(u,y).  This companion freezes G, its degree-twelve branch
discriminant, a squarefree specialization, the three distinct points over
weighted infinity, and the retained third flux.

The exceptional wall ``y=0`` is checked separately.  The first two fluxes
force u and the Faber parameters to be constant.  The Keller one-form then
reduces every possible rational trajectory to a monomial cusp.  Finally,
an independently generated Faber sidecar has the unique polar term
``-21*T^3/1024``; polynomiality of the reduced mate rules out that cusp.

The genus-four and rational-primitive arguments are mathematical parts of
THM-2262.  This file verifies all symbolic identities on which they rest.
"""

from __future__ import annotations

import hashlib

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def faber_coefficients(
    degree: int,
    p: sp.Expr,
    q: sp.Expr,
    r: sp.Expr,
    extra: int = 3,
) -> list[sp.Expr]:
    """Laurent coefficients of P^(degree/4) through degree+extra."""
    exponent = sp.Rational(degree, 4)
    coefficients = [sp.Integer(1)]
    quartic = {2: p, 3: q, 4: r}
    for index in range(1, degree + extra + 1):
        value = sum(
            quartic[step]
            * ((exponent + 1) * step - index)
            * coefficients[index - step]
            for step in range(2, min(4, index) + 1)
            if step in quartic
        ) / index
        coefficients.append(sp.factor(value))
    return coefficients


def faber_observables(
    degree: int,
    p: sp.Expr,
    q: sp.Expr,
    r: sp.Expr,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return the three quartic Laurent observables (Phi, Psi, R)."""
    coefficients = faber_coefficients(degree, p, q, r)
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


def replace_even_q(
    expression: sp.Expr,
    q: sp.Symbol,
    t: sp.Symbol,
) -> sp.Expr:
    """Replace every q^(2j) by T^j, rejecting an odd q power."""
    polynomial = sp.Poly(sp.expand(expression), q)
    result = sp.Integer(0)
    for (power,), coefficient in polynomial.terms():
        require(power % 2 == 0, "unexpected odd q power")
        result += coefficient * t ** (power // 2)
    return sp.factor(result)


def replace_t_square(
    expression: sp.Expr,
    t: sp.Symbol,
    z: sp.Symbol,
) -> sp.Expr:
    """Replace every even T power by the corresponding Z power."""
    polynomial = sp.Poly(sp.expand(expression), t)
    result = sp.Integer(0)
    for (power,), coefficient in polynomial.terms():
        require(power % 2 == 0, "unexpected odd T power")
        result += coefficient * z ** (power // 2)
    return sp.factor(result)


def polynomial_associate(
    found: sp.Expr,
    expected: sp.Expr,
    generators: tuple[sp.Symbol, ...],
) -> sp.Expr:
    """Return the nonzero scalar associating two polynomials."""
    found_poly = sp.Poly(sp.expand(found), *generators)
    expected_poly = sp.Poly(sp.expand(expected), *generators)
    require(not found_poly.is_zero and not expected_poly.is_zero, "zero polynomial")
    ratio = sp.factor(found_poly.LC() / expected_poly.LC())
    require(
        sp.expand(found_poly.as_expr() - ratio * expected_poly.as_expr()) == 0,
        "polynomials are not associates",
    )
    return ratio


def faber_sidecar(
    degree: int,
    lvar: sp.Symbol,
    d: sp.Symbol,
    q: sp.Symbol,
    t: sp.Symbol,
    s: sp.Symbol,
) -> sp.Expr:
    """Compute E_(4j-2)-R_j(P,H) in the linear L coordinate.

    This is independent of the three-flux calculation: it constructs the
    whole polynomial part E_m and subtracts the finite binomial truncation.
    """
    coefficients = faber_coefficients(degree, 2 * d, q, d**2 - s, extra=0)
    w = sp.symbols("w")
    e_m = sum(
        coefficients[index] * w ** (degree - index)
        for index in range(degree + 1)
    )
    in_l = sp.cancel(e_m.subs(w, (lvar + s) / q))
    numerator, denominator = sp.fraction(in_l)
    in_l = sp.factor(
        replace_even_q(numerator, q, t)
        / replace_even_q(denominator, q, t)
    )

    j = (degree + 2) // 4
    h = (lvar + s) ** 2 / t + d
    p_full = h**2 + lvar
    truncation = sum(
        sp.binomial(sp.Rational(2 * j - 1, 2), index)
        * h ** (2 * j - 1 - 2 * index)
        * lvar**index
        for index in range(j)
    )
    return sp.factor(in_l - truncation)


def canonical_hash(expression: sp.Expr) -> str:
    """Stable digest of an expanded SymPy expression."""
    payload = sp.srepr(sp.Poly(sp.expand(expression)).as_expr()).encode()
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    d, q, s, t = sp.symbols("d q s T")
    alpha, beta, gamma, delta, psi0 = sp.symbols(
        "alpha beta gamma delta Psi0"
    )
    y, u, z, v = sp.symbols("y u Z v")
    p = 2 * d
    r = d**2 - s

    bank = {
        degree: faber_observables(degree, p, q, r)
        for degree in (2, 6, 10, 14, 18)
    }
    expected_18 = (
        sp.Rational(63, 128)
        * q
        * (
            2 * q**4 * d**2
            + q**4 * s
            - 20 * q**2 * d * s**2
            + 10 * s**4
        ),
        -sp.Rational(63, 256)
        * (
            -q**6 * d
            + 20 * q**4 * d**2 * s
            + 5 * q**4 * s**2
            - 40 * q**2 * d * s**3
            + 4 * s**5
        ),
        -sp.Rational(3, 512)
        * q
        * (
            -3 * q**6
            + 84 * q**4 * d**3
            + 210 * q**4 * d * s
            - 840 * q**2 * d**2 * s**2
            - 280 * q**2 * s**3
            + 420 * d * s**4
        ),
    )
    require(
        all(
            sp.factor(found - wanted) == 0
            for found, wanted in zip(bank[18], expected_18)
        ),
        "degree-18 Laurent bank mismatch",
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

    weights = {
        18: sp.Integer(1),
        14: alpha,
        10: beta,
        6: gamma,
        2: delta,
    }
    phi = replace_even_q(
        sp.factor(sum(weights[n] * bank[n][0] for n in weights) / q),
        q,
        t,
    )
    psi = replace_even_q(
        sp.factor(sum(weights[n] * bank[n][1] for n in weights) - psi0),
        q,
        t,
    )
    third = replace_even_q(
        sp.factor(sum(weights[n] * bank[n][2] for n in weights) / q),
        q,
        t,
    )

    raw_phi = replace_t_square(sp.factor(phi.subs(d, u / t)), t, z)
    expected_raw_phi = (
        126 * u**2
        + 560 * u * alpha * s
        - 160 * u * beta
        - 1260 * u * s**2
        - 14 * z * alpha
        + 63 * z * s
        - 560 * alpha * s**3
        + 480 * beta * s**2
        - 384 * gamma * s
        + 256 * delta
        + 630 * s**4
    ) / 128
    require(
        sp.factor(raw_phi - expected_raw_phi) == 0,
        "raw first flux mismatch",
    )

    centered_s = (y + 2 * alpha) / 9
    phi_centered = sp.factor(raw_phi.subs(s, centered_s))
    solved_z = sp.factor(sp.solve(phi_centered, z)[0])
    n2 = (
        45927 * u**2
        + 22680 * u * alpha**2
        - 58320 * u * beta
        - 5670 * u * y**2
        - 1680 * alpha**4
        - 2240 * alpha**3 * y
        + 8640 * alpha**2 * beta
        - 840 * alpha**2 * y**2
        + 8640 * alpha * beta * y
        - 31104 * alpha * gamma
        + 2160 * beta * y**2
        + 93312 * delta
        - 15552 * gamma * y
        + 35 * y**4
    )
    expected_z = -2 * n2 / (5103 * y)
    require(
        sp.factor(solved_z - expected_z) == 0,
        "first-flux Z elimination mismatch",
    )

    psi_in_uz = replace_t_square(sp.factor(psi.subs(d, u / t)), t, z)
    eliminated_psi = sp.factor(
        psi_in_uz.subs(s, centered_s).subs(z, solved_z)
    )
    g = (
        -5878656 * psi0 * y
        - 26040609 * u**3
        - 19289340 * u**2 * alpha**2
        + 49601160 * u**2 * beta
        + 1607445 * u**2 * y**2
        - 2222640 * u * alpha**4
        + 11430720 * u * alpha**2 * beta
        + 1111320 * u * alpha**2 * y**2
        + 17635968 * u * alpha * gamma
        - 20995200 * u * beta**2
        - 2857680 * u * beta * y**2
        - 52907904 * u * delta
        - 138915 * u * y**4
        + 235200 * alpha**6
        + 326144 * alpha**5 * y
        - 1814400 * alpha**4 * beta
        + 82320 * alpha**4 * y**2
        - 2096640 * alpha**3 * beta * y
        + 4354560 * alpha**3 * gamma
        - 62720 * alpha**3 * y**3
        + 3110400 * alpha**2 * beta**2
        - 423360 * alpha**2 * beta * y**2
        - 13063680 * alpha**2 * delta
        + 2612736 * alpha**2 * gamma * y
        - 30380 * alpha**2 * y**4
        + 3110400 * alpha * beta**2 * y
        - 11197440 * alpha * beta * gamma
        + 241920 * alpha * beta * y**3
        - 2612736 * alpha * delta * y
        - 653184 * alpha * gamma * y**2
        + 777600 * beta**2 * y**2
        + 33592320 * beta * delta
        - 5598720 * beta * gamma * y
        + 78120 * beta * y**4
        + 1959552 * delta * y**2
        - 435456 * gamma * y**3
        + 1127 * y**6
    )
    psi_numerator = sp.together(eliminated_psi).as_numer_denom()[0]
    psi_associate = polynomial_associate(
        psi_numerator,
        g,
        (u, y, alpha, beta, gamma, delta, psi0),
    )

    discriminant = sp.factor(sp.discriminant(g, u))
    discriminant_poly = sp.Poly(discriminant, y)
    require(discriminant_poly.degree() == 12, "branch discriminant degree")

    specialization = {
        alpha: 0,
        beta: 1,
        gamma: 1,
        delta: 1,
        psi0: 1,
    }
    specialized_discriminant = sp.Poly(
        sp.expand(discriminant.subs(specialization)),
        y,
    )
    specialized_gcd = sp.gcd(
        specialized_discriminant,
        specialized_discriminant.diff(),
    )
    require(
        specialized_gcd.degree() == 0,
        "chosen branch-discriminant specialization is not squarefree",
    )

    infinity_form = sp.Poly(
        sp.expand(g.subs(u, v * y**2)),
        y,
    ).LC()
    expected_infinity = (
        1127
        - 138915 * v
        + 1607445 * v**2
        - 26040609 * v**3
    )
    require(
        sp.factor(infinity_form - expected_infinity) == 0,
        "weighted infinity cubic mismatch",
    )
    infinity_discriminant = sp.discriminant(expected_infinity, v)
    expected_infinity_discriminant = sp.Integer(
        -153384762202971019112448
    )
    require(
        infinity_discriminant == expected_infinity_discriminant,
        "weighted infinity discriminant mismatch",
    )

    # Retain the third flux before imposing the spectral cubic.
    third_in_uz = replace_t_square(
        sp.factor(t * third.subs(d, u / t)),
        t,
        z,
    ) / t
    third_centered = sp.factor(third_in_uz.subs(s, centered_s))
    n3 = (
        183708 * u**3
        + 90720 * u**2 * alpha**2
        - 233280 * u**2 * beta
        - 22680 * u**2 * y**2
        + 51030 * u * z * y
        - 6720 * u * alpha**4
        - 8960 * u * alpha**3 * y
        + 34560 * u * alpha**2 * beta
        - 3360 * u * alpha**2 * y**2
        + 34560 * u * alpha * beta * y
        - 124416 * u * alpha * gamma
        + 8640 * u * beta * y**2
        + 373248 * u * delta
        - 62208 * u * gamma * y
        + 140 * u * y**4
        - 6561 * z**2
        + 13440 * z * alpha**3
        + 10080 * z * alpha**2 * y
        - 51840 * z * alpha * beta
        - 25920 * z * beta * y
        + 93312 * z * gamma
        - 840 * z * y**3
    )
    expected_third = -n3 / (373248 * t)
    # The retained expression has one odd T in its denominator after Z=T^2.
    require(
        sp.factor(third_centered - expected_third) == 0,
        "retained third flux mismatch",
    )

    # Exceptional y=0 wall.
    exceptional_phi = sp.factor(phi_centered.subs(y, 0))
    expected_exceptional_phi = (
        15309 * u**2
        + 7560 * u * alpha**2
        - 19440 * u * beta
        - 560 * alpha**4
        + 2880 * alpha**2 * beta
        - 10368 * alpha * gamma
        + 31104 * delta
    ) / 15552
    require(
        sp.factor(exceptional_phi - expected_exceptional_phi) == 0,
        "exceptional first flux mismatch",
    )

    exceptional_psi = sp.factor(
        psi_in_uz.subs(s, 2 * alpha / 9)
    )
    z_coefficient = sp.factor(sp.diff(exceptional_psi, z))
    expected_z_coefficient = (
        567 * u + 140 * alpha**2 - 360 * beta
    ) / 2304
    require(
        sp.factor(z_coefficient - expected_z_coefficient) == 0,
        "exceptional second-flux Z coefficient mismatch",
    )
    u0 = -20 * (7 * alpha**2 - 18 * beta) / 567
    delta0 = (
        490 * alpha**4
        - 2520 * alpha**2 * beta
        + 3402 * alpha * gamma
        + 2025 * beta**2
    ) / 10206
    psi0_exceptional = (
        2
        * (14 * alpha**2 - 45 * beta)
        * (56 * alpha**3 - 225 * alpha * beta + 486 * gamma)
        / 45927
    )
    require(
        sp.factor(
            exceptional_phi.subs({u: u0, delta: delta0})
        )
        == 0,
        "exceptional delta relation mismatch",
    )
    require(
        sp.factor(
            exceptional_psi.subs({u: u0, delta: delta0})
            .subs(z, 0)
            .subs(psi0, psi0_exceptional)
        )
        == 0,
        "exceptional Psi0 relation mismatch",
    )

    k_parameter = (
        4480 * alpha**3
        - 17280 * alpha * beta
        + 31104 * gamma
    )
    exceptional_third = sp.factor(
        third.subs(s, 2 * alpha / 9)
        .subs(d, u0 / t)
        .subs(delta, delta0)
    )
    expected_exceptional_third = (
        t * (2187 * t**2 - k_parameter) / 124416
    )
    require(
        sp.factor(exceptional_third - expected_exceptional_third) == 0,
        "exceptional third flux mismatch",
    )

    # Whole-polynomial Faber sidecars, independent of the Laurent flux bank.
    lvar = sp.symbols("L")
    sidecars = {
        degree: faber_sidecar(degree, lvar, d, q, t, s)
        for degree in (2, 6, 10, 14, 18)
    }
    expected_sidecar_18 = (
        -sp.Rational(21, 1024)
        * t
        * (
            -12 * lvar**3
            + 24 * lvar**2 * s
            + lvar * (12 * t * d - 36 * s**2)
            + t**2
            - 48 * t * d * s
            + 48 * s**3
        )
    )
    require(
        sp.factor(sidecars[18] - expected_sidecar_18) == 0,
        "degree-18 Faber sidecar mismatch",
    )
    full_sidecar = sp.factor(
        sidecars[18]
        + alpha * sidecars[14]
        + beta * sidecars[10]
        + gamma * sidecars[6]
        + delta * sidecars[2]
    )
    # On the monomial cusp, T has a pole and d*T is constant.  Every term
    # besides the displayed T^3 has valuation at least -ord(T), because L
    # is polynomial; lower seeds cannot contribute T^3.
    t_polynomial = sp.Poly(sp.expand(full_sidecar), t)
    require(
        t_polynomial.degree() == 3
        and sp.factor(t_polynomial.coeff_monomial(t**3))
        == -sp.Rational(21, 1024),
        "unique exceptional T^3 sidecar pole mismatch",
    )

    print("faber_flux_bank=PASS degrees=2,6,10,14,18")
    print("independent_multinomial_bank=PASS degrees=2,6,10,14,18")
    print("centered_coordinate=y=9*s-2*alpha")
    print(f"first_flux_Z={sp.factor(expected_z)}")
    print(f"spectral_G_associate={psi_associate}")
    print(f"spectral_G_sha256={canonical_hash(g)}")
    print(f"branch_discriminant_degree={discriminant_poly.degree()}")
    print(f"branch_discriminant_lc={discriminant_poly.LC()}")
    print(f"branch_discriminant_sha256={canonical_hash(discriminant)}")
    print("squarefree_specialization=(alpha,beta,gamma,delta,Psi0)=(0,1,1,1,1)")
    print(f"squarefree_specialization_gcd_degree={specialized_gcd.degree()}")
    print(f"weighted_infinity_cubic={expected_infinity}")
    print(f"weighted_infinity_discriminant={infinity_discriminant}")
    print(f"exceptional_u={sp.factor(u0)}")
    print(f"exceptional_delta={sp.factor(delta0)}")
    print(f"exceptional_Psi0={sp.factor(psi0_exceptional)}")
    print(f"exceptional_third_over_q={sp.factor(expected_exceptional_third)}")
    print("exceptional_keller_primitive=(1701*T^3-K*T)'")
    print("exceptional_sidecar_T3_coefficient=-21/1024")
    print("squarefree_discriminant_implies_genus_four=MATHEMATICAL_PROOF_REQUIRED")
    print("monomial_cusp_sidecar_pole=MATHEMATICAL_PROOF_REQUIRED")
    print("status=THM2262_DEGREE18_TRIGONAL_SPECTRAL_EXACT_REFEREE")


if __name__ == "__main__":
    main()
