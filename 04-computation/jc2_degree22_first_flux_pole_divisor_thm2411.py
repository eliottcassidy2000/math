#!/usr/bin/env python3
"""Exact referee for the degree-22 first-flux pole divisor.

In the genuine nonsplit exact-square-prefix quartic coordinates

    P = w^4 + 2*d*w^2 + q*w + (d^2-s),      q^2=T,

the reduced degree-22 mate is normalized by the legal target translation
``P -> P + 2*alpha/11``.  Its remaining constant parameters are denoted
``B,C,D,E,W``.  With

    u=d*T,       Z=T^2,       y=11*s,

the first two constant Faber fluxes are the explicit polynomials ``N1,N2``.
On the first-flux pole divisor

    A = 616*B - 1089*u + 63*y^2 = 0,

the residual first equation solves for ``E`` and the second becomes a
quadratic ``F2(Z,y)``.  This file independently verifies the Faber bank,
the finite target-translation law, ``N1,N2``, the pole-divisor elimination,
the fixed-constant ``E`` obstruction, the sextic discriminant, the exact
perfect-square projected subbranch, and the retained full-mate
whole-polynomial degree-22 Faber sidecar.

The constant-field argument that this obstruction empties the genuine
``A=0`` trajectory locus is a mathematical part of THM-2411.  The
``H0*S3^2`` and ``H2*S2^2`` calculations are projected coefficient-cone
controls, not surviving Keller branches.
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
    """Laurent coefficients of ``(w^4+p*w^2+q*w+r)^(degree/4)``."""
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


def faber_seed(
    degree: int,
    w: sp.Symbol,
    d: sp.Symbol,
    q: sp.Symbol,
    s: sp.Expr,
) -> sp.Expr:
    """Polynomial part ``E_degree`` in depressed quartic coordinates."""
    coefficients = faber_coefficients(degree, 2 * d, q, d**2 - s, extra=0)
    return sp.expand(
        sum(
            coefficients[index] * w ** (degree - index)
            for index in range(degree + 1)
        )
    )


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
    """Independent finite multinomial coefficient of ``P^(degree/4)``."""
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
    """Replace every ``q^(2j)`` by ``T^j``, rejecting odd powers."""
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
    """Replace every even ``T`` power by the corresponding ``Z`` power."""
    polynomial = sp.Poly(sp.expand(expression), t)
    result = sp.Integer(0)
    for (power,), coefficient in polynomial.terms():
        require(power % 2 == 0, "unexpected odd T power")
        result += coefficient * z ** (power // 2)
    return sp.factor(result)


def faber_sidecar(
    degree: int,
    lvar: sp.Symbol,
    d: sp.Symbol,
    q: sp.Symbol,
    t: sp.Symbol,
    s: sp.Symbol,
) -> sp.Expr:
    """Compute ``E_(4j-2)-R_j(P,H)`` in the linear ``L`` coordinate."""
    coefficients = faber_coefficients(degree, 2 * d, q, d**2 - s, extra=0)
    w = sp.symbols("w_sidecar")
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
    require(
        sp.factor(p_full - (h**2 + lvar)) == 0,
        "sidecar coordinate definition failed",
    )
    return sp.factor(in_l - truncation)


def canonical_hash(expression: sp.Expr, *generators: sp.Symbol) -> str:
    """Stable digest of an expanded SymPy polynomial."""
    polynomial = sp.Poly(sp.expand(expression), *generators)
    payload = sp.srepr(polynomial.as_expr()).encode()
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    w, d, q, s, t = sp.symbols("w d q s T")
    alpha, beta, gamma, delta, epsilon, psi0 = sp.symbols(
        "alpha beta gamma delta epsilon Psi0"
    )
    b, cpar, dpar, epar, wpar = sp.symbols("B C D E W")
    y, u, z = sp.symbols("y u Z")
    p = 2 * d
    r = d**2 - s
    degrees = (2, 6, 10, 14, 18, 22)

    # Independent Laurent recurrence and finite multinomial reconstruction.
    bank = {
        degree: faber_observables(degree, p, q, r)
        for degree in degrees
    }
    for degree in degrees:
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

    # Exact finite target-translation normal form.
    shift = sp.Rational(2, 11) * alpha
    b_bar = sp.factor(beta - sp.Rational(9, 22) * alpha**2)
    c_bar = sp.factor(
        gamma
        - sp.Rational(7, 11) * alpha * beta
        + sp.Rational(21, 121) * alpha**3
    )
    d_bar = sp.factor(
        delta
        - sp.Rational(5, 11) * alpha * gamma
        + sp.Rational(35, 242) * alpha**2 * beta
        - sp.Rational(315, 10648) * alpha**4
    )
    e_bar = sp.factor(
        epsilon
        - sp.Rational(3, 11) * alpha * delta
        + sp.Rational(1815, 29282) * alpha**2 * gamma
        - sp.Rational(385, 29282) * alpha**3 * beta
        + sp.Rational(63, 29282) * alpha**5
    )
    w_bar = sp.factor(
        psi0
        + (
            105 * alpha**6
            - 770 * alpha**4 * beta
            + 4840 * alpha**3 * gamma
            - 31944 * alpha**2 * delta
            + 234256 * alpha * epsilon
        )
        / 644204
    )

    seeds = {
        degree: faber_seed(degree, w, d, q, s)
        for degree in degrees
    }
    shifted_seeds = {
        degree: faber_seed(degree, w, d, q, s - shift)
        for degree in degrees
    }
    for j, degree in enumerate(degrees, start=1):
        translated = sum(
            sp.binomial(sp.Rational(2 * j - 1, 2), index)
            * shift**index
            * seeds[4 * (j - index) - 2]
            for index in range(j)
        )
        require(
            sp.factor(shifted_seeds[degree] - translated) == 0,
            f"Faber translation failed at degree {degree}",
        )

    original_mate = (
        seeds[22]
        + alpha * seeds[18]
        + beta * seeds[14]
        + gamma * seeds[10]
        + delta * seeds[6]
        + epsilon * seeds[2]
    )
    normalized_mate = (
        shifted_seeds[22]
        + b_bar * shifted_seeds[14]
        + c_bar * shifted_seeds[10]
        + d_bar * shifted_seeds[6]
        + e_bar * shifted_seeds[2]
    )
    require(
        sp.factor(original_mate - normalized_mate) == 0,
        "whole degree-22 mate translation mismatch",
    )

    shifted_bank = {
        degree: faber_observables(
            degree, p, q, d**2 - (s - shift)
        )
        for degree in degrees
    }
    original_weights = {
        22: 1,
        18: alpha,
        14: beta,
        10: gamma,
        6: delta,
        2: epsilon,
    }
    normalized_weights = {
        22: 1,
        18: 0,
        14: b_bar,
        10: c_bar,
        6: d_bar,
        2: e_bar,
    }
    original_observables = tuple(
        sp.factor(
            sum(original_weights[n] * bank[n][index] for n in degrees)
        )
        for index in range(3)
    )
    normalized_observables = tuple(
        sp.factor(
            sum(
                normalized_weights[n] * shifted_bank[n][index]
                for n in degrees
            )
        )
        for index in range(3)
    )
    require(
        sp.factor(original_observables[0] - normalized_observables[0]) == 0,
        "first flux translation mismatch",
    )
    require(
        sp.factor(original_observables[1] - normalized_observables[1])
        == sp.factor(psi0 - w_bar),
        "second flux translation mismatch",
    )
    require(
        sp.factor(original_observables[2] - normalized_observables[2]) == 0,
        "third flux translation mismatch",
    )

    # Work from now on in the normalized alpha=0 chart, with y=11*s.
    weights = {22: 1, 14: b, 10: cpar, 6: dpar, 2: epar}
    phi = replace_even_q(
        sp.factor(sum(weights[n] * bank[n][0] for n in weights) / q),
        q,
        t,
    )
    psi = replace_even_q(
        sp.factor(
            sum(weights[n] * bank[n][1] for n in weights) - wpar
        ),
        q,
        t,
    )
    phi_in_uz = replace_t_square(
        sp.factor(phi.subs({d: u / t, s: y / 11})),
        t,
        z,
    )
    psi_in_uz = replace_t_square(
        sp.factor(psi.subs({d: u / t, s: y / 11})),
        t,
        z,
    )

    pole_a = 616 * b - 1089 * u + 63 * y**2
    pole_k = (
        -745360 * b * u * y
        + 6160 * b * y**3
        + 2342560 * cpar * u
        - 58080 * cpar * y**2
        + 511104 * dpar * y
        - 3748096 * epar
        + 922383 * u**2 * y
        - 25410 * u * y**3
        + 63 * y**5
    )
    n1 = 1331 * pole_a * z + 4 * pole_k
    n2 = (
        15944049 * z**2
        + 65591680 * b * z * y
        - 206145280 * cpar * z
        - 162339408 * z * u * y
        + 2236080 * z * y**3
        + 1443016960 * b * u**2
        - 71554560 * b * u * y**2
        + 98560 * b * y**4
        + 449771520 * cpar * u * y
        - 1239040 * cpar * y**3
        - 1978994688 * dpar * u
        + 16355328 * dpar * y**2
        - 239878144 * epar * y
        - 1319329792 * wpar
        - 1190488992 * u**3
        + 147581280 * u**2 * y**2
        - 1219680 * u * y**4
        + 672 * y**6
    )
    weighted_generators = (y, u, z, b, cpar, dpar, epar, wpar)
    weighted_values = (1, 2, 3, 2, 3, 4, 5, 6)

    def weighted_degrees(expression: sp.Expr) -> set[int]:
        polynomial = sp.Poly(sp.expand(expression), *weighted_generators)
        return {
            sum(power * weight for power, weight in zip(monomial, weighted_values))
            for monomial, _ in polynomial.terms()
        }

    require(
        weighted_degrees(pole_a) == {2}
        and weighted_degrees(pole_k) == {5}
        and weighted_degrees(n1) == {5}
        and weighted_degrees(n2) == {6},
        "normalized flux weight mismatch",
    )
    require(
        sp.Poly(n1, *weighted_generators).content() == 1
        and sp.Poly(n2, *weighted_generators).content() == 1,
        "normalized integer flux numerator is not primitive",
    )
    require(
        sp.factor(phi_in_uz + n1 / 7496192) == 0,
        "normalized first flux N1 mismatch",
    )
    require(
        sp.factor(psi_in_uz - n2 / 1319329792) == 0,
        "normalized second flux N2 mismatch",
    )

    # The exact A=0 pole divisor and the residual first-flux equation K=0.
    pole_u = sp.factor((616 * b + 63 * y**2) / 1089)
    pole_e = sp.factor(
        (
            -71148 * b**2 * y
            + 745360 * b * cpar
            + 5082 * b * y**3
            + 43560 * cpar * y**2
            + 287496 * dpar * y
            + 945 * y**5
        )
        / 2108304
    )
    require(sp.factor(pole_a.subs(u, pole_u)) == 0, "A=0 u mismatch")
    require(
        sp.factor(pole_k.subs({u: pole_u, epar: pole_e})) == 0,
        "A=0 residual K mismatch",
    )
    pole_p5 = (
        945 * y**5
        + 5082 * b * y**3
        + 43560 * cpar * y**2
        + (-71148 * b**2 + 287496 * dpar) * y
        + 745360 * b * cpar
    )
    require(
        sp.factor(
            pole_k.subs(u, pole_u)
            - sp.Rational(16, 9) * (pole_p5 - 2108304 * epar)
        )
        == 0,
        "A=0 fixed-E quintic associate mismatch",
    )
    pole_p5_poly = sp.Poly(pole_p5 - 2108304 * epar, y)
    require(
        pole_p5_poly.degree() == 5 and pole_p5_poly.LC() == 945,
        "A=0 fixed-E obstruction lost its degree-five leading term",
    )

    f2 = (
        3874403907 * z**2
        - 6375511296 * b * z * y
        - 50093303040 * cpar * z
        - 1738775808 * z * y**3
        + 59838693376 * b**3
        + 10491199488 * b**2 * y**2
        + 41215426560 * b * cpar * y
        - 272021815296 * b * dpar
        + 587575296 * b * y**4
        + 4817387520 * cpar * y**3
        - 31794757632 * dpar * y**2
        - 320597139456 * wpar
        + 20901888 * y**6
    )
    require(
        sp.factor(
            n2.subs({u: pole_u, epar: pole_e}) - f2 / 243
        )
        == 0,
        "A=0 second-flux quadratic associate mismatch",
    )
    require(weighted_degrees(f2) == {6}, "F2 weight mismatch")
    f2_associate = sp.Rational(1, 243)

    sextic = (
        42525 * y**6
        + 205821 * b * y**4
        + 1568160 * cpar * y**3
        + (-1920996 * b**2 + 7762392 * dpar) * y**2
        - 14609056 * b**3
        + 66411576 * b * dpar
        + 39530700 * cpar**2
        + 78270786 * wpar
    )
    discriminant = sp.factor(sp.discriminant(f2, z))
    discriminant_scalar = sp.Integer(63478233612288)
    require(
        sp.factor(discriminant - discriminant_scalar * sextic) == 0,
        "quadratic discriminant sextic mismatch",
    )
    sextic_poly = sp.Poly(sextic, y)
    require(
        sextic_poly.degree() == 6 and sextic_poly.LC() == 42525,
        "sextic degree or leading coefficient mismatch",
    )
    require(weighted_degrees(sextic) == {6}, "R6 weight mismatch")

    # Exact perfect-square branch.  Coefficient comparison is triangular:
    # y^5 kills the quadratic term, y^4/y^3 fix the remaining cubic,
    # y forces B*C=0, y^2 fixes D, and the constant fixes W.
    square_linear = sp.Rational(121, 50) * b
    square_constant = sp.Rational(1936, 105) * cpar
    square_model = 42525 * (
        y**3 + square_linear * y + square_constant
    ) ** 2
    square_d = sp.Rational(22141, 79200) * b**2
    square_w = -(
        2080981 * b**3 + 13186800 * cpar**2
    ) / 41164200
    y_coefficient_obstruction = sp.factor(
        sp.Poly(square_model - sextic, y).coeff_monomial(y)
    )
    require(
        y_coefficient_obstruction
        == sp.Rational(18974736, 5) * b * cpar,
        "perfect-square BC obstruction mismatch",
    )
    require(
        sp.factor(
            (sextic - square_model).subs(
                {dpar: square_d, wpar: square_w}
            )
        )
        == -sp.Rational(18974736, 5) * b * cpar * y,
        "perfect-square residual mismatch",
    )
    for square_control in (
        {b: 1, cpar: 0},
        {b: 0, cpar: 1},
        {b: 0, cpar: 0},
    ):
        require(
            sp.factor(
                (sextic - square_model)
                .subs({dpar: square_d, wpar: square_w})
                .subs(square_control)
            )
            == 0,
            "perfect-square positive control failed",
        )

    # The projected H2*S2^2 coefficient-cone control is nonvacuous.
    h2_d = sp.factor(-sp.Rational(3 * 42525, 7762392))
    h2_w = sp.factor(-sp.Rational(2 * 42525, 78270786))
    h2_squarefree = sp.Poly(y**2 - 2, y)
    h2_square = sp.Poly(y**2 + 1, y)
    h2_model = 42525 * h2_squarefree.as_expr() * h2_square.as_expr() ** 2
    require(
        sp.factor(
            sextic.subs(
                {b: 0, cpar: 0, dpar: h2_d, wpar: h2_w}
            )
            - h2_model
        )
        == 0,
        "H2*S2^2 positive control failed",
    )
    require(
        sp.gcd(h2_squarefree, h2_squarefree.diff()).degree() == 0
        and sp.gcd(h2_squarefree, h2_square).degree() == 0,
        "H2*S2^2 control is not squarefree and coprime",
    )
    generic_sextic = sp.Poly(
        sextic.subs({b: 1, cpar: 1, dpar: 1, wpar: 1}),
        y,
    )
    require(
        sp.gcd(generic_sextic, generic_sextic.diff()).degree() == 0,
        "generic squarefree hostile control failed",
    )

    # Whole-polynomial degree-22 full-mate sidecar, independent of N1/N2.
    # Individual Faber seeds need only lie in K[z]; polynomiality applies to
    # Q minus the full polynomial truncation, so lower-seed sidecars must be
    # retained.
    lvar = sp.symbols("L")
    sidecar = faber_sidecar(22, lvar, d, q, t, s)
    expected_sidecar = (
        sp.Rational(33, 2048)
        * t
        * (
            14 * lvar**4
            - 28 * lvar**3 * s
            - 14 * lvar**2 * t * d
            + 42 * lvar**2 * s**2
            - lvar * t**2
            + 56 * lvar * t * d * s
            - 56 * lvar * s**3
            + 14 * t**2 * d**2
            + 6 * t**2 * s
            - 140 * t * d * s**2
            + 70 * s**4
        )
    )
    require(
        sp.factor(sidecar - expected_sidecar) == 0,
        "degree-22 whole-polynomial Faber sidecar mismatch",
    )
    lower_sidecar_14 = faber_sidecar(14, lvar, d, q, t, s)
    lower_sidecar_10 = faber_sidecar(10, lvar, d, q, t, s)
    lower_sidecar_6 = faber_sidecar(6, lvar, d, q, t, s)
    lower_sidecar_2 = faber_sidecar(2, lvar, d, q, t, s)
    expected_sidecar_14 = (
        sp.Rational(35, 128)
        * t
        * (lvar**2 - 2 * lvar * s - t * d + 3 * s**2)
    )
    expected_sidecar_10 = -sp.Rational(5, 16) * t * (-lvar + 2 * s)
    expected_sidecar_6 = sp.Rational(3, 8) * t
    require(
        sp.factor(lower_sidecar_14 - expected_sidecar_14) == 0
        and sp.factor(lower_sidecar_10 - expected_sidecar_10) == 0
        and sp.factor(lower_sidecar_6 - expected_sidecar_6) == 0
        and sp.factor(lower_sidecar_2) == 0,
        "lower whole-polynomial Faber sidecar mismatch",
    )
    full_mate_sidecar = sp.factor(
        sidecar
        + b * lower_sidecar_14
        + cpar * lower_sidecar_10
        + dpar * lower_sidecar_6
        + epar * lower_sidecar_2
    )
    expected_full_mate_sidecar = sp.factor(
        expected_sidecar
        + b * expected_sidecar_14
        + cpar * expected_sidecar_10
        + dpar * expected_sidecar_6
    )
    require(
        sp.factor(full_mate_sidecar - expected_full_mate_sidecar) == 0,
        "degree-22 full-mate polynomial sidecar mismatch",
    )

    print("faber_flux_bank=PASS degrees=2,6,10,14,18,22")
    print("independent_multinomial_bank=PASS degrees=2,6,10,14,18,22")
    print("target_translation=PASS c=2*alpha/11")
    print(f"normalized_B={b_bar}")
    print(f"normalized_C={c_bar}")
    print(f"normalized_D={d_bar}")
    print(f"normalized_E={e_bar}")
    print(f"normalized_W={w_bar}")
    print("normalized_fluxes=N1/(-7496192),N2/1319329792")
    print(f"A0_u={pole_u}")
    print(f"A0_E={pole_e}")
    print("constant_E_obstruction=degree5,leading_coefficient=945")
    print("A0_trajectory_locus=EMPTY_BY_CONSTANT_FIELD_AND_NONSPLIT_DECK")
    print(f"F2_associate={f2_associate}")
    print(f"F2_discriminant_scalar={discriminant_scalar}")
    print(f"sextic_sha256={canonical_hash(sextic, y, b, cpar, dpar, wpar)}")
    print("projected_perfect_square_iff=BC=0,D=22141*B^2/79200,"
          "W=-(2080981*B^3+13186800*C^2)/41164200")
    print(
        "projected_H2_control="
        f"B=C=0,D={h2_d},W={h2_w},"
        "R6=42525*(y^2-2)*(y^2+1)^2"
    )
    print(
        "projected_generic_control=(B,C,D,W)=(1,1,1,1),"
        "squarefree_degree=6"
    )
    print("degree22_full_mate_whole_polynomial_sidecar=PASS")
    print("constant_field_emptiness=MATHEMATICAL_PROOF_REQUIRED")
    print("status=THM2411_DEGREE22_A0_EMPTY_EXACT_REFEREE")


if __name__ == "__main__":
    main()
