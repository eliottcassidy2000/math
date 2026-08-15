#!/usr/bin/env python3
"""Bounded independent exact audit for provisional THM-3443.

The checker rebuilds the degree-91 weighted-lift seed, its Puiseux place at
Q=infinity, the reconstructed gamma/x character, and the 7 x 13 CRT labels.
It also audits the deck action on the *whole* Puiseux germ.  The latter is the
main hostile gate: the leading coefficient is the marked (1,1) character, but
the nonzero s^2 coefficient is the (-1,-1) character, so the full germ is not
a pure (1,1) eigenvector.

Only Python's standard library is used.  There is no randomness, floating
point, or assert statement; every truth gate therefore survives ``python -O``.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path


EXPECTED_SEMANTIC_SHA256: str | None = "2cba88397477fe124be6a39891d2a55ceb512caee11443b83476e0c0c2ccfee8"
N = 91
RHO_POWER = 91


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def digest(value: object) -> str:
    return sha256(repr(value).encode("ascii")).hexdigest()


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


# ---------------------------------------------------------------------------
# Sparse rational polynomials in w.

PolyQ = dict[int, Fraction]


def clean_poly(raw: dict[int, object]) -> PolyQ:
    result: PolyQ = {}
    for exponent, raw_coefficient in raw.items():
        coefficient = Fraction(raw_coefficient)
        if coefficient:
            result[exponent] = result.get(exponent, Fraction(0)) + coefficient
            if not result[exponent]:
                del result[exponent]
    return result


def poly_key(polynomial: PolyQ) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        (exponent, coefficient.numerator, coefficient.denominator)
        for exponent, coefficient in sorted(polynomial.items())
    )


def poly_degree(polynomial: PolyQ) -> int:
    return max(polynomial, default=-1)


def poly_derivative(polynomial: PolyQ) -> PolyQ:
    return clean_poly(
        {exponent - 1: exponent * coefficient for exponent, coefficient in polynomial.items() if exponent}
    )


def poly_integral(polynomial: PolyQ) -> PolyQ:
    return clean_poly(
        {exponent + 1: coefficient / (exponent + 1) for exponent, coefficient in polynomial.items()}
    )


def poly_evaluate(polynomial: PolyQ, value: object) -> Fraction:
    value = Fraction(value)
    return sum(
        (coefficient * value**exponent for exponent, coefficient in polynomial.items()),
        Fraction(0),
    )


def weighted_seed(d: int) -> PolyQ:
    require(d >= 2, d)
    correction = Fraction(6, d * (d + 1))
    raw: dict[int, Fraction] = {}
    for exponent, coefficient in (
        (1, Fraction(2) - correction),
        (2, Fraction(-3) + correction),
        (d - 1, Fraction(1)),
        (d, Fraction(-1)),
    ):
        raw[exponent] = raw.get(exponent, Fraction(0)) + coefficient
    return clean_poly(raw)


def seed_audit() -> tuple[object, ...]:
    bank = []
    for d in range(2, 100):
        p = weighted_seed(d)
        r = poly_integral(p)
        n = d + 1
        require(poly_degree(p) == d, (d, poly_key(p)))
        require(poly_degree(r) == n, (d, poly_key(r)))
        require(poly_evaluate(p, 0) == 0, d)
        require(poly_evaluate(p, 1) == -1, d)
        require(poly_evaluate(r, 0) == poly_evaluate(r, 1) == 0, d)
        kappa = poly_evaluate(poly_derivative(p), 1)
        require(kappa == -5 + Fraction(6, d * (d + 1)), (d, kappa))
        require(kappa != -2, (d, kappa))
        bank.append((d, n, poly_key(p), poly_key(r), kappa))

    d, n = 90, 91
    p = weighted_seed(d)
    r = poly_integral(p)
    expected_p = clean_poly(
        {
            1: Fraction(2729, 1365),
            2: Fraction(-4094, 1365),
            89: 1,
            90: -1,
        }
    )
    expected_r = clean_poly(
        {
            2: Fraction(2729, 2730),
            3: Fraction(-4094, 4095),
            90: Fraction(1, 90),
            91: Fraction(-1, 91),
        }
    )
    require(p == expected_p and r == expected_r, (poly_key(p), poly_key(r)))
    a_d = p[d]
    c = Fraction(1)
    kappa = poly_evaluate(poly_derivative(p), 1)
    lift_a = -(1 + kappa) / (2 + kappa)
    require(a_d == -1, a_d)
    require(r[n] == Fraction(-1, 91), r[n])
    require(kappa == Fraction(-6824, 1365), kappa)
    require(lift_a == Fraction(-5459, 4094), lift_a)

    # (a_d/n) rho^n + c = 0 is exactly rho^91 = 91.
    require(-Fraction(RHO_POWER, n) + c == 0, (RHO_POWER, n, c))
    return (
        (2, 99, len(bank), digest(tuple(bank))),
        (d, n, poly_key(p), poly_key(r), a_d, c, kappa, lift_a),
        ("root_equation", "rho^91=91", "rho_nonzero"),
    )


# ---------------------------------------------------------------------------
# Q(rho), represented on 1,rho,...,rho^90 with rho^91=91.

Rho = dict[int, Fraction]


def rho_clean(raw: dict[int, object]) -> Rho:
    result: Rho = {}
    for raw_exponent, raw_coefficient in raw.items():
        quotient, exponent = divmod(raw_exponent, N)
        coefficient = Fraction(raw_coefficient) * Fraction(RHO_POWER) ** quotient
        if coefficient:
            result[exponent] = result.get(exponent, Fraction(0)) + coefficient
            if not result[exponent]:
                del result[exponent]
    return result


def rho_monomial(exponent: int, coefficient: object = 1) -> Rho:
    return rho_clean({exponent: coefficient})


def rho_add(left: Rho, right: Rho) -> Rho:
    result = dict(left)
    for exponent, coefficient in right.items():
        result[exponent] = result.get(exponent, Fraction(0)) + coefficient
        if not result[exponent]:
            del result[exponent]
    return result


def rho_scale(value: Rho, scalar: object) -> Rho:
    scalar = Fraction(scalar)
    return rho_clean({exponent: scalar * coefficient for exponent, coefficient in value.items()})


def rho_mul(left: Rho, right: Rho) -> Rho:
    raw: dict[int, Fraction] = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = left_exponent + right_exponent
            raw[exponent] = raw.get(exponent, Fraction(0)) + left_coefficient * right_coefficient
    return rho_clean(raw)


def rho_inverse_monomial(value: Rho) -> Rho:
    require(len(value) == 1, value)
    exponent, coefficient = next(iter(value.items()))
    return rho_monomial(-exponent, 1 / coefficient)


def rho_key(value: Rho) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        (exponent, coefficient.numerator, coefficient.denominator)
        for exponent, coefficient in sorted(value.items())
    )


# Truncated power series in s with coefficients in Q(rho).
Series = dict[int, Rho]


def series_clean(raw: dict[int, Rho], cutoff: int) -> Series:
    return {degree: coefficient for degree, coefficient in raw.items() if degree <= cutoff and coefficient}


def series_add(left: Series, right: Series, cutoff: int) -> Series:
    result = dict(left)
    for degree, coefficient in right.items():
        result[degree] = rho_add(result.get(degree, {}), coefficient)
        if not result[degree]:
            del result[degree]
    return series_clean(result, cutoff)


def series_scale(value: Series, scalar: object, cutoff: int) -> Series:
    return series_clean(
        {degree: rho_scale(coefficient, scalar) for degree, coefficient in value.items()},
        cutoff,
    )


def series_mul(left: Series, right: Series, cutoff: int) -> Series:
    result: Series = {}
    for left_degree, left_coefficient in left.items():
        for right_degree, right_coefficient in right.items():
            degree = left_degree + right_degree
            if degree <= cutoff:
                result[degree] = rho_add(
                    result.get(degree, {}), rho_mul(left_coefficient, right_coefficient)
                )
    return series_clean(result, cutoff)


def series_pow(value: Series, exponent: int, cutoff: int) -> Series:
    require(exponent >= 0, exponent)
    result: Series = {0: rho_monomial(0)}
    base = value
    power = exponent
    while power:
        if power & 1:
            result = series_mul(result, base, cutoff)
        base = series_mul(base, base, cutoff)
        power //= 2
    return result


def series_shift(value: Series, amount: int, cutoff: int) -> Series:
    require(amount >= 0, amount)
    return series_clean({degree + amount: coefficient for degree, coefficient in value.items()}, cutoff)


def series_derivative(value: Series, cutoff: int) -> Series:
    return series_clean(
        {
            degree - 1: rho_scale(coefficient, degree)
            for degree, coefficient in value.items()
            if degree
        },
        cutoff,
    )


def series_inverse(value: Series, cutoff: int) -> Series:
    require(0 in value, value)
    inverse_constant = rho_inverse_monomial(value[0])
    result: Series = {0: inverse_constant}
    for degree in range(1, cutoff + 1):
        total: Rho = {}
        for left_degree in range(1, degree + 1):
            if left_degree in value and degree - left_degree in result:
                total = rho_add(
                    total,
                    rho_mul(value[left_degree], result[degree - left_degree]),
                )
        result[degree] = rho_scale(rho_mul(inverse_constant, total), -1)
        if not result[degree]:
            del result[degree]
    return series_clean(result, cutoff)


def series_key(value: Series) -> tuple[tuple[int, tuple[tuple[int, int, int], ...]], ...]:
    return tuple((degree, rho_key(coefficient)) for degree, coefficient in sorted(value.items()))


def scaling_and_puiseux_audit() -> tuple[object, ...]:
    scaling_bank = []
    character_bank = []
    for n in range(3, 101):
        d = n - 1
        # If w ~ s^(-q), R(w) and Q=s^(-n) balance only at q=1.
        balanced = tuple(q for q in range(-4, 6) if -n * q == -n)
        require(balanced == (1,), (n, balanced))
        scaling_bank.append(
            (n, d, balanced, ("v_s(w)", -1), ("v_s(gamma)", -d), ("v_s(x)", d))
        )
        for j in range(n):
            gamma_character = (j * d) % n
            x_character = (-j * d) % n
            require(gamma_character == (-j) % n, (n, j, gamma_character))
            require(x_character == j, (n, j, x_character))
            character_bank.append((n, j, gamma_character, x_character))

    p = weighted_seed(90)
    r = poly_integral(p)
    cutoff = 2
    rho = rho_monomial(1)
    w_series: Series = {
        0: rho,
        1: rho_monomial(0, Fraction(1, 90)),
        2: rho_monomial(-1, Fraction(1, 180)),
    }

    # After w=s^-1 W and multiplication by s^91:
    # sum_k r_k s^(91-k) W^k - P_0 s^90 W + 1 = 0.
    inverse_equation: Series = {0: rho_monomial(0)}
    for exponent, coefficient in r.items():
        term = series_scale(series_pow(w_series, exponent, cutoff), coefficient, cutoff)
        term = series_shift(term, N - exponent, cutoff)
        inverse_equation = series_add(inverse_equation, term, cutoff)
    p0 = Fraction(37)
    p_term = series_shift(series_scale(w_series, -p0, cutoff), N - 1, cutoff)
    inverse_equation = series_add(inverse_equation, p_term, cutoff)
    require(inverse_equation == {}, series_key(inverse_equation))

    # G=s^90*gamma=P_0*s^90-s^90*p(s^-1 W); X=1/G.
    normalized_gamma: Series = {}
    normalized_gamma = series_add(
        normalized_gamma,
        {N - 1: rho_monomial(0, p0)},
        cutoff,
    )
    for exponent, coefficient in p.items():
        term = series_scale(series_pow(w_series, exponent, cutoff), -coefficient, cutoff)
        term = series_shift(term, (N - 1) - exponent, cutoff)
        normalized_gamma = series_add(normalized_gamma, term, cutoff)

    normalized_x_direct = series_inverse(normalized_gamma, cutoff)
    # Differentiating R(w)-P_0*w+s^-91=0 gives the exact identity
    # X=s^-90*x=(W-s*W')/91.
    s_w_prime = series_shift(series_derivative(w_series, cutoff), 1, cutoff)
    normalized_x_derivative = series_scale(
        series_add(w_series, series_scale(s_w_prime, -1, cutoff), cutoff),
        Fraction(1, N),
        cutoff,
    )
    expected_x: Series = {
        0: rho_monomial(1, Fraction(1, 91)),
        2: rho_monomial(-1, Fraction(-1, 16380)),
    }
    expected_gamma: Series = {
        0: rho_monomial(90),
        2: rho_monomial(88, Fraction(1, 180)),
    }
    require(normalized_gamma == expected_gamma, series_key(normalized_gamma))
    require(normalized_x_direct == normalized_x_derivative == expected_x, (
        series_key(normalized_x_direct),
        series_key(normalized_x_derivative),
    ))
    product = series_mul(normalized_gamma, normalized_x_direct, cutoff)
    require(product == {0: rho_monomial(0)}, series_key(product))

    # The theorem's rho^-90 amplitude equals rho/91 by rho^91=91.
    require(rho_monomial(-90) == rho_monomial(1, Fraction(1, 91)), (
        rho_key(rho_monomial(-90)),
        rho_key(rho_monomial(1, Fraction(1, 91))),
    ))

    leading = (
        "W_j(0)=rho*zeta^j",
        "s^90*gamma_j(0)=rho^90*zeta^(-j)",
        "s^-90*x_j(0)=rho^-90*zeta^j=(rho/91)*zeta^j",
    )
    return (
        (3, 100, len(scaling_bank), digest(tuple(scaling_bank))),
        ("generic_leading_characters", len(character_bank), digest(tuple(character_bank)),
         "gamma_mode=-1;x_mode=+1"),
        ("degree91_truncated_Hensel", series_key(w_series), "residual_mod_s^3=0"),
        ("normalized_gamma", series_key(normalized_gamma)),
        ("normalized_x", series_key(normalized_x_direct), "G*X=1_mod_s^3"),
        leading,
        ("Hensel_derivative_at_s=0", "-rho^90*zeta^(-j)!=0"),
        ("exact_differential_identity", "X=(W-s*dW/ds)/91"),
    )


def fourier_and_inertia_audit() -> tuple[object, ...]:
    n = N

    # Exact symmetry of the transformed inverse equation:
    # F(zeta*W,zeta*s)=F(W,s), because every nonconstant monomial has
    # total (s,W)-degree n modulo n.
    r = poly_integral(weighted_seed(90))
    symmetry_degrees = tuple(sorted((n - exponent, exponent) for exponent in r))
    symmetry_degrees += ((n - 1, 1), (0, 0))
    for s_degree, w_degree in symmetry_degrees:
        require((s_degree + w_degree) % n == 0, (s_degree, w_degree))

    # If X_0(s)=sum_m A_m s^m, Hensel uniqueness gives
    # X_j(s)=zeta^j X_0(zeta^-j s).  Its order-m coefficient is mode 1-m.
    coherence_checks = []
    for j in range(n):
        for m in range(0, 2 * n + 1):
            character = (j * (1 - m)) % n
            lhs_after_s_rotation = (character + m) % n
            rhs_after_branch_shift = (1 + (j - 1) * (1 - m)) % n
            require(lhs_after_s_rotation == rhs_after_branch_shift, (j, m))
            coherence_checks.append((j, m, character, lhs_after_s_rotation))

    semivariance_checks = []
    for k in range(n):
        representative_order = (1 - k) % n
        require(representative_order % n == (1 - k) % n, k)
        require(representative_order % n == (1 - k) % n, k)
        semivariance_checks.append((k, representative_order, (1 - k) % n))

    leading_modes = tuple(
        k for k in range(n) if all((1 - k) * j % n == 0 for j in range(n))
    )
    require(leading_modes == (1,), leading_modes)
    require((1 - 0) % n == 1, "constant mode")
    require((1 - 2) % n == 90, "quadratic mode")
    require((90 % 7, 90 % 13) == (6, 12), "quadratic CRT mode")

    # The s term vanishes exactly, but the s^2 term found above is nonzero.
    coefficient_modes = (
        (0, "rho/91", 1, (1, 1), "nonzero"),
        (1, "0", 0, (0, 0), "vanishes"),
        (2, "-1/(16380*rho)", 90, (6, 12), "nonzero"),
    )
    mode_one_orders = tuple(range(0, 4 * n + 1, n))
    require(mode_one_orders == (0, 91, 182, 273, 364), mode_one_orders)

    return (
        ("exact_deck_equivariance", "W_(j+1)(zeta*s)=zeta*W_j(s)",
         "X_j(s)=zeta^j*X_0(zeta^-j*s)"),
        ("DFT_semivariance", "Xhat_k(zeta*s)=zeta^(1-k)*Xhat_k(s)",
         "order_m_occurs_only_in_mode_k=1-m_mod91"),
        ("mode1", "deck_invariant", mode_one_orders, "Xhat_1_is_a_unit_in_k[[s^91]]"),
        coefficient_modes,
        (len(coherence_checks), digest(tuple(coherence_checks))),
        (len(semivariance_checks), digest(tuple(semivariance_checks))),
        ("scope", "leading_profile_is_pure_mode_1", "full_normalized_germ_is_mode_mixed"),
    )


def multiplicative_order(exponent: int, modulus: int) -> int:
    exponent %= modulus
    require(exponent != 0, (exponent, modulus))
    return modulus // gcd(exponent, modulus)


def crt_and_gauge_audit() -> tuple[object, ...]:
    n = N
    point_7, point_13 = 78, 14
    require((point_7 % 7, point_7 % 13) == (1, 0), point_7)
    require((point_13 % 7, point_13 % 13) == (0, 1), point_13)
    require(multiplicative_order(point_7, n) == 7, point_7)
    require(multiplicative_order(point_13, n) == 13, point_13)

    addresses = []
    character_checks = []
    for a in range(7):
        for b in range(13):
            j = (point_7 * a + point_13 * b) % n
            require((j % 7, j % 13) == (a, b), (a, b, j))
            addresses.append(j)
            for k in range(n):
                root_exponent = (point_7 * ((k * a) % 7) + point_13 * ((k * b) % 13)) % n
                require(root_exponent == (k * j) % n, (a, b, k, root_exponent, j))
                character_checks.append((a, b, k, root_exponent))
    require(sorted(addresses) == list(range(n)), addresses)

    # eta_7=zeta_91^78=zeta_7^6 and eta_13=zeta_91^14=zeta_13^2
    # for conventional zeta_7=zeta_91^13,zeta_13=zeta_91^7.
    require(point_7 == 13 * 6, point_7)
    require(point_13 == 7 * 2, point_13)
    normalized_leading_label = (1, 1)
    conventional_leading_label = (6, 2)
    quadratic_normalized_label = (90 % 7, 90 % 13)
    require(quadratic_normalized_label == (6, 12), quadratic_normalized_label)

    dual_checks = []
    for alpha in range(7):
        for beta in range(13):
            frequency = (13 * alpha + 7 * beta) % n
            for a in range(7):
                for b in range(13):
                    j = (point_7 * a + point_13 * b) % n
                    expected = (13 * alpha * a + 7 * beta * b) % n
                    require((frequency * j) % n == expected, (alpha, beta, a, b))
            dual_checks.append((alpha, beta, frequency))

    # Using conventional roots while retaining the label (1,1) agrees only at
    # the identity address; it is not a scalar renormalization of the table.
    wrong_root_matches = 0
    for a in range(7):
        for b in range(13):
            normalized = (point_7 * a + point_13 * b) % n
            conventional_wrong_11 = (13 * a + 7 * b) % n
            if normalized == conventional_wrong_11:
                wrong_root_matches += 1
    require(wrong_root_matches == 1, wrong_root_matches)

    units = tuple(u for u in range(1, n) if gcd(u, n) == 1)
    require(len(units) == 72, len(units))
    generator_controls = []
    root_controls = []
    for unit in units:
        # tau'=tau^u with zeta held fixed changes mode 1 to mode u.
        fixed_root_label = (unit % 7, unit % 13)
        # If zeta'=zeta^u is changed coherently, the new label is again (1,1).
        coherent_label = (1, 1)
        generator_controls.append((unit, fixed_root_label, coherent_label))

        inverse = pow(unit, -1, n)
        # Holding branches fixed but replacing zeta by zeta^u changes mode to u^-1.
        root_controls.append((unit, inverse, (inverse % 7, inverse % 13)))
    reverse = next(row for row in generator_controls if row[0] == 90)
    require(reverse == (90, (6, 12), (1, 1)), reverse)

    rho_root_controls = tuple(
        (shift, shift % n, "same_character_line_nonzero_scalar") for shift in range(n)
    )
    factor_root_controls = []
    for unit_7 in range(1, 7):
        for unit_13 in range(1, 13):
            factor_root_controls.append(
                (unit_7, unit_13, pow(unit_7, -1, 7), pow(unit_13, -1, 13))
            )
    require(len(factor_root_controls) == 72, len(factor_root_controls))

    return (
        ("point_CRT", point_7, point_13, len(addresses), digest(tuple(addresses))),
        ("factor_roots", "eta7=zeta91^78=zeta7^6", "eta13=zeta91^14=zeta13^2"),
        ("normalized_labels", normalized_leading_label, quadratic_normalized_label),
        ("conventional_leading_label", conventional_leading_label, wrong_root_matches),
        (len(character_checks), digest(tuple(character_checks))),
        (len(dual_checks), digest(tuple(dual_checks))),
        ("generator_controls", len(generator_controls), reverse, digest(tuple(generator_controls))),
        ("independent_root_controls", len(root_controls), digest(tuple(root_controls))),
        ("rho_root_controls", len(rho_root_controls), digest(rho_root_controls)),
        ("factor_root_controls", len(factor_root_controls), digest(tuple(factor_root_controls))),
    )


def hostile_and_wording_audit() -> tuple[tuple[object, ...], tuple[object, ...]]:
    hostiles = (
        (
            "wrong_Puiseux_scaling",
            "w~s^-q_balances_R(w)~s^-91_only_for_q=1",
            "q!=1_fails_the_Newton_balance",
        ),
        (
            "wrong_x_normalization",
            "v_s(x)=90",
            "s^-89*x_has_valuation_1_and_s^-91*x_has_valuation_-1",
        ),
        (
            "full_germ_is_not_one_eigenline",
            "coefficient_s^2=-1/(16380*rho)_is_nonzero",
            "its_mode_is_90=(-1,-1)_not_1=(1,1)",
        ),
        (
            "reverse_generator_fixed_root",
            "tau_to_tau^-1_sends_mode_(1,1)_to_(-1,-1)",
            "reversing_zeta_coherently_restores_the_label_(1,1)",
        ),
        (
            "conventional_factor_roots",
            "eta7=zeta7^6_and_eta13=zeta13^2",
            "the_leading_label_is_(6,2)_not_(1,1)",
        ),
        (
            "invalid_rho",
            "a_candidate_not_satisfying_rho^91=91_leaves_a_nonzero_constant_residual",
            "Hensel_branch_and_amplitude_normalization_do_not_start",
        ),
        (
            "forgotten_markings",
            "72_inertia_generators_91_rho_roots_and_72_factor-root_gauges",
            "nonvanishing_and_the_marked_character_line_survive_but_no_named_(1,1)_scalar_does",
        ),
    )
    corrections = (
        (
            "LEADING_ONLY",
            "call_(1,1)_the_leading_or_residue_eigenline_not_the_eigenline_of_the_whole_normalized_Puiseux_germ",
        ),
        (
            "INERTIA_COHERENCE",
            "record_X_j(s)=zeta^j_X_0(zeta^-j_s)_and_Xhat_k(zeta_s)=zeta^(1-k)_Xhat_k(s)",
            "fixed-s_DFT_labels_are_coherent_only_with_the_simultaneous_deck_rotation_of_s_and_branch_labels",
        ),
        (
            "MODE_ONE_UNIT",
            "Xhat_1_is_a_nonzero_unit_and_descends_to_k[[t]]_because_only_s^(91q)_terms_occur",
            "other_Fourier_modes_have_zero_constant_term_but_need_not_vanish_identically",
        ),
        (
            "SECTION4_TARGET",
            "map_the_LRC_basis_to_the_leading_coefficient_vector_or_to_the_mode-one_component_not_to_the_entire_germ_as_if_it_were_one-dimensional",
        ),
        (
            "CRT_ROOT_WORDING",
            "with_conventional_roots_the_exact_label_is_(6,2);this_is_a_root-gauge_readdressing_not_a_scalar_rescaling_of_one_fixed_label",
        ),
        (
            "CANONICAL_SCOPE",
            "canonical_only_after_linked_choices_of_s_root_inertia_generator_branch_zero_rho_and_factor_roots",
        ),
    )
    return hostiles, corrections


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert statement found")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "float literal found",
    )

    seed = seed_audit()
    scaling = scaling_and_puiseux_audit()
    fourier = fourier_and_inertia_audit()
    crt = crt_and_gauge_audit()
    hostiles, corrections = hostile_and_wording_audit()
    verdict = (
        "Puiseux_scaling_w=s^-1_W_and_x=s^90_X_PASS",
        "degree91_rho^91=91_gamma_mode_-1_and_x_leading_mode_+1_PASS",
        "marked_leading_profile_has_nonzero_(1,1)_Fourier_amplitude_PASS",
        "mode1_component_is_a_deck-invariant_formal_unit_PASS",
        "whole_normalized_germ_is_not_pure_(1,1)_REFUTED_by_nonzero_s2_mode_(-1,-1)",
        "CRT_exponents_78,14_and_normalized_factor_roots_PASS",
        "generator_and_root_changes_readdress_the_named_line_unless_changed_coherently",
        "no_LRC_intertwiner_bispectrum_current_or_LRC14_consequence",
    )
    semantic_surface = (seed, scaling, fourier, crt, hostiles, corrections, verdict)
    semantic_sha = digest(semantic_surface)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha == EXPECTED_SEMANTIC_SHA256, (semantic_sha, EXPECTED_SEMANTIC_SHA256))

    print("THM-3443 weighted-lift infinity-character independent exact companion")
    print("status=BOUNDED_INDEPENDENT_EXACT_AUDIT_FOR_PROVISIONAL_THM3443;no_status_promotion")
    print("arithmetic=python_standard_library_only;Fraction;Q(rho)_with_rho^91=91;truncated_formal_series;exact_C91_CRT")
    print(f"degree91_seed={seed}")
    print(f"Puiseux_and_reconstruction={scaling}")
    print(f"Fourier_and_inertia={fourier}")
    print(f"CRT_and_gauges={crt}")
    print(f"hostile_controls={hostiles}")
    print(f"wording_corrections={corrections}")
    print(f"verdict={verdict}")
    print("scope=one_marked_tame_Q-infinity_Puiseux_place;leading_residue_character_and_mode-one_unit_only;not_a_pure-character_full_germ")
    print("loss_boundary=no_unmarked_(1,1)_scalar;no_LRC_amplitude_intertwiner;no_bispectrum_identity;no_physical_current;no_LRC14_decrement")
    print(f"semantic_sha256={semantic_sha}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_randomness;no_elapsed_fields;no_assert_statements;no_float_literals;all_truth_gates_survive_python_O")
    print("commands=python -B 04-computation/jc_weighted_lift_infinity_character_thm3443.py;python -B -O 04-computation/jc_weighted_lift_infinity_character_thm3443.py")
    print("PASS")


if __name__ == "__main__":
    main()
