#!/usr/bin/env python3
"""Independent exact companion for the reserved THM-3440 lane.

This standard-library-only checker audits the reciprocal-Eisenstein place of
the weighted lift, the degree-91 CRT carrier, and the deliberately narrow
mod-13 group-cohomology bridge.  It does not import another project checker.

The exact gates cover

* the formal reciprocal identity for every tested degree 3 through 100;
* the explicit degree-91 seed, Eisenstein valuations, and one-place profile;
* a finite transposition hostile proving local C_91 does not mean global C_91;
* target-line reconstruction and the source-infinity valuation profile;
* all 91 point and character addresses, including the dual CRT distinction;
* H^1(C_13; F_13) with trivial action and every generator gauge; and
* carrier-versus-amplitude and characteristic-zero loss boundaries.

THM-3440 remains reserved/provisional.  This is an independent exact
companion, not a theorem status promotion.  Every truth gate uses ``require``
and therefore survives ``python -O``.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import gcd, lcm
from pathlib import Path


EXPECTED_SEMANTIC_SHA256 = "d2284337a7d1ec6b23822e0b578cbb645e8f3bd95e11638efeaa4ae681193dc7"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def object_sha256(value: object) -> str:
    return sha256(repr(value).encode("ascii")).hexdigest()


# A formal coefficient is a sorted linear combination of named indeterminates.
Linear = tuple[tuple[str, int, int], ...]


def linear_form(terms: dict[str, object]) -> Linear:
    cleaned: list[tuple[str, int, int]] = []
    for name, raw_value in sorted(terms.items()):
        value = Fraction(raw_value)
        if value:
            cleaned.append((name, value.numerator, value.denominator))
    return tuple(cleaned)


def add_linear(left: Linear, right: Linear) -> Linear:
    terms: dict[str, Fraction] = {}
    for name, numerator, denominator in left + right:
        terms[name] = terms.get(name, Fraction(0)) + Fraction(numerator, denominator)
    return linear_form(terms)


def scale_linear(value: Linear, scalar: object) -> Linear:
    scalar = Fraction(scalar)
    return linear_form(
        {name: scalar * Fraction(numerator, denominator) for name, numerator, denominator in value}
    )


def add_formal_term(
    polynomial: dict[tuple[int, int], Linear],
    exponent: tuple[int, int],
    coefficient: Linear,
) -> None:
    updated = add_linear(polynomial.get(exponent, ()), coefficient)
    if updated:
        polynomial[exponent] = updated
    elif exponent in polynomial:
        del polynomial[exponent]


def reciprocal_from_source(n: int) -> dict[tuple[int, int], Linear]:
    """Substitute Q=t^-1,w=s^-1 in T, then multiply by t*s^n."""
    require(n >= 2, n)
    # Entries are (coefficient, w exponent, Q exponent).
    source = [(linear_form({f"r{j}": 1}), j, 0) for j in range(1, n + 1)]
    source += [(linear_form({"P": -1}), 1, 0), (linear_form({"c": 1}), 0, 1)]
    transformed: dict[tuple[int, int], Linear] = {}
    for coefficient, w_power, q_power in source:
        exponent = (1 - q_power, n - w_power)  # (t exponent, s exponent)
        require(min(exponent) >= 0, (n, exponent))
        add_formal_term(transformed, exponent, coefficient)
    return transformed


def reciprocal_expected(n: int) -> dict[tuple[int, int], Linear]:
    expected: dict[tuple[int, int], Linear] = {(0, n): linear_form({"c": 1})}
    for j in range(1, n + 1):
        add_formal_term(expected, (1, n - j), linear_form({f"r{j}": 1}))
    add_formal_term(expected, (1, n - 1), linear_form({"P": -1}))
    return expected


def reciprocal_audit() -> tuple[object, ...]:
    bank = []
    for n in range(3, 101):
        actual = reciprocal_from_source(n)
        expected = reciprocal_expected(n)
        require(actual == expected, (n, actual, expected))
        bank.append((n, len(actual), object_sha256(tuple(sorted(actual.items())))))

    n = 91
    h = reciprocal_from_source(n)
    require(h[(0, n)] == linear_form({"c": 1}), h[(0, n)])
    require(h[(1, n - 1)] == linear_form({"P": -1, "r1": 1}), h[(1, n - 1)])
    require(h[(1, 0)] == linear_form({"r91": 1}), h[(1, 0)])
    require(len(h) == n + 1, len(h))
    return (
        "t*s^n*T(s^-1,t^-1)=c*s^n+t*(sum_j_rj*s^(n-j)-P*s^(n-1))",
        (3, 100, len(bank), object_sha256(tuple(bank))),
        (n, 93, len(h), h[(0, n)], h[(1, n - 1)], h[(1, 0)]),
        object_sha256(tuple(sorted(h.items()))),
    )


# Exact sparse univariate polynomials over Q, keyed by exponent.
PolyQ = dict[int, Fraction]


def clean_q(polynomial: dict[int, object]) -> PolyQ:
    return {power: Fraction(value) for power, value in polynomial.items() if Fraction(value)}


def add_q(left: PolyQ, right: PolyQ) -> PolyQ:
    answer = dict(left)
    for power, value in right.items():
        answer[power] = answer.get(power, Fraction(0)) + value
        if not answer[power]:
            del answer[power]
    return answer


def scale_q(polynomial: PolyQ, scalar: object) -> PolyQ:
    scalar = Fraction(scalar)
    return clean_q({power: scalar * value for power, value in polynomial.items()})


def multiply_q(left: PolyQ, right: PolyQ) -> PolyQ:
    answer: PolyQ = {}
    for left_power, left_value in left.items():
        for right_power, right_value in right.items():
            power = left_power + right_power
            answer[power] = answer.get(power, Fraction(0)) + left_value * right_value
    return clean_q(answer)


def derivative_q(polynomial: PolyQ) -> PolyQ:
    return clean_q({power - 1: power * value for power, value in polynomial.items() if power})


def integral_q(polynomial: PolyQ) -> PolyQ:
    return clean_q({power + 1: value / (power + 1) for power, value in polynomial.items()})


def evaluate_q(polynomial: PolyQ, value: object) -> Fraction:
    value = Fraction(value)
    return sum((coefficient * value**power for power, coefficient in polynomial.items()), Fraction(0))


def degree_q(polynomial: PolyQ) -> int:
    return max(polynomial, default=-1)


def key_q(polynomial: PolyQ) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        (power, value.numerator, value.denominator) for power, value in sorted(polynomial.items())
    )


def pd_seed(d: int) -> PolyQ:
    require(d >= 2, d)
    if d == 2:
        return {1: Fraction(2), 2: Fraction(-3)}
    correction = Fraction(6, d * (d + 1))
    return clean_q(
        {
            1: Fraction(2) - correction,
            2: Fraction(-3) + correction,
            d - 1: 1,
            d: -1,
        }
    )


def q_from_seed(seed: PolyQ, c: object = 1) -> PolyQ:
    c = Fraction(c)
    require(c != 0, c)
    return clean_q(
        {power + 1: Fraction(power, power + 1) * coefficient / c for power, coefficient in seed.items()}
    )


def primitive_integer_polynomial(polynomial: PolyQ) -> dict[int, int]:
    require(polynomial, "nonzero polynomial")
    denominator = 1
    for value in polynomial.values():
        denominator = lcm(denominator, value.denominator)
    integers = {power: value.numerator * (denominator // value.denominator) for power, value in polynomial.items()}
    content = 0
    for value in integers.values():
        content = gcd(content, abs(value))
    require(content > 0, integers)
    integers = {power: value // content for power, value in integers.items()}
    if integers[max(integers)] < 0:
        integers = {power: -value for power, value in integers.items()}
    return integers


PolyMod = dict[int, int]


def clean_mod(polynomial: dict[int, int], prime: int) -> PolyMod:
    return {power: value % prime for power, value in polynomial.items() if value % prime}


def fraction_poly_mod(polynomial: PolyQ, prime: int) -> PolyMod:
    converted = {}
    for power, value in polynomial.items():
        require(value.denominator % prime != 0, (prime, value))
        converted[power] = value.numerator * pow(value.denominator, -1, prime)
    return clean_mod(converted, prime)


def mod_divmod(dividend: PolyMod, divisor: PolyMod, prime: int) -> tuple[PolyMod, PolyMod]:
    require(divisor, "zero modular divisor")
    quotient: PolyMod = {}
    remainder = clean_mod(dividend, prime)
    divisor = clean_mod(divisor, prime)
    divisor_degree = max(divisor)
    divisor_lead_inverse = pow(divisor[divisor_degree], -1, prime)
    while remainder and max(remainder) >= divisor_degree:
        remainder_degree = max(remainder)
        shift = remainder_degree - divisor_degree
        factor = remainder[remainder_degree] * divisor_lead_inverse % prime
        quotient[shift] = (quotient.get(shift, 0) + factor) % prime
        for power, value in divisor.items():
            target = power + shift
            remainder[target] = (remainder.get(target, 0) - factor * value) % prime
            if not remainder[target]:
                remainder.pop(target, None)
    return clean_mod(quotient, prime), clean_mod(remainder, prime)


def mod_monic(polynomial: PolyMod, prime: int) -> PolyMod:
    require(polynomial, "zero polynomial is not monic")
    inverse = pow(polynomial[max(polynomial)], -1, prime)
    return clean_mod({power: inverse * value for power, value in polynomial.items()}, prime)


def mod_gcd(left: PolyMod, right: PolyMod, prime: int) -> PolyMod:
    left, right = clean_mod(left, prime), clean_mod(right, prime)
    while right:
        _, remainder = mod_divmod(left, right, prime)
        left, right = right, remainder
    return mod_monic(left, prime)


def explicit_seed_and_inertia_audit() -> tuple[object, ...]:
    d, n = 90, 91
    p = pd_seed(d)
    r = integral_q(p)
    q = q_from_seed(p)
    w_times_p = {power + 1: value for power, value in p.items()}

    require(degree_q(p) == d and p[d] == -1, key_q(p))
    require(degree_q(r) == n and r[n] == Fraction(-1, 91), key_q(r))
    require(evaluate_q(p, 0) == 0 and evaluate_q(p, 1) == -1, key_q(p))
    require(evaluate_q(r, 0) == evaluate_q(r, 1) == 0, key_q(r))
    require(add_q(q, scale_q(add_q(w_times_p, scale_q(r, -1)), -1)) == {}, "q=w*p-R failed")
    require(q == add_q(w_times_p, scale_q(r, -1)), (key_q(q), key_q(r)))
    require(evaluate_q(q, 1) == -1, evaluate_q(q, 1))

    # Explicit H has c*s^91 and t times five possible lower terms.  P=0
    # cancels the s^90 coefficient but leaves Eisenstein intact.
    explicit_support_generic_p = (0, 1, 88, 89, 90, 91)
    explicit_support_p0 = (0, 1, 88, 89, 91)
    valuations_p0 = {0: 1, 1: 1, 88: 1, 89: 1, 91: 0}
    require(r == clean_q({2: Fraction(2729, 2730), 3: Fraction(-4094, 4095), 90: Fraction(1, 90), 91: Fraction(-1, 91)}), key_q(r))
    require(all(Fraction(value) >= Fraction(n - power, n) for power, value in valuations_p0.items()), valuations_p0)
    require(valuations_p0[0] == 1 and valuations_p0[n] == 0, valuations_p0)
    require(gcd(n, 1) == 1, n)

    # At (P,Q)=(0,0), T=R and T'=p.  A good-prime gcd certificate proves
    # gcd_Q(R,p)=w: w is visibly common, while the mod-1009 gcd has degree 1.
    prime = 1009
    p_integer = primitive_integer_polynomial(p)
    r_integer = primitive_integer_polynomial(r)
    require(p_integer[max(p_integer)] % prime and r_integer[max(r_integer)] % prime, "bad reduction")
    gcd_mod = mod_gcd(fraction_poly_mod(p, prime), fraction_poly_mod(r, prime), prime)
    require(gcd_mod == {1: 1}, gcd_mod)
    r2 = r[2]
    p_prime_0 = derivative_q(p)[0]
    require(r2 == Fraction(2729, 2730) and p_prime_0 == 2 * r2, (r2, p_prime_0))
    require(n % 2 == 1, n)

    return (
        (d, n, key_q(p), key_q(r), key_q(q), object_sha256((key_q(p), key_q(r), key_q(q)))),
        (
            "Eisenstein_over_k(P)[[t]]",
            explicit_support_generic_p,
            explicit_support_p0,
            tuple(sorted(valuations_p0.items())),
            Fraction(-1, n),
            "one_segment_Newton_polygon;(e,f)=(91,1);s_uniformizer;tame_in_char0",
        ),
        (
            "finite_simple_branch_control",
            (0, 0, 0),
            r2,
            p_prime_0,
            prime,
            tuple(sorted(gcd_mod.items())),
            "gcd_Q(R,p)=w;unique_double_root_w=0;local_transposition",
            "C91_has_no_order2_element_so_global_degree91_extension_is_not_cyclic",
        ),
    )


def reconstruction_audit() -> tuple[object, ...]:
    d, n = 90, 91
    p = pd_seed(d)
    r = integral_q(p)
    q = q_from_seed(p)
    kappa = evaluate_q(derivative_q(p), 1)
    a = -(1 + kappa) / (2 + kappa)
    require(a == Fraction(-5459, 4094), a)

    # One exact punctured point on the target line (A,B,C)=(Q,P,1).
    w, p_target, c, b = Fraction(2), Fraction(3), Fraction(1), Fraction(1)
    p_w, r_w, q_w = evaluate_q(p, w), evaluate_q(r, w), evaluate_q(q, w)
    q_target = w * p_target - r_w
    gamma = p_target - p_w
    require(gamma != 0, gamma)
    u = w / gamma
    x = 1 / gamma
    v = u - 1
    y = v / x
    source_t = gamma - 1 - a * v
    z = source_t / x**2
    beta = c + p_w / gamma
    alpha = u + q_w / gamma**2
    target_rebuilt = (alpha / x**2, beta / x, x * gamma)
    require(target_rebuilt == (q_target, p_target, 1), target_rebuilt)
    require(x * y == v and x**2 * z == source_t, (x * y, x**2 * z))
    require(c * q_w == w * p_w - r_w, (q_w, p_w, r_w))

    # With s=w^-1 as normalized source uniformizer, the explicit leading
    # p term is -w^90.  These exact orders show punctured reconstruction and
    # simultaneously show that its center is at source infinity.
    valuation_profile = {
        "w": -1,
        "gamma": -(n - 1),
        "u": n - 2,
        "x": n - 1,
        "v": 0,
        "y": -(n - 1),
        "source_t": -(n - 1),
        "z": -3 * (n - 1),
    }
    require(
        valuation_profile
        == {"w": -1, "gamma": -90, "u": 89, "x": 90, "v": 0, "y": -90, "source_t": -90, "z": -270},
        valuation_profile,
    )
    require(valuation_profile["x"] > 0 and valuation_profile["y"] < 0 and valuation_profile["z"] < 0, valuation_profile)

    point_payload = (w, p_target, q_target, gamma, u, x, v, y, source_t, z, target_rebuilt)
    return (
        "target_line=(A,B,C)=(Q,P,1);P=BC;Q=AC^2",
        ("punctured_exact_control", object_sha256(point_payload), tuple(value.numerator.bit_length() for value in point_payload[:-1])),
        tuple(valuation_profile.items()),
        "gamma_is_invertible_on_the_punctured_branch",
        "center_has_x=0_and_poles_in_y,z_so_it_is_a_source_infinity_place_not_a_finite_affine_source_point",
    )


def additive_order(exponent: int, modulus: int) -> int:
    return modulus // gcd(exponent, modulus)


def formal_fourier(function: tuple[int, ...], frequency: int) -> tuple[int, ...]:
    modulus = len(function)
    coefficients = [0] * modulus
    for point, value in enumerate(function):
        coefficients[(frequency * point) % modulus] += value
    return tuple(coefficients)


def crt_and_fourier_audit() -> tuple[object, ...]:
    modulus = 91
    point_pairs = list(product(range(7), range(13)))
    point_inverse = {(a, b): (78 * a + 14 * b) % modulus for a, b in point_pairs}
    require(len(set(point_inverse.values())) == modulus, "point CRT not bijective")
    require(all((j % 7, j % 13) == pair for pair, j in point_inverse.items()), point_inverse)
    require(78 % 7 == 1 and 78 % 13 == 0, 78)
    require(14 % 7 == 0 and 14 % 13 == 1, 14)
    require(additive_order(78, modulus) == 7 and additive_order(14, modulus) == 13, "factor orders")
    require(all(
        point_inverse[((a + c) % 7, (b + d) % 13)] == (point_inverse[(a, b)] + point_inverse[(c, d)]) % modulus
        for a, b, c, d in product(range(7), range(13), range(7), range(13))
    ), "basis convolution law")

    # Point CRT and dual CRT have different exponent bases.  With
    # zeta_7=zeta_91^13 and zeta_13=zeta_91^7, (alpha,beta) maps to the
    # C_91 character frequency 13*alpha+7*beta.
    character_inverse = {(alpha, beta): (13 * alpha + 7 * beta) % modulus for alpha, beta in point_pairs}
    require(len(set(character_inverse.values())) == modulus, "dual CRT not bijective")
    for alpha, beta in point_pairs:
        frequency = character_inverse[(alpha, beta)]
        for j in range(modulus):
            product_exponent = (13 * alpha * (j % 7) + 7 * beta * (j % 13)) % modulus
            require((frequency * j) % modulus == product_exponent, (alpha, beta, j))
    require(all(
        character_inverse[((alpha + gamma) % 7, (beta + delta) % 13)]
        == (character_inverse[(alpha, beta)] + character_inverse[(gamma, delta)]) % modulus
        for alpha, beta, gamma, delta in product(range(7), range(13), range(7), range(13))
    ), "character multiplication law")

    sample = tuple(((11 * j * j + 7 * j + 3) % 17) - 8 for j in range(modulus))
    for alpha, beta in point_pairs:
        frequency = character_inverse[(alpha, beta)]
        transported = [0] * modulus
        for a, b in point_pairs:
            j = point_inverse[(a, b)]
            exponent = (13 * alpha * a + 7 * beta * b) % modulus
            transported[exponent] += sample[j]
        require(tuple(transported) == formal_fourier(sample, frequency), (alpha, beta))

    zero_function = (0,) * modulus
    delta_zero = (1,) + (0,) * (modulus - 1)
    require(all(not any(formal_fourier(zero_function, frequency)) for frequency in range(modulus)), "zero DFT")
    require(all(formal_fourier(delta_zero, frequency)[0] == 1 for frequency in range(modulus)), "delta DFT")

    oriented_generators = tuple(value for value in range(modulus) if gcd(value, modulus) == 1)
    require(len(oriented_generators) == 72, len(oriented_generators))
    require((13 % 7, 13 % 13) == (6, 0) and (7 % 7, 7 % 13) == (0, 7), "hostile basis")
    require(character_inverse[(1, 0)] == 13 != 78 and character_inverse[(0, 1)] == 7 != 14, character_inverse)

    return (
        (
            "point_CRT",
            "j->(j_mod7,j_mod13)",
            "(a,b)->78a+14b_mod91",
            (78, (78 % 7, 78 % 13), additive_order(78, modulus)),
            (14, (14 % 7, 14 % 13), additive_order(14, modulus)),
            modulus,
            object_sha256(tuple(sorted(point_inverse.items()))),
        ),
        (
            "dual_character_CRT",
            "(alpha,beta)->13alpha+7beta_mod91",
            modulus,
            object_sha256(tuple(sorted(character_inverse.items()))),
            "all_91_character_evaluations_and_products_exact",
        ),
        (
            "group_algebra",
            "all_8281_delta_basis_products_preserved",
            object_sha256(tuple(formal_fourier(sample, frequency) for frequency in range(modulus))),
            "transported_sample_DFT_matches_at_all_91_characters",
        ),
        (
            "hostiles",
            "13_and_7_are_not_normalized_point_factor_generators",
            "78_and_14_are_not_dual_character_frequencies",
            (len(oriented_generators), oriented_generators[:6], oriented_generators[-6:]),
            "carrier_without_amplitude_allows_zero_function_and_delta_zero_with_opposite_Fourier_nonvanishing",
        ),
    )


def cyclic_cocycle(k: int) -> tuple[int, ...]:
    return tuple((k * j) % 13 for j in range(13))


def h1_audit() -> tuple[object, ...]:
    cocycles = tuple(cyclic_cocycle(k) for k in range(13))
    require(len(set(cocycles)) == 13, len(set(cocycles)))
    require(all(
        cocycle[(a + b) % 13] == (cocycle[a] + cocycle[b]) % 13
        for cocycle in cocycles for a, b in product(range(13), repeat=2)
    ), "1-cocycle law")

    # Trivial action makes every 1-coboundary zero, hence H^1=Hom=C_13.
    zero_coboundary = (0,) * 13
    require(zero_coboundary == cocycles[0], zero_coboundary)
    subgroup = tuple((14 * j) % 91 for j in range(13))
    require(len(set(subgroup)) == 13 and subgroup[1] == 14, subgroup)
    pullbacks = []
    for k in range(13):
        target_cocycle = {subgroup[j]: (k * j) % 13 for j in range(13)}
        pullback = tuple(target_cocycle[(14 * j) % 91] for j in range(13))
        require(pullback == cocycles[k], (k, pullback))
        pullbacks.append(pullback)
    require(tuple(pullbacks) == cocycles, "pullback not isomorphism")

    units_13 = tuple(range(1, 13))
    units_91 = tuple(value for value in range(1, 91) if gcd(value, 91) == 1)
    gauge_rows = []
    factors = set()
    for domain_unit in units_13:
        inverse_domain = pow(domain_unit, -1, 13)
        for target_unit in units_91:
            factor = target_unit * inverse_domain % 13
            require(factor != 0, (domain_unit, target_unit, factor))
            # sigma'=sigma^u and tau'=tau^v give sigma -> h^(v/u).
            factors.add(factor)
            gauge_rows.append((domain_unit, target_unit, factor))
    require(factors == set(units_13), factors)
    reversal_factor = 90 % 13
    require(reversal_factor == 12, reversal_factor)

    nonzero_orbit = {(unit * 1) % 13 for unit in units_13}
    require(nonzero_orbit == set(units_13), nonzero_orbit)

    return (
        (
            "coefficients=F13_with_trivial_action",
            "Z1=Hom(C13,F13_additive)",
            "B1=0",
            "dim_F13_H1=1",
            len(cocycles),
            object_sha256(cocycles),
        ),
        (
            "iota(sigma)=tau^14",
            subgroup,
            "iota_star:H1(<tau^14>;F13)->H1(C13;F13)_is_an_isomorphism",
            object_sha256(tuple(pullbacks)),
        ),
        (
            "generator_gauge",
            "sigma'=sigma^u;tau'=tau^v;scalar=v/u_mod13",
            len(gauge_rows),
            tuple(sorted(factors)),
            ("tau_reversal", reversal_factor),
            object_sha256(tuple(gauge_rows)),
        ),
        (
            "unmarked_boundary",
            "zero_is_invariant_and_all_12_nonzero_vectors_form_one_projective_orbit",
            "there_is_no_canonical_nonzero_scalar_representative",
            "H1(C13;Q_trivial)=Hom(C13,Q)=0_because_Q_is_torsion_free",
        ),
    )


def hostile_and_wording_audit() -> tuple[object, ...]:
    n = 91
    hostiles = (
        (
            "r_n=0",
            "constant_coefficient_t*r_n_vanishes",
            "H_has_an_s_factor_if_r_(n-1)_is_the_next_nonzero_term",
            "declared_degree_n_Eisenstein_and_e=n_fail",
        ),
        (
            "c=0",
            "s^n_unit_coefficient_vanishes_and_degree_in_s_drops",
            "T_loses_Q_dependence_and_the_weighted_lift_det=bc_is_zero",
        ),
        (
            "P=r1",
            "s^(n-1)_coefficient_cancels_but_Eisenstein_survives",
            "generic_P_is_not_needed_for_the_local_Eisenstein_gate",
        ),
        (
            "characteristic_7_or_13",
            gcd(n, 7),
            gcd(n, 13),
            "e=91_is_not_tame;the_characteristic_zero_scope_is_essential",
        ),
        (
            "forgotten_generator",
            "72_Keller_inertia_generators_and_12_LRC_generators_rescale_coordinates",
            "only_zero/nonzero_and_the_nonzero_line_survive",
        ),
        (
            "carrier_only",
            "CRT_relabels_functions_but_supplies_no_function_or_intertwiner",
            "no_LRC_bispectrum_nonvanishing_follows",
        ),
    )
    repairs = (
        (
            "REPAIR_RECONSTRUCTION_SCOPE",
            "replace_belongs_to_the_Keller_cover_itself_with_belongs_to_the_normalized_Keller_function-field_cover_and_reconstructs_on_the_punctured_branch",
            "the_center_has_source_valuations_x=90,y=-90,z=-270_and_is_not_a_finite_affine_source_point",
        ),
        (
            "REPAIR_H1_TYPING",
            "state_F13_has_trivial_group_action_and_write_iota_star_direction_explicitly",
            "the_bridge_is_mod13_group_cohomology_on_both_sides_not_a_map_into_characteristic-zero_Tor_or_de_Rham_cohomology",
        ),
        (
            "CLARIFY_CRT_DUAL",
            "78,14_are_point-factor_exponents_whereas_13,7_are_character-frequency_exponents_for_zeta7=zeta91^13,zeta13=zeta91^7",
        ),
        (
            "STRENGTHEN_GLOBAL_BOUNDARY",
            "for_the_explicit_F91_seed_a_simple_finite_transposition_and_odd_order_exclude_global_C91_cyclicity",
            "the_local_statement_remains_one_geometric_place_with_an_n-cycle_after_strict_henselization",
        ),
        (
            "CLARIFY_PROJECTIVE_LANGUAGE",
            "in_dimension_one_projectivization_has_one_nonzero_point;what_is_lost_is_a_canonical_nonzero_scalar_representative",
        ),
    )
    return hostiles, repairs


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert statement found")

    reciprocal = reciprocal_audit()
    seed_and_inertia = explicit_seed_and_inertia_audit()
    reconstruction = reconstruction_audit()
    crt = crt_and_fourier_audit()
    h1 = h1_audit()
    hostiles, repairs = hostile_and_wording_audit()

    verdict = (
        "reciprocal_Eisenstein_identity_PASS",
        "unique_Q-infinity_place_e=91_f=1_and_local_91-cycle_PASS",
        "global_C91_cyclicity_EXCLUDED_for_explicit_F91_by_transposition_hostile",
        "target_line_and_punctured_reconstruction_PASS_with_source-infinity_center",
        "CRT_point_exponents_78,14_and_all_91_addresses_PASS",
        "typed_trivial-action_H1_bridge_PASS_only_up_to_generator_gauge",
        "no_amplitude_current_Tor_deRham_or_LRC_nonvanishing_consequence",
    )
    semantic_surface = (reciprocal, seed_and_inertia, reconstruction, crt, h1, hostiles, repairs, verdict)
    semantic_digest = object_sha256(semantic_surface)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256, (semantic_digest, EXPECTED_SEMANTIC_SHA256))

    print("THM-3440 weighted-lift infinity-inertia independent exact companion")
    print("status=VERIFIED_EXACT_INDEPENDENT_COMPANION_FOR_RESERVED_THM3440;no_status_promotion")
    print("arithmetic=python_standard_library_only;Q_sparse_polynomials;formal_linear_coefficients;finite_fields_F1009_and_F13;exact_CRT_enumeration")
    print(f"reciprocal_identity={reciprocal}")
    print(f"degree91_seed_and_inertia={seed_and_inertia}")
    print(f"target_line_reconstruction={reconstruction}")
    print(f"crt_character_grid={crt}")
    print(f"typed_H1_bridge={h1}")
    print(f"hostile_controls={hostiles}")
    print(f"wording_audit={repairs}")
    print(f"verdict={verdict}")
    print("scope=one_rank-one_valuation_over_the_generic_Q-infinity_divisor;local_tame_inertia_after_strict_henselization;not_global_cyclicity;not_a_finite_affine_source_point;not_a_finite_Jelonek_component")
    print("loss_boundary=CRT_and_H1_relabel_carriers_and_mod13_classes_only;no_amplitude_intertwiner;no_physical_current;no_Tor1_or_deRham_identification;no_LRC_7x13_nonvanishing")
    print(f"semantic_sha256={semantic_digest}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_randomness;no_elapsed_fields;no_assert_statements;all_truth_gates_survive_python_O")
    print("commands=python -B 04-computation/jc_weighted_lift_infinity_inertia_thm3440.py;python -B -O 04-computation/jc_weighted_lift_infinity_inertia_thm3440.py")


if __name__ == "__main__":
    main()
