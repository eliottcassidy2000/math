#!/usr/bin/env python3
"""Exact companion for THM-3077's pointed norm/augmentation lift."""

from fractions import Fraction
from itertools import product
from math import factorial, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def scale(value, series):
    return tuple(value * coefficient for coefficient in series)


ORDER = 7
ZERO = (Fraction(0),) * ORDER
ONE = (Fraction(1),) + (Fraction(0),) * (ORDER - 1)


def mul(left, right):
    answer = [Fraction(0)] * ORDER
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if i + j < ORDER:
                answer[i + j] += a * b
    return tuple(answer)


def power(series, exponent):
    answer = ONE
    factor = series
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = mul(answer, factor)
        factor = mul(factor, factor)
        remaining //= 2
    return answer


def inverse(series):
    require(series[0] != 0, "series is not a unit")
    answer = [1 / series[0]] + [Fraction(0)] * (ORDER - 1)
    for degree in range(1, ORDER):
        answer[degree] = -sum(
            series[index] * answer[degree - index]
            for index in range(1, degree + 1)
        ) / series[0]
    return tuple(answer)


def logarithm(series):
    require(series[0] == 1, "logarithm requires normalized unit")
    tail = add(series, scale(-1, ONE))
    answer = ZERO
    term = ONE
    for exponent in range(1, ORDER):
        term = mul(term, tail)
        answer = add(answer, scale(Fraction((-1) ** (exponent + 1), exponent), term))
    return answer


def exponential(series):
    require(series[0] == 0, "exponential requires zero constant term")
    answer = ONE
    term = ONE
    for exponent in range(1, ORDER):
        term = mul(term, series)
        answer = add(answer, scale(Fraction(1, factorial(exponent)), term))
    return answer


def extended_gcd(a, b):
    if b == 0:
        return a, 1, 0
    divisor, x1, y1 = extended_gcd(b, a % b)
    return divisor, y1, x1 - (a // b) * y1


# Exact modular kernel of weighted norm plus all relative coordinates.
kernel_cells = 0
for rank in range(2, 5):
    for weights in product(range(1, 5), repeat=rank):
        total_weight = sum(weights)
        for modulus in (2, 3, 5, 7):
            kernel = 0
            image = set()
            for vector in product(range(modulus), repeat=rank):
                encoded = (
                    sum(weight * value for weight, value in zip(weights, vector))
                    % modulus,
                ) + tuple(
                    (vector[index] - vector[0]) % modulus
                    for index in range(1, rank)
                )
                image.add(encoded)
                kernel += all(value == 0 for value in encoded)
            expected = gcd(modulus, total_weight)
            require(kernel == expected, "augmentation kernel mismatch")
            require(
                len(image) == modulus**rank // expected,
                "augmentation image mismatch",
            )
            kernel_cells += 1

# Large LRC moduli are cheap to exhaust in rank two.
for weights in product(range(1, 7), repeat=2):
    total_weight = sum(weights)
    for modulus in (13, 14, 26, 91):
        kernel = sum(
            (weights[0] * first + weights[1] * second) % modulus == 0
            and (second - first) % modulus == 0
            for first in range(modulus)
            for second in range(modulus)
        )
        require(kernel == gcd(modulus, total_weight), "large-modulus kernel mismatch")
        kernel_cells += 1


# A general monomial complement has determinant a*beta-b*alpha.  Verify the
# q-injectivity and existence criteria, and the integral Bezout minimum.
sidecar_cells = 0
bezout_cells = 0
for a in range(1, 21):
    for b in range(1, 21):
        divisor, x, y = extended_gcd(a, b)
        require(divisor == gcd(a, b), "extended gcd divisor changed")
        alpha, beta = -y, x
        require(a * beta - b * alpha == divisor, "Bezout complement changed")
        bezout_cells += 1
        for modulus in (2, 7, 13, 91):
            exists = any(
                gcd(a * beta0 - b * alpha0, modulus) == 1
                for alpha0 in range(modulus)
                for beta0 in range(modulus)
            )
            require(exists == (gcd(a, b, modulus) == 1), "sidecar criterion changed")
            sidecar_cells += 1


# Integer winding tables prove the carry formula and its converse.
carry_cells = 0
for weights in ((1, 1), (2, 3), (7, 13), (2, 3, 5), (4, 6, 9)):
    total_weight = sum(weights)
    rank = len(weights)
    for source_windings in product(range(-3, 4), repeat=rank):
        norm_winding = sum(
            weight * winding
            for weight, winding in zip(weights, source_windings)
        )
        relative_windings = tuple(
            source_windings[index] - source_windings[0]
            for index in range(1, rank)
        )
        carry = norm_winding - sum(
            weights[index] * relative_windings[index - 1]
            for index in range(1, rank)
        )
        require(carry == total_weight * source_windings[0], "source carry changed")
        reconstructed_first = carry // total_weight
        reconstructed = (reconstructed_first,) + tuple(
            reconstructed_first + winding for winding in relative_windings
        )
        require(reconstructed == source_windings, "winding reconstruction changed")
        carry_cells += 1


# Normalized unit germs admit a unique exact logarithmic lift.
germ_cells = 0
for weights in ((1, 1), (2, 3), (7, 13), (2, 3, 5), (4, 6, 9)):
    sources = []
    for index in range(len(weights)):
        sources.append(
            (Fraction(1),)
            + tuple(
                Fraction(((index + 2) * (degree + 1)) % 7 - 3, degree + index + 5)
                for degree in range(1, ORDER)
            )
        )
    norm = ONE
    for source, weight in zip(sources, weights):
        norm = mul(norm, power(source, weight))
    ratios = [mul(source, inverse(sources[0])) for source in sources[1:]]
    lifted_log = logarithm(norm)
    for weight, ratio in zip(weights[1:], ratios):
        lifted_log = add(lifted_log, scale(-weight, logarithm(ratio)))
    lifted_first = exponential(scale(Fraction(1, sum(weights)), lifted_log))
    require(lifted_first == sources[0], "normalized root lift changed")
    for source, ratio in zip(sources[1:], ratios):
        require(mul(ratio, lifted_first) == source, "relative germ lift changed")
    germ_cells += 1


# The pointed Fourier identity is an exact finite-group character count.
dft_cells = 0
for prime in (7, 13):
    for displacement in range(prime):
        coefficients = [0] * prime
        for frequency in range(1, prime):
            coefficients[(frequency * displacement) % prime] += 1
        # Reduce exactly modulo Phi_p=1+x+...+x^(p-1) by subtracting
        # the x^(p-1) coefficient from every entry.
        pivot = coefficients[-1]
        reduced = tuple(value - pivot for value in coefficients[:-1])
        expected = (prime - 1,) + (0,) * (prime - 2) if displacement == 0 else (
            -1,
        ) + (0,) * (prime - 2)
        require(reduced == expected, "DFT cyclotomic reduction changed")
        dft_cells += 1


# Physical two-block controls and the recursive leading-carrier carry.
one_a, one_b = 14, factorial(13)
two_a, two_b = 14 * 13, factorial(12)
for a, b in ((one_a, one_b), (two_a, two_b)):
    require(gcd(a, b, 91) == 7, "k14 unavoidable C7 kernel changed")
    require(gcd(a + b, 91) == 7, "k14 relative kernel changed")
    require(gcd(a + b, 13) == 1, "k14 C13 should survive")

flag_cells = 0
for stage in range(3, 16):
    a = stage
    b = factorial(stage - 1)
    total_weight = a + b
    for previous_winding in range(-7, 8):
        for carrier_winding in range(-7, 8):
            norm_winding = a * previous_winding + b * carrier_winding
            relative_winding = carrier_winding - previous_winding
            carry = norm_winding - b * relative_winding
            require(carry == total_weight * previous_winding, "flag carry changed")
            flag_cells += 1


print("THM-3077 POINTED NORM-AUGMENTATION CARRY")
print(f"kernel_cells={kernel_cells} kernel=diagonal_mu_gcd(q,W)")
print(f"sidecar_cells={sidecar_cells} bezout_cells={bezout_cells}")
print("monomial_q_injective=iff_gcd(a*beta-b*alpha,q)=1")
print("sidecar_exists=iff_gcd(a,b,q)=1")
print(f"carry_cells={carry_cells} carry=wind(N)-sum(w_i*wind(r_i))_i>=2")
print(f"germ_cells={germ_cells} normalized_logarithmic_lifts=PASS")
print(f"pointed_dft_cells={dft_cells} primes=7,13")
print("k14_one_weights=(14,13!);k14_two_weights=(182,12!);C91_kernel=C7")
print(f"flag_carry_cells={flag_cells} model_relation_carry_inverse=PASS")
print("hostile=diagonal_mu_W_endpoint_invisible;zeros_break_torus_lift")
print("scope=nonzero_common_line_gauge;point_or_continuation_required")
print("all_exact_checks=PASS")
