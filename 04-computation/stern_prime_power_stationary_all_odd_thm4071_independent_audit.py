#!/usr/bin/env python3
"""Independent exact audit for THM-4071.

This path reconstructs Stern signs from Euclidean continued-fraction depth,
computes tournament stars directly, and uses extended Euclid throughout.  Its
spectral audit reduces exponent polynomials by generic long division by the
prime-power cyclotomic polynomial and realizes critical classes through the
p-adic Morse coordinate z-z^{-1}.  It shares no packet, inverse, star,
normal-form, or stationary-histogram routine with the primary script.

There are no floats or Python assertions.
"""

from collections import Counter
from hashlib import sha256
from itertools import product
from math import gcd


PACKET_PRIME_MAX = 43
PACKET_Q_MAX = 50000
LIFT_PRIME_MAX = 19
LIFT_Q_MAX = 5000
SPECTRAL_UNIVERSES = (
    (3, 1), (3, 2), (3, 3), (3, 4),
    (5, 2), (5, 3), (7, 2), (11, 2),
)
MIXED_CASES = (
    (45, (9, 5)),
    (63, (9, 7)),
    (75, (3, 25)),
    (225, (9, 25)),
)
UNIT_CACHE = {}
PRINCIPAL_CACHE = {}


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def primes_through(limit):
    result = []
    for candidate in range(3, limit + 1, 2):
        prime = True
        divisor = 3
        while divisor * divisor <= candidate:
            if candidate % divisor == 0:
                prime = False
                break
            divisor += 2
        if prime:
            result.append(candidate)
    return result


def bezout_inverse(value, modulus):
    old_remainder, remainder = value, modulus
    old_coefficient, coefficient = 1, 0
    while remainder:
        quotient = old_remainder // remainder
        old_remainder, remainder = (
            remainder,
            old_remainder - quotient * remainder,
        )
        old_coefficient, coefficient = (
            coefficient,
            old_coefficient - quotient * coefficient,
        )
    check(old_remainder == 1, "extended inverse outside unit group")
    return old_coefficient % modulus


def unit_pairs(modulus):
    if modulus not in UNIT_CACHE:
        UNIT_CACHE[modulus] = tuple(
            (unit, bezout_inverse(unit, modulus))
            for unit in range(1, modulus) if gcd(unit, modulus) == 1
        )
    return UNIT_CACHE[modulus]


def depth_sign(numerator, denominator):
    check(0 < numerator < denominator, "depth fraction outside open interval")
    check(gcd(numerator, denominator) == 1, "depth fraction not reduced")
    digit_sum = 0
    left = numerator
    right = denominator
    while left:
        digit_sum += right // left
        right, left = left, right % left
    return 1 if (digit_sum - 1) % 2 == 0 else -1


def depth_packet(modulus):
    return sum(
        depth_sign(numerator, modulus)
        for numerator in range(1, modulus)
        if gcd(numerator, modulus) == 1
    )


def direct_lower_star(modulus):
    total = 0
    for numerator in range(1, modulus):
        common = gcd(numerator, modulus)
        total += depth_sign(numerator // common, modulus // common)
    return total


def half_box_extended(modulus):
    half = (modulus - 1) // 2
    count = 0
    for left in range(1, half + 1):
        if gcd(left, modulus) != 1:
            continue
        right = bezout_inverse(4 * left, modulus)
        if right <= half:
            count += 1
    return count


def sign_from_inverse(unit, modulus):
    inverse = bezout_inverse(unit, modulus)
    return 1 if (unit + inverse) % 2 == 0 else -1


def direct_lift_layer(prime, exponent):
    modulus = prime ** exponent
    lifted_modulus = prime * modulus
    correlations = Counter()
    examples = {}
    sign_words = {}
    total = 0
    near_full_control = None
    for unit, inverse in unit_pairs(modulus):
        base_sign = 1 if (unit + inverse) % 2 == 0 else -1
        signs = []
        for digit in range(prime):
            lifted_unit = unit + digit * modulus
            signs.append(sign_from_inverse(lifted_unit, lifted_modulus))
        fiber_sum = sum(signs)
        correlation = base_sign * fiber_sum
        correlations[correlation] += 1
        examples.setdefault(correlation, (unit, inverse, tuple(signs)))
        sign_words[unit] = tuple("+" if sign == 1 else "-" for sign in signs)
        total += fiber_sum
        if unit == 1:
            check(correlation == 2 - prime,
                  "x=1 near-full lift obstruction failed")
            near_full_control = correlation
    check(total == depth_packet(lifted_modulus),
          "direct lift partition did not recover depth packet")
    check(near_full_control is not None, "missing x=1 lift control")
    return correlations, examples, sign_words


def p_adic_valuation(value, prime, exponent):
    if value == 0:
        return exponent
    valuation = 0
    while valuation < exponent and value % prime == 0:
        value //= prime
        valuation += 1
    return valuation


def direct_kernel_histogram(prime, exponent, left_frequency, right_frequency):
    modulus = prime ** exponent
    histogram = [0] * modulus
    for unit, inverse in unit_pairs(modulus):
        phase = (left_frequency * unit + right_frequency * inverse) % modulus
        histogram[phase] += 1
    return histogram


def cyclotomic_long_remainder(histogram, prime):
    """Generic monic division by Phi_(p^a), not column normalization."""
    modulus = len(histogram)
    block = modulus // prime
    degree_phi = modulus - block
    coefficients = list(histogram)
    for degree in range(modulus - 1, degree_phi - 1, -1):
        leading = coefficients[degree]
        if leading == 0:
            continue
        shift = degree - degree_phi
        for digit in range(prime):
            coefficients[shift + digit * block] -= leading
    return tuple(coefficients[:degree_phi])


def cyclotomic_equal(left, right, prime):
    return (
        cyclotomic_long_remainder(left, prime)
        == cyclotomic_long_remainder(right, prime)
    )


def add_to_histogram(histogram, phase, count=1):
    histogram[phase % len(histogram)] += count


def hensel_square_root(ratio, prime, exponent):
    root = None
    for candidate in range(1, prime):
        if candidate * candidate % prime == ratio % prime:
            root = candidate
            break
    if root is None:
        return None
    modulus = prime
    while modulus < prime ** exponent:
        next_modulus = prime * modulus
        lifted = None
        for digit in range(prime):
            candidate = root + digit * modulus
            if (candidate * candidate - ratio) % next_modulus == 0:
                lifted = candidate
                break
        check(lifted is not None, "Hensel square-root lift failed")
        root = lifted
        modulus = next_modulus
    return root


def principal_morse_data(modulus, prime):
    if modulus not in PRINCIPAL_CACHE:
        data = []
        image = Counter()
        for digit_block in range(modulus // prime):
            principal = 1 + prime * digit_block
            inverse = bezout_inverse(principal, modulus)
            morse = (principal - inverse) % modulus
            check(morse % prime == 0, "Morse coordinate left pZ")
            image[morse] += 1
            data.append((principal, morse))
        check(image == Counter(range(0, modulus, prime)),
              "z-z^-1 was not a bijection onto pZ")
        PRINCIPAL_CACHE[modulus] = tuple(data)
    return PRINCIPAL_CACHE[modulus]


def morse_stationary_histogram(
        prime, exponent, left_frequency, right_frequency):
    modulus = prime ** exponent
    left_frequency %= modulus
    right_frequency %= modulus
    result = [0] * modulus
    common = min(
        p_adic_valuation(left_frequency, prime, exponent),
        p_adic_valuation(right_frequency, prime, exponent),
    )
    if common == exponent:
        result[0] = modulus - modulus // prime
        return result
    if common:
        divisor = prime ** common
        lower = morse_stationary_histogram(
            prime,
            exponent - common,
            left_frequency // divisor,
            right_frequency // divisor,
        )
        for phase, count in enumerate(lower):
            result[(divisor * phase) % modulus] += divisor * count
        return result
    if exponent == 1:
        return direct_kernel_histogram(
            prime, exponent, left_frequency, right_frequency
        )
    if left_frequency % prime == 0 or right_frequency % prime == 0:
        return result
    ratio = right_frequency * bezout_inverse(left_frequency, modulus) % modulus
    root = hensel_square_root(ratio, prime, exponent)
    if root is None:
        return result
    coefficient = left_frequency * root % modulus
    for sign in (1, -1):
        for principal, morse in principal_morse_data(modulus, prime):
            phase = sign * (2 * coefficient + coefficient * morse * morse)
            add_to_histogram(result, phase)
            represented_unit = sign * root * principal * principal % modulus
            direct_phase = (
                left_frequency * represented_unit
                + right_frequency * bezout_inverse(represented_unit, modulus)
            ) % modulus
            check(direct_phase == phase % modulus,
                  "p-adic Morse identity failed")
    return result


def independent_frequency_category(
        prime, exponent, left_frequency, right_frequency):
    common = min(
        p_adic_valuation(left_frequency, prime, exponent),
        p_adic_valuation(right_frequency, prime, exponent),
    )
    if common == exponent:
        return "origin"
    if exponent - common == 1:
        return "prime-residual"
    divisor = prime ** common
    reduced_left = left_frequency // divisor
    reduced_right = right_frequency // divisor
    if reduced_left % prime == 0 or reduced_right % prime == 0:
        return "derivative-unit-zero"
    ratio = (
        reduced_right
        * bezout_inverse(reduced_left % prime, prime)
    ) % prime
    soluble = any(
        candidate * candidate % prime == ratio
        for candidate in range(1, prime)
    )
    return "two-critical-classes" if soluble else "no-critical-class"


def audit_spectral_universe(prime, exponent):
    modulus = prime ** exponent
    categories = Counter()
    for left_frequency in range(modulus):
        for right_frequency in range(modulus):
            direct = direct_kernel_histogram(
                prime, exponent, left_frequency, right_frequency
            )
            morse = morse_stationary_histogram(
                prime, exponent, left_frequency, right_frequency
            )
            check(cyclotomic_equal(direct, morse, prime),
                  "generic cyclotomic audit failed at %s"
                  % ((prime, exponent, left_frequency, right_frequency),))
            category = independent_frequency_category(
                prime, exponent, left_frequency, right_frequency
            )
            categories[category] += 1
    return categories


def audit_q9_degenerates():
    zero = [0] * 9
    minus_three = [0] * 9
    minus_three[0] = -3
    check(direct_kernel_histogram(3, 2, 0, 0)[0] == 6,
          "K_9 origin was not phi(9)")
    check(cyclotomic_equal(
        direct_kernel_histogram(3, 2, 1, 0), zero, 3
    ), "primitive K_9 axis did not vanish")
    check(cyclotomic_equal(
        direct_kernel_histogram(3, 2, 3, 0), minus_three, 3
    ), "residual-prime K_9 axis was not -3")
    check(cyclotomic_equal(
        direct_kernel_histogram(3, 2, 1, 2), zero, 3
    ), "nonsquare K_9 frequency did not vanish")
    check(cyclotomic_equal(
        direct_kernel_histogram(3, 2, 1, 1),
        morse_stationary_histogram(3, 2, 1, 1),
        3,
    ), "stationary K_9 frequency failed")
    return "origin=6;primitive_axis=0;residual_axis=-3;nonsquare=0;stationary=PASS"


def audit_gauss_norm(prime, max_exponent):
    checks = 0
    for exponent in range(max_exponent + 1):
        modulus = prime ** exponent
        for coefficient in range(1, prime):
            gauss = [0] * modulus
            for value in range(modulus):
                add_to_histogram(gauss, coefficient * value * value)
            norm = [0] * modulus
            for left_phase, left_count in enumerate(gauss):
                if left_count == 0:
                    continue
                for right_phase, right_count in enumerate(gauss):
                    if right_count:
                        add_to_histogram(
                            norm,
                            left_phase - right_phase,
                            left_count * right_count,
                        )
            target = [0] * modulus
            target[0] = modulus
            check(cyclotomic_equal(norm, target, prime),
                  "generic Gauss norm audit failed")
            checks += 1
    return checks


def tuple_crt_kernel_audit(modulus, blocks):
    block_data = []
    for block in blocks:
        quotient = modulus // block
        crt_unit = bezout_inverse(quotient, block)
        block_data.append((block, quotient, crt_unit, unit_pairs(block)))
    global_units = unit_pairs(modulus)
    digest = sha256()
    for left_frequency in range(modulus):
        for right_frequency in range(modulus):
            direct = Counter(
                (left_frequency * unit + right_frequency * inverse) % modulus
                for unit, inverse in global_units
            )
            local_product = Counter()
            unit_lists = [data[3] for data in block_data]
            for local_tuple in product(*unit_lists):
                phase = 0
                for data, (unit, inverse) in zip(block_data, local_tuple):
                    block, quotient, crt_unit, unused = data
                    local_phase = (
                        crt_unit * left_frequency * unit
                        + crt_unit * right_frequency * inverse
                    ) % block
                    phase += quotient * local_phase
                local_product[phase % modulus] += 1
            check(direct == local_product,
                  "tuple CRT kernel factorization failed")
            digest.update(
                ("%d:%d:%s;" % (
                    left_frequency,
                    right_frequency,
                    tuple(sorted(direct.items())),
                )).encode()
            )
    return modulus * modulus, digest.hexdigest()


def parity_plaquette(left_modulus, right_modulus):
    modulus = left_modulus * right_modulus
    left_idempotent = (
        right_modulus * bezout_inverse(right_modulus, left_modulus)
    )
    right_idempotent = (
        left_modulus * bezout_inverse(left_modulus, right_modulus)
    )
    check(left_idempotent + right_idempotent == modulus + 1,
          "independent CRT carry failed")
    flux = (
        (-1) ** left_idempotent
        * (-1) ** right_idempotent
        * (-1)
    )
    check(flux == -1, "independent parity plaquette became rank one")
    return flux


def main():
    print("THM-4071 independent exact audit")
    print("method=Euclid_depth;direct_star;generic_Phi_division;p_adic_Morse;tuple_CRT")
    print("arithmetic=integer/cyclotomic_coefficients;floats=False;asserts=False")
    print()

    packet_digest = sha256()
    packet_count = 0
    nonsquarefree_count = 0
    hostile_9 = None
    for prime in primes_through(PACKET_PRIME_MAX):
        exponent = 0
        modulus = 1
        divisor_star = 0
        while True:
            exponent += 1
            modulus *= prime
            if modulus > PACKET_Q_MAX:
                break
            packet = depth_packet(modulus)
            divisor_star += packet
            direct_star = direct_lower_star(modulus)
            phi = modulus - modulus // prime
            half_box = half_box_extended(modulus)
            check(packet == 4 * half_box - phi,
                  "Euclidean depth disagreed with half box")
            check(direct_star == divisor_star,
                  "direct tournament star disagreed with divisor star")
            record = (
                prime, exponent, modulus, packet, direct_star, half_box
            )
            packet_digest.update((repr(record) + ";").encode())
            packet_count += 1
            nonsquarefree_count += int(exponent >= 2)
            if prime == 3 and exponent == 2:
                hostile_9 = record
    check(hostile_9 == (3, 2, 9, -2, 0, 1),
          "p=3,a=2 packet/star hostile changed")
    print("depth_packet_rows=%d nonsquarefree_rows=%d" % (
        packet_count, nonsquarefree_count
    ))
    print("depth_packet_digest=%s" % packet_digest.hexdigest())
    print("p3_a2=(q,S,B,N)=(9,-2,0,1)")
    print()

    lift_layers = 0
    lift_fibers = 0
    hostile_words = None
    first_full = None
    for prime in primes_through(LIFT_PRIME_MAX):
        exponent = 1
        while prime ** (exponent + 1) <= LIFT_Q_MAX:
            correlations, examples, words = direct_lift_layer(prime, exponent)
            lift_layers += 1
            lift_fibers += sum(correlations.values())
            if prime == 3 and exponent == 1:
                hostile_words = (words[1], words[2], dict(correlations))
            if prime in correlations and first_full is None:
                first_full = (prime, exponent, examples[prime])
            exponent += 1
    check(hostile_words == (("+", "-", "-"), ("-", "-", "+"), {-1: 2}),
          "direct p=3,a=2 lift words changed")
    check(first_full is not None, "no full lift correlation found")
    print("direct_lift_layers=%d direct_lift_fibers=%d" % (
        lift_layers, lift_fibers
    ))
    print("p3_a2_lift_words=%s" % (hostile_words,))
    print("first_full_direct_correlation=%s" % (first_full,))
    print("noncontractive_x1=C_p(0,-1)=2-p verified_on_every_layer")
    print()

    spectral_pairs = 0
    print("generic cyclotomic stationary audit")
    for prime, exponent in SPECTRAL_UNIVERSES:
        modulus = prime ** exponent
        categories = audit_spectral_universe(prime, exponent)
        spectral_pairs += modulus * modulus
        print("  p=%d a=%d q=%d categories=%s" % (
            prime, exponent, modulus, dict(sorted(categories.items()))
        ))
    print("q9_degenerates=%s" % audit_q9_degenerates())
    gauss_checks = 0
    for prime, max_exponent in ((3, 4), (5, 3), (7, 3), (11, 2)):
        gauss_checks += audit_gauss_norm(prime, max_exponent)
    print("spectral_pairs=%d generic_Gauss_checks=%d" % (
        spectral_pairs, gauss_checks
    ))
    print()

    mixed_pairs = 0
    for modulus, blocks in MIXED_CASES:
        pairs, digest = tuple_crt_kernel_audit(modulus, blocks)
        mixed_pairs += pairs
        flux = parity_plaquette(blocks[0], blocks[1])
        print("mixed_tuple_CRT q=%d blocks=%s pairs=%d flux=%d digest=%s" % (
            modulus, blocks, pairs, flux, digest
        ))
    print("mixed_tuple_CRT_pairs=%d" % mixed_pairs)
    print()
    print("PASS: independent THM-4071 construction agrees")


if __name__ == "__main__":
    main()
