#!/usr/bin/env python3
"""Primary exact audit for THM-4071.

The audit uses canonical modular inverses for Stern packets, the exact affine
inverse-lift law, and a specialized prime-power cyclotomic normal form for the
p-adic stationary formula.  Mixed prime-power CRT kernels are compared as
full exponent histograms.  There are no floats or Python assertions.
"""

from collections import Counter
from fractions import Fraction
from math import gcd


PACKET_PRIME_MAX = 101
PACKET_Q_MAX = 400000
LIFT_PRIME_MAX = 29
LIFT_Q_MAX = 100000
KLOOSTERMAN_UNIVERSES = (
    (3, 1), (3, 2), (3, 3), (3, 4),
    (5, 2), (5, 3), (7, 2), (11, 2),
)
CRT_CASES = (
    (45, (9, 5)),
    (63, (9, 7)),
    (75, (3, 25)),
    (225, (9, 25)),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def primes_through(limit):
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for prime in range(2, limit + 1):
        if sieve[prime]:
            for composite in range(prime * prime, limit + 1, prime):
                sieve[composite] = False
    return [prime for prime in range(3, limit + 1, 2) if sieve[prime]]


def valuation_mod(value, prime, exponent):
    if value == 0:
        return exponent
    valuation = 0
    while valuation < exponent and value % prime == 0:
        value //= prime
        valuation += 1
    return valuation


def stern_packet(modulus, prime):
    total = 0
    for unit in range(1, modulus):
        if unit % prime:
            inverse = pow(unit, -1, modulus)
            total += 1 if (unit + inverse) % 2 == 0 else -1
    return total


def half_box_count(modulus, prime):
    half = (modulus - 1) // 2
    count = 0
    for left in range(1, half + 1):
        if left % prime and pow(4 * left, -1, modulus) <= half:
            count += 1
    return count


def affine_correlation(prime, alpha, beta):
    total = 0
    for top_digit in range(prime):
        inverse_digit = (alpha + beta * top_digit) % prime
        total += 1 if (top_digit + inverse_digit) % 2 == 0 else -1
    return total


def verify_lift_layer(prime, exponent):
    modulus = prime ** exponent
    lifted_modulus = prime * modulus
    correlations = Counter()
    examples = {}
    signed_total = 0
    fibers = 0
    for unit in range(1, modulus):
        if unit % prime == 0:
            continue
        inverse = pow(unit, -1, modulus)
        carry = (unit * inverse - 1) // modulus
        alpha = (-inverse * carry) % prime
        beta = (-inverse * inverse) % prime
        correlation = affine_correlation(prime, alpha, beta)
        base_sign = 1 if (unit + inverse) % 2 == 0 else -1
        direct_fiber = 0
        for top_digit in range(prime):
            lifted_unit = unit + top_digit * modulus
            lifted_inverse = pow(lifted_unit, -1, lifted_modulus)
            require((lifted_inverse - inverse) % modulus == 0,
                    "inverse failed to reduce")
            inverse_digit = (lifted_inverse - inverse) // modulus
            require(inverse_digit == (alpha + beta * top_digit) % prime,
                    "affine inverse-lift law failed")
            direct_fiber += (
                1 if (lifted_unit + lifted_inverse) % 2 == 0 else -1
            )
        require(direct_fiber == base_sign * correlation,
                "affine fiber correlation failed")
        signed_total += direct_fiber
        correlations[correlation] += 1
        examples.setdefault(correlation,
                            (unit, inverse, carry, alpha, beta))
        fibers += 1
    require(signed_total == stern_packet(lifted_modulus, prime),
            "lift recursion did not recover the packet")
    return fibers, correlations, examples


def add_phase(histogram, phase, multiplicity=1):
    histogram[phase % len(histogram)] += multiplicity


def kloosterman_histogram(prime, exponent, left_frequency, right_frequency):
    modulus = prime ** exponent
    histogram = [0] * modulus
    for unit in range(1, modulus):
        if unit % prime:
            phase = (
                left_frequency * unit
                + right_frequency * pow(unit, -1, modulus)
            ) % modulus
            histogram[phase] += 1
    return histogram


def cyclotomic_normal_form(histogram, prime):
    """Reduce by Phi_(p^a) using its disjoint residue columns."""
    modulus = len(histogram)
    column_count = modulus // prime
    normal = []
    for residue in range(column_count):
        pivot = histogram[residue + (prime - 1) * column_count]
        for digit in range(prime - 1):
            normal.append(
                histogram[residue + digit * column_count] - pivot
            )
    return tuple(normal)


def stationary_histogram(prime, exponent, left_frequency, right_frequency):
    modulus = prime ** exponent
    left_frequency %= modulus
    right_frequency %= modulus
    result = [0] * modulus
    common = min(
        valuation_mod(left_frequency, prime, exponent),
        valuation_mod(right_frequency, prime, exponent),
    )
    if common == exponent:
        result[0] = modulus - modulus // prime
        return result
    if common:
        divisor = prime ** common
        lower = stationary_histogram(
            prime,
            exponent - common,
            left_frequency // divisor,
            right_frequency // divisor,
        )
        for phase, count in enumerate(lower):
            result[(divisor * phase) % modulus] += divisor * count
        return result
    if exponent == 1:
        return kloosterman_histogram(
            prime, exponent, left_frequency, right_frequency
        )
    if left_frequency % prime == 0 or right_frequency % prime == 0:
        return result
    ratio = right_frequency * pow(left_frequency, -1, modulus) % modulus
    roots = [
        root for root in range(1, modulus)
        if root % prime and root * root % modulus == ratio
    ]
    if not roots:
        return result
    require(len(roots) == 2, "unit square did not have two roots")
    root = min(roots)
    coefficient = left_frequency * root % modulus
    for multiple in range(0, modulus, prime):
        add_phase(result, 2 * coefficient + coefficient * multiple * multiple)
        add_phase(result, -2 * coefficient - coefficient * multiple * multiple)
    return result


def frequency_category(prime, exponent, left_frequency, right_frequency):
    common = min(
        valuation_mod(left_frequency, prime, exponent),
        valuation_mod(right_frequency, prime, exponent),
    )
    if common == exponent:
        return "origin"
    residual_exponent = exponent - common
    if residual_exponent == 1:
        return "prime-residual"
    divisor = prime ** common
    reduced_left = left_frequency // divisor
    reduced_right = right_frequency // divisor
    if reduced_left % prime == 0 or reduced_right % prime == 0:
        return "one-unit-vanishing"
    ratio = reduced_right * pow(reduced_left, -1, prime) % prime
    square = any(root * root % prime == ratio for root in range(1, prime))
    return "stationary-square" if square else "nonresidue-vanishing"


def verify_kloosterman_universe(prime, exponent):
    modulus = prime ** exponent
    categories = Counter()
    for left_frequency in range(modulus):
        for right_frequency in range(modulus):
            direct = kloosterman_histogram(
                prime, exponent, left_frequency, right_frequency
            )
            stationary = stationary_histogram(
                prime, exponent, left_frequency, right_frequency
            )
            require(
                cyclotomic_normal_form(direct, prime)
                == cyclotomic_normal_form(stationary, prime),
                "stationary formula failed at %s"
                % ((prime, exponent, left_frequency, right_frequency),),
            )
            categories[frequency_category(
                prime, exponent, left_frequency, right_frequency
            )] += 1
    return categories


def verify_gauss_norm(prime, max_exponent):
    checks = 0
    for exponent in range(max_exponent + 1):
        modulus = prime ** exponent
        for coefficient in range(1, prime):
            gauss = [0] * modulus
            for value in range(modulus):
                add_phase(gauss, coefficient * value * value)
            norm_square = [0] * modulus
            for left_phase, left_count in enumerate(gauss):
                if left_count == 0:
                    continue
                for right_phase, right_count in enumerate(gauss):
                    if right_count:
                        add_phase(
                            norm_square,
                            left_phase - right_phase,
                            left_count * right_count,
                        )
            target = [0] * modulus
            target[0] = modulus
            require(
                cyclotomic_normal_form(norm_square, prime)
                == cyclotomic_normal_form(target, prime),
                "quadratic Gauss norm failed",
            )
            checks += 1
    return checks


def verify_parity_plaquettes():
    count = 0
    for left_modulus in range(3, 26, 2):
        for right_modulus in range(left_modulus + 2, 26, 2):
            if gcd(left_modulus, right_modulus) != 1:
                continue
            modulus = left_modulus * right_modulus
            left_idempotent = (
                right_modulus * pow(right_modulus, -1, left_modulus)
            )
            right_idempotent = (
                left_modulus * pow(left_modulus, -1, right_modulus)
            )
            require(left_idempotent + right_idempotent == modulus + 1,
                    "CRT carry identity failed")
            flux = (
                (-1) ** left_idempotent
                * (-1) ** right_idempotent
                * (-1)
            )
            require(flux == -1, "canonical parity became rank one")
            count += 1
    return count


def verify_prime_power_block_crt(modulus, blocks):
    running_product = 1
    for block in blocks:
        require(gcd(running_product, block) == 1,
                "CRT blocks are not coprime")
        running_product *= block
    require(running_product == modulus, "CRT block product mismatch")
    global_units = [
        (unit, pow(unit, -1, modulus))
        for unit in range(1, modulus) if gcd(unit, modulus) == 1
    ]
    local_data = {}
    for block in blocks:
        quotient = modulus // block
        crt_unit = pow(quotient, -1, block)
        units = [
            (unit, pow(unit, -1, block))
            for unit in range(1, block) if gcd(unit, block) == 1
        ]
        local_data[block] = (quotient, crt_unit, units)
    for left_frequency in range(modulus):
        for right_frequency in range(modulus):
            direct = Counter(
                (left_frequency * unit + right_frequency * inverse) % modulus
                for unit, inverse in global_units
            )
            convolved = Counter({0: 1})
            for block in blocks:
                quotient, crt_unit, units = local_data[block]
                local = Counter(
                    quotient * (
                        crt_unit * left_frequency * unit
                        + crt_unit * right_frequency * inverse
                    ) % modulus
                    for unit, inverse in units
                )
                next_convolution = Counter()
                for old_phase, old_count in convolved.items():
                    for new_phase, new_count in local.items():
                        next_convolution[(old_phase + new_phase) % modulus] += (
                            old_count * new_count
                        )
                convolved = next_convolution
            require(direct == convolved,
                    "prime-power block CRT factorization failed")
    return modulus * modulus


def main():
    print("THM-4071 primary exact audit")
    print("method=inverse_parity;affine_lift;stationary_cyclotomic;block_CRT")
    print("arithmetic=integer/GF(2)/cyclotomic_coefficients;floats=False;asserts=False")
    print()

    print("prime-power packets (a:q:S:B)")
    packet_rows = 0
    nonsquarefree_rows = 0
    largest_modulus = 0
    strongest_packet = (Fraction(0), None)
    for prime in primes_through(PACKET_PRIME_MAX):
        modulus = 1
        exponent = 0
        star = 0
        cells = []
        while True:
            exponent += 1
            modulus *= prime
            if modulus > PACKET_Q_MAX:
                break
            packet = stern_packet(modulus, prime)
            star += packet
            phi = modulus - modulus // prime
            half_box = half_box_count(modulus, prime)
            require(packet == 4 * half_box - phi,
                    "half-box packet identity failed")
            cells.append("%d:%d:%+d:%+d" % (
                exponent, modulus, packet, star
            ))
            ratio = Fraction(abs(packet) ** 2, modulus)
            if ratio > strongest_packet[0]:
                strongest_packet = (
                    ratio, (prime, exponent, modulus, packet)
                )
            packet_rows += 1
            nonsquarefree_rows += int(exponent >= 2)
            largest_modulus = max(largest_modulus, modulus)
        print("  p=%d %s" % (prime, ",".join(cells)))
    print("packet_rows=%d nonsquarefree_rows=%d largest_q=%d" % (
        packet_rows, nonsquarefree_rows, largest_modulus
    ))
    print("max_absS_squared_over_q=%s" % (strongest_packet,))
    print()

    print("affine lift audit")
    lift_layers = 0
    lift_fibers = 0
    hostile_3_2 = None
    first_full_correlation = None
    for prime in primes_through(LIFT_PRIME_MAX):
        exponent = 1
        while prime ** (exponent + 1) <= LIFT_Q_MAX:
            fibers, correlations, examples = verify_lift_layer(prime, exponent)
            lift_layers += 1
            lift_fibers += fibers
            if prime == 3 and exponent == 1:
                hostile_3_2 = (dict(sorted(correlations.items())), examples)
            if prime in (3, 5, 7) and exponent <= 3:
                print("  p=%d %d->%d C_hist=%s" % (
                    prime, exponent, exponent + 1,
                    dict(sorted(correlations.items())),
                ))
            if prime in correlations and first_full_correlation is None:
                first_full_correlation = (
                    prime, exponent, examples[prime]
                )
            exponent += 1
    require(hostile_3_2 is not None, "missing p=3,a=2 hostile")
    print("lift_layers=%d lift_fibers=%d" % (lift_layers, lift_fibers))
    print("p3_a2_hostile=%s" % (hostile_3_2,))
    print("first_actual_full_affine_correlation=%s" % (
        first_full_correlation,
    ))
    print()

    print("prime-power Kloosterman stationary audit")
    frequency_pairs = 0
    for prime, exponent in KLOOSTERMAN_UNIVERSES:
        modulus = prime ** exponent
        categories = verify_kloosterman_universe(prime, exponent)
        frequency_pairs += modulus * modulus
        print("  p=%d a=%d q=%d categories=%s" % (
            prime, exponent, modulus, dict(sorted(categories.items()))
        ))
    gauss_checks = 0
    for prime, max_exponent in ((3, 4), (5, 3), (7, 3), (11, 3)):
        gauss_checks += verify_gauss_norm(prime, max_exponent)
    print("kloosterman_pairs=%d gauss_norm_checks=%d" % (
        frequency_pairs, gauss_checks
    ))
    print()

    plaquettes = verify_parity_plaquettes()
    packet_3 = stern_packet(3, 3)
    packet_5 = stern_packet(5, 5)
    packet_15 = sum(
        1 if (unit + pow(unit, -1, 15)) % 2 == 0 else -1
        for unit in range(1, 15) if gcd(unit, 15) == 1
    )
    require(packet_3 * packet_5 == 0 and packet_15 == 8,
            "nonmultiplicative packet hostile failed")
    print("CRT parity plaquettes=%d all_flux=-1" % plaquettes)
    print("nonmultiplicative_hostile=S3*S5=%d S15=%d" % (
        packet_3 * packet_5, packet_15
    ))
    block_pairs = 0
    for modulus, blocks in CRT_CASES:
        pairs = verify_prime_power_block_crt(modulus, blocks)
        block_pairs += pairs
        print("  mixed_block_CRT q=%d blocks=%s pairs=%d" % (
            modulus, blocks, pairs
        ))
    print("mixed_block_CRT_pairs=%d" % block_pairs)
    print()
    print("PASS: THM-4071 exact controls verified")


if __name__ == "__main__":
    main()
