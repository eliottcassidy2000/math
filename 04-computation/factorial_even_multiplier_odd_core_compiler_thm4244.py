#!/usr/bin/env python3
"""Primary exact audit for the general complementary-block compiler.

This checker uses native integer bit operations for the barcode and a
coefficient-formula / hand-written finite-field Euclidean algorithm for the
factorial hostile.  Every truth gate raises RuntimeError explicitly, so the
normal and optimized interpreters execute the same checks.
"""

from hashlib import sha256
from math import comb


MAX_A_EXCLUSIVE = 162
PRIME_BOUND = 320
EXPONENT_BOUND = 12
SUFFIX_MAX_A = 256
SUFFIX_HEIGHT_MULTIPLIERS = 128
A14_HEIGHT_BOUND = 1 << 18
MODULAR_PRIMES = (1000003, 1000033)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def is_prime(integer):
    if integer < 2:
        return False
    if integer % 2 == 0:
        return integer == 2
    divisor = 3
    while divisor * divisor <= integer:
        if integer % divisor == 0:
            return False
        divisor += 2
    return True


def primes_below(bound):
    return tuple(integer for integer in range(2, bound) if is_prime(integer))


def valuation_two(integer):
    require(integer > 0, integer)
    value = 0
    while integer % 2 == 0:
        integer //= 2
        value += 1
    return value


def pair_representatives(multiplier):
    return tuple(range(2, multiplier // 2, 2))


def direct_candidates(multiplier, height, reset_limit):
    target = multiplier * height
    return tuple(
        test
        for test in range(1, reset_limit + 1)
        if (test * height) & target == test * height
    )


def predicted_candidates(multiplier, height, reset_limit):
    answer = []
    for test in pair_representatives(multiplier):
        complement = multiplier - test
        if (test * height) & (complement * height) == 0:
            answer.extend(
                visible
                for visible in (test, complement)
                if visible <= reset_limit
            )
    return tuple(sorted(answer))


def odd_core_candidates(multiplier, height, reset_limit):
    power = valuation_two(multiplier)
    scale = 1 << power
    odd_core = multiplier // scale
    answer = []
    for reduced in range(1, (odd_core - 1) // 2 + 1):
        reduced_complement = odd_core - reduced
        if (reduced * height) & (reduced_complement * height) == 0:
            answer.extend(
                visible
                for visible in (scale * reduced, scale * reduced_complement)
                if visible <= reset_limit
            )
    return tuple(sorted(answer))


def barcode_audit():
    digest = sha256()
    cells = 0
    truncated_cells = 0
    fixed_points = 0
    automatic_pairs = 0
    for multiplier in range(2, MAX_A_EXCLUSIVE, 2):
        multiplier_valuation = valuation_two(multiplier)
        for prime in primes_below(PRIME_BOUND):
            if prime <= multiplier or prime < 5:
                continue
            reset_limit = min(multiplier - 1, (prime - 1) // 2)
            require(reset_limit >= multiplier // 2, (multiplier, prime, reset_limit))
            if reset_limit < multiplier - 2:
                truncated_cells += EXPONENT_BOUND
            for exponent in range(1, EXPONENT_BOUND + 1):
                height = prime**exponent
                direct = direct_candidates(multiplier, height, reset_limit)
                predicted = predicted_candidates(multiplier, height, reset_limit)
                odd_core = odd_core_candidates(multiplier, height, reset_limit)
                require(direct == predicted, (multiplier, prime, exponent, direct, predicted))
                require(direct == odd_core, (
                    multiplier,
                    prime,
                    exponent,
                    direct,
                    odd_core,
                ))

                for odd_test in range(1, reset_limit + 1, 2):
                    require(odd_test not in direct, (multiplier, prime, exponent, odd_test))

                fixed = multiplier // 2
                require(fixed <= reset_limit, (multiplier, prime, reset_limit, fixed))
                require(fixed not in direct, (multiplier, prime, exponent, fixed, direct))
                require((fixed * height) & (fixed * height) == fixed * height, fixed)
                fixed_points += 1

                overlaps = []
                for test in pair_representatives(multiplier):
                    complement = multiplier - test
                    overlap = (test * height) & (complement * height)
                    overlaps.append(bool(overlap))
                    require(test <= reset_limit, (multiplier, prime, test, reset_limit))
                    upper_visible = complement <= reset_limit
                    visibility_formula = prime >= 2 * complement + 1
                    require(upper_visible == visibility_formula, (
                        multiplier,
                        prime,
                        test,
                        complement,
                        reset_limit,
                    ))
                    if valuation_two(test) < multiplier_valuation:
                        require(valuation_two(complement) == valuation_two(test), (
                            multiplier,
                            test,
                            complement,
                        ))
                        require(overlap != 0, (multiplier, height, test, complement))
                        automatic_pairs += 1

                require((direct == ()) == all(overlaps), (
                    multiplier,
                    prime,
                    exponent,
                    direct,
                    overlaps,
                ))
                record = f"{multiplier}:{prime}:{exponent}:{reset_limit}:{direct}\n"
                digest.update(record.encode())
                cells += 1
    return cells, truncated_cells, fixed_points, automatic_pairs, digest.hexdigest()


def universal_suffix_audit():
    digest = sha256()
    heights = 0
    pairs = 0
    for multiplier in range(2, SUFFIX_MAX_A + 1, 2):
        power = valuation_two(multiplier)
        scale = 1 << power
        odd_core = multiplier // scale
        exponent = (odd_core - 1).bit_length() + 1
        modulus = 1 << exponent
        require((1 << (exponent - 1)) >= odd_core, (
            multiplier,
            odd_core,
            exponent,
            modulus,
        ))
        for height_multiplier in range(1, SUFFIX_HEIGHT_MULTIPLIERS + 1):
            height = height_multiplier * modulus - 1
            require(height % 2 == 1 and height % modulus == modulus - 1, height)
            for reduced in range(1, (odd_core - 1) // 2 + 1):
                reduced_complement = odd_core - reduced
                left_residue = (reduced * height) % modulus
                right_residue = (reduced_complement * height) % modulus
                require(left_residue == modulus - reduced, (
                    multiplier,
                    height,
                    reduced,
                    left_residue,
                ))
                require(right_residue == modulus - reduced_complement, (
                    multiplier,
                    height,
                    reduced_complement,
                    right_residue,
                ))
                common_bit = 1 << (exponent - 1)
                require(left_residue & common_bit, (
                    multiplier,
                    height,
                    reduced,
                    left_residue,
                ))
                require(right_residue & common_bit, (
                    multiplier,
                    height,
                    reduced_complement,
                    right_residue,
                ))
                reduced_overlap = (reduced * height) & (reduced_complement * height)
                require(reduced_overlap & common_bit, (
                    multiplier,
                    height,
                    reduced,
                    reduced_complement,
                    reduced_overlap,
                ))
                original_overlap = (
                    (scale * reduced * height)
                    & (scale * reduced_complement * height)
                )
                require(original_overlap == reduced_overlap << power, (
                    multiplier,
                    height,
                    reduced,
                    reduced_complement,
                    original_overlap,
                    reduced_overlap,
                ))
                pairs += 1
            require(direct_candidates(multiplier, height, multiplier - 1) == (), (
                multiplier,
                odd_core,
                exponent,
                height,
            ))
            digest.update(
                f"{multiplier}:{power}:{odd_core}:{exponent}:{height}\n".encode()
            )
            heights += 1
    return heights, pairs, digest.hexdigest()


def power_two_invariance_audit():
    digest = sha256()
    cells = 0
    shift_checks = 0
    for odd_core in range(1, 34, 2):
        for height in range(1, 512, 2):
            reference = None
            for power in range(1, 7):
                scale = 1 << power
                multiplier = scale * odd_core
                direct = direct_candidates(multiplier, height, multiplier - 1)
                compiled = odd_core_candidates(multiplier, height, multiplier - 1)
                require(direct == compiled, (
                    odd_core,
                    height,
                    power,
                    direct,
                    compiled,
                ))
                normalized = tuple(candidate // scale for candidate in direct)
                if reference is None:
                    reference = normalized
                require(normalized == reference, (
                    odd_core,
                    height,
                    power,
                    normalized,
                    reference,
                ))
                for reduced in range(1, (odd_core - 1) // 2 + 1):
                    reduced_complement = odd_core - reduced
                    reduced_overlap = (
                        (reduced * height) & (reduced_complement * height)
                    )
                    original_overlap = (
                        (scale * reduced * height)
                        & (scale * reduced_complement * height)
                    )
                    require(original_overlap == reduced_overlap << power, (
                        odd_core,
                        height,
                        power,
                        reduced,
                        original_overlap,
                        reduced_overlap,
                    ))
                    shift_checks += 1
                digest.update(
                    f"{odd_core}:{height}:{power}:{normalized}\n".encode()
                )
                cells += 1
    return cells, shift_checks, digest.hexdigest()


def a14_suffix_and_prime_power_audit():
    digest = sha256()
    height_cells = 0
    for height in range(7, A14_HEIGHT_BOUND, 8):
        overlaps = tuple((test * height) & ((14 - test) * height) for test in (2, 4, 6))
        require(all(overlaps), (height, overlaps))
        require(direct_candidates(14, height, 13) == (), (height, overlaps))
        digest.update(f"H:{height}:{overlaps}\n".encode())
        height_cells += 1

    prime_power_cells = 0
    for prime in primes_below(2000):
        if prime <= 14 or prime % 8 != 7:
            continue
        reset_limit = min(13, (prime - 1) // 2)
        for exponent in range(1, 20, 2):
            height = prime**exponent
            require(height % 8 == 7, (prime, exponent, height % 8))
            require(direct_candidates(14, height, reset_limit) == (), (
                prime,
                exponent,
                reset_limit,
            ))
            digest.update(f"P:{prime}:{exponent}:{reset_limit}\n".encode())
            prime_power_cells += 1
    return height_cells, prime_power_cells, digest.hexdigest()


def hostile_records():
    requested = (
        (2, 5, 1, ()),
        (4, 5, 1, ()),
        (12, 17, 1, (4, 8)),
        (14, 17, 1, (2, 4, 6, 8)),
        (14, 19, 1, ()),
        (14, 23, 1, ()),
        (14, 23, 2, (2, 4, 6, 8, 10)),
        (14, 29, 2, (4, 10)),
        (14, 73, 1, (2, 4, 6, 8, 10, 12)),
    )
    answer = []
    for multiplier, prime, exponent, expected in requested:
        height = prime**exponent
        reset_limit = min(multiplier - 1, (prime - 1) // 2)
        direct = direct_candidates(multiplier, height, reset_limit)
        predicted = predicted_candidates(multiplier, height, reset_limit)
        require(direct == expected == predicted, (
            multiplier,
            prime,
            exponent,
            expected,
            direct,
            predicted,
        ))
        overlaps = tuple(
            (test * height) & ((multiplier - test) * height)
            for test in pair_representatives(multiplier)
        )
        answer.append((multiplier, prime, exponent, height, reset_limit, direct, overlaps))
    return tuple(answer)


def trim(poly, modulus):
    result = [coefficient % modulus for coefficient in poly]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def polynomial_remainder(dividend, divisor, modulus):
    remainder = trim(dividend, modulus)
    divisor = trim(divisor, modulus)
    require(divisor != [0], "division by zero polynomial")
    inverse_lead = pow(divisor[-1], -1, modulus)
    while len(remainder) >= len(divisor) and remainder != [0]:
        scale = remainder[-1] * inverse_lead % modulus
        shift = len(remainder) - len(divisor)
        for index, coefficient in enumerate(divisor):
            target = index + shift
            remainder[target] = (remainder[target] - scale * coefficient) % modulus
        remainder = trim(remainder, modulus)
    return remainder


def polynomial_gcd(left, right, modulus):
    left = trim(left, modulus)
    right = trim(right, modulus)
    while right != [0]:
        left, right = right, polynomial_remainder(left, right, modulus)
    inverse_lead = pow(left[-1], -1, modulus)
    return tuple(coefficient * inverse_lead % modulus for coefficient in left)


def factorial_table(bound, modulus):
    table = [1] * (bound + 1)
    for integer in range(1, bound + 1):
        table[integer] = table[integer - 1] * integer % modulus
    return table


def moment_polynomial(order, parameter_d, modulus):
    factorials = factorial_table(2 * order, modulus)
    coefficients = []
    for v_degree in range(order + 1):
        remaining = order - v_degree
        inner = 0
        for linear_count in range(remaining + 1):
            term = comb(remaining, linear_count)
            term *= pow(parameter_d, remaining - linear_count, modulus)
            term *= -1 if linear_count % 2 else 1
            term *= factorials[2 * v_degree + linear_count]
            inner = (inner + term) % modulus
        coefficients.append(comb(order, v_degree) * inner % modulus)
    return trim(coefficients, modulus)


def multinomial_control_polynomial(order, parameter_d, modulus):
    factorials = factorial_table(2 * order, modulus)
    coefficients = [0] * (order + 1)
    for v_count in range(order + 1):
        for linear_count in range(order - v_count + 1):
            constant_count = order - v_count - linear_count
            term = comb(order, v_count) * comb(order - v_count, linear_count)
            term *= pow(parameter_d, constant_count, modulus)
            term *= -1 if linear_count % 2 else 1
            term *= factorials[2 * v_count + linear_count]
            coefficients[v_count] = (coefficients[v_count] + term) % modulus
    return trim(coefficients, modulus)


def modular_nonconverse_audit():
    records = []
    multiplier = 14
    prime = 17
    exponent = 1
    height = prime**exponent
    total_order = multiplier * height
    parameter_d = total_order + 1
    reset_limit = min(multiplier - 1, (prime - 1) // 2)
    candidates = direct_candidates(multiplier, height, reset_limit)
    require(candidates == (2, 4, 6, 8), candidates)
    for modulus in MODULAR_PRIMES:
        require(is_prime(modulus), modulus)
        require(modulus > 2 * total_order, (modulus, total_order))
        positive = polynomial_gcd((1, 2, 1), (2, 3, 1), modulus)
        require(positive == (1, 1), (modulus, positive))
        for small_order in range(9):
            formula = moment_polynomial(small_order, small_order + 1, modulus)
            multinomial = multinomial_control_polynomial(
                small_order,
                small_order + 1,
                modulus,
            )
            require(formula == multinomial, (modulus, small_order, formula, multinomial))

        left = moment_polynomial(total_order - 1, parameter_d, modulus)
        right = moment_polynomial(total_order, parameter_d, modulus)
        require(len(left) - 1 == total_order - 1, (modulus, "left degree", len(left)))
        require(len(right) - 1 == total_order, (modulus, "right degree", len(right)))
        common = polynomial_gcd(left, right, modulus)
        require(common == (1,), (modulus, common))
        records.append((modulus, total_order, len(left) - 1, len(right) - 1, 0))
    return candidates, tuple(records)


def main():
    overall = sha256()
    cells, truncated, fixed, automatic, barcode_digest = barcode_audit()
    print(
        "barcode_universe="
        f"even_2<=a<{MAX_A_EXCLUSIVE},primes_a<p<{PRIME_BOUND},"
        f"1<=k<={EXPONENT_BOUND}"
    )
    print(
        f"barcode_cells={cells} truncated_cells={truncated} "
        f"fixed_point_checks={fixed} automatic_pair_checks={automatic}"
    )
    print(f"barcode_semantic_sha256={barcode_digest}")
    overall.update(barcode_digest.encode())

    suffix_heights, suffix_pairs, suffix_digest = universal_suffix_audit()
    print(
        "universal_suffix_universe="
        f"even_2<=a<={SUFFIX_MAX_A},a=2^q*b,2^(ell-1)>=b,"
        f"H=m*2^ell-1,"
        f"1<=m<={SUFFIX_HEIGHT_MULTIPLIERS}"
    )
    print(f"universal_suffix_heights={suffix_heights} pair_checks={suffix_pairs}")
    print(f"universal_suffix_semantic_sha256={suffix_digest}")
    overall.update(suffix_digest.encode())

    invariant_cells, invariant_shifts, invariant_digest = power_two_invariance_audit()
    print(
        "power_two_invariance_universe="
        "odd_1<=b<=33,odd_1<=H<512,1<=q<=6"
    )
    print(
        f"power_two_invariance_cells={invariant_cells} "
        f"overlap_shift_checks={invariant_shifts}"
    )
    print(f"power_two_invariance_semantic_sha256={invariant_digest}")
    overall.update(invariant_digest.encode())

    a14_heights, a14_prime_powers, a14_digest = a14_suffix_and_prime_power_audit()
    print(
        f"a14_suffix_heights={a14_heights} prime_odd_exponent_cells={a14_prime_powers}"
    )
    print(f"a14_semantic_sha256={a14_digest}")
    overall.update(a14_digest.encode())

    hostile = hostile_records()
    for record in hostile:
        print(f"hostile={record}")
        overall.update(f"{record}\n".encode())

    candidates, modular_records = modular_nonconverse_audit()
    print(f"nonconverse_candidate_row=a:14,p:17,k:1,C:{candidates}")
    for record in modular_records:
        print(
            "modular_nonconverse="
            f"q:{record[0]},N:{record[1]},degrees:{record[2]}/{record[3]},"
            f"gcd_degree:{record[4]}"
        )
        overall.update(f"{record}\n".encode())

    print(f"overall_semantic_sha256={overall.hexdigest()}")
    print("status=PASS runtime_require_gates=True optimization_invariant=True")


if __name__ == "__main__":
    main()
