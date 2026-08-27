#!/usr/bin/env python3
"""Structurally independent audit for THM-4243.

The barcode audit represents binary expansions as explicit sets of occupied
positions.  The hostile factorial check separately builds the two actual
moment polynomials modulo large primes and runs a polynomial Euclidean
algorithm.  No code is imported from the primary companion.
"""

from hashlib import sha256
from math import comb


ODD_HEIGHT_BOUND = 1 << 15
PRIME_BOUND = 350
EXPONENT_BOUND = 19
RESET_LIMITS = (5, 6, 8, 9)
MODULAR_PRIMES = (1000003, 1000033)
GCD_ROWS = ((11, 1), (13, 1), (17, 1), (23, 1), (43, 1))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def occupied_positions(integer):
    positions = set()
    position = 0
    while integer:
        integer, digit = divmod(integer, 2)
        if digit:
            positions.add(position)
        position += 1
    return positions


def set_submask_candidates(height, limit):
    target = occupied_positions(10 * height)
    return tuple(
        multiplier
        for multiplier in range(1, limit + 1)
        if occupied_positions(multiplier * height) <= target
    )


def support_pair_prediction(height, limit):
    first_disjoint = not (
        occupied_positions(2 * height) & occupied_positions(8 * height)
    )
    second_disjoint = not (
        occupied_positions(4 * height) & occupied_positions(6 * height)
    )
    answer = []
    for multiplier in range(1, limit + 1):
        if multiplier in (2, 8) and first_disjoint:
            answer.append(multiplier)
        if multiplier in (4, 6) and second_disjoint:
            answer.append(multiplier)
    return tuple(answer)


def support_closes(height):
    first = occupied_positions(2 * height) & occupied_positions(8 * height)
    second = occupied_positions(4 * height) & occupied_positions(6 * height)
    return bool(first) and bool(second)


def independent_height_audit():
    digest = sha256()
    cells = 0
    for height in range(1, ODD_HEIGHT_BOUND, 2):
        for limit in RESET_LIMITS:
            direct = set_submask_candidates(height, limit)
            predicted = support_pair_prediction(height, limit)
            require(direct == predicted, (height, limit, direct, predicted))
            require((direct == ()) == support_closes(height), (height, limit))
            digest.update(f"{height}:{limit}:{direct}\n".encode())
            cells += 1
    return cells, digest.hexdigest()


def is_prime(integer):
    if integer < 2:
        return False
    divisor = 2
    while divisor * divisor <= integer:
        if integer % divisor == 0:
            return False
        divisor += 1
    return True


def primes_below(bound):
    return tuple(integer for integer in range(2, bound) if is_prime(integer))


def independent_prime_power_audit():
    digest = sha256()
    rows = 0
    for prime in primes_below(PRIME_BOUND):
        if prime < 11:
            continue
        limit = min(9, (prime - 1) // 2)
        for exponent in range(1, EXPONENT_BOUND + 1):
            height = prime ** exponent
            direct = set_submask_candidates(height, limit)
            predicted = support_pair_prediction(height, limit)
            require(direct == predicted, (prime, exponent, direct, predicted))
            if prime % 8 in (5, 7) and exponent % 2:
                require(direct == (), (prime, exponent, direct))
            digest.update(f"{prime}:{exponent}:{direct}\n".encode())
            rows += 1
    return rows, digest.hexdigest()


def independent_modulus_64_audit():
    forcing = []
    witnesses = []
    low_mask = (1 << 7) - 1
    for residue in range(1, 64, 2):
        first = occupied_positions((2 * residue) & low_mask) & occupied_positions(
            (8 * residue) & low_mask
        )
        second = occupied_positions((4 * residue) & low_mask) & occupied_positions(
            (6 * residue) & low_mask
        )
        if first and second:
            forcing.append(residue)
        else:
            require(not support_closes(residue), residue)
            witnesses.append((residue, set_submask_candidates(residue, 9)))
    expected = tuple(
        residue
        for residue in range(1, 64, 2)
        if residue % 8 in (5, 7)
        or residue % 32 == 27
        or residue in (41, 57)
    )
    require(tuple(forcing) == expected, (forcing, expected))
    return tuple(forcing), tuple(witnesses)


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
            remainder[target] = (
                remainder[target] - scale * coefficient
            ) % modulus
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
    """Return A_order^(parameter_d)(v), low coefficient first."""
    factorials = factorial_table(2 * order, modulus)
    coefficients = []
    for v_degree in range(order + 1):
        remaining = order - v_degree
        inner = 0
        for linear_count in range(remaining + 1):
            sign = -1 if linear_count % 2 else 1
            term = comb(remaining, linear_count)
            term *= pow(parameter_d, remaining - linear_count, modulus)
            term *= sign
            term *= factorials[2 * v_degree + linear_count]
            inner = (inner + term) % modulus
        coefficients.append(comb(order, v_degree) * inner % modulus)
    return trim(coefficients, modulus)


def multinomial_moment_polynomial(order, parameter_d, modulus):
    """Small-order control using explicit counts of constant/linear/v slots."""
    factorials = factorial_table(2 * order, modulus)
    coefficients = [0] * (order + 1)
    for v_count in range(order + 1):
        for linear_count in range(order - v_count + 1):
            constant_count = order - v_count - linear_count
            multinomial = comb(order, v_count) * comb(order - v_count, linear_count)
            term = multinomial * pow(parameter_d, constant_count, modulus)
            term *= -1 if linear_count % 2 else 1
            term *= factorials[linear_count + 2 * v_count]
            coefficients[v_count] = (coefficients[v_count] + term) % modulus
    return trim(coefficients, modulus)


def modular_gcd_audit():
    records = []
    for modulus in MODULAR_PRIMES:
        require(is_prime(modulus), modulus)
        positive = polynomial_gcd((1, 2, 1), (2, 3, 1), modulus)
        require(positive == (1, 1), (modulus, positive))
        for small_order in range(9):
            direct = moment_polynomial(small_order, small_order + 1, modulus)
            multinomial = multinomial_moment_polynomial(
                small_order, small_order + 1, modulus
            )
            require(direct == multinomial, (modulus, small_order))

        for prime, exponent in GCD_ROWS:
            height = prime ** exponent
            total_order = 10 * height
            parameter_d = total_order + 1
            require(modulus > 2 * total_order, (modulus, total_order))
            left = moment_polynomial(total_order - 1, parameter_d, modulus)
            right = moment_polynomial(total_order, parameter_d, modulus)
            require(len(left) - 1 == total_order - 1, (modulus, prime, "left degree"))
            require(len(right) - 1 == total_order, (modulus, prime, "right degree"))
            gcd = polynomial_gcd(left, right, modulus)
            require(len(gcd) == 1, (modulus, prime, exponent, gcd))
            limit = min(9, (prime - 1) // 2)
            record = (
                modulus,
                prime,
                exponent,
                total_order,
                set_submask_candidates(height, limit),
                len(gcd) - 1,
            )
            records.append(record)
    return tuple(records)


def main():
    cells, height_digest = independent_height_audit()
    print(
        f"explicit_support_universe=1<=H<{ODD_HEIGHT_BOUND},H_odd "
        f"limits={RESET_LIMITS} cells={cells} classification=True"
    )
    print(f"explicit_support_semantic_sha256={height_digest}")

    rows, prime_digest = independent_prime_power_audit()
    print(
        f"prime_power_universe=11<=p<{PRIME_BOUND},p_prime,"
        f"1<=k<={EXPONENT_BOUND} cells={rows} classification=True"
    )
    print(f"prime_power_semantic_sha256={prime_digest}")

    forcing, witnesses = independent_modulus_64_audit()
    print(f"independent_mod64_uniform_forcing_residues={forcing}")
    print(f"independent_mod64_nonforcing_witnesses={witnesses}")
    print("independent_mod64_maximality=True")

    records = modular_gcd_audit()
    print("polynomial_positive_control=gcd((x+1)^2,(x+1)(x+2))=x+1")
    for record in records:
        modulus, prime, exponent, total_order, candidates, gcd_degree = record
        print(
            f"modular_gcd=q:{modulus},p:{prime},k:{exponent},N:{total_order},"
            f"C:{candidates},gcd_degree:{gcd_degree}"
        )
    print("modular_gcd_hostiles=True degree_preservation=True small_multinomial_controls=True")


if __name__ == "__main__":
    main()
