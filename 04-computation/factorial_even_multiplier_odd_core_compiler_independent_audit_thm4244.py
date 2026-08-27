#!/usr/bin/env python3
"""Structurally independent exact audit for the general pair compiler.

Binary data are explicit sets of occupied positions, complement blocks are
built as orbits rather than smaller representatives, factorial coefficients
come from an explicit three-slot multinomial sum, and SymPy supplies the
finite-field polynomial gcd.  RuntimeError gates remain active under -O.
"""

from hashlib import sha256
from math import comb, factorial

from sympy import Poly, symbols


MAX_A_EXCLUSIVE = 130
PRIME_BOUND = 280
EXPONENT_BOUND = 9
SUFFIX_MAX_A = 192
SUFFIX_HEIGHT_MULTIPLIERS = 96
A14_HEIGHT_BOUND = 1 << 17
MODULAR_PRIMES = (1000003, 1000033)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def occupied_positions(integer):
    positions = set()
    place = 0
    while integer:
        integer, digit = divmod(integer, 2)
        if digit:
            positions.add(place)
        place += 1
    return frozenset(positions)


def odd_primes_below(bound):
    primes = []
    for candidate in range(3, bound, 2):
        composite = False
        for divisor in primes:
            if divisor * divisor > candidate:
                break
            if candidate % divisor == 0:
                composite = True
                break
        if not composite:
            primes.append(candidate)
    return tuple(primes)


def direct_support_candidates(multiplier, height, reset_limit):
    target = occupied_positions(multiplier * height)
    return tuple(
        test
        for test in range(1, reset_limit + 1)
        if occupied_positions(test * height) <= target
    )


def complement_orbits(multiplier):
    unseen = set(range(2, multiplier, 2))
    orbits = []
    while unseen:
        test = min(unseen)
        complement = multiplier - test
        orbit = tuple(sorted({test, complement}))
        unseen.difference_update(orbit)
        orbits.append(orbit)
    return tuple(orbits)


def orbit_prediction(multiplier, height, reset_limit):
    answer = []
    for orbit in complement_orbits(multiplier):
        if len(orbit) == 1:
            continue
        left, right = orbit
        disjoint = not (
            occupied_positions(left * height) & occupied_positions(right * height)
        )
        if disjoint:
            answer.extend(test for test in orbit if test <= reset_limit)
    return tuple(sorted(answer))


def two_adic_decomposition(integer):
    power = 0
    while integer % 2 == 0:
        integer //= 2
        power += 1
    return power, integer


def odd_core_orbit_prediction(multiplier, height, reset_limit):
    power, odd_core = two_adic_decomposition(multiplier)
    scale = 1 << power
    answer = []
    for reduced in range(1, (odd_core - 1) // 2 + 1):
        reduced_complement = odd_core - reduced
        overlap = (
            occupied_positions(reduced * height)
            & occupied_positions(reduced_complement * height)
        )
        if not overlap:
            answer.extend(
                visible
                for visible in (scale * reduced, scale * reduced_complement)
                if visible <= reset_limit
            )
    return tuple(sorted(answer))


def support_barcode_audit():
    digest = sha256()
    cells = 0
    fixed_points = 0
    one_visible_blocks = 0
    two_visible_blocks = 0
    for multiplier in range(2, MAX_A_EXCLUSIVE, 2):
        orbits = complement_orbits(multiplier)
        expected_even_singletons = 1 if multiplier % 4 == 0 else 0
        require(sum(1 for orbit in orbits if len(orbit) == 1) == expected_even_singletons, (
            multiplier,
            orbits,
        ))
        for prime in odd_primes_below(PRIME_BOUND):
            if prime <= multiplier or prime < 5:
                continue
            reset_limit = min(multiplier - 1, (prime - 1) // 2)
            require(2 * reset_limit >= multiplier, (multiplier, prime, reset_limit))
            for exponent in range(1, EXPONENT_BOUND + 1):
                height = prime**exponent
                direct = direct_support_candidates(multiplier, height, reset_limit)
                predicted = orbit_prediction(multiplier, height, reset_limit)
                core_predicted = odd_core_orbit_prediction(
                    multiplier,
                    height,
                    reset_limit,
                )
                require(direct == predicted, (multiplier, prime, exponent, direct, predicted))
                require(direct == core_predicted, (
                    multiplier,
                    prime,
                    exponent,
                    direct,
                    core_predicted,
                ))

                fixed = multiplier // 2
                require(fixed <= reset_limit, (multiplier, prime, fixed, reset_limit))
                require(fixed not in direct, (multiplier, prime, exponent, fixed))
                fixed_support = occupied_positions(fixed * height)
                require(bool(fixed_support), (multiplier, height, fixed))
                require(fixed_support & fixed_support == fixed_support, fixed_support)
                fixed_points += 1

                nontrivial_overlap_flags = []
                for orbit in orbits:
                    if len(orbit) == 1:
                        require(orbit[0] == fixed, (multiplier, orbit, fixed))
                        continue
                    left, right = orbit
                    require(left < multiplier // 2 < right, (multiplier, orbit))
                    require(left <= reset_limit, (multiplier, prime, orbit, reset_limit))
                    overlap = occupied_positions(left * height) & occupied_positions(right * height)
                    nontrivial_overlap_flags.append(bool(overlap))
                    visible = tuple(test for test in orbit if test <= reset_limit)
                    require(len(visible) in (1, 2), (multiplier, prime, orbit, visible))
                    if len(visible) == 1:
                        one_visible_blocks += 1
                    else:
                        two_visible_blocks += 1
                require((direct == ()) == all(nontrivial_overlap_flags), (
                    multiplier,
                    prime,
                    exponent,
                    direct,
                    nontrivial_overlap_flags,
                ))
                digest.update(f"{multiplier}:{prime}:{exponent}:{reset_limit}:{direct}\n".encode())
                cells += 1
    return (
        cells,
        fixed_points,
        one_visible_blocks,
        two_visible_blocks,
        digest.hexdigest(),
    )


def independent_suffix_audit():
    digest = sha256()
    heights = 0
    orbit_checks = 0
    for multiplier in range(2, SUFFIX_MAX_A + 1, 2):
        power, odd_core = two_adic_decomposition(multiplier)
        exponent = 1
        half_modulus = 1
        while half_modulus < odd_core:
            half_modulus *= 2
            exponent += 1
        modulus = 2 * half_modulus
        for height_multiplier in range(1, SUFFIX_HEIGHT_MULTIPLIERS + 1):
            height = height_multiplier * modulus - 1
            for reduced in range(1, (odd_core - 1) // 2 + 1):
                reduced_complement = odd_core - reduced
                left_word = occupied_positions((reduced * height) % modulus)
                right_word = occupied_positions((reduced_complement * height) % modulus)
                common_position = exponent - 1
                require(common_position in left_word and common_position in right_word, (
                    multiplier,
                    height,
                    reduced,
                    reduced_complement,
                    left_word,
                    right_word,
                ))
                reduced_overlap = (
                    occupied_positions(reduced * height)
                    & occupied_positions(reduced_complement * height)
                )
                require(common_position in reduced_overlap, (
                    multiplier,
                    height,
                    reduced,
                    reduced_complement,
                    reduced_overlap,
                ))
                original_overlap = (
                    occupied_positions((1 << power) * reduced * height)
                    & occupied_positions((1 << power) * reduced_complement * height)
                )
                require(
                    original_overlap
                    == frozenset(position + power for position in reduced_overlap),
                    (
                        multiplier,
                        height,
                        reduced,
                        reduced_complement,
                        original_overlap,
                        reduced_overlap,
                    ),
                )
                orbit_checks += 1
            require(
                direct_support_candidates(multiplier, height, multiplier - 1) == (),
                (multiplier, odd_core, exponent, height),
            )
            digest.update(
                f"{multiplier}:{power}:{odd_core}:{exponent}:{height}\n".encode()
            )
            heights += 1
    return heights, orbit_checks, digest.hexdigest()


def independent_power_two_invariance():
    digest = sha256()
    cells = 0
    shift_checks = 0
    for odd_core in range(1, 26, 2):
        for height in range(1, 384, 2):
            reference = None
            for power in range(1, 6):
                scale = 1 << power
                multiplier = scale * odd_core
                direct = direct_support_candidates(multiplier, height, multiplier - 1)
                compiled = odd_core_orbit_prediction(multiplier, height, multiplier - 1)
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
                        occupied_positions(reduced * height)
                        & occupied_positions(reduced_complement * height)
                    )
                    original_overlap = (
                        occupied_positions(scale * reduced * height)
                        & occupied_positions(scale * reduced_complement * height)
                    )
                    shifted = frozenset(position + power for position in reduced_overlap)
                    require(original_overlap == shifted, (
                        odd_core,
                        height,
                        power,
                        reduced,
                        original_overlap,
                        shifted,
                    ))
                    shift_checks += 1
                digest.update(
                    f"{odd_core}:{height}:{power}:{normalized}\n".encode()
                )
                cells += 1
    return cells, shift_checks, digest.hexdigest()


def independent_a14_audit():
    digest = sha256()
    cells = 0
    for height in range(7, A14_HEIGHT_BOUND, 8):
        actual = []
        for orbit in complement_orbits(14):
            if len(orbit) == 1:
                continue
            left, right = orbit
            overlap = occupied_positions(left * height) & occupied_positions(right * height)
            require(bool(overlap), (height, orbit, overlap))
            actual.append(tuple(sorted(overlap)))
        require(direct_support_candidates(14, height, 13) == (), (height, actual))
        digest.update(f"{height}:{actual}\n".encode())
        cells += 1

    prime_cells = 0
    for prime in odd_primes_below(1400):
        if prime <= 14 or prime % 8 != 7:
            continue
        reset_limit = min(13, (prime - 1) // 2)
        for exponent in range(1, 18, 2):
            height = prime**exponent
            require(direct_support_candidates(14, height, reset_limit) == (), (
                prime,
                exponent,
            ))
            digest.update(f"P:{prime}:{exponent}:{reset_limit}\n".encode())
            prime_cells += 1
    return cells, prime_cells, digest.hexdigest()


def factorial_table(bound, modulus):
    table = [1]
    for integer in range(1, bound + 1):
        table.append(table[-1] * integer % modulus)
    return tuple(table)


def multinomial_moment_coefficients(order, parameter_d, modulus):
    factorials = factorial_table(2 * order, modulus)
    coefficients = [0] * (order + 1)
    for quadratic_count in range(order + 1):
        for linear_count in range(order - quadratic_count + 1):
            constant_count = order - quadratic_count - linear_count
            ways = comb(order, quadratic_count) * comb(
                order - quadratic_count,
                linear_count,
            )
            contribution = ways
            contribution *= pow(parameter_d, constant_count, modulus)
            contribution *= -1 if linear_count % 2 else 1
            contribution *= factorials[2 * quadratic_count + linear_count]
            coefficients[quadratic_count] = (
                coefficients[quadratic_count] + contribution
            ) % modulus
    return tuple(coefficients)


def word_expansion_control(order, parameter_d, modulus):
    states = {(0, 0): 1}
    for _ in range(order):
        next_states = {}
        for (x_degree, v_degree), coefficient in states.items():
            additions = (
                ((x_degree, v_degree), coefficient * parameter_d),
                ((x_degree + 1, v_degree), -coefficient),
                ((x_degree + 2, v_degree + 1), coefficient),
            )
            for address, value in additions:
                next_states[address] = (next_states.get(address, 0) + value) % modulus
        states = next_states
    coefficients = [0] * (order + 1)
    for (x_degree, v_degree), coefficient in states.items():
        coefficients[v_degree] = (
            coefficients[v_degree] + coefficient * factorial(x_degree)
        ) % modulus
    return tuple(coefficients)


def independent_modular_nonconverse():
    variable = symbols("v")
    multiplier = 14
    prime = 17
    exponent = 1
    height = prime**exponent
    total_order = multiplier * height
    parameter_d = total_order + 1
    reset_limit = min(multiplier - 1, (prime - 1) // 2)
    candidates = direct_support_candidates(multiplier, height, reset_limit)
    require(candidates == (2, 4, 6, 8), candidates)
    records = []
    for modulus in MODULAR_PRIMES:
        require(modulus > 2 * total_order, (modulus, total_order))
        positive_left = Poly(variable**2 + 2 * variable + 1, variable, modulus=modulus)
        positive_right = Poly(variable**2 + 3 * variable + 2, variable, modulus=modulus)
        require(positive_left.gcd(positive_right).degree() == 1, modulus)
        for small_order in range(9):
            multi = multinomial_moment_coefficients(
                small_order,
                small_order + 1,
                modulus,
            )
            word = word_expansion_control(small_order, small_order + 1, modulus)
            require(multi == word, (modulus, small_order, multi, word))

        left_coefficients = multinomial_moment_coefficients(
            total_order - 1,
            parameter_d,
            modulus,
        )
        right_coefficients = multinomial_moment_coefficients(
            total_order,
            parameter_d,
            modulus,
        )
        left = Poly.from_list(list(reversed(left_coefficients)), gens=variable, modulus=modulus)
        right = Poly.from_list(list(reversed(right_coefficients)), gens=variable, modulus=modulus)
        require(left.degree() == total_order - 1, (modulus, left.degree()))
        require(right.degree() == total_order, (modulus, right.degree()))
        common = left.gcd(right)
        require(common.degree() == 0, (modulus, common))
        records.append((modulus, total_order, left.degree(), right.degree(), common.degree()))
    return candidates, tuple(records)


def hostile_records():
    cells = (
        (2, 5, 1),
        (4, 5, 1),
        (12, 17, 1),
        (14, 17, 1),
        (14, 19, 1),
        (14, 23, 1),
        (14, 23, 2),
        (14, 29, 2),
        (14, 73, 1),
    )
    records = []
    for multiplier, prime, exponent in cells:
        height = prime**exponent
        reset_limit = min(multiplier - 1, (prime - 1) // 2)
        direct = direct_support_candidates(multiplier, height, reset_limit)
        predicted = orbit_prediction(multiplier, height, reset_limit)
        require(direct == predicted, (multiplier, prime, exponent, direct, predicted))
        records.append((multiplier, prime, exponent, height, reset_limit, direct))
    return tuple(records)


def main():
    overall = sha256()
    cells, fixed, one_visible, two_visible, barcode_digest = support_barcode_audit()
    print(
        "support_barcode_universe="
        f"even_2<=a<{MAX_A_EXCLUSIVE},primes_a<p<{PRIME_BOUND},"
        f"1<=k<={EXPONENT_BOUND}"
    )
    print(
        f"support_barcode_cells={cells} fixed_orbit_checks={fixed} "
        f"one_visible_blocks={one_visible} two_visible_blocks={two_visible}"
    )
    print(f"support_barcode_semantic_sha256={barcode_digest}")
    overall.update(barcode_digest.encode())

    suffix_heights, suffix_orbits, suffix_digest = independent_suffix_audit()
    print(
        "support_suffix_universe="
        f"even_2<=a<={SUFFIX_MAX_A},a=2^q*b,2^(ell-1)>=b,"
        f"H=m*2^ell-1,1<=m<={SUFFIX_HEIGHT_MULTIPLIERS}"
    )
    print(f"support_suffix_heights={suffix_heights} orbit_checks={suffix_orbits}")
    print(f"support_suffix_semantic_sha256={suffix_digest}")
    overall.update(suffix_digest.encode())

    invariant_cells, invariant_shifts, invariant_digest = independent_power_two_invariance()
    print(
        "support_power_two_invariance_universe="
        "odd_1<=b<=25,odd_1<=H<384,1<=q<=5"
    )
    print(
        f"support_power_two_invariance_cells={invariant_cells} "
        f"overlap_shift_checks={invariant_shifts}"
    )
    print(f"support_power_two_invariance_semantic_sha256={invariant_digest}")
    overall.update(invariant_digest.encode())

    a14_cells, prime_cells, a14_digest = independent_a14_audit()
    print(f"support_a14_cells={a14_cells} prime_odd_exponent_cells={prime_cells}")
    print(f"support_a14_semantic_sha256={a14_digest}")
    overall.update(a14_digest.encode())

    for record in hostile_records():
        print(f"support_hostile={record}")
        overall.update(f"{record}\n".encode())

    candidates, modular_records = independent_modular_nonconverse()
    print(f"sympy_nonconverse_candidate_row=a:14,p:17,k:1,C:{candidates}")
    for record in modular_records:
        print(
            "sympy_modular_nonconverse="
            f"q:{record[0]},N:{record[1]},degrees:{record[2]}/{record[3]},"
            f"gcd_degree:{record[4]}"
        )
        overall.update(f"{record}\n".encode())

    print(f"overall_semantic_sha256={overall.hexdigest()}")
    print("status=PASS runtime_require_gates=True optimization_invariant=True")


if __name__ == "__main__":
    main()
