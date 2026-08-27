#!/usr/bin/env python3
"""Primary exact audit for THM-4243.

This companion checks the candidate barcode directly with integer bit
operations, verifies the closed two-pair formula on a declared exhaustive
universe, audits prime powers, and computes the maximal uniform suffix
certificate modulo 64.  It does not construct factorial-moment coefficients.
"""

from hashlib import sha256


ODD_HEIGHT_BOUND = 1 << 18
PRIME_BOUND = 1000
EXPONENT_BOUND = 30
RESET_LIMITS = (5, 6, 8, 9)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def primes_below(bound):
    result = []
    for candidate in range(2, bound):
        prime = True
        divisor = 2
        while divisor * divisor <= candidate:
            if candidate % divisor == 0:
                prime = False
                break
            divisor += 1
        if prime:
            result.append(candidate)
    return result


def is_submask(left, right):
    return left & ~right == 0


def direct_candidates(height, limit):
    return tuple(
        multiplier
        for multiplier in range(1, limit + 1)
        if is_submask(multiplier * height, 10 * height)
    )


def first_pair_overlaps(height):
    return (2 * height) & (8 * height) != 0


def second_pair_overlaps(height):
    return (4 * height) & (6 * height) != 0


def predicted_candidates(height, limit):
    first_is_carry_free = not first_pair_overlaps(height)
    second_is_carry_free = not second_pair_overlaps(height)
    answer = []
    for multiplier in range(1, limit + 1):
        if multiplier in (2, 8) and first_is_carry_free:
            answer.append(multiplier)
        if multiplier in (4, 6) and second_is_carry_free:
            answer.append(multiplier)
    return tuple(answer)


def closes(height):
    return first_pair_overlaps(height) and second_pair_overlaps(height)


def exhaustive_height_audit():
    digest = sha256()
    cells = 0
    empty_by_limit = {}
    candidate_shapes = {}
    for limit in RESET_LIMITS:
        empty = 0
        shapes = {}
        for height in range(1, ODD_HEIGHT_BOUND, 2):
            direct = direct_candidates(height, limit)
            predicted = predicted_candidates(height, limit)
            require(direct == predicted, (height, limit, direct, predicted))
            require(all(multiplier % 2 == 0 for multiplier in direct), (height, direct))
            require((direct == ()) == closes(height), (height, limit, direct))
            if height % 8 in (5, 7):
                require(direct == (), (height, limit, direct))
            empty += direct == ()
            shapes[direct] = shapes.get(direct, 0) + 1
            digest.update(f"{height}:{limit}:{direct}\n".encode())
            cells += 1
        empty_by_limit[limit] = empty
        candidate_shapes[limit] = tuple(sorted(shapes.items()))
    return cells, empty_by_limit, candidate_shapes, digest.hexdigest()


def maximal_modulus_64_suffix_audit():
    # If H changes by 64q, each of 2H,4H,6H,8H changes by a multiple
    # of 128.  Therefore their bits below position seven depend only on
    # H modulo 64.
    require(all((multiplier * 64) % 128 == 0 for multiplier in (2, 4, 6, 8)),
            "low-product invariance failed")

    forcing = []
    nonforcing_witnesses = []
    for residue in range(1, 64, 2):
        first_low_overlap = ((2 * residue) % 128) & ((8 * residue) % 128)
        second_low_overlap = ((4 * residue) % 128) & ((6 * residue) % 128)
        if first_low_overlap and second_low_overlap:
            forcing.append(residue)
            for quotient in range(128):
                require(closes(residue + 64 * quotient), (residue, quotient))
        else:
            # The least positive representative is an exact witness that
            # this whole residue class cannot be a uniform forcing class.
            require(not closes(residue), residue)
            nonforcing_witnesses.append((residue, direct_candidates(residue, 9)))

    expected = tuple(
        residue
        for residue in range(1, 64, 2)
        if residue % 8 in (5, 7)
        or residue % 32 == 27
        or residue in (41, 57)
    )
    require(tuple(forcing) == expected, (forcing, expected))
    require(len(forcing) + len(nonforcing_witnesses) == 32, "residue partition")
    return tuple(forcing), tuple(nonforcing_witnesses)


def prime_power_audit():
    digest = sha256()
    rows = 0
    root_family_rows = 0
    closed_rows = 0
    for prime in primes_below(PRIME_BOUND):
        if prime < 11:
            continue
        limit = min(9, (prime - 1) // 2)
        for exponent in range(1, EXPONENT_BOUND + 1):
            height = prime ** exponent
            direct = direct_candidates(height, limit)
            predicted = predicted_candidates(height, limit)
            require(direct == predicted, (prime, exponent, direct, predicted))
            if prime % 8 in (5, 7) and exponent % 2 == 1:
                require(height % 8 in (5, 7), (prime, exponent, height % 8))
                require(direct == (), (prime, exponent, direct))
                root_family_rows += 1
            closed_rows += direct == ()
            digest.update(
                f"{prime}:{exponent}:{height & ((1 << 64) - 1)}:{direct}\n".encode()
            )
            rows += 1
    return rows, root_family_rows, closed_rows, digest.hexdigest()


def congruence_family_audit():
    power_27_table = {}
    for exponent_residue in (1, 3, 5, 7):
        bases = tuple(
            residue
            for residue in range(1, 32, 2)
            if pow(residue, exponent_residue, 32) == 27
        )
        power_27_table[exponent_residue] = bases
    expected_table = {1: (27,), 3: (3,), 5: (11,), 7: (19,)}
    require(power_27_table == expected_table, (power_27_table, expected_table))

    square_bases = tuple(
        residue
        for residue in range(1, 32, 2)
        if pow(residue, 2, 64) in (41, 57)
    )
    require(square_bases == (11, 13, 19, 21), square_bases)

    # Strictness of the modulus-64 suffix observer: this admissible
    # prime-power height closes by higher bits although its residue class
    # is not uniformly forcing modulo 64.
    outside_height = 11 ** 3
    require(outside_height % 64 == 51, outside_height % 64)
    require(closes(outside_height), outside_height)
    require(not closes(51), "least residue must remain the hostile witness")
    return power_27_table, square_bases, outside_height


def main():
    cells, empty_counts, shapes, height_digest = exhaustive_height_audit()
    print(
        f"odd_height_universe=1<=H<{ODD_HEIGHT_BOUND},H_odd "
        f"limits={RESET_LIMITS} cells={cells} exact_classification=True"
    )
    print(f"empty_counts_by_limit={empty_counts}")
    print(f"candidate_shape_histograms={shapes}")
    print(f"height_semantic_sha256={height_digest}")

    forcing, witnesses = maximal_modulus_64_suffix_audit()
    print(f"mod64_uniform_forcing_residues={forcing}")
    print(f"mod64_nonforcing_witnesses={witnesses}")
    print("mod64_maximality=True distinction=uniform_suffix_not_unbounded_classification")

    rows, root_rows, closed_rows, pp_digest = prime_power_audit()
    print(
        f"prime_power_universe=11<=p<{PRIME_BOUND},p_prime,"
        f"1<=k<={EXPONENT_BOUND} cells={rows} root_family_rows={root_rows} "
        f"closed_rows={closed_rows} classification=True"
    )
    print(f"prime_power_semantic_sha256={pp_digest}")

    table, squares, outside = congruence_family_audit()
    print(f"p_power_27_mod32_table={table}")
    print(f"square_family_p_mod32={squares}")
    print(
        f"strict_suffix_hostile=H:{outside},H_mod64:{outside % 64},"
        f"C:{direct_candidates(outside, 9)},residue_51_C:{direct_candidates(51, 9)}"
    )


if __name__ == "__main__":
    main()
