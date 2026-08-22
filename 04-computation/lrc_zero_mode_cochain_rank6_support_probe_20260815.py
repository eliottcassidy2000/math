#!/usr/bin/env python3
"""Exact companion for THM-3416, the global zero-mode-cochain rank-six support.

The all-q proof has two inputs.  THM-3414 classifies fixed-zero covers by at
most six owners.  At half twist, capacity leaves four target-free exceptional
orders 3, 5, 17, 29.  Successively anchoring one exceptional block and pricing
the other five only on its complement rules out all four.  The 17/29 mixed
case is decided at lcm 493, not by a bounded-modulus conjecture.

An independent rare-coordinate cover solver audits both twists for every
2 <= Q <= 300.  The finite primitive census is a hostile control only; the
all-q theorem comes from the analytic tail bounds, exact finite order lists,
THM-3414, and divisor ancestry.  Runtime gates survive python -O.
"""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCAN_LIMIT = 300
LOWER_BASES = (8, 9, 10, 12)
RANK_SIX_BASES = (11, 15, 23, 25)
ALL_BASES = tuple(sorted(LOWER_BASES + RANK_SIX_BASES))
PINNED = (
    (
        "THM-3414-fixed-zero-six-owner-bases",
        ROOT / "01-canon/theorems/THM-3414-fixed-zero-six-owner-base-classification.md",
        "5568a4e93bc4478566335e2722c488c999797462eeb7c95af364b20dba41e998",
    ),
    (
        "THM-3415-global-rank-five",
        ROOT / "01-canon/theorems/THM-3415-zero-mode-cochain-global-rank-five-support.md",
        "de8d2f615695070dc75cad38ad4a9c22308d8bf900cd6f9a69cd2003f815a14d",
    ),
)

EXPECTED_EVENT_DIGEST = "354654b682bb9f9796e11e6f67cc4511bcb4d403a5ba3604e2471abbe14e0706"
EXPECTED_SEMANTIC_DIGEST = "18abbdc82a40b0299d9dc59cd745d52a09e9706d3cf0cb0379ab5ed2f064df22"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def prime_factors(value):
    factors = []
    trial = 2
    while trial * trial <= value:
        if value % trial == 0:
            factors.append(trial)
            while value % trial == 0:
                value //= trial
        trial += 1
    if value > 1:
        factors.append(value)
    return tuple(factors)


def divisors(value):
    return tuple(divisor for divisor in range(2, value + 1) if value % divisor == 0)


def danger_mask(q, residue, epsilon):
    modulus = 2 * q
    mask = 0
    for sheet in range(q):
        phase_word = residue * (2 * sheet + epsilon) % modulus
        centred = min(phase_word, modulus - phase_word)
        if 14 * centred < modulus:
            mask |= 1 << sheet
    return mask


def augmented_bank(q, epsilon):
    modulus = q if epsilon == 0 else 2 * q
    primes = prime_factors(modulus)
    grouped = {}
    raw = 0
    for residue in range(1, modulus):
        if residue % q == 0:
            continue
        sheet_mask = danger_mask(q, residue, epsilon)
        if not sheet_mask:
            continue
        raw += 1
        mask = sheet_mask
        for offset, prime in enumerate(primes):
            if residue % prime:
                mask |= 1 << (q + offset)
        old = grouped.get(mask)
        if old is None or residue < old:
            grouped[mask] = residue

    unique = tuple(sorted(((mask, residue) for mask, residue in grouped.items()), key=lambda x: x[1]))
    maximal = tuple(
        item
        for item in unique
        if not any(item[0] != other[0] and item[0] | other[0] == other[0] for other in unique)
    )
    full = (1 << (q + len(primes))) - 1
    return modulus, primes, raw, unique, maximal, full


def rare_coordinate_cover(full, items, cap):
    masks = tuple(mask for mask, _ in items)
    residues = tuple(residue for _, residue in items)
    width = full.bit_length()
    by_bit = tuple(
        tuple(index for index, mask in enumerate(masks) if mask & (1 << bit))
        for bit in range(width)
    )
    union = 0
    for mask in masks:
        union |= mask
    if union != full:
        return (), 0, 0

    nodes = 0
    branches = 0

    @lru_cache(maxsize=None)
    def solve(state, slots):
        nonlocal nodes, branches
        nodes += 1
        if state == full:
            return ()
        if slots == 0:
            return None
        missing = full ^ state
        gains = sorted(((mask & missing).bit_count() for mask in masks), reverse=True)
        if sum(gains[:slots]) < missing.bit_count():
            return None
        missing_bits = tuple(bit for bit in range(width) if missing & (1 << bit))
        pivot = min(
            missing_bits,
            key=lambda bit: (
                sum(1 for index in by_bit[bit] if masks[index] | state != state),
                bit,
            ),
        )
        candidates = sorted(
            (index for index in by_bit[pivot] if masks[index] | state != state),
            key=lambda index: (-(masks[index] & missing).bit_count(), residues[index]),
        )
        for index in candidates:
            branches += 1
            suffix = solve(state | masks[index], slots - 1)
            if suffix is not None:
                return (index,) + suffix
        return None

    chosen = solve(0, cap)
    if chosen is None:
        return (), nodes, branches
    witness = tuple(sorted(residues[index] for index in chosen))
    joined = 0
    for index in chosen:
        joined |= masks[index]
    require(joined == full and len(witness) <= cap, (cap, witness, joined, full))
    return witness, nodes, branches


def scan_record(q, epsilon):
    modulus, primes, raw, unique, maximal, full = augmented_bank(q, epsilon)
    witness5, nodes5, branches5 = rare_coordinate_cover(full, maximal, 5)
    if witness5:
        witness = witness5
        nodes6 = branches6 = 0
    else:
        witness, nodes6, branches6 = rare_coordinate_cover(full, maximal, 6)
    rank = len(witness) if witness else None
    orders = tuple(sorted(q // gcd(q, residue) for residue in witness))
    if witness:
        require(gcd(modulus, *witness) == 1, (q, epsilon, witness))
        require(len(set(witness)) == len(witness), (q, epsilon, witness))
    return (
        q,
        epsilon,
        raw,
        len(unique),
        len(maximal),
        rank,
        witness,
        orders,
        nodes5 + nodes6,
        branches5 + branches6,
        len(primes),
    )


def zero_block_count(order):
    return 1 + 2 * ((order - 1) // 14)


def half_block_count(order):
    odd_word_count = 2 * ((((order - 1) // 7) + 1) // 2)
    if order % 2 == 0:
        return odd_word_count
    return max(odd_word_count, zero_block_count(order))


def exact_order_masks(q, order, epsilon=1):
    modulus = q if epsilon == 0 else 2 * q
    return tuple(
        sorted(
            {
                danger_mask(q, residue, epsilon)
                for residue in range(1, modulus)
                if residue % q
                and q // gcd(q, residue) == order
                and danger_mask(q, residue, epsilon)
            }
        )
    )


def anchor_complement_record(anchor, order):
    q = lcm(anchor, order)
    full = (1 << q) - 1
    anchors = tuple(
        mask
        for mask in exact_order_masks(q, anchor)
        if mask.bit_count() * anchor == half_block_count(anchor) * q
    )
    others = exact_order_masks(q, order)
    require(anchors and others, (anchor, order, q, len(anchors), len(others)))
    maximum = max((other & (full ^ base)).bit_count() for base in anchors for other in others)
    return anchor, order, q, maximum, Fraction(maximum, q), len(anchors), len(others)


def critical_density_orders(threshold, last):
    require(Fraction(last + 1 + 6, 7 * (last + 1)) < threshold, (threshold, last))
    return tuple(
        order
        for order in range(2, last + 1)
        if Fraction(half_block_count(order), order) >= threshold
    )


def target_free(order):
    return not any(order % base == 0 for base in ALL_BASES)


def half_twist_gate_audit():
    # h(m) <= (m+6)/7 makes each list all-order after its displayed cutoff.
    six_critical = critical_density_orders(Fraction(1, 6), 36)
    require(
        six_critical == (3, 5, 8, 9, 10, 11, 12, 15, 17, 22, 23, 24, 29, 36),
        six_critical,
    )
    five_critical = critical_density_orders(Fraction(1, 5), 15)
    require(five_critical == (3, 5, 8, 9, 10, 15), five_critical)
    exceptional = tuple(order for order in six_critical if target_free(order))
    require(exceptional == (3, 5, 17, 29), exceptional)

    quotas = {
        3: Fraction(2, 15),
        5: Fraction(4, 25),
        17: Fraction(14, 85),
        29: Fraction(24, 145),
    }
    exceptional_matrix = tuple(
        anchor_complement_record(anchor, order)
        for anchor in exceptional
        for order in exceptional
    )
    by_pair = {(row[0], row[1]): row[4] for row in exceptional_matrix}
    strict_arcs = tuple(
        (anchor, order)
        for anchor in exceptional
        for order in exceptional
        if anchor != order and by_pair[(anchor, order)] < quotas[anchor]
    )
    require(
        strict_arcs == ((3, 17), (3, 29), (5, 17), (5, 29), (17, 29), (29, 17)),
        strict_arcs,
    )
    require(by_pair[(3, 5)] == quotas[3] and lcm(3, 5) == 15, by_pair)

    # When 3|m, every danger word runs through residue classes mod 3 in a
    # balanced interval, so at most ceil(2 h(m)/3) points avoid the anchor.
    # The analytic inequality is strict for m>=33; only 3,6,21 survive the
    # target-free finite check below 33.
    require(15 * (2 * 33 + 26) < 2 * 21 * 33, "order-three tail")
    divisible_three_small = tuple(
        order for order in range(3, 33, 3) if target_free(order)
    )
    require(divisible_three_small == (3, 6, 21), divisible_three_small)
    require(half_block_count(6) == 0, half_block_count(6))
    order_three_small = tuple(
        anchor_complement_record(3, order)
        for order in divisible_three_small
        if half_block_count(order)
    )
    require(all(row[4] < quotas[3] for row in order_three_small), order_three_small)
    coprime_three_five_density = tuple(
        order for order in five_critical if target_free(order) and order % 3
    )
    require(coprime_three_five_density == (5,), coprime_three_five_density)

    critical_5 = critical_density_orders(quotas[5], 50)
    candidates_5 = tuple(
        order for order in critical_5 if target_free(order) and order != 3
    )
    require(candidates_5 == (5, 17, 29, 31, 37, 43), candidates_5)
    controls_5 = tuple(anchor_complement_record(5, order) for order in candidates_5)
    require(all(row[4] < quotas[5] for row in controls_5), controls_5)

    critical_17 = critical_density_orders(quotas[17], 39)
    candidates_17 = tuple(
        order for order in critical_17 if target_free(order) and order not in (3, 5)
    )
    require(candidates_17 == (17, 29), candidates_17)
    controls_17 = tuple(anchor_complement_record(17, order) for order in candidates_17)
    require(all(row[4] < quotas[17] for row in controls_17), controls_17)
    require(by_pair[(17, 29)] == Fraction(70, 493), by_pair[(17, 29)])

    critical_29 = critical_density_orders(quotas[29], 37)
    candidates_29 = tuple(
        order
        for order in critical_29
        if target_free(order) and order not in (3, 5, 17)
    )
    require(candidates_29 == (29,), candidates_29)
    controls_29 = tuple(anchor_complement_record(29, order) for order in candidates_29)
    require(all(row[4] < quotas[29] for row in controls_29), controls_29)

    return (
        six_critical,
        five_critical,
        exceptional,
        tuple(sorted(quotas.items())),
        exceptional_matrix,
        strict_arcs,
        order_three_small,
        critical_5,
        controls_5,
        critical_17,
        controls_17,
        critical_29,
        controls_29,
    )


def augmented_witness_record(q, epsilon, residues):
    modulus = q if epsilon == 0 else 2 * q
    primes = prime_factors(modulus)
    full = (1 << (q + len(primes))) - 1
    joined = 0
    masks = []
    for residue in residues:
        mask = danger_mask(q, residue, epsilon)
        require(mask, (q, epsilon, residue))
        masks.append(mask)
        augmented = mask
        for offset, prime in enumerate(primes):
            if residue % prime:
                augmented |= 1 << (q + offset)
        joined |= augmented
    require(joined == full, (q, epsilon, residues, joined, full))
    require(gcd(modulus, *residues) == 1, (q, epsilon, residues))
    multiplicities = tuple(sum(mask >> sheet & 1 for mask in masks) for sheet in range(q))
    xor = 0
    for mask in masks:
        xor ^= mask
    return (
        q,
        epsilon,
        residues,
        tuple(q // gcd(q, residue) for residue in residues),
        tuple(mask.bit_count() for mask in masks),
        tuple(sorted(Counter(multiplicities).items())),
        xor.bit_count(),
    )


def construction_audit():
    atoms = (
        augmented_witness_record(11, 1, (1, 2, 3, 5, 7, 9)),
        augmented_witness_record(15, 0, (1, 2, 3, 4, 5, 7)),
        augmented_witness_record(15, 1, (1, 4, 6, 7, 8, 10)),
        augmented_witness_record(23, 1, (1, 4, 5, 7, 9, 11)),
        augmented_witness_record(25, 1, (1, 9, 10, 11, 19, 21)),
    )
    require(atoms[0][5:] == (((1, 11),), 11), atoms[0])
    require(atoms[1][5:] == (((1, 14), (6, 1)), 14), atoms[1])
    require(atoms[2][5:] == (((1, 14), (4, 1)), 14), atoms[2])
    require(atoms[3][5:] == (((1, 23),), 23), atoms[3])
    require(atoms[4][5:] == (((1, 25),), 25), atoms[4])
    return atoms


def in_rank_six_support(value):
    return (
        any(value % base == 0 for base in RANK_SIX_BASES)
        and not any(value % base == 0 for base in LOWER_BASES)
    )


def fibonacci_period(modulus):
    values = []
    previous, current = 0, 1
    while True:
        values.append(previous)
        previous, current = current, (previous + current) % modulus
        if (previous, current) == (0, 1):
            return tuple(values)


BERGGREN_MATRICES = (
    ((1, -2, 2), (2, -1, 2), (2, -2, 3)),
    ((1, 2, 2), (2, 1, 2), (2, 2, 3)),
    ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3)),
)


def matrix_child(vector, matrix, modulus):
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(3)) % modulus
        for row in range(3)
    )


def arithmetic_transport_audit():
    modulus = lcm(*ALL_BASES)
    require(modulus == 455400, modulus)
    support = tuple(residue for residue in range(modulus) if in_rank_six_support(residue))
    density = Fraction(len(support), modulus)
    require((len(support), density) == (55000, Fraction(25, 207)), (len(support), density))

    primes = prime_factors(modulus)
    jordan_two = modulus * modulus
    for prime in primes:
        jordan_two = jordan_two // (prime * prime) * (prime * prime - 1)
    farey_accepting = 0
    for denominator in support:
        numerator_count = modulus
        for prime in primes:
            if denominator % prime == 0:
                numerator_count = numerator_count // prime * (prime - 1)
        farey_accepting += numerator_count
    farey_density = Fraction(farey_accepting, jordan_two)
    require(
        (jordan_two, farey_accepting, farey_density)
        == (131383296000, 16682793600, Fraction(157981, 1244160)),
        (jordan_two, farey_accepting, farey_density),
    )

    fib = fibonacci_period(modulus)
    require(len(fib) == 1200, len(fib))
    fib_indices = tuple(index for index, value in enumerate(fib) if in_rank_six_support(value))
    fib_residues = tuple(index for index in range(150) if in_rank_six_support(fib[index]))
    expected_fib_residues = (10, 20, 25, 40, 50, 70, 80, 100, 110, 125, 130, 140)
    require(fib_residues == expected_fib_residues, fib_residues)
    require(
        all(
            in_rank_six_support(fib[index])
            == (
                (index % 10 == 0 or index % 25 == 0)
                and index % 6 != 0
                and index % 15 != 0
            )
            for index in range(len(fib))
        ),
        "Fibonacci index law",
    )
    require((len(fib_indices), Fraction(len(fib_indices), len(fib))) == (96, Fraction(2, 25)), len(fib_indices))

    u_residues = tuple(
        index
        for index in range(99)
        if in_rank_six_support(4 * index * index + 12 * index + 11)
    )
    expected_u = (0, 8, 11, 22, 30, 33, 44, 52, 63, 66, 74, 85, 88, 96)
    require(u_residues == expected_u, u_residues)
    require(
        all(
            in_rank_six_support(4 * index * index + 12 * index + 11)
            == (index % 11 in (0, 8) and index % 9 not in (1, 5))
            for index in range(99)
        ),
        "U-spine law",
    )

    odd_modulus = lcm(9, 11, 15, 23, 25)
    require(odd_modulus == 56925, odd_modulus)
    level = ((3, 4, 5),)
    level_counts = []
    for depth in range(11):
        accepting = sum(in_rank_six_support(2 * vector[2] + 1) for vector in level)
        level_counts.append((depth, len(level), len(set(level)), accepting))
        level = tuple(
            matrix_child(vector, matrix, odd_modulus)
            for vector in level
            for matrix in BERGGREN_MATRICES
        )
    require(
        tuple(level_counts)
        == (
            (0, 1, 1, 1),
            (1, 3, 3, 0),
            (2, 9, 9, 2),
            (3, 27, 27, 5),
            (4, 81, 81, 10),
            (5, 243, 243, 41),
            (6, 729, 729, 133),
            (7, 2187, 2187, 378),
            (8, 6561, 6561, 1210),
            (9, 19683, 19683, 3519),
            (10, 59049, 59045, 10634),
        ),
        level_counts,
    )

    return (
        modulus,
        len(support),
        density,
        jordan_two,
        farey_accepting,
        farey_density,
        len(fib),
        fib_residues,
        len(fib_indices),
        u_residues,
        tuple(level_counts),
    )


def main():
    dependencies = tuple((name, lf_hash(path)) for name, path, _ in PINNED)
    for name, path, expected in PINNED:
        require(lf_hash(path) == expected, ("dependency changed", name, lf_hash(path), expected))

    half_gate = half_twist_gate_audit()
    atoms = construction_audit()
    arithmetic = arithmetic_transport_audit()

    records = tuple(
        scan_record(q, epsilon)
        for q in range(2, SCAN_LIMIT + 1)
        for epsilon in (0, 1)
    )
    by_key = {(row[0], row[1]): row for row in records}
    primitive_rank_six = tuple(
        (row[0], row[1], row[6], row[7]) for row in records if row[5] == 6
    )
    for q in range(2, SCAN_LIMIT + 1):
        observed_at_most_six = any(
            by_key[(divisor, epsilon)][5] is not None
            for divisor in divisors(q)
            for epsilon in (0, 1)
        )
        theoretical_at_most_six = any(q % base == 0 for base in ALL_BASES)
        require(observed_at_most_six == theoretical_at_most_six, (q, observed_at_most_six))

    exact_rank_six = tuple(q for q in range(2, SCAN_LIMIT + 1) if in_rank_six_support(q))
    require(exact_rank_six[:4] == (11, 15, 22, 23), exact_rank_six[:10])
    require(25 in exact_rank_six and 17 not in exact_rank_six and 29 not in exact_rank_six, exact_rank_six[:20])

    event_surface = tuple((row[0], row[1], row[5], row[6], row[7]) for row in records)
    event_digest = sha256(repr(event_surface).encode("ascii")).hexdigest()
    semantic_surface = (
        half_gate,
        atoms,
        arithmetic,
        exact_rank_six,
        primitive_rank_six,
    )
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_EVENT_DIGEST is not None:
        require(event_digest == EXPECTED_EVENT_DIGEST, (event_digest, EXPECTED_EVENT_DIGEST))
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST, (semantic_digest, EXPECTED_SEMANTIC_DIGEST))

    total_raw = sum(row[2] for row in records)
    total_unique = sum(row[3] for row in records)
    total_maximal = sum(row[4] for row in records)
    total_nodes = sum(row[8] for row in records)
    total_branches = sum(row[9] for row in records)
    primitive_keys = tuple((row[0], row[1]) for row in primitive_rank_six)

    print("THM-3416 zero-mode-cochain global rank-six support exact companion")
    print(f"dependency_sha256_lf={dependencies}")
    print("status=PROVED all_q_global_rank6_downstream_of_THM3414_computer_assisted_fixed_zero_atlas_plus_elementary_half_twist_anchor_complements;FINITE_EXACT_Q2_Q300_hostile_census;no_LRC14_decrement")
    print("global_rank6=rho_ZMC(q)=6_iff_(11|q_or_15|q_or_23|q_or_25|q)_and_8,9,10,12_all_not_divide_q")
    print(f"half_six_capacity=(critical_orders,target_free_exceptional)={(half_gate[0], half_gate[2])}")
    print(f"exceptional_anchor_quotas={half_gate[3]}")
    print(f"exceptional_complement_matrix={half_gate[4]}")
    print(f"generalized_tournament=(vertices,strict_arcs,missing_target_tie,bidirected_pair)={(half_gate[2], half_gate[5], ((3,5,'equality',15),), ((17,29),))}")
    print(f"order3_divisible_small_controls={half_gate[6]}")
    print(f"order5_controls=(critical,candidates)={(half_gate[7], half_gate[8])}")
    print(f"order17_controls=(critical,candidates)={(half_gate[9], half_gate[10])}")
    print(f"order29_controls=(critical,candidates)={(half_gate[11], half_gate[12])}")
    print(f"rank6_atoms=(q,epsilon,residues,orders,sizes,multiplicity_hist,xor_size)={atoms}")
    print(f"scan_Q2_Q300=(raw,unique,maximal,nodes,branches)={(total_raw,total_unique,total_maximal,total_nodes,total_branches)}")
    print(f"primitive_rank6_Q2_Q300=(count,keys,digest)={(len(primitive_keys),primitive_keys,sha256(repr(primitive_rank_six).encode('ascii')).hexdigest())}")
    print(f"global_rank6_Q2_Q300=(count,first,last,digest)={(len(exact_rank_six),exact_rank_six[:20],exact_rank_six[-10:],sha256(repr(exact_rank_six).encode('ascii')).hexdigest())}")
    print(f"support_mod455400=(classes,density,harmonic_coefficient)={(arithmetic[1],arithmetic[2],arithmetic[2])}")
    print(f"Farey_Stern_Brocot_mod455400=(primitive_states,accepting,density)={(arithmetic[3],arithmetic[4],arithmetic[5])}")
    print(f"Fibonacci_transport=(period,index_residues_mod150,count_per_period,density)={(arithmetic[6],arithmetic[7],arithmetic[8],Fraction(arithmetic[8],arithmetic[6]))}")
    print(f"Berggren_U_spine=(rank6_n_mod99,density)={(arithmetic[9],Fraction(len(arithmetic[9]),99))}")
    print(f"Berggren_ternary_levels_mod56925=(depth,nodes,states,rank6)={arithmetic[10]}")
    print(f"event_sha256={event_digest}")
    print(f"semantic_sha256={semantic_digest}")
    print("scope=unrestricted_distinct_positive_transverse_zero_complete_mode_cochains;OR_equals_XOR_only_for_q11,q23,q25_sheet_atoms;generalized_tournament_is_an_exceptional_order_proof_sidecar_not_the_cover_carrier;LRC14_open")


if __name__ == "__main__":
    main()
