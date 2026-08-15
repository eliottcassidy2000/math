#!/usr/bin/env python3
"""Exact companion for the zero-mode-cochain global rank-five support.

This is deliberately a different search route from the union-state BFS and
combination census in lrc_unrestricted_zero_mode_cochain_rank_probe_20260815.
It branches on a rare uncovered sheet/prime-breaker coordinate, memoizes only
the covered state and remaining budget, and tracks no quotient-order profile
assumption.  The exact Q<=500 census is a hostile audit, not the all-q proof.

The all-q proof uses exact block-density formulas.  At half twist, either an
order-3 block is absent and five-cover capacity forces orders 8,9,10, or an
order-3 block is present and its complement forces orders 8,9,12.  At zero
twist, THM-3408 excludes base-free five-covers; the only residual target-free
case is odd with 15|Q and v_3(Q)=1, where the order-3 complement is too large
for four blocks.  Together with the rank-four theorem and Q=10,12 partitions,
this proves global rank five exactly on

    (10|q or 12|q) and 8 not|q and 9 not|q.

No LRC(14) ledger consequence is claimed.  Runtime gates survive python -O.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCAN_LIMIT = 500
PINNED = (
    (
        "THM-3405-zero-mode-gauge",
        ROOT / "01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md",
        "d3e7dbeeb85c6f897bd9e31270bd0b6602ae4feac3b46a45eb5ce23ae5d24fe0",
    ),
    (
        "THM-3408-fixed-zero-duality",
        ROOT / "01-canon/theorems/THM-3408-fixed-zero-additive-order-duality-and-six-core-corridor.md",
        "895d77c09fc3c19b3d06799b0109432d10116b8e4a07c9d53d75c9e4e3b11883",
    ),
    (
        "rank-four-source",
        ROOT / "04-computation/lrc_unrestricted_zero_mode_cochain_rank_probe_20260815.py",
        "70c176e1a056d285471a07d1d011a26070ff288c3cee5ec39971d349d416de31",
    ),
    (
        "rank-four-output-MISTAKE-391-repaired",
        ROOT / "05-knowledge/results/lrc_unrestricted_zero_mode_cochain_rank_probe_20260815.out",
        "52684d84bba6076c760285937e29cfd4a81c998324d6a9019e2919d4f764ab5d",
    ),
)

EXPECTED_EVENT_DIGEST = "28f6ed1fa9322355bbde1b302f27536ac9ceef62cdb7c532cde4e3dd10369798"
EXPECTED_SEMANTIC_DIGEST = "8c4a9bc6fee718a2c6f6ec282b14bb7f7968aec40ec37b9c8d12a3e13e8528dc"


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
    all_union = 0
    for mask in masks:
        all_union |= mask
    if all_union != full:
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
        gains = sorted(
            ((mask & missing).bit_count() for mask in masks),
            reverse=True,
        )
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
        candidates = tuple(
            sorted(
                (index for index in by_bit[pivot] if masks[index] | state != state),
                key=lambda index: (-(masks[index] & missing).bit_count(), residues[index]),
            )
        )
        for index in candidates:
            branches += 1
            answer = solve(state | masks[index], slots - 1)
            if answer is not None:
                return (index,) + answer
        return None

    chosen = solve(0, cap)
    if chosen is None:
        return (), nodes, branches
    witness = tuple(sorted(residues[index] for index in chosen))
    joined = 0
    for index in chosen:
        joined |= masks[index]
    require(joined == full, (cap, witness, joined, full))
    require(len(witness) <= cap, (cap, witness))
    return witness, nodes, branches


def record(q, epsilon):
    modulus, primes, raw, unique, maximal, full = augmented_bank(q, epsilon)
    witness4, nodes4, branches4 = rare_coordinate_cover(full, maximal, 4)
    if witness4:
        rank = len(witness4)
        witness = witness4
        nodes5 = branches5 = 0
    else:
        witness5, nodes5, branches5 = rare_coordinate_cover(full, maximal, 5)
        rank = len(witness5) if witness5 else None
        witness = witness5
    orders = tuple(q // gcd(q, residue) for residue in witness)
    augmented_gcd = gcd(modulus, *witness) if witness else None
    if witness:
        require(augmented_gcd == 1, (q, epsilon, witness, augmented_gcd))
        require(len(set(witness)) == len(witness), (q, epsilon, witness))
    return (
        q,
        epsilon,
        raw,
        len(unique),
        len(maximal),
        rank,
        witness,
        tuple(sorted(orders)),
        nodes4 + nodes5,
        branches4 + branches5,
    )


def zero_block_count(order):
    return 1 + 2 * ((order - 1) // 14)


def half_block_count(order):
    odd_word_count = 2 * ((((order - 1) // 7) + 1) // 2)
    if order % 2 == 0:
        return odd_word_count
    return max(odd_word_count, zero_block_count(order))


def exact_order_masks(q, order, epsilon):
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


def order_three_complement_record(order, epsilon):
    q = lcm(3, order)
    order_three = exact_order_masks(q, 3, epsilon)
    require(len(order_three) == 1, (q, order, epsilon, order_three))
    full = (1 << q) - 1
    complement = full ^ order_three[0]
    masks = exact_order_masks(q, order, epsilon)
    maximum = max(((mask & complement).bit_count() for mask in masks), default=0)
    return order, q, maximum, Fraction(maximum, q), len(masks)


def density_and_complement_audit():
    # Both z(m) and h(m) are at most (m+6)/7.  Therefore m>=16 is
    # automatically below density 1/5, and m>=37 is automatically below
    # density 1/6.  The displayed finite lists are consequently all-order
    # classifications, not cutoff conjectures.
    require(5 * (16 + 6) < 7 * 16, "one-fifth tail")
    require(6 * (37 + 6) < 7 * 37, "one-sixth tail")

    half_five = tuple(order for order in range(2, 16) if 5 * half_block_count(order) >= order)
    require(half_five == (3, 5, 8, 9, 10, 15), half_five)
    half_six = tuple(order for order in range(2, 37) if 6 * half_block_count(order) >= order)
    require(
        half_six == (3, 5, 8, 9, 10, 11, 12, 15, 17, 22, 23, 24, 29, 36),
        half_six,
    )
    half_complements = tuple(order_three_complement_record(order, 1) for order in half_six)
    require(
        tuple(order for order, q, count, _, _ in half_complements if 6 * count > q) == (9,),
        half_complements,
    )
    require(
        tuple(order for order, q, count, _, _ in half_complements if 6 * count == q)
        == (8, 12),
        half_complements,
    )

    zero_five = tuple(order for order in range(2, 16) if 5 * zero_block_count(order) >= order)
    require(zero_five == (2, 3, 4, 5, 15), zero_five)
    zero_six = tuple(order for order in range(2, 37) if 6 * zero_block_count(order) >= order)
    require(zero_six == (2, 3, 4, 5, 6, 15, 16, 17, 18, 29, 30), zero_six)
    odd_single_three = tuple(order for order in zero_six if order % 2 and order % 9)
    require(odd_single_three == (3, 5, 15, 17, 29), odd_single_three)
    zero_complements = tuple(order_three_complement_record(order, 0) for order in odd_single_three)
    require(
        all(6 * count < q for _, q, count, _, _ in zero_complements),
        zero_complements,
    )

    return half_five, half_six, half_complements, zero_five, zero_six, zero_complements


def totient(value):
    answer = value
    for prime in prime_factors(value):
        answer = answer // prime * (prime - 1)
    return answer


def coprime_prefix_count(value, bound):
    return sum(1 for residue in range(1, bound + 1) if gcd(value, residue) == 1)


def fixed_zero_alpha(order):
    cutoff = (order - 1) // 14
    return Fraction(2 * coprime_prefix_count(order, cutoff), totient(order))


def fixed_zero_base_gate_audit():
    # THM-3408 proves for every base-free order that alpha<=1/5, with
    # equality exactly on these four orders.  The small calculations below
    # verify the prime-loss stratum that makes fivefold equality impossible.
    equality_orders = (22, 33, 44, 50)
    require(
        tuple((order, fixed_zero_alpha(order)) for order in equality_orders)
        == tuple((order, Fraction(1, 5)) for order in equality_orders),
        equality_orders,
    )
    loss_eleven = tuple(
        (order, order // gcd(order, 11), fixed_zero_alpha(order // gcd(order, 11)))
        for order in equality_orders[:3]
    )
    require(loss_eleven == ((22, 2, 0), (33, 3, 0), (44, 4, 0)), loss_eleven)
    loss_five = (50, 50 // gcd(50, 5), fixed_zero_alpha(50 // gcd(50, 5)))
    require(loss_five == (50, 10, 0), loss_five)
    return equality_orders, loss_eleven, loss_five


def partition_record(q, residues):
    masks = tuple(danger_mask(q, residue, 1) for residue in residues)
    full = (1 << q) - 1
    union = 0
    xor = 0
    for mask in masks:
        union |= mask
        xor ^= mask
    require(union == full and xor == full, (q, residues, masks, union, xor))
    require(sum(mask.bit_count() for mask in masks) == q, (q, residues, masks))
    require(gcd(2 * q, *residues) == 1, (q, residues))
    sheets = tuple(tuple(sheet for sheet in range(q) if mask & (1 << sheet)) for mask in masks)
    return q, residues, sheets


def augmented_witness_record(q, epsilon, residues):
    modulus, _, _, unique, _, full = augmented_bank(q, epsilon)
    by_residue = {residue: mask for mask, residue in unique}
    joined = 0
    for residue in residues:
        require(residue in by_residue, (q, epsilon, residue))
        joined |= by_residue[residue]
    require(joined == full, (q, epsilon, residues, joined, full))
    require(gcd(modulus, *residues) == 1, (q, epsilon, residues))
    return q, epsilon, residues


def construction_audit():
    q10 = partition_record(10, (1, 3, 4, 7, 9))
    q12 = partition_record(12, (1, 5, 7, 8, 11))
    require(
        q10[2] == ((0, 9), (3, 6), (2, 7), (1, 8), (4, 5)),
        q10,
    )
    require(
        q12[2] == ((0, 11), (2, 9), (3, 8), (1, 4, 7, 10), (5, 6)),
        q12,
    )

    breaker_lifts = []
    for q in range(16, SCAN_LIMIT + 1, 8):
        scale = q // 8
        breaker_lifts.append(augmented_witness_record(q, 1, (1, scale, 3 * scale, 5 * scale, 7 * scale)))
    for q in range(18, SCAN_LIMIT + 1, 9):
        scale = q // 9
        breaker_lifts.append(augmented_witness_record(q, 1, (1, scale, 5 * scale, 6 * scale, 7 * scale)))
    zero_exceptions = (
        augmented_witness_record(16, 0, (1, 3, 5, 7, 8)),
        augmented_witness_record(18, 0, (1, 5, 6, 7, 9)),
    )
    return q10, q12, zero_exceptions, tuple(breaker_lifts)


def fibonacci_period(modulus):
    values = []
    previous, current = 0, 1
    while True:
        values.append(previous)
        previous, current = current, (previous + current) % modulus
        if (previous, current) == (0, 1):
            return tuple(values)


def in_rank_five_support(value):
    return (value % 10 == 0 or value % 12 == 0) and value % 8 and value % 9


def arithmetic_transport_audit():
    support_residues = tuple(residue for residue in range(360) if in_rank_five_support(residue))
    require(len(support_residues) == 32, support_residues)
    require(Fraction(len(support_residues), 360) == Fraction(4, 45), support_residues)

    primitive_states = 0
    accepting_states = 0
    for numerator in range(360):
        for denominator in range(360):
            if gcd(numerator, denominator, 360) != 1:
                continue
            primitive_states += 1
            if in_rank_five_support(denominator):
                accepting_states += 1
    require((primitive_states, accepting_states) == (82944, 4128), (primitive_states, accepting_states))
    require(Fraction(accepting_states, primitive_states) == Fraction(43, 864), accepting_states)

    fib = fibonacci_period(360)
    require(len(fib) == 120, len(fib))
    fib_rank_four = tuple(index for index, value in enumerate(fib) if value % 8 == 0 or value % 9 == 0)
    fib_rank_five = tuple(index for index, value in enumerate(fib) if in_rank_five_support(value))
    require(fib_rank_four == tuple(range(0, 120, 6)), fib_rank_four)
    require(fib_rank_five == (15, 45, 75, 105), fib_rank_five)

    berggren_rank_four = tuple(
        index
        for index in range(9)
        if (4 * index * index + 12 * index + 11) % 9 == 0
    )
    berggren_rank_five = tuple(
        index
        for index in range(360)
        if in_rank_five_support(4 * index * index + 12 * index + 11)
    )
    require(berggren_rank_four == (1, 5), berggren_rank_four)
    require(berggren_rank_five == (), berggren_rank_five)
    return (
        support_residues,
        (primitive_states, accepting_states, Fraction(accepting_states, primitive_states)),
        fib_rank_four,
        fib_rank_five,
        berggren_rank_four,
        berggren_rank_five,
    )


def main():
    dependencies = tuple((name, lf_hash(path)) for name, path, _ in PINNED)
    for name, path, expected in PINNED:
        require(lf_hash(path) == expected, ("dependency changed", name, path, lf_hash(path), expected))

    density_audit = density_and_complement_audit()
    fixed_zero_audit = fixed_zero_base_gate_audit()
    constructions = construction_audit()
    arithmetic = arithmetic_transport_audit()

    records = tuple(record(q, epsilon) for q in range(2, SCAN_LIMIT + 1) for epsilon in (0, 1))
    primitive = tuple(row for row in records if row[5] is not None)
    primitive_rank4 = tuple((row[0], row[1], row[6], row[7]) for row in primitive if row[5] <= 4)
    primitive_rank5 = tuple((row[0], row[1], row[6], row[7]) for row in primitive if row[5] == 5)
    require(
        primitive_rank4
        == (
            (8, 1, (1, 3, 5, 7), (8, 8, 8, 8)),
            (9, 1, (1, 5, 6, 7), (3, 9, 9, 9)),
        ),
        primitive_rank4,
    )
    primitive_rank5_keys = tuple((q, epsilon) for q, epsilon, _, _ in primitive_rank5)
    expected_zero_rank5 = ((16, 0), (18, 0))
    expected_half_rank5 = tuple(
        (q, 1)
        for q in range(2, SCAN_LIMIT + 1)
        if q in (10, 12) or (q > 8 and q % 8 == 0) or (q > 9 and q % 9 == 0)
    )
    require(
        tuple(key for key in primitive_rank5_keys if key[1] == 0) == expected_zero_rank5,
        primitive_rank5_keys,
    )
    require(
        tuple(key for key in primitive_rank5_keys if key[1] == 1) == expected_half_rank5,
        primitive_rank5_keys,
    )

    by_key = {(row[0], row[1]): row for row in records}
    global_rows = []
    for q in range(2, SCAN_LIMIT + 1):
        candidates = []
        for quotient in divisors(q):
            for epsilon in (0, 1):
                rank = by_key[(quotient, epsilon)][5]
                if rank is not None:
                    candidates.append((rank, quotient, epsilon))
        best = min((candidate[0] for candidate in candidates), default=None)
        minimizers = tuple(candidate[1:] for candidate in candidates if candidate[0] == best)
        global_rows.append((q, best, minimizers))
    global_rows = tuple(global_rows)
    global_rank5 = tuple(row for row in global_rows if row[1] == 5)
    predicted_global_rank5 = tuple(
        q
        for q in range(2, SCAN_LIMIT + 1)
        if (q % 10 == 0 or q % 12 == 0) and q % 8 and q % 9
    )
    require(tuple(row[0] for row in global_rank5) == predicted_global_rank5, global_rank5)
    global_rank4 = tuple(row[0] for row in global_rows if row[1] == 4)
    require(
        global_rank4
        == tuple(q for q in range(2, SCAN_LIMIT + 1) if q % 8 == 0 or q % 9 == 0),
        global_rank4,
    )

    event_digest = sha256(repr(records).encode("ascii")).hexdigest()
    require(event_digest == EXPECTED_EVENT_DIGEST, event_digest)
    search_totals = (
        sum(row[2] for row in records),
        sum(row[3] for row in records),
        sum(row[4] for row in records),
        sum(row[8] for row in records),
        sum(row[9] for row in records),
    )
    require(search_totals == (371760, 184338, 183424, 381717, 386458), search_totals)

    semantic = sha256(
        repr(
            (
                density_audit,
                fixed_zero_audit,
                constructions[:3],
                len(constructions[3]),
                arithmetic,
                primitive_rank4,
                primitive_rank5_keys,
                global_rank4,
                global_rank5,
                search_totals,
                event_digest,
            )
        ).encode("ascii")
    ).hexdigest()
    require(semantic == EXPECTED_SEMANTIC_DIGEST, semantic)

    half_complements = density_audit[2]
    zero_complements = density_audit[5]
    support_residues, farey_states, fib_rank4, fib_rank5, berggren_rank4, berggren_rank5 = arithmetic

    print("LRC ZERO-MODE-COCHAIN RANK-FIVE SUPPORT PROBE")
    print(f"source_sha256_lf={lf_hash(Path(__file__))}")
    print(f"dependency_sha256_lf={dependencies}")
    print("status=PROVED-ELEMENTARY global_rank5_iff_(10|q_or_12|q)_and_8_not|q_and_9_not|q;FINITE-EXACT independent_rare-coordinate_Q2_Q500_primitive_audit;no_LRC14_decrement")
    print("all_q_reduction=THM3405_divisor_minimum_plus_rank4_iff_8|q_or_9|q;every_rank_at_most_5_primitive_quotient_is_divisible_by_8_or_9_or_10_or_12")
    print(f"half_density_gate=(orders_with_5h>=m,orders_with_6h>=m)={(density_audit[0], density_audit[1])}")
    print(f"half_order3_complement_gate=(greater_than_1/6,equal_1/6)={(tuple(row[0] for row in half_complements if 6 * row[2] > row[1]), tuple(row[0] for row in half_complements if 6 * row[2] == row[1]))}")
    print(f"zero_basefree_gate=(alpha_equal_1/5,prime11_losses,prime5_loss)={fixed_zero_audit}")
    print(f"zero_odd_15_order3_complement_controls={zero_complements}")
    print(f"rank5_atoms_Q10_Q12_disjoint_OR_equals_XOR={(constructions[0], constructions[1])}")
    print(f"primitive_breaker_constructions_through_Q500=(zero_exceptions,count_half_lifts)={(constructions[2], len(constructions[3]))}")
    print(f"primitive_rank4=(Q,epsilon,residues,orders)={primitive_rank4}")
    print(f"primitive_rank5_Q2_Q500=(zero_keys,half_count,all_keys_sha256)={(expected_zero_rank5, len(expected_half_rank5), sha256(repr(primitive_rank5_keys).encode('ascii')).hexdigest())}")
    print(f"global_rank5_support_mod360=(residue_count,density,harmonic_coefficient,residues)={(len(support_residues), Fraction(len(support_residues), 360), Fraction(4,45), support_residues)}")
    print(f"global_rank5_Q2_Q500=(count,first,last,rows_sha256)={(len(global_rank5), global_rank5[0], global_rank5[-1], sha256(repr(global_rank5).encode('ascii')).hexdigest())}")
    print(f"Farey_Stern_Brocot_mod360=(primitive_states,rank5_accepting_states,density)={farey_states}")
    print(f"Fibonacci_period360_rank_support=(rank4_indices,rank5_indices)={(fib_rank4, fib_rank5)}")
    print(f"Berggren_U_spine_support=(rank4_n_mod9,rank5_n_mod360)={(berggren_rank4, berggren_rank5)}")
    print(f"search_totals=(raw,unique,maximal,nodes,branches)={search_totals}")
    print(f"event_sha256={event_digest}")
    print(f"semantic_sha256={semantic}")
    print("scope=global_minimum_over_unrestricted_distinct_positive_transverse_zero_complete_cochains;primitive_Q500_pattern_beyond_support_necessity_is_finite_only;no_tournament_equivalence")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
