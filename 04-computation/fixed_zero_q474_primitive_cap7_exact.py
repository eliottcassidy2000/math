#!/usr/bin/env python3
"""Exact fixed-zero Q=474 primitive-augmented cap-seven decision.

The target result is stronger than primitive augmented infeasibility: no
literal union of at most seven distinct fixed-zero owner blocks covers all
474 sheets.  The proof is a finite reflection-orbit fresh-gain certificate.

The only exceptional least order is two.  Its block is exactly the even
sheets, and restriction of every remaining block to the odd sheets is a
half-twist block on 237 sheets.  A second reflection certificate excludes a
six-block half-twist cover there.  Order six is literally dominated by the
order-two block.  For every other possible least quotient order, unit
multiplication normalizes one owner and the six largest fresh orbit gains do
not fill the remaining reflection orbits.

The primitive gate is encoded by one breaker bit for every prime dividing
the owner modulus.  All truth gates use ``require`` and therefore survive
``python -O``.  This file uses only the Python standard library and imports
no repository computation module.
"""

from __future__ import annotations

import ast
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from math import gcd
from pathlib import Path


TARGET_Q = 474
TARGET_CAP = 7
EXPECTED_SEMANTIC_SHA256 = "c1f1ace0f12e3e74396997c62e0ab0aa18a7d4bf1ad24f60ec3dbfdd74368ecd"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def prime_factors(value):
    factors = []
    trial = 2
    remaining = value
    while trial * trial <= remaining:
        if remaining % trial == 0:
            factors.append(trial)
            while remaining % trial == 0:
                remaining //= trial
        trial += 1
    if remaining > 1:
        factors.append(remaining)
    return tuple(factors)


def fixed_sheet_mask(q, residue):
    mask = 0
    for sheet in range(q):
        phase = residue * sheet % q
        if 14 * min(phase, q - phase) < q:
            mask |= 1 << sheet
    return mask


def half_sheet_mask(q, residue):
    modulus = 2 * q
    mask = 0
    for sheet in range(q):
        phase = residue * (2 * sheet + 1) % modulus
        if 14 * min(phase, modulus - phase) < modulus:
            mask |= 1 << sheet
    return mask


def reflection_orbits(q, epsilon):
    require(epsilon in (0, 1), (q, epsilon))
    unseen = set(range(q))
    answer = []
    while unseen:
        sheet = min(unseen)
        mate = (-sheet if epsilon == 0 else -1 - sheet) % q
        orbit = tuple(sorted({sheet, mate}))
        answer.append(orbit)
        unseen.difference_update(orbit)
    return tuple(answer)


def compress_reflection(mask, orbits):
    compressed = 0
    for index, orbit in enumerate(orbits):
        bits = tuple((mask >> sheet) & 1 for sheet in orbit)
        require(len(set(bits)) == 1, (orbit, bits))
        if bits[0]:
            compressed |= 1 << index
    return compressed


def build_bank(q, epsilon, augmented):
    require(epsilon in (0, 1), (q, epsilon, augmented))
    modulus = q if epsilon == 0 else 2 * q
    primes = prime_factors(modulus)
    grouped = {}
    raw = 0
    for residue in range(1, modulus):
        if residue % q == 0:
            continue
        sheet_mask = (
            fixed_sheet_mask(q, residue)
            if epsilon == 0
            else half_sheet_mask(q, residue)
        )
        if not sheet_mask:
            continue
        raw += 1
        mask = sheet_mask
        if augmented:
            for offset, prime in enumerate(primes):
                if residue % prime:
                    mask |= 1 << (q + offset)
        grouped.setdefault(mask, []).append(residue)

    items = tuple(
        (mask, min(aliases), tuple(aliases))
        for mask, aliases in sorted(grouped.items(), key=lambda row: min(row[1]))
    )
    maximal = tuple(
        item
        for item in items
        if not any(
            item[0] != other[0] and item[0] | other[0] == other[0]
            for other in items
        )
    )
    width = q + (len(primes) if augmented else 0)
    full = (1 << width) - 1
    return {
        "q": q,
        "epsilon": epsilon,
        "augmented": augmented,
        "modulus": modulus,
        "primes": primes,
        "raw": raw,
        "items": items,
        "maximal": maximal,
        "full": full,
    }


def quotient_order(q, residue):
    return q // gcd(q, residue)


def additive_order(modulus, residue):
    return modulus // gcd(modulus, residue)


def unit_normalizer(modulus, residue, canonical):
    for unit in range(1, modulus):
        if gcd(unit, modulus) == 1 and unit * residue % modulus == canonical:
            return unit
    return None


def histogram(values):
    return tuple(sorted(Counter(values).items()))


def exact_union_cover(full, items, cap):
    """Exact rare-coordinate union search used only for hostile controls."""

    masks = tuple(item[0] for item in items)
    residues = tuple(item[1] for item in items)
    width = full.bit_length()
    by_bit = tuple(
        tuple(index for index, mask in enumerate(masks) if mask & (1 << bit))
        for bit in range(width)
    )
    require(all(by_bit), (width, tuple(i for i, row in enumerate(by_bit) if not row)))

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
            (
                (mask & missing).bit_count()
                for mask in masks
                if mask | state != state
            ),
            reverse=True,
        )
        if len(gains) < 1 or sum(gains[:slots]) < missing.bit_count():
            return None

        missing_bits = tuple(
            bit for bit in range(width) if missing & (1 << bit)
        )
        pivot = min(
            missing_bits,
            key=lambda bit: (
                sum(
                    1
                    for index in by_bit[bit]
                    if masks[index] | state != state
                ),
                bit,
            ),
        )
        candidates = sorted(
            (
                index
                for index in by_bit[pivot]
                if masks[index] | state != state
            ),
            key=lambda index: (
                -(masks[index] & missing).bit_count(),
                residues[index],
            ),
        )
        for index in candidates:
            branches += 1
            suffix = solve(state | masks[index], slots - 1)
            if suffix is not None:
                return (index,) + suffix
        return None

    chosen = solve(0, cap)
    witness = () if chosen is None else tuple(sorted(residues[i] for i in chosen))
    if witness:
        joined = 0
        for index in chosen:
            joined |= masks[index]
        require(joined == full, (cap, witness))
        require(len(witness) == len(set(witness)) <= cap, witness)
    return witness, nodes, branches, solve.cache_info().hits


def bank_summary(bank, order_function):
    return (
        bank["raw"],
        len(bank["items"]),
        len(bank["maximal"]),
        tuple(
            sorted(
                Counter(order_function(item[1]) for item in bank["maximal"]).items()
            )
        ),
    )


def target_fixed_certificate():
    q = TARGET_Q
    literal = build_bank(q, 0, False)
    augmented = build_bank(q, 0, True)
    require(literal["primes"] == (2, 3, 79), literal["primes"])
    require(augmented["primes"] == (2, 3, 79), augmented["primes"])

    literal_summary = bank_summary(literal, lambda r: quotient_order(q, r))
    augmented_summary = bank_summary(augmented, lambda r: quotient_order(q, r))
    require(
        literal_summary
        == (473, 237, 236, ((2, 1), (3, 1), (79, 39), (158, 39), (237, 78), (474, 78))),
        literal_summary,
    )
    require(
        augmented_summary
        == (473, 237, 237, ((2, 1), (3, 1), (6, 1), (79, 39), (158, 39), (237, 78), (474, 78))),
        augmented_summary,
    )

    literal_order_sizes = histogram(
        (quotient_order(q, item[1]), item[0].bit_count())
        for item in literal["items"]
    )
    require(
        literal_order_sizes
        == (
            ((2, 237), 1),
            ((3, 158), 1),
            ((6, 79), 1),
            ((79, 66), 39),
            ((158, 69), 39),
            ((237, 66), 78),
            ((474, 67), 78),
        ),
        literal_order_sizes,
    )

    by_order = {
        quotient_order(q, item[1]): item
        for item in literal["items"]
        if quotient_order(q, item[1]) in (2, 3, 6)
    }
    require(set(by_order) == {2, 3, 6}, by_order)
    order_two = by_order[2][0]
    order_six = by_order[6][0]
    require(order_six | order_two == order_two and order_six != order_two, "order-6 dominance")

    fixed_orbits = reflection_orbits(q, 0)
    require(len(fixed_orbits) == 238, len(fixed_orbits))
    compressed = tuple(
        (
            compress_reflection(item[0], fixed_orbits),
            item[1],
            quotient_order(q, item[1]),
        )
        for item in literal["maximal"]
    )
    full_orbits = (1 << len(fixed_orbits)) - 1

    normalizer_checks = 0
    for residue in range(1, q):
        order = quotient_order(q, residue)
        canonical = q // order
        require(unit_normalizer(q, residue, canonical) is not None, (residue, order))
        normalizer_checks += 1
    require(normalizer_checks == 473, normalizer_checks)

    rows = []
    gain_histograms = []
    for least_order in (3, 79, 158, 237, 474):
        canonical = q // least_order
        first = compress_reflection(fixed_sheet_mask(q, canonical), fixed_orbits)
        require(quotient_order(q, canonical) == least_order, (least_order, canonical))
        candidates = tuple(
            mask
            for mask, _, order in compressed
            if order >= least_order and mask != first
        )
        gains = tuple(((mask & ~first) & full_orbits).bit_count() for mask in candidates)
        top = tuple(sorted(gains, reverse=True)[: TARGET_CAP - 1])
        missing = len(fixed_orbits) - first.bit_count()
        total = sum(top)
        require(len(top) == TARGET_CAP - 1 and total < missing, (least_order, missing, top))
        rows.append(
            (
                least_order,
                canonical,
                first.bit_count(),
                missing,
                top,
                total,
                missing - total,
            )
        )
        gain_histograms.append((least_order, histogram(gains)))

    expected_rows = (
        (3, 158, 80, 158, (23, 23, 23, 23, 23, 23), 138, 20),
        (79, 6, 34, 204, (33, 33, 31, 31, 30, 30), 188, 16),
        (158, 3, 35, 203, (33, 33, 31, 31, 31, 31), 190, 13),
        (237, 2, 34, 204, (31, 30, 30, 30, 30, 30), 181, 23),
        (474, 1, 34, 204, (31, 31, 30, 30, 30, 30), 182, 22),
    )
    require(tuple(rows) == expected_rows, rows)

    order_two_even_cells = 0
    for sheet in range(q):
        observed = bool(order_two & (1 << sheet))
        require(observed == (sheet % 2 == 0), (sheet, observed))
        order_two_even_cells += 1

    restriction_cells = 0
    for residue in range(1, q):
        if residue == q // 2:
            continue
        fixed = fixed_sheet_mask(q, residue)
        half = half_sheet_mask(q // 2, residue)
        for sheet in range(q // 2):
            fixed_value = bool(fixed & (1 << (2 * sheet + 1)))
            half_value = bool(half & (1 << sheet))
            require(fixed_value == half_value, (residue, sheet))
            restriction_cells += 1

    return (
        literal_summary,
        augmented_summary,
        literal_order_sizes,
        len(fixed_orbits),
        tuple(rows),
        tuple(gain_histograms),
        normalizer_checks,
        order_two_even_cells,
        restriction_cells,
    )


def half_237_certificate():
    q = 237
    modulus = 2 * q
    bank = build_bank(q, 1, False)
    require(bank["raw"] == 470, bank["raw"])
    require(len(bank["items"]) == len(bank["maximal"]) == 235, len(bank["items"]))

    order_histogram = histogram(additive_order(modulus, item[1]) for item in bank["items"])
    require(
        order_histogram == ((3, 1), (79, 39), (158, 39), (237, 78), (474, 78)),
        order_histogram,
    )

    orbits = reflection_orbits(q, 1)
    require(len(orbits) == 119, len(orbits))
    compressed = tuple(
        (
            compress_reflection(item[0], orbits),
            item[1],
            additive_order(modulus, item[1]),
        )
        for item in bank["maximal"]
    )
    full_orbits = (1 << len(orbits)) - 1

    normalizer_checks = 0
    for residue in range(1, modulus):
        if residue % q == 0 or not half_sheet_mask(q, residue):
            continue
        order = additive_order(modulus, residue)
        canonical = modulus // order
        require(unit_normalizer(modulus, residue, canonical) is not None, (residue, order))
        normalizer_checks += 1
    require(normalizer_checks == 470, normalizer_checks)

    rows = []
    gain_histograms = []
    for least_order in (3, 79, 158, 237, 474):
        canonical = modulus // least_order
        first = compress_reflection(half_sheet_mask(q, canonical), orbits)
        candidates = tuple(
            mask
            for mask, _, order in compressed
            if order >= least_order and mask != first
        )
        gains = tuple(((mask & ~first) & full_orbits).bit_count() for mask in candidates)
        top = tuple(sorted(gains, reverse=True)[: TARGET_CAP - 2])
        missing = len(orbits) - first.bit_count()
        total = sum(top)
        require(len(top) == TARGET_CAP - 2 and total < missing, (least_order, missing, top))
        rows.append(
            (
                least_order,
                canonical,
                first.bit_count(),
                missing,
                top,
                total,
                missing - total,
            )
        )
        gain_histograms.append((least_order, histogram(gains)))

    expected_rows = (
        (3, 158, 40, 79, (12, 12, 12, 12, 12), 60, 19),
        (79, 6, 17, 102, (18, 18, 18, 18, 18), 90, 12),
        (158, 3, 18, 101, (18, 18, 18, 18, 18), 90, 11),
        (237, 2, 17, 102, (17, 17, 17, 17, 17), 85, 17),
        (474, 1, 17, 102, (17, 17, 17, 17, 17), 85, 17),
    )
    require(tuple(rows) == expected_rows, rows)
    return (
        (bank["raw"], len(bank["items"]), len(bank["maximal"])),
        order_histogram,
        len(orbits),
        tuple(rows),
        tuple(gain_histograms),
        normalizer_checks,
    )


def control_record(q, augmented, cap):
    bank = build_bank(q, 0, augmented)
    witness, nodes, branches, hits = exact_union_cover(
        bank["full"], bank["maximal"], cap
    )
    if witness:
        joined = 0
        by_residue = {item[1]: item[0] for item in bank["items"]}
        for residue in witness:
            joined |= by_residue[residue]
        require(joined == bank["full"], (q, augmented, cap, witness))
        require(len(witness) == len(set(witness)), witness)
        if augmented:
            require(gcd(bank["modulus"], *witness) == 1, witness)
    return (
        q,
        "augmented" if augmented else "literal",
        cap,
        witness if witness else None,
        nodes,
        branches,
        hits,
        bank["raw"],
        len(bank["items"]),
        len(bank["maximal"]),
    )


def hostile_controls():
    q29_cap6 = control_record(29, True, 6)
    q29_cap7 = control_record(29, True, 7)
    q21_literal = control_record(21, False, 7)
    q58_lit6 = control_record(58, False, 6)
    q58_lit7 = control_record(58, False, 7)
    q58_aug7 = control_record(58, True, 7)

    require(q29_cap6[3] is None, q29_cap6)
    require(q29_cap7[3] == (1, 4, 5, 6, 7, 9, 13), q29_cap7)
    require(q21_literal[3] is None, q21_literal)
    require(q58_lit6[3] is None, q58_lit6)
    require(q58_lit7[3] == (2, 8, 10, 12, 14, 18, 26), q58_lit7)
    require(gcd(58, *q58_lit7[3]) == 2, q58_lit7)
    require(q58_aug7[3] is None, q58_aug7)
    return (
        q29_cap6,
        q29_cap7,
        q21_literal,
        q58_lit6,
        q58_lit7,
        q58_aug7,
    )


def main():
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert found")

    fixed = target_fixed_certificate()
    half = half_237_certificate()
    controls = hostile_controls()
    semantic_surface = (fixed, half, controls, "Q474_LITERAL_UNSAT_IMPLIES_AUGMENTED_UNSAT")
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256, semantic_digest)

    print("Fixed-zero Q=474 primitive-augmented cap-seven exact decision")
    print("status=FINITE_EXACT literal_cap7_NEGATIVE;primitive_augmented_cap7_NEGATIVE")
    print("universe=fixed_zero_residues_1_to_Q-1;strict_1_over_14_sheet_predicate;literal_sheet_union;augmented_prime_breakers_2,3,79")
    print("augmented_gate=all_three_breaker_bits_covered_iff_gcd(474,r_1,...,r_k)=1;literal_NEGATIVE_is_the_stronger_result")
    print("distinct_owner_semantics=residues_are_distinct;identical_literal_or_augmented_masks_are_collapsed_because_duplicates_never_enlarge_an_at_most_cap_union")
    print("solver_normalization=reflection_orbits;least_order_owner;unit_multiplier_to_canonical_Q_over_m;order6_dominated_by_order2;order2_odd_coset_transplant_to_half_Q237")
    print(f"target_banks=(literal,augmented)={(fixed[0], fixed[1])}")
    print(f"target_literal_order_size_hist={fixed[2]}")
    print(f"fixed_reflection_orbits={fixed[3]}")
    print("fixed_fresh_gain_rows=(least_order,canonical_residue,first_orbits,missing_orbits,top_six_fresh,sum,deficit)=" + repr(fixed[4]))
    print(f"fixed_unit_normalizer_checks={fixed[6]}")
    print(f"order2_transplant_cells=(even_sheet_cells,odd_restriction_cells)={(fixed[7], fixed[8])}")
    print(f"half_Q237_bank_and_orders={(half[0], half[1], half[2])}")
    print("half_Q237_fresh_gain_rows=(least_additive_order,canonical_residue,first_orbits,missing_orbits,top_five_fresh,sum,deficit)=" + repr(half[3]))
    print(f"half_Q237_unit_normalizer_checks={half[5]}")
    print("structural_obstruction=every_non_order2_least_order_has_positive_reflection_fresh_gain_deficit;order2_would_force_the_independently_excluded_half_Q237_cap6_cover")
    print(f"hostile_controls={controls}")
    print("control_meaning=Q29_primitive_rank7_positive;Q21_literal_cap7_negative;Q58_literal_rank7_pullback_positive_but_augmented_cap7_negative")
    print(f"semantic_sha256={semantic_digest}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=standard_library_only;no_elapsed_fields;truth_gates_survive_python_O")
    print("commands=python -B 04-computation/fixed_zero_q474_primitive_cap7_exact.py;python -B -O 04-computation/fixed_zero_q474_primitive_cap7_exact.py")
    print("scope=Q474_fixed_zero_literal_and_primitive_augmented_only;no_all_q_antichain_completeness;no_LRC14_consequence")


if __name__ == "__main__":
    main()
