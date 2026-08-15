#!/usr/bin/env python3
"""Exact companion for the odd target-free p=7 half-twist obstruction.

The p=7 activity/partition clause needed here is proved directly from the
proved rank-six support theorem THM-3416, independently of the broader proved
THM-3429 classification.  Suppose at most seven distinct blocks cover at joint
quotient period Q=7M.  Joint period supplies a 7-active owner.  If an
inactive owner exists, at most six inactive blocks descend exactly to M.  As M
is target-free, THM-3416 says that they miss a base sheet.  At most six active
blocks remain, and each hits at most one point of the missed seven-fibre, an
impossibility.  Hence all owners are active; fibre capacity then forces exactly
seven owners and one hit by each owner on every fibre.

Write Q=7M and d=gcd(Q,r).  Activity means 7 does not divide r, so the
reduced order m=Q/d is 7k.  Because Q is odd, d,m,k are odd.  Direct counting
of the strict half-twist block gives

    |B_(Q,r)| = M       if r is even,
                 M-d   if r is odd.

One hit in every fibre forces size M, hence every owner residue is even.
All even half-twist blocks contain the common reflection-fixed sheet
(Q-1)/2.  Seven such blocks cannot be an exact partition.  Independently,
their augmented gcd with 2Q is even.  Thus the p=7 lane is empty even in the
literal joint-period universe; the augmented primitive contradiction is a
second typing obstruction.

The all-Q argument above is algebraic.  The finite computation is a hostile
referee: it checks the exact size/fibre identity including repeated factors
of seven, replays positive and near-positive boundaries, and uses a separate
exact union solver on small controls.  Only the Python standard library is
used, and every truth gate survives ``python -O``.
"""

from __future__ import annotations

import ast
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM3416 = ROOT / "01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md"
THM3416_SHA256_LF = "42a9309145de51d1bb6fca0b7c1945302ff37a63a3183e1dfed838c07118e8bf"
TARGET_FREE_BASES = (8, 9, 10, 11, 12, 15, 23, 25)
EXPECTED_SEMANTIC_SHA256 = "8599e1b8ced3588d68d4cdd176de5932b72cb1ac4bc434b2dab14114f3914d89"


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


def target_free(q):
    return not any(q % base == 0 for base in TARGET_FREE_BASES)


def half_mask(q, residue):
    modulus = 2 * q
    mask = 0
    for sheet in range(q):
        phase = residue * (2 * sheet + 1) % modulus
        if 14 * min(phase, modulus - phase) < modulus:
            mask |= 1 << sheet
    return mask


def quotient_order(q, residue):
    return q // gcd(q, residue)


def joined_period(q, residues):
    return lcm(*(quotient_order(q, residue) for residue in residues))


def multiplicities(q, residues):
    masks = tuple(half_mask(q, residue) for residue in residues)
    values = tuple(
        sum((mask >> sheet) & 1 for mask in masks)
        for sheet in range(q)
    )
    return masks, values


def cover_certificate(q, residues):
    require(len(residues) == len(set(residues)), (q, residues))
    masks, values = multiplicities(q, residues)
    union_size = sum(value > 0 for value in values)
    order_histogram = tuple(
        sorted(Counter(quotient_order(q, residue) for residue in residues).items())
    )
    multiplicity_histogram = tuple(sorted(Counter(values).items()))
    return (
        q,
        residues,
        tuple(mask.bit_count() for mask in masks),
        order_histogram,
        joined_period(q, residues),
        gcd(2 * q, *residues),
        sum(residue % 7 != 0 for residue in residues) if q % 7 == 0 else None,
        union_size,
        multiplicity_histogram,
        all(value == 1 for value in values),
    )


def p7_residue_record(q, residue):
    require(q % 2 == 1 and q % 7 == 0, (q, residue))
    require(0 < residue < 2 * q and residue % q and residue % 7, (q, residue))
    d = gcd(q, residue)
    order = q // d
    require(order % 7 == 0 and d % 7 != 0, (q, residue, d, order))
    k = order // 7
    require(k % 2 == 1, (q, residue, order, k))

    mask = half_mask(q, residue)
    base_size = q // 7
    expected_size = base_size if residue % 2 == 0 else base_size - d
    require(mask.bit_count() == expected_size, (q, residue, mask.bit_count(), expected_size))

    fibre_counts = tuple(
        sum((mask >> (base + step * base_size)) & 1 for step in range(7))
        for base in range(base_size)
    )
    require(max(fibre_counts) <= 1, (q, residue, fibre_counts))
    expected_zero_fibres = 0 if residue % 2 == 0 else d
    require(fibre_counts.count(0) == expected_zero_fibres, (q, residue, fibre_counts))
    require(fibre_counts.count(1) == expected_size, (q, residue, fibre_counts))

    fixed_sheet = (q - 1) // 2
    fixed_hit = bool(mask & (1 << fixed_sheet))
    require(fixed_hit == (residue % 2 == 0), (q, residue, fixed_sheet, fixed_hit))
    return (
        "even" if residue % 2 == 0 else "odd",
        d,
        order,
        k,
        expected_size,
        expected_zero_fibres,
        int(fixed_hit),
    )


def build_bank(q, augmented):
    modulus = 2 * q
    primes = prime_factors(modulus)
    grouped = {}
    raw = 0
    for residue in range(1, modulus):
        if residue % q == 0:
            continue
        sheet_mask = half_mask(q, residue)
        if not sheet_mask:
            continue
        raw += 1
        mask = sheet_mask
        if augmented:
            for offset, prime in enumerate(primes):
                if residue % prime:
                    mask |= 1 << (q + offset)
        grouped.setdefault(mask, residue)

    items = tuple(sorted(((mask, residue) for mask, residue in grouped.items()), key=lambda row: row[1]))
    maximal = tuple(
        item
        for item in items
        if not any(
            item[0] != other[0] and item[0] | other[0] == other[0]
            for other in items
        )
    )
    full = (1 << (q + (len(primes) if augmented else 0))) - 1
    return raw, items, maximal, full


def exact_union_cover(q, augmented, cap):
    raw, items, maximal, full = build_bank(q, augmented)
    masks = tuple(mask for mask, _ in maximal)
    residues = tuple(residue for _, residue in maximal)
    width = full.bit_length()
    by_bit = tuple(
        tuple(index for index, mask in enumerate(masks) if mask & (1 << bit))
        for bit in range(width)
    )
    require(all(by_bit), (q, augmented, tuple(i for i, row in enumerate(by_bit) if not row)))
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
        if not gains or sum(gains[:slots]) < missing.bit_count():
            return None
        missing_bits = tuple(bit for bit in range(width) if missing & (1 << bit))
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
    witness = None if chosen is None else tuple(sorted(residues[index] for index in chosen))
    if witness is not None:
        require(len(witness) == len(set(witness)) <= cap, witness)
        masks_by_residue = {residue: mask for mask, residue in items}
        joined = 0
        for residue in witness:
            joined |= masks_by_residue[residue]
        require(joined == full, (q, augmented, cap, witness))
        if augmented:
            require(gcd(2 * q, *witness) == 1, witness)
    return (
        q,
        "augmented" if augmented else "literal",
        cap,
        witness,
        nodes,
        branches,
        solve.cache_info().hits,
        raw,
        len(items),
        len(maximal),
    )


def identity_census():
    moduli = tuple(
        q
        for q in range(21, 316, 14)
        if target_free(q)
    )
    expected_moduli = (21, 35, 49, 91, 119, 133, 147, 203, 217, 245, 259, 273, 287, 301)
    require(moduli == expected_moduli, moduli)

    active_residues = 0
    phase_cells = 0
    fibre_cells = 0
    fixed_sheet_checks = 0
    records = []
    for q in moduli:
        profiles = []
        even = 0
        odd = 0
        for residue in range(1, 2 * q):
            if residue % q == 0 or residue % 7 == 0:
                continue
            profile = p7_residue_record(q, residue)
            profiles.append(profile)
            if residue % 2 == 0:
                even += 1
            else:
                odd += 1
            active_residues += 1
            phase_cells += q
            fibre_cells += q
            fixed_sheet_checks += 1
        records.append(
            (
                q,
                len(profiles),
                even,
                odd,
                tuple(sorted(Counter(profiles).items())),
            )
        )
    require(active_residues == 4080, active_residues)
    require(phase_cells == fibre_cells == 906360, (phase_cells, fibre_cells))
    require(fixed_sheet_checks == active_residues, fixed_sheet_checks)
    return moduli, active_residues, phase_cells, fibre_cells, fixed_sheet_checks, tuple(records)


def sharp_controls():
    q13 = cover_certificate(13, (1, 2, 3, 5, 7, 9, 11))
    q51 = cover_certificate(51, (1, 11, 12, 18, 23, 34, 35))
    q91 = cover_certificate(91, (7, 14, 21, 35, 49, 63, 77))
    q21_near = cover_certificate(21, (2, 4, 6, 8, 10, 16, 20))

    require(q13[4:] == (13, 1, None, 13, ((1, 13),), True), q13)
    require(q51[3:] == (((3, 1), (17, 2), (51, 4)), 51, 1, None, 51, ((1, 42), (2, 4), (3, 3), (4, 2)), False), q51)
    require(q91[3:] == (((13, 7),), 13, 7, 0, 91, ((1, 91),), True), q91)
    require(q21_near[3:] == (((7, 1), (21, 6)), 21, 2, 7, 15, ((0, 6), (1, 14), (7, 1)), False), q21_near)
    require(all(size == 3 for size in q21_near[2]), q21_near)

    searches = (
        exact_union_cover(13, True, 6),
        exact_union_cover(13, True, 7),
        exact_union_cover(21, False, 7),
        exact_union_cover(35, False, 7),
        exact_union_cover(49, False, 7),
        exact_union_cover(51, True, 6),
        exact_union_cover(51, True, 7),
        exact_union_cover(91, False, 7),
        exact_union_cover(91, True, 7),
    )
    expected_witnesses = (
        None,
        (1, 2, 3, 5, 7, 9, 11),
        None,
        None,
        None,
        None,
        (1, 11, 12, 18, 23, 34, 35),
        (7, 14, 21, 35, 49, 63, 77),
        None,
    )
    require(tuple(row[3] for row in searches) == expected_witnesses, searches)
    return q13, q51, q91, q21_near, searches


def main():
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert found")
    require(lf_sha256(THM3416) == THM3416_SHA256_LF, lf_sha256(THM3416))

    census = identity_census()
    controls = sharp_controls()
    semantic_surface = (
        THM3416_SHA256_LF,
        census,
        controls,
        "THM3416_DESCENT_DERIVES_P7_PARTITION_THEN_PARITY_FORCES_FIXED_SHEET_COLLISION_AND_EVEN_AUGMENTED_GCD",
    )
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256, semantic_digest)

    print("Odd target-free p=7 half-twist parity obstruction")
    print("status=PROVED_FROM_THM3416_PLUS_ELEMENTARY_P7_FIBRE_AND_PARITY_LEMMAS all_Q_literal_joint_period_p7_lane_NEGATIVE;primitive_augmented_NEGATIVE")
    print(f"dependency_THM3416_sha256_lf={THM3416_SHA256_LF}")
    print("typed_literal_universe=Q_odd_composite_target_free;7_divides_Q;at_most_seven_distinct_transverse_half_twist_residues_mod_2Q(Q_not_divide_r);literal_sheet_cover;joint_quotient_period_Q")
    print("typed_augmented_universe=literal_universe_plus_gcd(2Q,R)=1;literal_contradiction_precedes_augmented_gate")
    print("distinct_owner_semantics=owners_are_distinct_residues_mod_2Q;equal_sheet_masks_if_any_still_count_as_separate_owners;descent_uses_at_most_six_inactive_owners")
    print("direct_descent=joint_period_supplies_a_7_active_owner;if_any_owner_is_7_inactive_then_at_most_six_inactive_blocks_pull_back_from_target_free_M=Q/7_and_THM3416_makes_them_miss_a_base_sheet;at_most_six_active_blocks_each_hit_at_most_one_point_of_that_7_fibre")
    print("derived_partition=no_inactive_owner;cover_and_cap_seven_force_exactly_seven_active_owners_and_one_hit_per_owner_per_7_fibre;global_OR_equals_XOR")
    print("all_Q_block_identity=Q=7M,d=gcd(Q,r),7_not_divide_r:even_r_gives_size_M_and_zero_missed_fibres;odd_r_gives_size_M-d_and_d_missed_fibres")
    print("parity_conclusion=one_hit_on_every_fibre_forces_all_seven_residues_even")
    print("literal_contradiction=all_even_blocks_contain_the_common_sheet_(Q-1)/2;total_mass_Q_has_overlap_defect_at_least_6")
    print("augmented_contradiction=all_residues_even_implies_gcd(2Q,R)>=2")
    print(f"finite_identity_census=(moduli,active_residues,phase_cells,fibre_cells,fixed_sheet_checks)={census[:5]}")
    print(f"finite_profile_records={census[5]}")
    print(f"sharp_certificates=(Q13_partition,Q51_mixed_positive,Q91_pullback_partition,Q21_joint_period_near_miss)={controls[:4]}")
    print(f"exact_search_controls={controls[4]}")
    print("hostile_meaning=Q91_shows_partition_without_7_activity_or_joint_period;Q21_shows_7_activity_fibre_saturation_and_joint_period_without_cover_or_augmented_primitivity")
    print(f"semantic_sha256={semantic_digest}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=standard_library_only;no_elapsed_fields;truth_gates_survive_python_O")
    print("commands=python -B 04-computation/odd_target_free_p7_half_twist_parity_obstruction.py;python -B -O 04-computation/odd_target_free_p7_half_twist_parity_obstruction.py")
    print("scope=standalone_closure_of_only_the_p7_lane;proof_independently_based_on_THM3416;remaining_p3_p5_p17_p29_lanes_fixed_zero_arbitrary_time_and_LRC14_remain_open")


if __name__ == "__main__":
    main()
