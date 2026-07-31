#!/usr/bin/env python3
"""Exact hostile probe for the THM-2818/THM-2819 wrap coincidence.

The script works on the common physical circle of length T.  It reconstructs
the three target-twelve right-cofiber rows from the constructors below
THM-2818, compares the THM-2819 odometer displacement 13*tau with the
exceptional half-step, and tests the induced translation on interval copies.

The comparison retains:

* the six native/pulled factor predicates;
* all thirteen source/target carrier-twist masks;
* all 169 endpoint-address intersection profiles;
* the complete U/V contributor ancestry chamber.

No assertion statement is used.
"""

from bisect import bisect_right
from collections import Counter
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "lrc14_right_cofiber_positive_copy_stratification_thm2818.py":
        "85edac9bb03f1fef198343268f4fb1226bec122d45ded79a049f8fa9a73882a8",
}
for name, expected in PINNED.items():
    payload = (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    require(
        sha256(payload).hexdigest() == expected,
        f"canonical input changed: {name}",
    )


import lrc14_right_cofiber_positive_copy_stratification_thm2818 as copies


P = 13
T = copies.T
H = copies.HALF_STEP
TAU = copies.SHIFT
W = P * TAU
R = copies.physical.R
FACTOR_NAMES = copies.FACTOR_NAMES


def split_blocks(pieces):
    blocks = []
    current = []
    for index, piece in enumerate(pieces):
        if index and piece[0] - pieces[index - 1][0] != H:
            blocks.append(tuple(current))
            current = []
        current.append(piece)
    if current:
        blocks.append(tuple(current))
    return tuple(blocks)


def cell_factors(full_module, e3, clocks, clock, target):
    universe = ((0, full_module.T),)
    return (
        e3,
        clocks[clock],
        full_module.subtract_comb(
            universe, full_module.W[1], 182, -13, 13
        ),
        full_module.subtract_comb(
            universe,
            full_module.W[2],
            182,
            -14 * target - 13,
            -14 * target + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.C2, 182, -13, 13
        ),
        full_module.subtract_comb(
            universe,
            full_module.C3,
            182,
            14 * target - 13,
            14 * target + 13,
        ),
    )


def index_supports(present_sets):
    return {
        address: (
            support,
            tuple(right for _left, right in support),
        )
        for address, support in present_sets.items()
    }


def index_one(support):
    return support, tuple(right for _left, right in support)


def local_intersection_profile(interval, indexed_support):
    """Intersection with one interval in coordinates relative to its left."""
    support, rights = indexed_support
    left, right = interval
    index = bisect_right(rights, left)
    rows = []
    while index < len(support) and support[index][0] < right:
        a = max(left, support[index][0])
        b = min(right, support[index][1])
        if a < b:
            rows.append((a - left, b - left))
        index += 1
    return tuple(rows)


def contains_fast(interval, indexed_support):
    return local_intersection_profile(interval, indexed_support) == (
        (0, interval[1] - interval[0]),
    )


def meets_fast(interval, indexed_support):
    return bool(local_intersection_profile(interval, indexed_support))


def endpoint_profile(interval, indexed_present):
    return tuple(
        local_intersection_profile(interval, indexed_present[address])
        for address in copies.endpoint_base.KEYS
    )


def target_endpoint_profile(interval, indexed_present):
    translated = (interval[0] + TAU, interval[1] + TAU)
    return endpoint_profile(translated, indexed_present)


def address_translation_matches(before, after):
    before_bank = dict(zip(copies.endpoint_base.KEYS, before))
    after_bank = dict(zip(copies.endpoint_base.KEYS, after))
    return tuple(
        delta
        for delta in copies.endpoint_base.KEYS
        if all(
            after_bank[address]
            == before_bank[
                (
                    (address[0] + delta[0]) % P,
                    (address[1] + delta[1]) % P,
                )
            ]
            for address in copies.endpoint_base.KEYS
        )
    )


def index_factor_faces(factors):
    return (
        tuple(index_one(factor) for factor in factors),
        tuple(index_one(copies.physical.overlap.shift_union(
            factor, -TAU
        )) for factor in factors),
        tuple(index_one(factor) for factor in factors),
        tuple(index_one(copies.physical.overlap.shift_union(
            factor, TAU
        )) for factor in factors),
    )


def factor_signature(interval, indexed_faces):
    target_interval = (interval[0] + TAU, interval[1] + TAU)
    return (
        tuple(contains_fast(interval, factor)
              for factor in indexed_faces[0]),
        tuple(contains_fast(interval, factor)
              for factor in indexed_faces[1]),
        tuple(contains_fast(target_interval, factor)
              for factor in indexed_faces[2]),
        tuple(contains_fast(target_interval, factor)
              for factor in indexed_faces[3]),
    )


def index_carrier_twists(source, target):
    unit = T // P
    source_twists = tuple(
        index_one(copies.support_of(
            copies.physical.overlap.shift_weighted(
                source, -twist * unit
            )
        ))
        for twist in range(P)
    )
    target_twists = tuple(
        index_one(copies.support_of(
            copies.physical.overlap.shift_weighted(
                target, twist * unit
            )
        ))
        for twist in range(P)
    )
    return source_twists, target_twists


def carrier_signature(interval, indexed_twists):
    target_interval = (interval[0] + TAU, interval[1] + TAU)
    return (
        tuple(meets_fast(interval, support)
              for support in indexed_twists[0]),
        tuple(meets_fast(target_interval, support)
              for support in indexed_twists[1]),
    )


def count_profile_types(profile):
    return (
        sum(not row for row in profile),
        sum(row == ((0, copies.LENGTH),) for row in profile),
        sum(bool(row) and row != ((0, copies.LENGTH),) for row in profile),
    )


def shifted_interval(piece, shift):
    return piece[0] + shift, piece[1] + shift, piece[2]


def main():
    require(
        (T, R, TAU, W)
        == (
            297_836_897_838_480,
            4_826_809,
            431_933_040,
            5_615_129_520,
        ),
        "THM-2818 physical constants stopped matching THM-2819",
    )
    require(
        H * 2 * P**5 == T
        and W == 14 * H
        and R * H * 2 == 13 * T
        and R * W == 91 * T,
        "scale identity changed",
    )

    (
        _module,
        _rails,
        _present,
        details,
        full_module,
        e3,
        clocks,
        q_pairs,
        _delayed,
        _source_weight,
        _target_weight,
        _rail_common,
    ) = copies.physical_setup()
    ancestry_arguments, u_events, v_events = copies.ancestry_setup()
    present_sets = copies.endpoint_base.present_cache()
    indexed_present = index_supports(present_sets)

    rows = []
    all_interior_factor_equal = True
    all_interior_carrier_equal = True
    all_interior_source_endpoint_equal = True
    all_interior_target_endpoint_equal = True
    all_interior_fixed_address_equal = True
    global_source_address_deltas = None
    global_target_address_deltas = None

    for clock in (1, 2, 3):
        (
            source,
            target_physical,
            _target_pullback,
            _common,
            pieces,
        ) = copies.physical_objects(
            details, full_module, e3, clocks, clock, 12
        )
        live, dead = copies.semantic_partition(pieces, q_pairs[clock])
        live_set = {piece[:2] for piece in live}
        dead_set = {piece[:2] for piece in dead}
        blocks = split_blocks(pieces)
        factors = cell_factors(full_module, e3, clocks, clock, 12)
        indexed_factor_faces = index_factor_faces(factors)
        indexed_carrier_twists = index_carrier_twists(
            source, target_physical
        )

        ancestry_report = copies.label_report(
            pieces,
            ancestry_arguments,
            u_events,
            v_events,
        )
        block_rows = []
        factor_pair_histogram = Counter()
        carrier_pair_histogram = Counter()
        source_endpoint_pair_histogram = Counter()
        target_endpoint_pair_histogram = Counter()
        exit_first_factor_failures = Counter()
        exit_first_carrier_failures = Counter()
        exit_source_endpoint_delta_histogram = Counter()
        exit_target_endpoint_delta_histogram = Counter()
        row_source_deltas = None
        row_target_deltas = None

        for block_index, block in enumerate(blocks):
            require(len(block) > 14, "block too short for wrap probe")
            for index in range(len(block) - 14):
                before = block[index]
                after = block[index + 14]
                require(
                    shifted_interval(before, W) == after,
                    "14-half-step translation missed its interior copy",
                )
                require(
                    (before[:2] in live_set) == (after[:2] in live_set)
                    and (before[:2] in dead_set)
                    == (after[:2] in dead_set),
                    "wrap changed live/dead parity",
                )

                factor_before = factor_signature(
                    before[:2], indexed_factor_faces
                )
                factor_after = factor_signature(
                    after[:2], indexed_factor_faces
                )
                factor_pair_histogram[
                    factor_before == factor_after
                ] += 1

                carrier_before = carrier_signature(
                    before[:2], indexed_carrier_twists
                )
                carrier_after = carrier_signature(
                    after[:2], indexed_carrier_twists
                )
                carrier_pair_histogram[
                    carrier_before == carrier_after
                ] += 1

                source_endpoint_before = endpoint_profile(
                    before[:2], indexed_present
                )
                source_endpoint_after = endpoint_profile(
                    after[:2], indexed_present
                )
                target_endpoint_before = target_endpoint_profile(
                    before[:2], indexed_present
                )
                target_endpoint_after = target_endpoint_profile(
                    after[:2], indexed_present
                )
                source_equal = (
                    source_endpoint_before == source_endpoint_after
                )
                target_equal = (
                    target_endpoint_before == target_endpoint_after
                )
                source_endpoint_pair_histogram[source_equal] += 1
                target_endpoint_pair_histogram[target_equal] += 1

                source_deltas = set(address_translation_matches(
                    source_endpoint_before, source_endpoint_after
                ))
                target_deltas = set(address_translation_matches(
                    target_endpoint_before, target_endpoint_after
                ))
                row_source_deltas = (
                    source_deltas
                    if row_source_deltas is None
                    else row_source_deltas & source_deltas
                )
                row_target_deltas = (
                    target_deltas
                    if row_target_deltas is None
                    else row_target_deltas & target_deltas
                )

            left_boundary = block[:14]
            right_boundary = block[-14:]
            left_live = sum(piece[:2] in live_set for piece in left_boundary)
            right_live = sum(piece[:2] in live_set for piece in right_boundary)
            require(
                (left_live, right_live) == (7, 7),
                "finite-chain live boundary stopped being balanced",
            )
            overlap = block[14:]
            require(
                len(overlap) == len(block) - 14
                and sum(piece[:2] in live_set for piece in overlap)
                == sum(piece[:2] in live_set for piece in block) - 7,
                "wrap overlap flux count changed",
            )

            for piece in right_boundary:
                exited = shifted_interval(piece, W)
                factor_before = factor_signature(
                    piece[:2], indexed_factor_faces
                )
                factor_after = factor_signature(
                    exited[:2], indexed_factor_faces
                )
                failed_factor_bits = tuple(
                    (face, FACTOR_NAMES[index])
                    for face in range(4)
                    for index in range(6)
                    if factor_before[face][index]
                    != factor_after[face][index]
                )
                exit_first_factor_failures[
                    failed_factor_bits[0] if failed_factor_bits else None
                ] += 1

                carrier_before = carrier_signature(
                    piece[:2], indexed_carrier_twists
                )
                carrier_after = carrier_signature(
                    exited[:2], indexed_carrier_twists
                )
                failed_carrier_bits = tuple(
                    (side, index)
                    for side in range(2)
                    for index in range(P)
                    if carrier_before[side][index]
                    != carrier_after[side][index]
                )
                exit_first_carrier_failures[
                    failed_carrier_bits[0]
                    if failed_carrier_bits else None
                ] += 1

                source_deltas = address_translation_matches(
                    endpoint_profile(piece[:2], indexed_present),
                    endpoint_profile(exited[:2], indexed_present),
                )
                target_deltas = address_translation_matches(
                    target_endpoint_profile(piece[:2], indexed_present),
                    target_endpoint_profile(exited[:2], indexed_present),
                )
                exit_source_endpoint_delta_histogram[source_deltas] += 1
                exit_target_endpoint_delta_histogram[target_deltas] += 1

            block_rows.append((
                block_index,
                len(block),
                sum(piece[:2] in live_set for piece in block),
                sum(piece[:2] in dead_set for piece in block),
                len(block) - 14,
                sum(piece[:2] in live_set for piece in block[14:]),
                sum(piece[:2] in dead_set for piece in block[14:]),
                (left_live, 14 - left_live),
                (right_live, 14 - right_live),
            ))

        all_interior_factor_equal &= set(factor_pair_histogram) == {True}
        all_interior_carrier_equal &= set(carrier_pair_histogram) == {True}
        all_interior_source_endpoint_equal &= (
            set(source_endpoint_pair_histogram) == {True}
        )
        all_interior_target_endpoint_equal &= (
            set(target_endpoint_pair_histogram) == {True}
        )
        all_interior_fixed_address_equal &= (
            set(source_endpoint_pair_histogram) == {True}
            and set(target_endpoint_pair_histogram) == {True}
        )
        global_source_address_deltas = (
            row_source_deltas
            if global_source_address_deltas is None
            else global_source_address_deltas & row_source_deltas
        )
        global_target_address_deltas = (
            row_target_deltas
            if global_target_address_deltas is None
            else global_target_address_deltas & row_target_deltas
        )

        sample = blocks[0][0]
        rows.append({
            "clock": clock,
            "blocks": tuple(block_rows),
            "factor_pairs": tuple(sorted(factor_pair_histogram.items())),
            "carrier_pairs": tuple(sorted(carrier_pair_histogram.items())),
            "source_endpoint_pairs":
                tuple(sorted(source_endpoint_pair_histogram.items())),
            "target_endpoint_pairs":
                tuple(sorted(target_endpoint_pair_histogram.items())),
            "source_address_deltas":
                tuple(sorted(row_source_deltas or set())),
            "target_address_deltas":
                tuple(sorted(row_target_deltas or set())),
            "sample_source_endpoint_profile_types":
                count_profile_types(endpoint_profile(
                    sample[:2], indexed_present
                )),
            "sample_target_endpoint_profile_types":
                count_profile_types(target_endpoint_profile(
                    sample[:2], indexed_present
                )),
            "exit_first_factor_failures":
                tuple(sorted(
                    exit_first_factor_failures.items(),
                    key=lambda row: repr(row[0]),
                )),
            "exit_first_carrier_failures":
                tuple(sorted(
                    exit_first_carrier_failures.items(),
                    key=lambda row: repr(row[0]),
                )),
            "exit_source_endpoint_delta_histogram":
                tuple(sorted(
                    exit_source_endpoint_delta_histogram.items(),
                    key=lambda row: repr(row[0]),
                )),
            "exit_target_endpoint_delta_histogram":
                tuple(sorted(
                    exit_target_endpoint_delta_histogram.items(),
                    key=lambda row: repr(row[0]),
                )),
            "ancestry": (
                ancestry_report[0],
                ancestry_report[1],
                ancestry_report[3],
                ancestry_report[4],
                ancestry_report[5],
            ),
        })

    print("THM2818/2819 COFIBER WRAP FLUX HOSTILE PROBE")
    print(
        f"common_circle_T={T};tau_grid={TAU};"
        f"wrap=13tau={W};half_step={H};wrap_over_half_step={W // H}"
    )
    print(
        f"phase_increments=R*half/T=13/2;"
        f"R*wrap/T={R * W // T};"
        f"C3*half/T={full_module.C3 * H // T};"
        f"C3*wrap/T={full_module.C3 * W // T}"
    )
    for row in rows:
        print(f"clock_row={row}")
    print(
        "global_interior_preservation="
        f"factor:{all_interior_factor_equal},"
        f"carrier:{all_interior_carrier_equal},"
        f"source_endpoint_fixed:{all_interior_source_endpoint_equal},"
        f"target_endpoint_fixed:{all_interior_target_endpoint_equal},"
        f"all_fixed_address:{all_interior_fixed_address_equal}"
    )
    print(
        "global_endpoint_address_translation="
        f"source:{tuple(sorted(global_source_address_deltas or set()))},"
        f"target:{tuple(sorted(global_target_address_deltas or set()))}"
    )
    print(
        "flux_law=each block loses 14 overlap sites under the partial "
        "forward wrap, exactly 7 live and 7 dead; the incoming and outgoing "
        "boundary collars have the same (7,7) parity census"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
