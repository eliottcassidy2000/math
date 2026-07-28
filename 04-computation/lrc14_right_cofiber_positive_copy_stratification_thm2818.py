#!/usr/bin/env python3
"""Exact all-cell companion for THM-2818.

The computation reconstructs every nonzero THM-2771 right-cofiber cell,
separates raw pieces from delayed-coefficient-live copies, audits literal
THM-2584 ancestry chambers, and explains the multipliers 2,121,265,254 as
positive copy counts.  A selected-cell side audit retains native factors,
carrier masks, and endpoint support.  No Python ``assert`` is used.
"""

from bisect import bisect_right
from collections import Counter
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "lrc14_fully_marked_root_zero_target_profile_thm2749.py":
        "d67c852c52f88feaadb2fcaa0a9a07a212f2e47018040b455855df886200595e",
    "lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.py":
        "25cbed38026d61891173c687006250a69fe38aea56d67439406bd8bb60fa2552",
    "lrc14_root_zero_full_target_semantic_clutch_20260728.py":
        "208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8",
    "lrc14_full_arm_orbit_path_sheet_audit_thm2791.py":
        "1e00b6711db0d878fca70047f5f1532518084dbf6928551cd28fe51283fde543",
    "lrc14_extended_carrier_endpoint_lib.py":
        "4b3f9f195b1634e1e84a1bc8bccb878a1c8c44aec13f24d197f92547c9e36c57",
}
for name, expected in PINNED.items():
    actual = sha256(
        (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    require(actual == expected, f"pinned dependency changed: {name}")


import lrc14_extended_carrier_endpoint_lib as endpoint_base
import lrc14_fully_marked_root_zero_target_profile_thm2749 as marked
import lrc14_full_arm_orbit_path_sheet_audit_thm2791 as ancestry
import lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751 as wing
import lrc14_root_zero_full_target_semantic_clutch_20260728 as physical


P = 13
T = physical.T
SHIFT = physical.SHIFT
DEPTH = P**5
W1 = 27581135604
W2 = 27580222516
C = 103478815440
G0 = 5905329039529920
K1 = 483303
K2 = 483287
LENGTH = 26444880
HALF_STEP = T // (2 * DEPTH)
PATH = (59162, 26, 56658)
SELECTED_I = (142004992589460, 142005019034340)
SELECTED_MINUS = (142004190428100, 142004216872980)
SELECTED_PLUS = (142082000080020, 142082026524900)
FACTOR_NAMES = ("E3", "clock", "q1", "q2", "c2", "c3")


EXPECTED_MULTIPLIER = {}
for target in range(3, 12):
    EXPECTED_MULTIPLIER[1, target] = 2
EXPECTED_MULTIPLIER[1, 12] = 121
for target in range(2, 10):
    EXPECTED_MULTIPLIER[2, target] = 2
EXPECTED_MULTIPLIER[2, 12] = 265
for target in (2, 3, 4, 5, 6, 7, 10, 11):
    EXPECTED_MULTIPLIER[3, target] = 2
EXPECTED_MULTIPLIER[3, 12] = 254

EXPECTED_RAW_COUNT = {
    **{
        cell: 2
        for cell, multiplier in EXPECTED_MULTIPLIER.items()
        if multiplier == 2
    },
    (1, 12): 241,
    (2, 12): 528,
    (3, 12): 506,
}

EXPECTED_EXCEPTIONAL_BLOCKS = {
    1: ((145, 73, 72), (96, 48, 48)),
    2: ((143, 72, 71), (289, 145, 144), (96, 48, 48)),
    3: ((143, 72, 71), (289, 145, 144), (74, 37, 37)),
}

EXPECTED_LABELS = {
    1: (
        966606,
        28534,
        W1,
        "15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd",
    ),
    2: (
        966574,
        28534,
        W2,
        "bc32242480e738421f7f9d7c1f300c09077d9c76ef4d1156083a87673b4d598c",
    ),
    3: (
        966574,
        28534,
        W2,
        "bc32242480e738421f7f9d7c1f300c09077d9c76ef4d1156083a87673b4d598c",
    ),
}


def weighted_intersection(pieces, intervals):
    return tuple(
        physical.relative.private.old.intersect_weighted_union(
            pieces, intervals
        )
    )


def intersect_sorted(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def contained(interval, intervals):
    return intersect_sorted((interval,), intervals) == (interval,)


def support_of(weighted):
    return physical.overlap.merge_intervals(
        (left, right) for left, right, weight in weighted if weight
    )


def physical_setup():
    (
        module, _prefixes, _whole, _masses, rails, present, _starts
    ) = physical.relative.lift.m.core.build_carrier_data()
    pair_prefixes = physical.relative.private.build_pair_prefixes(module)
    _sv, _tv, rail_pairs, details = physical.overlap.overlap_vectors(
        module, pair_prefixes, rails, present, rail_index=8
    )
    require(
        tuple(rails[8][:3]) == (1, 4, 12)
        and all(a == b for _left, _right, a, b in rail_pairs),
        "rail-eight metadata changed",
    )
    full_module = physical.target.load_present_module()
    e3 = physical.target.exclusive_source(full_module, 3)
    fork = physical.target.deepest_fork(full_module)
    clocks = tuple(
        full_module.make_comb(
            full_module.C1, 182, 26 * clock - 13, 26 * clock + 13
        )
        for clock in range(7)
    )
    q_pairs = physical.q_restricted_pair_prefixes(
        full_module, pair_prefixes, fork
    )
    delayed = marked.marked_prefixes(
        module,
        marked.private.build_pair_prefixes(module),
        marked.two.deepest_fork(module),
    )
    source_weight, target_weight, rail_common = marked.rail_data(
        rails, marked.RAIL
    )
    return (
        module, rails, present, details, full_module, e3, clocks, q_pairs,
        delayed, source_weight, target_weight, rail_common,
    )


def physical_objects(details, full_module, e3, clocks, clock, target):
    section = physical.target.source_present_section(
        full_module, e3, clock, 0, target, clocks
    )
    source_base, target_base = details[clock]
    source = weighted_intersection(source_base, section)
    target_physical = weighted_intersection(target_base, section)
    target_pullback = physical.overlap.shift_weighted(target_physical, -SHIFT)
    aligned = physical.overlap.intersect_weighted_profiles(
        source, target_pullback
    )
    common = tuple(
        (left, right, a)
        for left, right, a, b in aligned
        if a == b
    )
    right = physical.subtract_weighted(target_pullback, common)
    return source, target_physical, target_pullback, common, right


def raw_right(module, rail_common, clocks, clock, target):
    source = marked.two.exclusive_source(module, 3)
    raw = tuple(marked.two.intersect_sorted(source, clocks[clock]))
    raw = tuple(module.subtract_comb(raw, module.W[1], 182, -13, 13))
    raw = tuple(module.subtract_comb(raw, module.C2, 182, -13, 13))
    raw = tuple(module.subtract_comb(
        raw, module.W[2], 182, -14 * target - 13, -14 * target + 13
    ))
    raw = tuple(module.subtract_comb(
        raw, module.C3, 182, 14 * target - 13, 14 * target + 13
    ))
    left = marked.intersect(rail_common, raw)
    pulled = marked.intersect(rail_common, marked.shift_union(raw, SHIFT))
    return wing.difference(pulled, left)


def ancestry_setup():
    e_set = tuple(ancestry.base.build_set(
        ancestry.base.PAT_E, ancestry.base.ZELL
    ))
    q_set = tuple(ancestry.base.build_set(
        ancestry.host.PAT_QB, ancestry.base.ZELL
    ))
    arguments = (
        *ancestry.scaled_intervals(q_set, DEPTH),
        *ancestry.scaled_intervals(e_set, DEPTH * P**2),
        *ancestry.scaled_intervals(e_set, DEPTH),
    )
    u_events = tuple(sorted(set(
        ancestry.mapped_events(q_set, DEPTH)
        + ancestry.mapped_events(e_set, DEPTH * P**2)
    )))
    v_events = ancestry.mapped_events(e_set, DEPTH, T // P)
    return arguments, u_events, v_events


def chamber(events, left, right):
    index = bisect_right(events, left)
    require(index > 0 and index < len(events), "ancestry chamber wrapped")
    require(events[index] > right, "ancestry wall splits a cell")
    return events[index - 1], events[index]


def label_report(pieces, arguments, u_events, v_events):
    hull = (pieces[0][0], pieces[-1][1])
    chambers = chamber(u_events, *hull), chamber(v_events, *hull)
    coordinate = (pieces[0][0] + pieces[0][1]) // 2
    u_labels, v_labels = ancestry.contributor_sets(coordinate, *arguments)
    supplied = (
        PATH[0] * P**2 + PATH[1] in u_labels
        and PATH[2] in v_labels
    )
    return (
        len(u_labels),
        len(v_labels),
        len(u_labels) * len(v_labels),
        ancestry.path_digest(u_labels, v_labels),
        supplied,
        chambers,
    )


def semantic_partition(pieces, q_pair):
    source_cache = {}
    target_cache = {}
    live = []
    dead = []
    hostile = []
    for piece in pieces:
        left, right, _weight = piece
        unit = ((left, right, 1),)
        target_unit = physical.overlap.shift_weighted(unit, SHIFT)
        source_value = physical.relative.private.delayed_carry_pair(
            unit, q_pair, source_cache
        )[12][1]
        target_value = physical.relative.private.delayed_carry_pair(
            target_unit, q_pair, target_cache
        )[6][1]
        if (source_value, target_value) == (C, C):
            live.append(piece)
        elif (source_value, target_value) == (0, 0):
            dead.append(piece)
        else:
            hostile.append((piece, source_value, target_value))
    require(not hostile, "a piece acquired signed/nonuniform content")
    return tuple(live), tuple(dead)


def alternating_blocks(pieces, live):
    live_support = {piece[:2] for piece in live}
    blocks = []
    current = []
    for index, piece in enumerate(pieces):
        if index and piece[0] - pieces[index - 1][0] != HALF_STEP:
            blocks.append(tuple(current))
            current = []
        current.append(piece[:2] in live_support)
    blocks.append(tuple(current))
    require(
        all(
            bit == (index % 2 == 0)
            for block in blocks
            for index, bit in enumerate(block)
        ),
        "exceptional chain stopped alternating from a live head",
    )
    return tuple(
        (len(block), sum(block), len(block) - sum(block))
        for block in blocks
    )


def section_factors(full_module, e3, clocks):
    universe = ((0, full_module.T),)
    return (
        e3,
        clocks[1],
        full_module.subtract_comb(
            universe, full_module.W[1], 182, -13, 13
        ),
        full_module.subtract_comb(
            universe, full_module.W[2], 182, -14 * 4 - 13, -14 * 4 + 13
        ),
        full_module.subtract_comb(
            universe, full_module.C2, 182, -13, 13
        ),
        full_module.subtract_comb(
            universe, full_module.C3, 182, 14 * 4 - 13, 14 * 4 + 13
        ),
    )


def factor_masks(interval, factors):
    target_interval = tuple(endpoint + SHIFT for endpoint in interval)
    return (
        tuple(contained(interval, factor) for factor in factors),
        tuple(contained(
            interval, physical.overlap.shift_union(factor, -SHIFT)
        ) for factor in factors),
        tuple(contained(target_interval, factor) for factor in factors),
        tuple(contained(
            target_interval, physical.overlap.shift_union(factor, SHIFT)
        ) for factor in factors),
    )


def carrier_masks(interval, source, target):
    target_interval = tuple(endpoint + SHIFT for endpoint in interval)
    unit = T // P
    source_mask = tuple(bool(intersect_sorted(
        (interval,),
        support_of(physical.overlap.shift_weighted(source, -twist * unit)),
    )) for twist in range(P))
    target_mask = tuple(bool(intersect_sorted(
        (target_interval,),
        support_of(physical.overlap.shift_weighted(target, twist * unit)),
    )) for twist in range(P))
    return source_mask, target_mask


def endpoint_mask(interval, present_sets):
    values = []
    partial = []
    for address in endpoint_base.KEYS:
        overlap = intersect_sorted((interval,), present_sets[address])
        if overlap == (interval,):
            values.append(1)
        elif not overlap:
            values.append(0)
        else:
            values.append(2)
            partial.append((address, overlap))
    require(not partial, "endpoint mask cuts a selected copy")
    return tuple(values)


def translated_masks(source, target):
    source_bank = dict(zip(endpoint_base.KEYS, source))
    target_bank = dict(zip(endpoint_base.KEYS, target))
    return tuple(
        delta
        for delta in endpoint_base.KEYS
        if all(
            target_bank[address]
            == source_bank[((address[0] + delta[0]) % P,
                            (address[1] + delta[1]) % P)]
            for address in endpoint_base.KEYS
        )
    )


def selected_sidecars(full_module, e3, clocks, source, target):
    intervals = (SELECTED_I, SELECTED_MINUS, SELECTED_PLUS)
    factors = section_factors(full_module, e3, clocks)
    masks = tuple(factor_masks(interval, factors) for interval in intervals)
    expected_masks = (
        ((True,) * 6,) * 4,
        (
            (False, True, True, True, True, True),
            (True,) * 6,
            (True,) * 6,
            (False, True, True, True, True, True),
        ),
        (
            (False, True, True, True, False, True),
            (True,) * 6,
            (True,) * 6,
            (False, True, True, True, False, True),
        ),
    )
    require(masks == expected_masks, "selected native/pulled masks changed")

    delta_zero = (True,) + (False,) * 12
    zero = (False,) * 13
    twist_masks = tuple(carrier_masks(interval, source, target)
                        for interval in intervals)
    require(
        twist_masks == (
            (delta_zero, delta_zero),
            (zero, delta_zero),
            (zero, delta_zero),
        ),
        "selected carrier masks changed",
    )

    present_sets = endpoint_base.present_cache()
    endpoint_rows = []
    for interval in intervals:
        source_mask = endpoint_mask(interval, present_sets)
        target_interval = tuple(endpoint + SHIFT for endpoint in interval)
        target_mask = endpoint_mask(target_interval, present_sets)
        endpoint_rows.append((source_mask, target_mask))
    endpoint_counts = tuple(
        ((row[0].count(0), row[0].count(1)),
         (row[1].count(0), row[1].count(1)))
        for row in endpoint_rows
    )
    require(
        endpoint_counts == (
            ((88, 81), (88, 81)),
            ((169, 0), (88, 81)),
            ((79, 90), (70, 99)),
        ),
        "selected endpoint support changed",
    )
    comparisons = (
        sum(a != b for a, b in zip(endpoint_rows[0][0], endpoint_rows[1][0])),
        sum(a != b for a, b in zip(endpoint_rows[0][0], endpoint_rows[2][0])),
        sum(a != b for a, b in zip(endpoint_rows[0][1], endpoint_rows[1][1])),
        sum(a != b for a, b in zip(endpoint_rows[0][1], endpoint_rows[2][1])),
        translated_masks(endpoint_rows[0][0], endpoint_rows[1][0]),
        translated_masks(endpoint_rows[0][0], endpoint_rows[2][0]),
        translated_masks(endpoint_rows[0][1], endpoint_rows[1][1]),
        translated_masks(endpoint_rows[0][1], endpoint_rows[2][1]),
    )
    require(
        comparisons == (81, 27, 0, 18, (), (), ((0, 0),), ()),
        "selected endpoint comparison changed",
    )
    return masks, twist_masks, endpoint_counts, comparisons


def main():
    (
        module, _rails, present, details, full_module, e3, clocks, q_pairs,
        delayed, _source_weight, target_weight, rail_common,
    ) = physical_setup()
    arguments, u_events, v_events = ancestry_setup()

    reports = []
    exceptional_blocks = {}
    first_raw_nonuniform = None
    selected_source = selected_target = None
    for cell in sorted(EXPECTED_MULTIPLIER):
        clock, target = cell
        multiplier = EXPECTED_MULTIPLIER[cell]
        source, target_physical, _target_pullback, _common, pieces = (
            physical_objects(details, full_module, e3, clocks, clock, target)
        )
        if cell == (1, 4):
            selected_source = source
            selected_target = target_physical
        expected_weight = W1 if clock == 1 else W2
        require(
            len(pieces) == EXPECTED_RAW_COUNT[cell]
            and Counter(weight for _left, _right, weight in pieces)
            == Counter({expected_weight: len(pieces)})
            and set(right - left for left, right, _weight in pieces) == {LENGTH},
            f"raw equal-copy geometry changed at {cell}",
        )
        live, dead = semantic_partition(pieces, q_pairs[clock])
        require(
            len(live) == multiplier and len(dead) == len(pieces) - multiplier,
            f"live/dead copy count changed at {cell}",
        )
        if live and dead and first_raw_nonuniform is None:
            first_raw_nonuniform = cell
        if target == 12:
            exceptional_blocks[clock] = alternating_blocks(pieces, live)
            require(
                exceptional_blocks[clock] == EXPECTED_EXCEPTIONAL_BLOCKS[clock],
                f"exceptional chain profile changed at {cell}",
            )

        labels = label_report(pieces, arguments, u_events, v_events)
        require(
            labels[:4] == EXPECTED_LABELS[clock] and labels[4],
            f"ancestry labels changed at {cell}",
        )

        raw_support = raw_right(module, rail_common, clocks, clock, target)
        total_vector, _masses = wing.coefficient_vector(
            module, delayed, present, target_weight, raw_support, 12,
            [{} for _ in range(7)],
        )
        copy_value = expected_weight * C
        expected_copy_vector = (0,) + (copy_value,) * 6
        expected_total_vector = (0,) + (multiplier * copy_value,) * 6
        copy_vector, _copy_masses = wing.coefficient_vector(
            module, delayed, present, target_weight, (live[0][:2],), 12,
            [{} for _ in range(7)],
        )
        require(
            copy_vector == expected_copy_vector
            and total_vector == expected_total_vector,
            f"coefficient copy sum changed at {cell}",
        )
        if dead:
            dead_vector, _dead_masses = wing.coefficient_vector(
                module, delayed, present, target_weight, (dead[0][:2],), 12,
                [{} for _ in range(7)],
            )
            require(dead_vector == (0,) * 7, f"dead copy contributes at {cell}")

        k = K1 if clock == 1 else K2
        copy_beta = (copy_value // P) % P
        total_beta = (multiplier * copy_value // P) % P
        require(
            copy_value == G0 * k
            and copy_value % P == 0
            and copy_value % P**2 != 0
            and total_beta == multiplier * copy_beta % P,
            f"copywise Bockstein changed at {cell}",
        )
        reports.append((
            cell, len(pieces), len(live), len(dead), expected_weight,
            labels[0], labels[1], labels[3], labels[5], copy_value,
            copy_beta, total_beta,
        ))

    require(
        len(reports) == 28 and first_raw_nonuniform == (1, 12),
        "cell census or first raw failure changed",
    )
    first_global_label_change = next(
        row[0] for row in reports if row[7] != EXPECTED_LABELS[1][3]
    )
    require(first_global_label_change == (2, 2), "label prototype moved")
    require(selected_source is not None and selected_target is not None,
            "selected cell missing")
    sidecars = selected_sidecars(
        full_module, e3, clocks, selected_source, selected_target
    )

    report_digest = sha256(repr(tuple(reports)).encode()).hexdigest()
    count_rows = tuple((row[0], row[1], row[2], row[3]) for row in reports)
    print("THM-2818 RIGHT-COFIBER POSITIVE COPY STRATIFICATION")
    print("status=VERIFIED-EXACT; no row exclusion or LRC14 conclusion")
    print(f"nonzero_cells={tuple(sorted(EXPECTED_MULTIPLIER))}")
    print(f"raw_live_dead_by_cell={count_rows}")
    print(f"full_report_sha256={report_digest}")
    print(
        f"copy_chambers=(clock1=(U,V,w)={EXPECTED_LABELS[1][:3]},"
        f"clock2_or3=(U,V,w)={EXPECTED_LABELS[2][:3]})"
    )
    print(
        f"copy_values=(clock1={W1 * C},clock2_or3={W2 * C});"
        "clock_vector=(0,value,value,value,value,value,value)"
    )
    print(
        "primitive_multipliers=(2,121,265,254)="
        "literal_positive_live_copy_counts;signed_cancellation=none"
    )
    print(
        f"exceptional_blocks={exceptional_blocks};half_step={HALF_STEP};"
        "word=1010...;live=sum_ceiling(length/2)"
    )
    print(
        f"first_raw_uniform_failure={first_raw_nonuniform};"
        f"first_global_label_change={first_global_label_change};"
        "U_loss=32;V_and_supplied_path=fixed"
    )
    print(
        "copywise_Bockstein=(clock1=9,clock2_or3=2);"
        "additive_in_all_28_cells=yes"
    )
    print(
        f"selected_factor_order={FACTOR_NAMES};masks={sidecars[0]}"
    )
    print(
        f"selected_carrier_masks={sidecars[1]};"
        f"endpoint_zero_one_counts={sidecars[2]};"
        f"endpoint_comparison={sidecars[3]}"
    )
    print(
        "boundary=positive copy counting explains every raw coefficient and "
        "local Bockstein, but native factors, source carrier support, and "
        "endpoint masks block a common physical cospan; target convolution "
        "and root-deck transport remain absent"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
