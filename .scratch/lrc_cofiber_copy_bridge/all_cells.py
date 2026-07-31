#!/usr/bin/env python3
"""All-cell geometric-copy audit for the THM-2771 right cofiber.

This follows the rail-eight experiment in ``audit.py`` through every one of
the 28 nonzero physical (present-clock,target-label) cells.  It distinguishes
raw right-wing pieces from the pieces that survive the delayed coefficient
filter, then tests whether the primitive multipliers 2, 121, 265, and 254
are literal positive copy counts.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import Counter
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BRIDGE_PATH = Path(__file__).with_name("audit.py")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


require(
    sha256(BRIDGE_PATH.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    == "a154a4fba6e4901264e2228b29b980eadff3702d69b04d0061c77c695c2462cf",
    "pinned rail-eight bridge audit changed",
)
spec = spec_from_file_location("lrc_cofiber_copy_bridge_audit", BRIDGE_PATH)
require(spec is not None and spec.loader is not None, "bridge import failed")
bridge = module_from_spec(spec)
spec.loader.exec_module(bridge)


P = bridge.P
T = bridge.T
SHIFT = bridge.SHIFT
W1 = bridge.W
W2 = 27580222516
C = bridge.C
G0 = bridge.G0
K1 = bridge.K1
K2 = 483287
LENGTH = bridge.I[1] - bridge.I[0]
HALF_STEP = T // (2 * bridge.DEPTH)
PATH = bridge.PATH

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


def physical_setup():
    (
        module, _prefixes, _whole, _masses, rails, present, _starts
    ) = bridge.physical.relative.lift.m.core.build_carrier_data()
    pair_prefixes = bridge.physical.relative.private.build_pair_prefixes(
        module
    )
    _sv, _tv, _rail_pairs, details = (
        bridge.physical.overlap.overlap_vectors(
            module, pair_prefixes, rails, present, rail_index=8
        )
    )
    full_module = bridge.physical.target.load_present_module()
    e3 = bridge.physical.target.exclusive_source(full_module, 3)
    fork = bridge.physical.target.deepest_fork(full_module)
    clocks = tuple(
        full_module.make_comb(
            full_module.C1, 182, 26 * clock - 13, 26 * clock + 13
        )
        for clock in range(7)
    )
    q_pairs = bridge.physical.q_restricted_pair_prefixes(
        full_module, pair_prefixes, fork
    )
    delayed = bridge.marked.marked_prefixes(
        module,
        bridge.marked.private.build_pair_prefixes(module),
        bridge.marked.two.deepest_fork(module),
    )
    source_weight, target_weight, rail_common = bridge.marked.rail_data(
        rails, bridge.marked.RAIL
    )
    return (
        module, rails, present, details, full_module, e3, clocks, q_pairs,
        delayed, source_weight, target_weight, rail_common,
    )


def physical_right(details, full_module, e3, clocks, clock, target):
    section = bridge.physical.target.source_present_section(
        full_module, e3, clock, 0, target, clocks
    )
    source_base, target_base = details[clock]
    source = bridge.weighted_intersection(source_base, section)
    target_physical = bridge.weighted_intersection(target_base, section)
    target_pullback = bridge.physical.overlap.shift_weighted(
        target_physical, -SHIFT
    )
    aligned = bridge.physical.overlap.intersect_weighted_profiles(
        source, target_pullback
    )
    common = tuple(
        (left, right, a)
        for left, right, a, b in aligned
        if a == b
    )
    return bridge.physical.subtract_weighted(target_pullback, common)


def raw_right(module, rail_common, clocks, clock, target):
    source = bridge.marked.two.exclusive_source(module, 3)
    raw = tuple(
        bridge.marked.two.intersect_sorted(source, clocks[clock])
    )
    raw = tuple(module.subtract_comb(
        raw, module.W[1], 182, -13, 13
    ))
    raw = tuple(module.subtract_comb(
        raw, module.C2, 182, -13, 13
    ))
    raw = tuple(module.subtract_comb(
        raw, module.W[2], 182,
        -14 * target - 13, -14 * target + 13,
    ))
    raw = tuple(module.subtract_comb(
        raw, module.C3, 182,
        14 * target - 13, 14 * target + 13,
    ))
    a = bridge.marked.intersect(rail_common, raw)
    b = bridge.marked.intersect(
        rail_common, bridge.marked.shift_union(raw, SHIFT)
    )
    return bridge.wing.difference(b, a)


def ancestry_setup():
    e_set = tuple(
        bridge.ancestry.base.build_set(
            bridge.ancestry.base.PAT_E,
            bridge.ancestry.base.ZELL,
        )
    )
    q_set = tuple(
        bridge.ancestry.base.build_set(
            bridge.ancestry.host.PAT_QB,
            bridge.ancestry.base.ZELL,
        )
    )
    arguments = (
        *bridge.ancestry.scaled_intervals(q_set, bridge.DEPTH),
        *bridge.ancestry.scaled_intervals(
            e_set, bridge.DEPTH * P**2
        ),
        *bridge.ancestry.scaled_intervals(e_set, bridge.DEPTH),
    )
    u_events = tuple(sorted(set(
        bridge.ancestry.mapped_events(q_set, bridge.DEPTH)
        + bridge.ancestry.mapped_events(
            e_set, bridge.DEPTH * P**2
        )
    )))
    v_events = bridge.ancestry.mapped_events(
        e_set, bridge.DEPTH, T // P
    )
    return arguments, u_events, v_events


def chamber(events, left, right):
    index = bisect_right(events, left)
    require(index > 0 and index < len(events), "chamber wrapped")
    require(events[index] > right, "an ancestry wall splits the cell")
    return events[index - 1], events[index]


def label_report(pieces, arguments, u_events, v_events):
    hull = (pieces[0][0], pieces[-1][1])
    chambers = (
        chamber(u_events, *hull),
        chamber(v_events, *hull),
    )
    coordinate = (pieces[0][0] + pieces[0][1]) // 2
    u_labels, v_labels = bridge.ancestry.contributor_sets(
        coordinate, *arguments
    )
    supplied = (
        PATH[0] * P**2 + PATH[1] in u_labels
        and PATH[2] in v_labels
    )
    return (
        len(u_labels),
        len(v_labels),
        len(u_labels) * len(v_labels),
        bridge.ancestry.path_digest(u_labels, v_labels),
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
        target_unit = bridge.physical.overlap.shift_weighted(
            unit, SHIFT
        )
        source_value = (
            bridge.physical.relative.private.delayed_carry_pair(
                unit, q_pair, source_cache
            )[12][1]
        )
        target_value = (
            bridge.physical.relative.private.delayed_carry_pair(
                target_unit, q_pair, target_cache
            )[6][1]
        )
        if (source_value, target_value) == (C, C):
            live.append(piece)
        elif (source_value, target_value) == (0, 0):
            dead.append(piece)
        else:
            hostile.append((piece, source_value, target_value))
    require(not hostile, "a copy acquired nonuniform/signed content")
    return tuple(live), tuple(dead)


def alternating_blocks(pieces, live):
    """Split exceptional t=12 pieces at large gaps and audit 1,0 alternation."""
    live_support = {piece[:2] for piece in live}
    blocks = []
    current = []
    for index, piece in enumerate(pieces):
        if (
            index
            and piece[0] - pieces[index - 1][0] != HALF_STEP
        ):
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
        "exceptional blocks stopped alternating live/dead from live",
    )
    return tuple(
        (len(block), sum(block), len(block) - sum(block))
        for block in blocks
    )


def main():
    (
        module, rails, present, details, full_module, e3, clocks, q_pairs,
        delayed, _source_weight, target_weight, rail_common,
    ) = physical_setup()
    arguments, u_events, v_events = ancestry_setup()

    reports = []
    exceptional_blocks = {}
    first_raw_nonuniform = None
    for cell in sorted(EXPECTED_MULTIPLIER):
        clock, target = cell
        multiplier = EXPECTED_MULTIPLIER[cell]
        pieces = physical_right(
            details, full_module, e3, clocks, clock, target
        )
        expected_weight = W1 if clock == 1 else W2
        require(
            len(pieces) == EXPECTED_RAW_COUNT[cell]
            and Counter(weight for _left, _right, weight in pieces)
            == Counter({expected_weight: len(pieces)})
            and set(right - left for left, right, _weight in pieces)
            == {LENGTH},
            f"raw equal-weight copy geometry changed at {cell}",
        )
        live, dead = semantic_partition(pieces, q_pairs[clock])
        require(
            len(live) == multiplier
            and len(dead) == len(pieces) - multiplier,
            f"live/dead copy count changed at {cell}",
        )
        if live and dead and first_raw_nonuniform is None:
            first_raw_nonuniform = cell
        if target == 12:
            exceptional_blocks[clock] = alternating_blocks(
                pieces, live
            )
            require(
                exceptional_blocks[clock]
                == EXPECTED_EXCEPTIONAL_BLOCKS[clock],
                f"exceptional block profile changed at {cell}",
            )

        labels = label_report(
            pieces, arguments, u_events, v_events
        )
        expected_labels = EXPECTED_LABELS[clock]
        require(
            labels[:4] == expected_labels
            and labels[4],
            f"ancestry labels changed at {cell}",
        )

        raw_support = raw_right(
            module, rail_common, clocks, clock, target
        )
        total_vector, _masses = bridge.wing.coefficient_vector(
            module, delayed, present, target_weight, raw_support, 12,
            [{} for _ in range(7)],
        )
        copy_value = expected_weight * C
        expected_copy_vector = (0,) + (copy_value,) * 6
        expected_total_vector = (0,) + (
            multiplier * copy_value,
        ) * 6
        copy_vector, _copy_masses = bridge.wing.coefficient_vector(
            module, delayed, present, target_weight,
            (live[0][:2],), 12, [{} for _ in range(7)],
        )
        require(
            copy_vector == expected_copy_vector
            and total_vector == expected_total_vector,
            f"clock-vector copy sum changed at {cell}",
        )
        if dead:
            dead_vector, _dead_masses = bridge.wing.coefficient_vector(
                module, delayed, present, target_weight,
                (dead[0][:2],), 12, [{} for _ in range(7)],
            )
            require(
                dead_vector == (0,) * 7,
                f"dead copy contributes at {cell}",
            )

        k = K1 if clock == 1 else K2
        copy_beta = (copy_value // P) % P
        total_beta = (multiplier * copy_value // P) % P
        require(
            copy_value == G0 * k
            and copy_value % P == 0
            and copy_value % P**2 != 0
            and total_beta == multiplier * copy_beta % P,
            f"copywise Bockstein additivity changed at {cell}",
        )
        reports.append((
            cell,
            len(pieces),
            len(live),
            len(dead),
            expected_weight,
            labels[0],
            labels[1],
            labels[3],
            labels[5],
            copy_value,
            copy_beta,
            total_beta,
        ))

    require(
        len(reports) == 28
        and first_raw_nonuniform == (1, 12),
        "support census or first raw nonuniform cell changed",
    )
    first_global_label_change = next(
        row[0]
        for row in reports
        if row[7] != EXPECTED_LABELS[1][3]
    )
    require(
        first_global_label_change == (2, 2),
        "first global label-chamber change moved",
    )

    print("THM-2771 ALL-CELL GEOMETRIC COPY AUDIT")
    print("status=FINITE-EXACT scratch; no LRC conclusion")
    print(
        "columns=(cell,raw_pieces,live_copies,dead_zero_copies,"
        "weight,U_count,V_count,label_digest,U/V_chambers,"
        "copy_value,copy_beta,total_beta)"
    )
    for row in reports:
        print(row)
    print(
        f"individual_copy_clock_vectors="
        f"(clock1={(0,) + (W1 * C,) * 6},"
        f"clock2_or3={(0,) + (W2 * C,) * 6})"
    )
    print(
        "primitive_multipliers=(2,121,265,254)="
        "literal_positive_live_copy_counts; signed_cancellation=none"
    )
    print(
        "exceptional_raw_partitions="
        "(e1,t12):241=121_live+120_zero;"
        "(e2,t12):528=265_live+263_zero;"
        "(e3,t12):506=254_live+252_zero"
    )
    print(
        f"exceptional_alternating_blocks={exceptional_blocks};"
        f"within_block_left_endpoint_step={HALF_STEP};"
        "each_block_word=1010...;"
        "live_count=sum_ceiling(block_length/2)"
    )
    print(
        f"first_raw_uniform_content_failure={first_raw_nonuniform};"
        "mechanism=delayed_terminal_filter_splits_c_from_zero"
    )
    print(
        "cellwise_live_uniformity=PASS_28/28;"
        "cellwise_ancestry_identity=PASS_by_wall_free_chambers"
    )
    print(
        f"first_global_label_prototype_change={first_global_label_change};"
        "mechanism=U loses 32 labels and weight changes "
        f"{W1}->{W2}; V and supplied path stay fixed"
    )
    print(
        "Bockstein=copywise_additive_in_every_cell;"
        "clock1_copy_beta=9;clock2_or3_copy_beta=2"
    )
    print(
        "boundary=literal positive copy counts explain the raw coefficient "
        "and local mod13 Bockstein, but no endpoint/root-deck/cospan map is "
        "created by counting"
    )


if __name__ == "__main__":
    main()
