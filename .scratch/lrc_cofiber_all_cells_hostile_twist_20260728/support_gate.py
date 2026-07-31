#!/usr/bin/env python3
"""Locate the first support gate for all 91 physical cofiber cells.

This is a small independent sidecar to the all-cell hostile audit.  It
compares the exact target support with THM-2809's eleven-label residual
without identifying their failure mechanisms.
"""

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


PHYSICAL_NAME = "lrc14_root_zero_full_target_semantic_clutch_20260728.py"
PHYSICAL_HASH = (
    "208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8"
)
actual_hash = sha256(
    (COMP / PHYSICAL_NAME).read_bytes().replace(b"\r\n", b"\n")
).hexdigest()
require(actual_hash == PHYSICAL_HASH, "pinned physical dependency changed")


import lrc14_root_zero_full_target_semantic_clutch_20260728 as physical


P = physical.P
SHIFT = physical.SHIFT


def weighted_intersection(weighted, intervals):
    return tuple(
        physical.relative.private.old.intersect_weighted_union(
            weighted, intervals
        )
    )


def setup():
    (
        module, _prefixes, _whole, _masses, rails, present, _starts
    ) = physical.relative.lift.m.core.build_carrier_data()
    pair_prefixes = physical.relative.private.build_pair_prefixes(module)
    _sv, _tv, _rail_pairs, details = (
        physical.overlap.overlap_vectors(
            module, pair_prefixes, rails, present, rail_index=8
        )
    )
    full_module = physical.target.load_present_module()
    e3 = physical.target.exclusive_source(full_module, 3)
    clocks = tuple(
        full_module.make_comb(
            full_module.C1,
            182,
            26 * clock - 13,
            26 * clock + 13,
        )
        for clock in range(7)
    )
    return details, full_module, e3, clocks


def stage_record(details, full_module, e3, clocks, clock, target):
    section = physical.target.source_present_section(
        full_module, e3, clock, 0, target, clocks
    )
    if not section:
        return "section_empty", (0, 0, 0, 0)

    source_base, target_base = details[clock]
    source = weighted_intersection(source_base, section)
    target_native = weighted_intersection(target_base, section)
    if not source and not target_native:
        return "both_bases_empty", (
            len(section), 0, 0, 0
        )
    if source and not target_native:
        return "source_only", (
            len(section), len(source), 0, 0
        )
    if target_native and not source:
        target_pullback = physical.overlap.shift_weighted(
            target_native, -SHIFT
        )
        return "target_only", (
            len(section), 0, len(target_native), len(target_pullback)
        )

    target_pullback = physical.overlap.shift_weighted(
        target_native, -SHIFT
    )
    aligned = physical.overlap.intersect_weighted_profiles(
        source, target_pullback
    )
    common = tuple(
        (left, right, left_weight)
        for left, right, left_weight, right_weight in aligned
        if left_weight == right_weight
    )
    require(
        all(
            left_weight == right_weight
            for _left, _right, left_weight, right_weight in aligned
        ),
        f"unequal common weights at {(clock, target)}",
    )
    right = physical.subtract_weighted(target_pullback, common)
    stage = "right_nonempty" if right else "fully_common"
    return stage, (
        len(section), len(source), len(target_native), len(right)
    )


def main():
    details, full_module, e3, clocks = setup()
    records = {}
    histogram = Counter()
    supports = {clock: [] for clock in range(7)}
    for clock in range(7):
        for target in range(P):
            stage, counts = stage_record(
                details, full_module, e3, clocks, clock, target
            )
            records[clock, target] = stage, counts
            histogram[stage] += 1
            if stage in {"right_nonempty", "target_only"}:
                supports[clock].append(target)

    support_rows = tuple(
        (clock, tuple(targets))
        for clock, targets in supports.items()
        if targets
    )
    support_union = tuple(
        sorted({target for _clock, targets in support_rows
                for target in targets})
    )
    require(
        len(records) == 91
        and sum(histogram.values()) == 91,
        "support-gate universe changed",
    )
    require(
        support_union == tuple(range(2, P)),
        "physical support union left the eleven-label residual",
    )

    label_zero = tuple(
        (clock, *records[clock, 0]) for clock in range(7)
    )
    label_one = tuple(
        (clock, *records[clock, 1]) for clock in range(7)
    )
    require(
        all(row[1] == "section_empty" for row in label_zero),
        "label zero stopped dying at the section gate",
    )
    require(
        records[2, 1][0] == records[3, 1][0] == "source_only",
        "label-one source-only hostile changed",
    )

    print("THM-2771 PHYSICAL SUPPORT-GATE SIDECAR")
    print("status=FINITE-EXACT scratch; no common-mechanism claim")
    print(f"physical_dependency_sha256={actual_hash}")
    print(f"stage_histogram={dict(histogram)}")
    print(f"right_support_by_clock={support_rows}")
    print(f"right_support_union={support_union}")
    print(f"label_zero_stage_rows={label_zero}")
    print(f"label_one_stage_rows={label_one}")
    print(
        "connection=the support union equals THM-2809's "
        "L_star={2,...,12}"
    )
    print(
        "hostile_boundary=label0 dies at the source-present section; "
        "label1 has a nonempty section and source-only weight at clocks "
        "2,3 but no target weight, so equality of support complements does "
        "not yet identify THM-2809's delayed-digit mechanism"
    )


if __name__ == "__main__":
    main()
