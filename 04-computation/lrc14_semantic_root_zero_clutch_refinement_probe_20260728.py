#!/usr/bin/env python3
"""Exact semantic-sheet refinement of the relative-present root-zero clutch.

This probe inserts the full source condition E3, the delayed word
Q_(3,{1,2}), and one common lawful target sheet U_(0,3) into the partial
adjacent-chart overlap of the root-zero clutch.  It then recomputes the two
seven-clock coefficient vectors independently.  It is exploratory until its
output is audited and routed; no LRC(14) conclusion is asserted here.
"""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import hashlib
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


DEPENDENCY = COMP / "lrc14_root_zero_overlap_clutch_20260728.py"
require(
    hashlib.sha256(DEPENDENCY.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    == "e10fa7c9a5a238461ef422ea314dc334f7e65ec1787cf65d4e4bea12b96aefb8",
    "root-zero clutch dependency changed",
)

import lrc14_root_zero_overlap_clutch_20260728 as clutch


relative = clutch.relative
old = relative.private.old
P = relative.P
T = relative.lift.m.T
SHEET = (0, 3)


def intersect_danger(pieces, module, speed, pd, lo, hi):
    danger = module.make_comb(speed, pd, lo, hi)
    return tuple(old.intersect_weighted_union(pieces, danger))


def intersect_safe(pieces, module, speed, pd, lo, hi):
    danger = module.make_comb(speed, pd, lo, hi)
    safe = clutch.complement(danger)
    return tuple(old.intersect_weighted_union(pieces, safe))


def restrict_e3_and_sheet(pieces, module, s, t):
    """Insert E3 and the four shifted factors defining U_(s,t)."""
    pieces = intersect_safe(pieces, module, 1, 91, -13, 13)
    for speed in relative.semantic.UNITS:
        pieces = intersect_safe(pieces, module, speed, 182, -13, 13)
    pieces = intersect_safe(
        pieces, module, relative.semantic.BLOCKERS[0], 182, -13, 13
    )
    pieces = intersect_safe(
        pieces, module, relative.semantic.BLOCKERS[1], 182, -13, 13
    )
    pieces = intersect_danger(
        pieces, module, relative.semantic.BLOCKERS[2], 182, -13, 13
    )

    # The lawful two-target sheet uses q1/c2 at colour s and q2/c3 at t.
    pieces = intersect_safe(
        pieces, module, relative.semantic.UNITS[0], 182,
        -14 * s - 13, -14 * s + 13,
    )
    pieces = intersect_safe(
        pieces, module, relative.semantic.BLOCKERS[1], 182,
        14 * s - 13, 14 * s + 13,
    )
    pieces = intersect_safe(
        pieces, module, relative.semantic.UNITS[1], 182,
        -14 * t - 13, -14 * t + 13,
    )
    pieces = intersect_safe(
        pieces, module, relative.semantic.BLOCKERS[2], 182,
        14 * t - 13, 14 * t + 13,
    )
    return pieces


def build_q3_pair_prefixes(module):
    """Delayed sector-zero prefixes with the full Q_(3,{1,2}) word."""
    word = relative.private.prior.sector_words(module)[0]
    by_clock = []
    for ell in range(7):
        qell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        # The inherited word already has guard/unit/C3 safe and C2 danger.
        # Add the missing unshifted C1 danger to obtain Q_(3,{1,2}).
        qell = module.intersect_comb(qell, module.C1, 182, -13, 13)
        by_h = []
        for h in range(P):
            pair = []
            for kappa in range(2):
                digit = old.sat.intersect_interval(
                    qell,
                    (2 * h + kappa) * T // (2 * P),
                    (2 * h + kappa + 1) * T // (2 * P),
                )
                pair.append(module.make_prefix(digit))
            by_h.append(tuple(pair))
        by_clock.append(tuple(by_h))
    return tuple(by_clock)


def main():
    module, _prefixes, _, _, rails, present, _starts = (
        relative.lift.m.core.build_carrier_data()
    )
    pair_prefixes = build_q3_pair_prefixes(module)
    target_pullback = clutch.shift_weighted(rails[8][3], -clutch.SHIFT)
    rail_pairs = clutch.intersect_weighted_profiles(
        rails[8][3], target_pullback
    )
    require(rail_pairs and all(row[2] == row[3] for row in rail_pairs),
            "rail-eight source/target weights stopped agreeing")

    source_vector = []
    target_vector = []
    source_masses = []
    target_masses = []
    for ell in range(7):
        source, target = clutch.restrict_to_relative_overlap(
            module, present, rail_pairs, ell
        )
        source = restrict_e3_and_sheet(source, module, *SHEET)
        target = restrict_e3_and_sheet(target, module, *SHEET)
        source_masses.append(sum((b - a) * w for a, b, w in source))
        target_masses.append(sum((b - a) * w for a, b, w in target))
        source_values = relative.private.delayed_carry_pair(
            source, pair_prefixes[ell][6], {}
        )
        target_values = relative.private.delayed_carry_pair(
            target, pair_prefixes[ell][6], {}
        )
        source_vector.append(source_values[12][1])
        target_vector.append(target_values[6][1])

    source_vector = tuple(source_vector)
    target_vector = tuple(target_vector)
    source_amplitude = 1812281403506324508838080
    target_amplitude = 1826551436254490256030720
    require(
        tuple(source_masses) == (929934280541992017600,) * 7
        and tuple(target_masses) == (929934304688494607040,) * 7
        and source_vector == (0,) + (source_amplitude,) * 6
        and target_vector == (0,) + (target_amplitude,) * 6
        and source_vector != target_vector,
        "semantic-sheet overlap coefficient table changed",
    )
    source_unit = relative.private.is_unit(
        source_vector, clutch.SOURCE_STATE[3], clutch.CONTENT
    ) if all(value % clutch.CONTENT == 0 for value in source_vector) else False
    target_unit = relative.private.is_unit(
        target_vector, clutch.TARGET_STATE[3], clutch.CONTENT
    ) if all(value % clutch.CONTENT == 0 for value in target_vector) else False
    source_profile = clutch.normalized_profile(
        source_vector, clutch.SOURCE_STATE[3]
    )[1]
    target_profile = clutch.normalized_profile(
        target_vector, clutch.TARGET_STATE[3]
    )[1]
    require(
        source_unit and target_unit
        and source_profile == (5, 0, 0, 0, 0, 0)
        and target_profile == (9, 0, 0, 0, 0, 0),
        "semantic-sheet private-unit profiles changed",
    )

    print("LRC14 SEMANTIC ROOT-ZERO CLUTCH REFINEMENT PROBE")
    print("status=EXPLORATORY EXACT recomputation")
    print(f"sheet={SHEET} rail=8 roots=(12,1) carries=(12,6)")
    print(f"source_weighted_masses={tuple(source_masses)}")
    print(f"target_weighted_masses={tuple(target_masses)}")
    print(f"source_vector={source_vector}")
    print(f"target_vector={target_vector}")
    print(f"vectors_equal={source_vector == target_vector}")
    print(f"content26=(source:{all(v % clutch.CONTENT == 0 for v in source_vector)},target:{all(v % clutch.CONTENT == 0 for v in target_vector)})")
    print(f"private_units=(source_root12:{source_unit},target_root1:{target_unit})")
    print(f"normalized_profiles=(source:{source_profile},target:{target_profile})")
    print("SCOPE: full E3, delayed Q_(3,{1,2}), and U_(0,3) inside the partial chart overlap; no endpoint current or LRC14 conclusion")


if __name__ == "__main__":
    main()
