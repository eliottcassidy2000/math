#!/usr/bin/env python3
"""Exploratory exact rows for the LRC address/chamber pullback."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

import lrc14_root_zero_full_target_semantic_clutch_20260728 as marked


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def selected_row(s, t):
    module, _prefixes, _, _, rails, present, _starts = (
        marked.relative.lift.m.core.build_carrier_data()
    )
    pair_prefixes = marked.relative.private.build_pair_prefixes(module)
    _raw_source, _raw_target, _rail_pairs, overlap_details = (
        marked.overlap.overlap_vectors(
            module, pair_prefixes, rails, present, rail_index=8
        )
    )
    full_module = marked.target.load_present_module()
    e3 = marked.target.exclusive_source(full_module, 3)
    fork = marked.target.deepest_fork(full_module)
    clock_comb = tuple(
        full_module.make_comb(
            full_module.C1, 182, 26 * clock - 13, 26 * clock + 13
        )
        for clock in range(7)
    )
    q_pairs = marked.q_restricted_pair_prefixes(
        full_module, pair_prefixes, fork
    )

    result = []
    for clock, (source_pieces, target_pieces) in enumerate(overlap_details):
        section = marked.target.source_present_section(
            full_module, e3, clock, s, t, clock_comb
        )
        common_source = marked.common_physical_cut(
            source_pieces, section, -1
        )
        common_target = marked.common_physical_cut(
            target_pieces, section, +1
        )
        one_source = tuple(
            marked.relative.private.old.intersect_weighted_union(
                source_pieces, section
            )
        )
        one_target = tuple(
            marked.relative.private.old.intersect_weighted_union(
                target_pieces, section
            )
        )
        wing_source = marked.subtract_weighted(one_source, common_source)
        wing_target = marked.subtract_weighted(one_target, common_target)
        pair = q_pairs[clock]

        universe = ((0, full_module.T),)
        factor_intervals = (
            e3,
            clock_comb[clock],
            full_module.subtract_comb(
                universe, full_module.W[1], 182,
                -14 * s - 13, -14 * s + 13,
            ),
            full_module.subtract_comb(
                universe, full_module.W[2], 182,
                -14 * t - 13, -14 * t + 13,
            ),
            full_module.subtract_comb(
                universe, full_module.C2, 182,
                14 * s - 13, 14 * s + 13,
            ),
            full_module.subtract_comb(
                universe, full_module.C3, 182,
                14 * t - 13, 14 * t + 13,
            ),
        )
        reconstructed = universe
        for factor in factor_intervals:
            reconstructed = marked.target.intersect_sorted(
                reconstructed, factor
            )
        require(tuple(reconstructed) == tuple(section),
                "section factorization changed")

        def coefficient(pieces, carry):
            return marked.relative.private.delayed_carry_pair(
                pieces, pair, {}
            )[carry][1]

        def mass(pieces):
            return sum((b - a) * weight for a, b, weight in pieces)

        def failure_partition(one, direction, carry):
            remaining = one
            partition = []
            for factor in factor_intervals:
                shifted = marked.overlap.shift_union(
                    factor, direction * marked.SHIFT
                )
                failed = tuple(
                    marked.relative.private.old.intersect_weighted_union(
                        remaining, marked.overlap.complement(shifted)
                    )
                )
                partition.append(
                    (
                        coefficient(failed, carry),
                        mass(failed),
                        len(failed),
                    )
                )
                remaining = tuple(
                    marked.relative.private.old.intersect_weighted_union(
                        remaining, shifted
                    )
                )
            return tuple(partition), remaining

        source_origins, source_remaining = failure_partition(
            one_source, -1, 12
        )
        target_origins, target_remaining = failure_partition(
            one_target, +1, 6
        )
        require(source_remaining == common_source,
                "source first-failure remainder changed")
        require(target_remaining == common_target,
                "target first-failure remainder changed")
        require(
            sum(row[0] for row in source_origins)
            == coefficient(wing_source, 12)
            and sum(row[0] for row in target_origins)
            == coefficient(wing_target, 6),
            "first-failure coefficient did not sum to the wing",
        )
        target_wing_pullback = marked.overlap.shift_weighted(
            wing_target, -marked.SHIFT
        )
        common_wing_pairs = marked.overlap.intersect_weighted_profiles(
            wing_source, target_wing_pullback
        )

        result.append(
            {
                "clock": clock,
                "one_source": coefficient(one_source, 12),
                "one_target": coefficient(one_target, 6),
                "common_source": coefficient(common_source, 12),
                "common_target": coefficient(common_target, 6),
                "wing_source": coefficient(wing_source, 12),
                "wing_target": coefficient(wing_target, 6),
                "wing_source_mass": mass(wing_source),
                "wing_target_mass": mass(wing_target),
                "wing_source_pieces": len(wing_source),
                "wing_target_pieces": len(wing_target),
                "source_origins": source_origins,
                "target_origins": target_origins,
                "common_wing_pair_count": len(common_wing_pairs),
                "common_wing_support_mass": sum(
                    right - left
                    for left, right, _sw, _tw in common_wing_pairs
                ),
            }
        )
    return tuple(result)


def main():
    for label in ((8, 3), (2, 3), (0, 10), (0, 3)):
        rows = selected_row(*label)
        print(f"label={label}")
        for key in (
            "one_source",
            "one_target",
            "common_source",
            "common_target",
            "wing_source",
            "wing_target",
            "wing_source_mass",
            "wing_target_mass",
            "wing_source_pieces",
            "wing_target_pieces",
            "common_wing_pair_count",
            "common_wing_support_mass",
        ):
            print(f"  {key}={tuple(row[key] for row in rows)}")
        print("  origin_order=(E3,clock,q1,q2,c2,c3)")
        print(
            "  source_origin_coefficients="
            f"{tuple(tuple(cell[0] for cell in row['source_origins']) for row in rows)}"
        )
        print(
            "  target_origin_coefficients="
            f"{tuple(tuple(cell[0] for cell in row['target_origins']) for row in rows)}"
        )
        require(
            all(
                row["one_source"] - row["common_source"]
                == row["wing_source"]
                and row["one_target"] - row["common_target"]
                == row["wing_target"]
                for row in rows
            ),
            "literal decomposition failed",
        )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
