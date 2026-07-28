#!/usr/bin/env python3
"""Scratch hostile census for the full physical-clock wing over 81 labels."""

from collections import Counter
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_root_zero_full_target_semantic_clutch_20260728 as marked


P = marked.P
CONTENT = marked.CONTENT


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def reduce_row(values, root):
    require(all(value % CONTENT == 0 for value in values), "content failure")
    inverse = pow(root, -1, P)
    residues = tuple(value // CONTENT * inverse % P for value in values)
    profile = tuple((residues[index] - residues[-1]) % P for index in range(6))
    det = marked.relative.private.old.sat.multiplication_determinant_7(profile)
    return residues, profile, det


def phi7_inverse(value):
    columns = []
    for index in range(6):
        basis = tuple(int(position == index) for position in range(6))
        columns.append(marked.multiply_phi_seven(value, basis))
    matrix = [
        [columns[column][row] % P for column in range(6)]
        + [int(row == 0)]
        for row in range(6)
    ]
    for column in range(6):
        pivot = next((row for row in range(column, 6) if matrix[row][column]), None)
        require(pivot is not None, "nonunit requested")
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        scale = pow(matrix[column][column], -1, P)
        matrix[column] = [(entry * scale) % P for entry in matrix[column]]
        for row in range(6):
            if row == column:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                (matrix[row][entry] - scale * matrix[column][entry]) % P
                for entry in range(7)
            ]
    inverse = tuple(matrix[row][-1] for row in range(6))
    require(
        marked.multiply_phi_seven(value, inverse) == (1, 0, 0, 0, 0, 0),
        "bad inverse",
    )
    return inverse


def main():
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
    q_pairs = marked.q_restricted_pair_prefixes(full_module, pair_prefixes, fork)
    common_s = (0, 1, 2, 3, 8, 9, 10, 11, 12)
    common_t = (3, 4, 5, 6, 7, 8, 9, 10, 11)

    determinant_pairs = Counter()
    augmentation_pairs = Counter()
    ratio_histogram = Counter()
    support_pairs = Counter()
    label_rows = []
    literal_checks = 0
    chosen = None
    for s in common_s:
        for t in common_t:
            source_wing = []
            target_wing = []
            for clock, (source_pieces, target_pieces) in enumerate(overlap_details):
                section = marked.target.source_present_section(
                    full_module, e3, clock, s, t, clock_comb
                )
                common_source = marked.common_physical_cut(source_pieces, section, -1)
                common_target = marked.common_physical_cut(target_pieces, section, +1)
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
                pair = q_pairs[clock]
                source_common_value = marked.relative.private.delayed_carry_pair(
                    common_source, pair, {}
                )[12][1]
                target_common_value = marked.relative.private.delayed_carry_pair(
                    common_target, pair, {}
                )[6][1]
                require(source_common_value == target_common_value, "common mismatch")
                source_one_value = marked.relative.private.delayed_carry_pair(
                    one_source, pair, {}
                )[12][1]
                target_one_value = marked.relative.private.delayed_carry_pair(
                    one_target, pair, {}
                )[6][1]
                source_difference = source_one_value - source_common_value
                target_difference = target_one_value - target_common_value
                literal_source = marked.subtract_weighted(one_source, common_source)
                literal_target = marked.subtract_weighted(one_target, common_target)
                literal_source_value = marked.relative.private.delayed_carry_pair(
                    literal_source, pair, {}
                )[12][1]
                literal_target_value = marked.relative.private.delayed_carry_pair(
                    literal_target, pair, {}
                )[6][1]
                literal_checks += int(
                    source_difference == literal_source_value
                    and target_difference == literal_target_value
                )
                source_wing.append(source_difference)
                target_wing.append(target_difference)

            source_wing = tuple(source_wing)
            target_wing = tuple(target_wing)
            source_residues, source_profile, source_det = reduce_row(source_wing, 12)
            target_residues, target_profile, target_det = reduce_row(target_wing, 1)
            determinant_pairs[(source_det, target_det)] += 1
            augmentation_pairs[(sum(source_residues) % P, sum(target_residues) % P)] += 1
            source_support = tuple(i for i, value in enumerate(source_profile) if value)
            target_support = tuple(i for i, value in enumerate(target_profile) if value)
            support_pairs[(source_support, target_support)] += 1
            ratio = None
            if source_det and target_det:
                ratio = marked.multiply_phi_seven(
                    target_profile, phi7_inverse(source_profile)
                )
                ratio_histogram[ratio] += 1
                require(
                    marked.multiply_phi_seven(ratio, source_profile) == target_profile,
                    "ratio failure",
                )
            label_rows.append(
                (s, t, source_profile, target_profile, source_det, target_det, ratio)
            )
            if (s, t) == (0, 3):
                chosen = label_rows[-1]

    require(literal_checks == 81 * 7, "literal wing check failure")
    print("all-81 physical-present wing census")
    print(f"labels={len(label_rows)} literal_object_subtractions={literal_checks}/567")
    print(f"determinant_pairs={dict(sorted(determinant_pairs.items()))}")
    print(f"augmentation_pairs={dict(sorted(augmentation_pairs.items()))}")
    print(f"support_pair_count={len(support_pairs)} support_pairs={dict(support_pairs)}")
    print(f"unit_ratio_labels={sum(ratio_histogram.values())}/81 distinct_ratios={len(ratio_histogram)}")
    for ratio, count in sorted(ratio_histogram.items(), key=lambda item: (-item[1], item[0])):
        determinant = marked.relative.private.old.sat.multiplication_determinant_7(ratio)
        print(f"ratio={ratio} count={count} det={determinant}")
    print(f"chosen_0_3={chosen}")
    exceptional = [row for row in label_rows if not row[4] or not row[5]]
    print(f"nonunit_labels={len(exceptional)}")
    for row in exceptional:
        print(f"nonunit={row}")


if __name__ == "__main__":
    main()
