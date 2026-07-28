#!/usr/bin/env python3
"""Exact all-label physical-clock wing chamber census.

This companion refines the MISTAKE-313 repair on the canonical rail-8
root-zero bank.  It constructs the one-sided and two-sided carriers before
subtraction, retains the physical-present clock, and verifies the resulting
source rectangle, target chamber law, unit ratios, and label characters.
"""

from collections import Counter
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

def require(condition, message):
    if not condition:
        raise RuntimeError(message)


DEPENDENCY = ROOT / "04-computation" / (
    "lrc14_root_zero_full_target_semantic_clutch_20260728.py"
)
EXPECTED_DEPENDENCY_SHA256 = (
    "208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8"
)
dependency_bytes = DEPENDENCY.read_bytes().replace(b"\r\n", b"\n")
require(
    sha256(dependency_bytes).hexdigest() == EXPECTED_DEPENDENCY_SHA256,
    "audited THM-2754 dependency changed",
)

import lrc14_root_zero_full_target_semantic_clutch_20260728 as marked


P = marked.P
CONTENT = marked.CONTENT


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


def nonconstant_character_buckets(rows, frequency_s, frequency_t):
    buckets = [[0] * P for _ in range(6)]
    for s, t, profile in rows:
        bucket = (frequency_s * s + frequency_t * t) % P
        for clock, value in enumerate(profile):
            buckets[clock][bucket] = (buckets[clock][bucket] + value) % P
    return tuple(
        len(set(clock_buckets)) > 1 for clock_buckets in buckets
    )


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
    source_table = []
    target_table = []
    s_coefficients = {
        0: (9, 2, 2), 1: (9, 2, 0), 2: (7, 0, 4),
        3: (7, 2, 3), 8: (7, 4, 4), 9: (7, 4, 4),
        10: (7, 4, 4), 11: (7, 4, 4), 12: (9, 2, 2),
    }
    for s, t, source_profile, target_profile, *_rest in label_rows:
        source_active = s not in (1, 3) and t not in (8, 9)
        a_s, b_s, c_s = s_coefficients[s]
        expected_source = (0, 0, 0, 12 * source_active % P, 0, 0)
        expected_target = (
            0,
            a_s,
            b_s * (t <= 9) % P,
            c_s * (t not in (8, 9)) % P,
            0,
            0,
        )
        require(source_profile == expected_source, "source rectangle law")
        require(target_profile == expected_target, "target chamber law")
        source_table.append((s, t, source_profile))
        target_table.append((s, t, target_profile))
    source_character_census = Counter()
    target_character_census = Counter()
    for frequency_s in range(P):
        for frequency_t in range(P):
            if (frequency_s, frequency_t) == (0, 0):
                continue
            source_character_census[
                nonconstant_character_buckets(
                    source_table, frequency_s, frequency_t
                )
            ] += 1
            target_character_census[
                nonconstant_character_buckets(
                    target_table, frequency_s, frequency_t
                )
            ] += 1
    require(
        all(any(pattern) for pattern in source_character_census)
        and all(any(pattern) for pattern in target_character_census),
        "a nontrivial label character vanished modulo the residue certificate",
    )
    require(
        determinant_pairs == Counter({
            (0, 1): 2, (0, 3): 11, (0, 4): 5, (0, 8): 8,
            (0, 10): 2, (0, 12): 4,
            (1, 3): 4, (1, 5): 20, (1, 8): 15, (1, 11): 10,
        })
        and ratio_histogram == Counter({
            (0, 4, 4, 4, 4, 10): 20,
            (9, 0, 0, 0, 0, 6): 15,
            (0, 2, 2, 2, 2, 6): 10,
            (11, 0, 0, 0, 0, 4): 4,
        })
        and source_character_census
        == Counter({(False, False, False, True, False, False): 168})
        and target_character_census
        == Counter({(False, True, True, True, False, False): 168}),
        "wing census or character law changed",
    )
    augmented_gain_histogram = Counter()
    for row in label_rows:
        source_det = row[4]
        if not source_det:
            continue
        source_augmentation = sum(row[2]) % P
        target_augmentation = sum(row[3]) % P
        augmented_gain_histogram[
            target_augmentation * pow(source_augmentation, -1, P) % P
        ] += 1
    require(
        augmented_gain_histogram == Counter({11: 20, 2: 19, 0: 10}),
        "augmented gain chamber law changed",
    )
    print("all-81 physical-present wing census")
    print(f"labels={len(label_rows)} literal_object_subtractions={literal_checks}/567")
    print(f"determinant_pairs={dict(sorted(determinant_pairs.items()))}")
    print(f"augmentation_pairs={dict(sorted(augmentation_pairs.items()))}")
    print(f"support_pair_count={len(support_pairs)} support_pairs={dict(support_pairs)}")
    print(f"unit_ratio_labels={sum(ratio_histogram.values())}/81 distinct_ratios={len(ratio_histogram)}")
    print("source_rectangle=S* x T*, S*=(0,2,8,9,10,11,12), T*=(3,4,5,6,7,10,11); profile=12z^3")
    print(f"target_s_coefficients_(a,b,c)={s_coefficients}")
    print("target_chamber_profile=a_s*z+b_s*1_(t<=9)*z^2+c_s*1_(t_not_8,9)*z^3")
    print(f"source_nontrivial_label_character_patterns={dict(source_character_census)}")
    print(f"target_nontrivial_label_character_patterns={dict(target_character_census)}")
    print(f"augmented_gain_on_49_source_unit_labels={dict(sorted(augmented_gain_histogram.items()))}")
    for ratio, count in sorted(ratio_histogram.items(), key=lambda item: (-item[1], item[0])):
        determinant = marked.relative.private.old.sat.multiplication_determinant_7(ratio)
        print(f"ratio={ratio} count={count} det={determinant}")
    ratio_symbols = {
        ratio: symbol
        for symbol, (ratio, _count) in zip(
            "ABCD",
            sorted(ratio_histogram.items(), key=lambda item: (-item[1], item[0])),
        )
    }
    print(f"ratio_symbols={ratio_symbols}; X=source-null/target-unit")
    for s in common_s:
        row = [entry for entry in label_rows if entry[0] == s]
        print(
            "ratio_grid_s=" + str(s) + " "
            + " ".join(ratio_symbols.get(entry[6], "X") for entry in row)
        )
    print(f"chosen_0_3={chosen}")
    exceptional = [row for row in label_rows if not row[4] or not row[5]]
    print(f"nonunit_labels={len(exceptional)}")


if __name__ == "__main__":
    main()
