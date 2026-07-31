#!/usr/bin/env python3
"""Exact semantic/full-target restriction of the R0/L1 overlap clutch."""

from collections import Counter
from math import gcd
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_root_zero_overlap_clutch_20260728 as clutch
import lrc14_two_target_present_semantic_attachment_probe_20260728 as semantic


P = clutch.P
H = 6
KAPPA = 1
BROADCAST_CONTENT = 57_297_240


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def build_semantic_prefixes(module, fork):
    result = []
    for word in clutch.relative.private.prior.sector_words(module):
        by_clock = []
        for ell in range(7):
            qell = module.subtract_comb(
                word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
            )
            qell = semantic.intersect_sorted(qell, fork)
            pair = []
            for kappa in range(2):
                digit = clutch.relative.private.old.sat.intersect_interval(
                    qell,
                    (2 * H + kappa) * module.T // (2 * P),
                    (2 * H + kappa + 1) * module.T // (2 * P),
                )
                pair.append(module.make_prefix(digit))
            by_clock.append(tuple(pair))
        result.append(tuple(by_clock))
    return tuple(result)


def semantic_sections(module, source, s, t):
    clocks = tuple(
        module.make_comb(
            module.C1, 182, 26 * ell - 13, 26 * ell + 13
        ) for ell in range(7)
    )
    return tuple(
        tuple(semantic.source_present_section(
            module, source, ell, s, t, clocks
        )) for ell in range(7)
    )


def restricted_vectors(module, pair_prefixes, rails, present,
                       semantic_prefixes, sections, rail_index,
                       equal_weight_only=False):
    target_pullback = clutch.shift_weighted(
        rails[rail_index][3], -clutch.SHIFT
    )
    rail_pairs = clutch.intersect_weighted_profiles(
        rails[rail_index][3], target_pullback
    )
    if equal_weight_only:
        rail_pairs = tuple(row for row in rail_pairs if row[2] == row[3])
    source_vector = []
    target_vector = []
    support_masses = []
    details = []
    for ell in range(7):
        source, target = clutch.restrict_to_relative_overlap(
            module, present, rail_pairs, ell
        )
        target_section_pullback = clutch.shift_union(
            sections[ell], -clutch.SHIFT
        )
        common_source_section = clutch.intersect_sorted(
            sections[ell], target_section_pullback
        )
        source = clutch.relative.private.old.intersect_weighted_union(
            source, common_source_section
        )

        source_section_forward = clutch.shift_union(
            sections[ell], clutch.SHIFT
        )
        common_target_section = clutch.intersect_sorted(
            sections[ell], source_section_forward
        )
        target = clutch.relative.private.old.intersect_weighted_union(
            target, common_target_section
        )

        source_values = clutch.relative.private.delayed_carry_pair(
            source,
            semantic_prefixes[0][ell],
            {},
        )
        target_values = clutch.relative.private.delayed_carry_pair(
            target,
            semantic_prefixes[0][ell],
            {},
        )
        source_vector.append(source_values[12][KAPPA])
        target_vector.append(target_values[6][KAPPA])

        source_carry = clutch.relative.private.old.intersect_weighted_comb(
            source, clutch.R // P, P, 12, 13
        )
        target_carry = clutch.relative.private.old.intersect_weighted_comb(
            target, clutch.R // P, P, 6, 7
        )
        direct_source = (
            clutch.relative.private.old.delayed_weighted_numerator(
                source_carry, semantic_prefixes[0][ell][KAPPA]
            )
        )
        direct_target = (
            clutch.relative.private.old.delayed_weighted_numerator(
                target_carry, semantic_prefixes[0][ell][KAPPA]
            )
        )
        require(
            direct_source == source_values[12][KAPPA]
            and direct_target == target_values[6][KAPPA],
            f"direct/carry-descent disagreement at clock {ell}",
        )
        require(
            clutch.shift_weighted(source, clutch.SHIFT) == tuple(target),
            f"source/target physical carriers do not translate at clock {ell}",
        )
        support_masses.append((
            clutch.interval_mass(tuple((a, b) for a, b, _ in source)),
            clutch.interval_mass(tuple((a, b) for a, b, _ in target)),
            len(source), len(target),
        ))
        details.append((tuple(source), tuple(target),
                        tuple(source_carry), tuple(target_carry)))
    return (tuple(source_vector), tuple(target_vector),
            tuple(support_masses), rail_pairs, tuple(details))


def primitive(values, root, content):
    if not values or not content or any(value % content for value in values):
        return False, None
    inverse = pow(root, -1, P)
    normalized = tuple((value // content) * inverse % P for value in values)
    reduced = tuple((normalized[index] - normalized[-1]) % P
                    for index in range(6))
    determinant = clutch.relative.private.old.sat.multiplication_determinant_7(
        reduced
    )
    return bool(determinant), (normalized, reduced, determinant)


def main():
    module, _prefixes, _, _, rails, present, _starts = (
        clutch.relative.lift.m.core.build_carrier_data()
    )
    pair_prefixes = clutch.relative.private.build_pair_prefixes(module)
    fork = semantic.deepest_fork(module)
    semantic_prefixes = build_semantic_prefixes(module, fork)
    source_e3 = semantic.exclusive_source(module, 3)
    chosen_sections = semantic_sections(module, source_e3, 0, 3)
    (chosen_source, chosen_target, chosen_masses, _pairs,
     chosen_details) = restricted_vectors(
        module, pair_prefixes, rails, present, semantic_prefixes,
        chosen_sections, 8,
    )
    print(f"chosen_0_3_source={chosen_source}")
    print(f"chosen_0_3_target={chosen_target}")
    print(f"chosen_0_3_masses={chosen_masses}")

    q_source = clutch.relative.semantic.frac(
        clutch.relative.Z + clutch.Fraction(7 * 6715, clutch.R)
    )
    q_target = clutch.relative.semantic.frac(
        clutch.relative.Z + clutch.Fraction(7 * 6716, clutch.R)
    )
    delayed = clutch.relative.semantic.frac(clutch.R * q_source)
    prefix = semantic_prefixes[0][1][KAPPA]
    prefix_intervals = tuple(
        (left, left + length)
        for left, length in zip(prefix[0], prefix[1])
    )
    require(
        q_target == q_source + clutch.Fraction(7, clutch.R)
        and semantic.strict_member(q_source, chosen_sections[1], module.T)
        and semantic.strict_member(q_target, chosen_sections[1], module.T)
        and semantic.strict_member(delayed, prefix_intervals, module.T)
        and semantic.strict_member(delayed, fork, module.T)
        and clutch.relative.semantic.semantic_record(q_source) == (3, (1, 2))
        and clutch.relative.semantic.semantic_record(q_target) == (3, (1, 2))
        and clutch.containing_weighted_piece(
            q_source * module.T, chosen_details[1][0]
        )
        and clutch.containing_weighted_piece(
            q_target * module.T, chosen_details[1][1]
        ),
        "original adjacent semantic witness left the strengthened carrier",
    )
    print(f"strengthened_witness={q_source}->{q_target} delayed={delayed}")
    for content in (BROADCAST_CONTENT, clutch.CONTENT):
        print(f"content={content} source_primitive="
              f"{primitive(chosen_source, 12, content)}")
        print(f"content={content} target_primitive="
              f"{primitive(chosen_target, 1, content)}")

    hostile_sections = semantic_sections(module, source_e3, 0, 0)
    hostile_source, hostile_target, *_ = restricted_vectors(
        module, pair_prefixes, rails, present, semantic_prefixes,
        hostile_sections, 8,
    )
    require(
        hostile_source == hostile_target == (0,) * 7,
        "t=0 E3/full-target hostile unexpectedly survived",
    )
    print("hostile_t0_E3_full_target=(0,0,0,0,0,0,0)")

    chosen_rail_rows = []
    chosen_content = 0
    for rail_index in range(14):
        source, target, masses, _pairs, _details = restricted_vectors(
            module, pair_prefixes, rails, present, semantic_prefixes,
            chosen_sections, rail_index, equal_weight_only=True,
        )
        for value in source + target:
            chosen_content = gcd(chosen_content, value)
        source_unit, source_profile = primitive(
            source, 12, BROADCAST_CONTENT
        )
        target_unit, target_profile = primitive(
            target, 1, BROADCAST_CONTENT
        )
        chosen_rail_rows.append((
            rail_index,
            tuple(index for index, value in enumerate(source) if value),
            source == target,
            source_unit,
            target_unit,
            source_profile[2] if source_profile else None,
            target_profile[2] if target_profile else None,
        ))
    require(
        all(row[2] for row in chosen_rail_rows),
        "an equal-weight semantic rail failed coefficient equality",
    )
    print(f"chosen_0_3_equal_weight_rail_rows={tuple(chosen_rail_rows)}")
    print(f"chosen_0_3_equal_weight_content={chosen_content} "
          f"content_over_broadcast={chosen_content // BROADCAST_CONTENT} "
          f"residual_mod13={(chosen_content // BROADCAST_CONTENT) % P}")

    return

    all_rows = []
    content = 0
    rail_stats = []
    for rail_index in range(14):
        rail_positive = rail_equal = 0
        clock_support = Counter()
        for s in range(P):
            for t in range(1, P):
                sections = semantic_sections(module, source_e3, s, t)
                source, target, masses, _pairs, _details = restricted_vectors(
                    module, pair_prefixes, rails, present,
                    semantic_prefixes, sections, rail_index,
                    equal_weight_only=True,
                )
                if any(source) or any(target):
                    rail_positive += 1
                if source == target and any(source):
                    rail_equal += 1
                for ell, value in enumerate(source):
                    if value:
                        clock_support[ell] += 1
                        content = gcd(content, value)
                for value in target:
                    if value:
                        content = gcd(content, value)
                all_rows.append((rail_index, s, t, source, target))
        rail_stats.append((rail_index, rail_positive, rail_equal,
                           tuple(sorted(clock_support.items()))))
    units = []
    for rail_index, s, t, source, target in all_rows:
        source_unit, _ = primitive(source, 12, content)
        target_unit, _ = primitive(target, 1, content)
        if source == target and source_unit and target_unit:
            units.append((rail_index, s, t))
    print(f"equal_weight_global_content={content}")
    print(f"equal_weight_rail_stats={tuple(rail_stats)}")
    print(f"equal_weight_both_unit_equal={len(units)} "
          f"first={tuple(units[:30])}")


if __name__ == "__main__":
    main()
