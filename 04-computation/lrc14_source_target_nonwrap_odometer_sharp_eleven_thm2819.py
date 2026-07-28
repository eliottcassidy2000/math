#!/usr/bin/env python3
"""Exact finite referee for THM-2819.

Translate the THM-2809 source-coordinate label bank by tau=7/13^6 into
THM-2749's marked target coordinate.  Labels 1,...,12 conjugate without
wrapping to target labels 0,...,11.  The canonical target label 12 instead
pulls back to the integer lift 13, not lift 0; its residual displacement is
91/13^6.  The exact C3 gate kills that wrapped chart.  Target labels
1,...,11 retain the positive source eleven-face with identical exact masses,
so the marked target maximum is sharply eleven.
"""

import hashlib
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


DEPENDENCIES = {
    "lrc14_slope7_fixed_configuration_carry_nerve_thm2672.py":
        "83ccf3a38660a92cc990bdf304fd4ea4475339731c3e7e92ad35383ef097f361",
    "lrc14_fully_marked_root_zero_clutch_thm2749.py":
        "93b46b2701db8f72d00fa2ae131f9d9afd3200f32998959af3bb2e1fa2f56841",
    "lrc14_sharp_labelwise_eleven_label_maximum_thm2809.py":
        "d74ea1db38238dcd95a32598d1f728d2e37a08be7cc8c951c87bc3de15af6551",
}

for dependency, expected_hash in DEPENDENCIES.items():
    payload = (COMP / dependency).read_bytes().replace(b"\r\n", b"\n")
    actual_hash = hashlib.sha256(payload).hexdigest()
    require(actual_hash == expected_hash,
            f"audited dependency changed: {dependency}")


import lrc14_predecessor_carry_private_root_atlas_thm2640 as atlas
import lrc14_slope7_fixed_configuration_carry_nerve_thm2672 as fixed
import lrc14_fully_marked_root_zero_clutch_thm2749 as marked
import lrc14_sharp_labelwise_eleven_label_maximum_thm2809 as source_face


P = 13
SOURCE_CARRY = 12
TARGET_CARRY = 6
RAIL_COUNT = 14
SHIFT = 7 * atlas.T // atlas.R
WRAP_SHIFT = 13 * SHIFT
EXPECTED_PREDEEP_ROWS = (
    (0, 784_513_810_080, 728_825_300_280, 0),
    (1, 356_159_643_840, 330_902_579_700, 0),
    (2, 473_353_040_160, 439_916_713_920, 0),
    (3, 625_685_860_800, 581_278_296_060, 0),
    (4, 620_815_595_400, 578_179_396_872, 0),
    (5, 438_008_061_491, 407_334_158_691, 0),
    (6, 394_663_389_120, 366_653_853_720, 0),
)


def interval_mass(intervals):
    return sum(right - left for left, right in intervals)


def weighted_mass(pieces):
    return sum((right - left) * weight
               for left, right, weight in pieces)


def translate(intervals, shift):
    """Literal image translation, as opposed to fixed.shift_union preimage."""
    return tuple(marked.clutch.shift_union(intervals, shift))


def preimage(intervals, label):
    return tuple(fixed.shift_union(
        intervals, label * SHIFT, atlas.T
    ))


def weighted_intersection(weighted, support):
    if not weighted or not support:
        return tuple()
    return tuple(atlas.old.intersect_weighted_union(
        weighted, support
    ))


def weighted_support(weighted):
    return tuple((left, right) for left, right, *_rest in weighted)


def target_face_row(details, rail, present, prefixes, labels):
    """Direct target-coordinate delayed row for one rail and label set."""
    values = []
    post_label_support = []
    for clock in range(7):
        common = tuple(
            row for row in details[clock][1] if row[2]
        )
        for label in labels:
            common = weighted_intersection(
                common, preimage(weighted_support(rail), label)
            )
            common = weighted_intersection(
                common, preimage(present[clock, 7], label)
            )
            if not common:
                break
        post_label_support.append(bool(common))
        value = (
            atlas.delayed_carry_pair(
                common, prefixes[0][clock], {}
            )[TARGET_CARRY][1]
            if common else 0
        )
        values.append(value)
    values = tuple(values)
    return (
        tuple(clock for clock, value in enumerate(values) if value),
        sum(values),
        tuple(post_label_support),
    )


def main():
    require(
        (atlas.P, atlas.R, atlas.T, SHIFT, WRAP_SHIFT)
        == (
            13,
            4_826_809,
            297_836_897_838_480,
            431_933_040,
            5_615_129_520,
        ),
        "canonical endpoint scales changed",
    )
    require(
        atlas.R * SHIFT == 7 * atlas.T
        and atlas.R * WRAP_SHIFT == 91 * atlas.T,
        "slope-seven or odometer lift ceased to be integral",
    )

    shard = atlas.shard((0, RAIL_COUNT))
    content = shard[1]
    metadata = shard[5]
    raw_rows = shard[6]
    require(
        content == 26
        and len(metadata) == len(raw_rows) == RAIL_COUNT,
        "first-fourteen rail shard changed",
    )
    flags = fixed.build_flags(raw_rows, content)

    (
        module,
        _prefixes,
        _whole,
        _masses,
        rails,
        present,
        _starts,
    ) = atlas.core.build_carrier_data()
    source = marked.semantic.exclusive_source(module, 3)
    terminal_fork = marked.semantic.deepest_fork(module)
    prefixes = marked.build_semantic_prefixes(module, terminal_fork)
    sections = marked.semantic_sections(module, source, 0, 4)

    # The carry and configuration rows match under delta -> epsilon=delta-1.
    carry_rows = []
    for delta in range(1, P):
        epsilon = delta - 1
        source_carry = (SOURCE_CARRY + 7 * delta) % P
        target_carry = (TARGET_CARRY + 7 * epsilon) % P
        require(
            source_carry == target_carry,
            f"carry conjugacy failed at source label {delta}",
        )
        carry_rows.append((delta, epsilon, source_carry))
    carry_rows = tuple(carry_rows)

    # The physical chart identity is checked on every rail and every present
    # bank.  fixed.shift_union is the preimage under x -> x+shift; translate
    # is the literal image.  Thus T_tau U_delta=U_(delta-1).
    rail_conjugacy_checks = 0
    present_conjugacy_checks = 0
    for rail_index in range(RAIL_COUNT):
        rail_support = weighted_support(rails[rail_index][3])
        for delta in range(1, P):
            source_chart = preimage(rail_support, delta)
            target_chart = preimage(rail_support, delta - 1)
            require(
                translate(source_chart, SHIFT) == target_chart,
                f"rail chart conjugacy failed at "
                f"{(rail_index, delta)}",
            )
            rail_conjugacy_checks += 1
    for clock in range(7):
        for delta in range(1, P):
            source_chart = preimage(present[clock, 7], delta)
            target_chart = preimage(present[clock, 7], delta - 1)
            require(
                translate(source_chart, SHIFT) == target_chart,
                f"present chart conjugacy failed at {(clock, delta)}",
            )
            present_conjugacy_checks += 1

    # THM-2749's complete marked source carrier and carry cell translate
    # literally into the marked target carrier and carry cell.
    endpoint_transport_checks = 0
    marked_details = []
    for rail_index in range(RAIL_COUNT):
        source_vector, target_vector, details = (
            marked.fully_marked_vectors(
                module,
                rails,
                present,
                prefixes,
                sections,
                rail_index,
            )
        )
        require(
            source_vector == target_vector,
            f"marked endpoint vector changed on rail {rail_index}",
        )
        for clock, row in enumerate(details):
            source_carrier, target_carrier = row[0], row[1]
            source_cell, target_cell = row[2], row[3]
            require(
                marked.clutch.shift_weighted(
                    source_carrier, SHIFT
                ) == target_carrier
                and marked.clutch.shift_weighted(
                    source_cell, SHIFT
                ) == target_cell,
                f"marked endpoint transport failed at "
                f"{(rail_index, clock)}",
            )
            endpoint_transport_checks += 1
        marked_details.append(details)
    marked_details = tuple(marked_details)

    # The target label zero is the image of source label one.  Its unique
    # compatible row has h=12,kappa=1, and its future half digit is disjoint
    # from the marked h=6,kappa=1 prefix.  Since R*tau=7 is integral, the
    # future coordinate is literally unchanged by endpoint translation.
    target_zero_rows = []
    prefix_disjoint_rows = []
    ordinary_prefixes = atlas.build_pair_prefixes(module)
    for rail_index in range(RAIL_COUNT):
        source_label_one = source_face.compatible_rows(
            flags, rail_index, 1
        )
        require(
            source_label_one
            == ((
                *source_face.LABEL_ONE_CONFIG,
                6,
                1,
                0,
                (169, 182),
            ),),
            f"target-zero/source-one row changed on rail {rail_index}",
        )
        target_zero_rows.append(source_label_one[0])
    for clock in range(7):
        label_zero_prefix = source_face.prefix_intervals(
            ordinary_prefixes[1][clock][12][1]
        )
        marked_prefix = source_face.prefix_intervals(
            prefixes[0][clock][1]
        )
        meet = source_face.interval_intersection(
            label_zero_prefix, marked_prefix
        )
        require(
            not meet,
            f"target label zero met the marked prefix at clock {clock}",
        )
        prefix_disjoint_rows.append((
            clock,
            len(label_zero_prefix),
            len(marked_prefix),
            len(meet),
        ))
    target_zero_rows = tuple(target_zero_rows)
    prefix_disjoint_rows = tuple(prefix_disjoint_rows)

    # The only canonical wrap is target epsilon=12.  Pulling it back through
    # the endpoint translation gives integer source lift 13, not lift 0.
    # Its residual displacement 13*tau=91/R is nontrivial physically.
    for clock in range(7):
        target_twelve = preimage(present[clock, 7], 12)
        source_pullback = translate(target_twelve, -SHIFT)
        source_lift_thirteen = preimage(present[clock, 7], 13)
        require(
            source_pullback == source_lift_thirteen,
            f"canonical wrap did not pull back to lift 13 at clock {clock}",
        )

    require(
        module.C3 * 91 % atlas.R == 0
        and module.C3 * 91 // atlas.R == 14,
        "odometer shift stopped stabilizing the C3 phase",
    )

    # Isolate the first universal failed factor.  The wrapped present packet
    # survives both present-complement gates with positive mass.  Its C3-safe
    # factor is invariant under 91/R, however, while the marked source band
    # (169,181)/182 is strictly C3-danger.  The deep intersection is empty
    # before any rail, delayed prefix, or other target label is imposed.
    predeep_rows = []
    for clock in range(7):
        base_present = tuple(present[clock, 7])
        wrapped_present = preimage(base_present, 13)
        source_safe = marked.clutch.complement(base_present)
        translated_target_safe = marked.clutch.shift_union(
            source_safe, -SHIFT
        )
        after_source_safe = marked.clutch.intersect_sorted(
            wrapped_present, source_safe
        )
        after_both_safe = marked.clutch.intersect_sorted(
            after_source_safe, translated_target_safe
        )
        deep_weighted = atlas.old.intersect_weighted_comb(
            tuple((left, right, 1)
                  for left, right in after_both_safe),
            module.C3,
            182,
            169,
            181,
        )
        row = (
            clock,
            interval_mass(after_source_safe),
            interval_mass(after_both_safe),
            weighted_mass(deep_weighted),
        )
        require(
            row == EXPECTED_PREDEEP_ROWS[clock],
            f"wrapped factor-stage row changed at clock {clock}",
        )
        require(
            row[1] > 0 and row[2] > 0 and row[3] == 0,
            f"wrap did not first die at the C3 gate at clock {clock}",
        )
        predeep_rows.append(row)
    predeep_rows = tuple(predeep_rows)

    # Direct target-coordinate reconstruction: labels 1,...,11 retain the
    # source positive eleven-face with identical exact clock supports/masses.
    # Adding the wrapped target label 12 kills it before delayed integration.
    target_eleven_rows = []
    target_twelve_omit_zero_rows = []
    for rail_index in range(RAIL_COUNT):
        eleven = target_face_row(
            marked_details[rail_index],
            rails[rail_index][3],
            present,
            prefixes,
            tuple(range(1, 12)),
        )
        twelve = target_face_row(
            marked_details[rail_index],
            rails[rail_index][3],
            present,
            prefixes,
            tuple(range(1, 13)),
        )
        expected = source_face.EXPECTED_ELEVEN_FACE_ROWS[rail_index]
        require(
            eleven[:2] == expected,
            f"positive target eleven-face changed on rail {rail_index}",
        )
        require(
            twelve[:2] == (tuple(), 0)
            and not any(twelve[2]),
            f"target twelve-face omitting zero survived on rail "
            f"{rail_index}",
        )
        target_eleven_rows.append(eleven[:2])
        target_twelve_omit_zero_rows.append(twelve[:2])
    target_eleven_rows = tuple(target_eleven_rows)
    target_twelve_omit_zero_rows = tuple(
        target_twelve_omit_zero_rows
    )

    require(
        len(set(target_zero_rows)) == 1
        and all(row[1] > 0 for row in target_eleven_rows)
        and set(target_twelve_omit_zero_rows) == {(tuple(), 0)},
        "target face classification lost sharpness",
    )

    print("THM2819 SOURCE-TARGET NONWRAP/ODOMETER SHARP ELEVEN EXACT REFEREE")
    print(f"first14_rail_metadata={tuple(metadata)}")
    print(f"tau_numerator_over_T={SHIFT}; R*tau=7")
    print(f"wrap_numerator_over_T={WRAP_SHIFT}; R*wrap=91")
    print(f"nonwrap_carry_rows={carry_rows}")
    print(
        f"chart_conjugacy_checks=rails:{rail_conjugacy_checks},"
        f"present:{present_conjugacy_checks}"
    )
    print(f"marked_endpoint_transport_checks={endpoint_transport_checks}")
    print(f"unique_target0_source1_row={target_zero_rows[0]}")
    print(f"target0_prefix_disjoint_rows={prefix_disjoint_rows}")
    print("canonical_target12_pullback=integer_source_lift13")
    print("odometer_cocycle=91/R=7/13^5")
    print("C3*91/R=14; wrapped_present_C3_safe_is_invariant")
    print(f"wrapped_factor_stage_rows={predeep_rows}")
    print("first_universal_wrap_failure=strict_C3_safe_vs_"
          "marked_danger_(169,181)/182")
    print("positive_target_eleven_labels=(1,2,3,4,5,6,7,8,9,10,11)")
    print(f"positive_target_eleven_rows={target_eleven_rows}")
    print(
        "target_twelve_face_census=12_containing_label0_die_at_"
        "future_digit+1_omitting_label0_dies_at_wrapped_C3_gate"
    )
    print("sharp_target_maximum_label_cardinality=11; simplex_dimension=10")
    print("scope=first14 THM2749 marked endpoint rails; no outside-rail, "
          "row, or LRC14 exclusion")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
