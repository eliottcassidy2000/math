#!/usr/bin/env python3
"""Exact root-zero chart-overlap clutch on the relative-present E3 fibre.

The relative-present carrier reaches root zero if a scalar THM-2657 lift is
forced to retain its source half-tooth edge.  Adjacent translated danger
teeth overlap, however.  This companion checks that the canonical
residue-7 -> residue-8 lift lands in the open overlap H^R_0 cap H^L_1,
recharts there, and retains an exact recomputed primitive-unit vector after
restricting to the overlap carrier.

This is a partial chart isomorphism on one canonical fibre.  It is not a
global edge-preserving action, endpoint current, row exclusion, or LRC(14)
proof.  All arithmetic is exact; there are no Python assertions.
"""

from __future__ import annotations

from bisect import bisect_right
from fractions import Fraction
from pathlib import Path
import hashlib
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

DEPENDENCIES = {
    "lrc14_relative_present_semantic_lift_probe_20260728.py":
        "f16754bd38ae0dfa0d7d91cc404b4447dbf359635101aa7b4223363f8064352f",
    "lrc14_central_half_odometer_full_local_cycle_thm2698.py":
        "45cc393a856c00342fdf84875a0bc5a6d4c3df196ab35bb9ac2aad3cfc966c25",
    "lrc14_full_physical_lift_fibre_thm2707.py":
        "f05a07b2fb22cb5b39ed7d14e66d26154ecc50fc214861dc6576c3bcfaed2412",
    "lrc14_half_c221_semantic_source_fibre_census_20260728.py":
        "0fbeb041007fea1b9e14f0ff6e82fc97ebf724b26c2c10ef85732b4c994b94cd",
    "lrc14_predecessor_carry_private_root_atlas_thm2640.py":
        "a28b03a5903256c1c1c294ea5af389c7991fc0a5ad6908f0f25a5b0cc6e71abf",
    "lrc14_replica_dichotomy_typed_row_opus_20260727.py":
        "6ba64a68a9fd008d2e06949b1f1cf75012f1f4e734f75f55ce0af58ae20ad7b9",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


for dependency, expected_hash in DEPENDENCIES.items():
    actual_hash = hashlib.sha256((COMP / dependency).read_bytes()).hexdigest()
    require(actual_hash == expected_hash,
            f"audited dependency changed: {dependency}")

import lrc14_relative_present_semantic_lift_probe_20260728 as relative


P = relative.P
R = relative.R
T = relative.lift.m.T
SHIFT = 7 * T // R
CONTENT = relative.CANONICAL_CONTENT
SOURCE_STATE = (7, 1, 12, 12)  # residue, script-edge, carry, root: right
TARGET_STATE = (8, 0, 6, 1)    # residue, script-edge, carry, root: left


def merge_intervals(intervals):
    """Merge a finite interval union, including touching endpoints."""
    out = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if out and left <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return tuple(out)


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


def interval_mass(intervals):
    return sum(right - left for left, right in intervals)


def cyclic_interval(left, right, period):
    """One open cyclic interval, represented by half-open support pieces."""
    length = right - left
    require(0 < length < period, "invalid cyclic interval length")
    start = left % period
    end = start + length
    if end <= period:
        return ((start, end),)
    return ((0, end - period), (start, period))


def shift_union(intervals, shift, period=T):
    out = []
    for left, right in intervals:
        length = right - left
        start = (left + shift) % period
        end = start + length
        if end <= period:
            out.append((start, end))
        else:
            out.extend(((start, period), (0, end - period)))
    return merge_intervals(out)


def shift_weighted(pieces, shift, period=T):
    out = []
    for left, right, weight in pieces:
        length = right - left
        start = (left + shift) % period
        end = start + length
        if end <= period:
            out.append((start, end, weight))
        else:
            out.extend(((start, period, weight),
                        (0, end - period, weight)))
    return tuple(sorted(out))


def complement(intervals, period=T):
    intervals = merge_intervals(intervals)
    out = []
    cursor = 0
    for left, right in intervals:
        if cursor < left:
            out.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < period:
        out.append((cursor, period))
    return tuple(out)


def intersect_weighted_profiles(source, target):
    """Common support with both step-profile weights retained."""
    out = []
    i = j = 0
    while i < len(source) and j < len(target):
        left = max(source[i][0], target[j][0])
        right = min(source[i][1], target[j][1])
        if left < right:
            out.append((left, right, source[i][2], target[j][2]))
        if source[i][1] <= target[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def restrict_to_relative_overlap(module, present, rail_pairs, ell):
    """Return source pieces and forward-translated target pieces.

    The input lives in the source coordinate.  It already lies on both the
    source rail and the pullback of the same target rail.  We impose both
    present complements and the exact deep overlap

        H^R_12 cap T_tau^(-1) H^L_1 = (169,181)/182.
    """
    source = tuple((a, b, source_weight)
                   for a, b, source_weight, _target_weight in rail_pairs)
    target = tuple((a, b, target_weight)
                   for a, b, _source_weight, target_weight in rail_pairs)
    present_complement = complement(present[ell, 7])
    target_complement_pullback = shift_union(
        present_complement, -SHIFT
    )
    old = relative.private.old
    source = old.intersect_weighted_union(source, present_complement)
    source = old.intersect_weighted_union(
        source, target_complement_pullback
    )
    source = old.intersect_weighted_comb(
        source, module.C3, 182, 169, 181
    )
    target = old.intersect_weighted_union(target, present_complement)
    target = old.intersect_weighted_union(
        target, target_complement_pullback
    )
    target = old.intersect_weighted_comb(
        target, module.C3, 182, 169, 181
    )
    return tuple(source), shift_weighted(target, SHIFT)


def overlap_vectors(module, pair_prefixes, rails, present, rail_index,
                    equal_weight_only=False):
    target_pullback = shift_weighted(rails[rail_index][3], -SHIFT)
    rail_pairs = intersect_weighted_profiles(
        rails[rail_index][3], target_pullback
    )
    if equal_weight_only:
        rail_pairs = tuple(
            row for row in rail_pairs if row[2] == row[3]
        )
    require(rail_pairs, "empty source/target rail overlap")

    source_vector = []
    target_vector = []
    details = []
    for ell in range(7):
        source, target = restrict_to_relative_overlap(
            module, present, rail_pairs, ell
        )
        source_values = relative.private.delayed_carry_pair(
            source, pair_prefixes[0][ell][6], {}
        )
        target_values = relative.private.delayed_carry_pair(
            target, pair_prefixes[0][ell][6], {}
        )
        source_vector.append(source_values[12][1])
        target_vector.append(target_values[6][1])
        details.append((source, target))
    return (tuple(source_vector), tuple(target_vector),
            rail_pairs, tuple(details))


def normalized_profile(vector, root):
    inverse = pow(root, -1, P)
    values = tuple((value // CONTENT) * inverse % P for value in vector)
    reduced = tuple((values[index] - values[-1]) % P
                    for index in range(6))
    return values, reduced


def containing_weighted_piece(value, pieces):
    return tuple(
        (left, right, weight)
        for left, right, weight in pieces
        if left < value < right
    )


def weighted_slack(value, pieces):
    containing = containing_weighted_piece(value, pieces)
    require(containing, "point left the weighted carrier")
    return min(
        min(value - left, right - value)
        for left, right, _weight in containing
    ) / T


def main():
    require((P, R, T, SHIFT) == (
        13, 4826809, 297836897838480, 431933040
    ), "canonical scales changed")

    module, _prefixes, _, _, rails, present, starts = (
        relative.lift.m.core.build_carrier_data()
    )
    pair_prefixes = relative.private.build_pair_prefixes(module)
    require(module.C3 == 742586 and len(rails) == 162,
            "canonical deep blocker or rail bank changed")

    # Every right half-tooth overlaps the next left half-tooth in a strip
    # of width 12/182.  This is a cover overlap, not an edge identity.
    adjacent_overlaps = []
    for root in range(P):
        right = cyclic_interval(14 * root, 14 * root + 13, 182)
        next_left = cyclic_interval(
            14 * (root + 1) - 13, 14 * (root + 1), 182
        )
        expected = cyclic_interval(14 * root + 1, 14 * root + 13, 182)
        overlap = intersect_sorted(right, next_left)
        require(overlap == expected and interval_mass(overlap) == 12,
                f"adjacent chart overlap changed at root {root}")
        adjacent_overlaps.append(overlap)

    source_overlap = ((169, 181),)
    target_overlap = ((1, 13),)
    require(
        intersect_sorted(
            cyclic_interval(168, 181, 182),
            shift_union(cyclic_interval(1, 14, 182), -14, 182),
        ) == source_overlap
        and shift_union(source_overlap, 14, 182) == target_overlap
        and shift_union(target_overlap, -14, 182) == source_overlap,
        "forward/reverse root-zero overlap law changed",
    )

    delayed = relative.semantic.frac(R * relative.Z)
    eta = 2 * delayed - 1
    require(delayed == Fraction(11195237, 20792408)
            and eta == Fraction(799033, 10396204)
            and eta - Fraction(1, 14)
            == Fraction(56447, 14 * 742586),
            "overlap-threshold parameter changed")

    n_source, n_target = 6715, 6716
    q_source = relative.semantic.frac(
        relative.Z + Fraction(7 * n_source, R)
    )
    q_target = relative.semantic.frac(
        relative.Z + Fraction(7 * n_target, R)
    )
    require(
        q_source == Fraction(47850889647341, 100360982066072)
        and q_target == Fraction(47851035194197, 100360982066072)
        and q_target == q_source + Fraction(7, R),
        "adjacent semantic lift witness changed",
    )
    require(
        relative.semantic.frac(R * q_source)
        == relative.semantic.frac(R * q_target) == delayed,
        "translation did not preserve the delayed coordinate",
    )

    deep_source = relative.semantic.frac(module.C3 * q_source) * 182
    deep_target = relative.semantic.frac(module.C3 * q_target) * 182
    require(
        deep_source == 168 + 14 * eta
        == Fraction(125553481, 742586)
        and deep_target == 14 * eta
        == Fraction(799033, 742586)
        and relative.semantic.frac(
            module.C3 * Fraction(7, R)
        ) * 182 == 14,
        "deep-coordinate transport changed",
    )
    overlap_grid_margin = min(
        deep_source - 169, 181 - deep_source,
        deep_target - 1, 13 - deep_target,
    )
    overlap_q_margin = overlap_grid_margin / (182 * module.C3)
    require(
        overlap_grid_margin == Fraction(56447, 742586)
        and overlap_q_margin == Fraction(56447, 100360982066072)
        and overlap_q_margin == 56447 * relative.Q_RADIUS,
        "root-zero chart-overlap margin changed",
    )

    require(
        relative.semantic.semantic_record(q_source) == (3, (1, 2))
        and relative.semantic.semantic_record(q_target) == (3, (1, 2)),
        "semantic E3 fork witness changed",
    )
    source_labels = relative.full_target_labels(q_source)
    target_labels = relative.full_target_labels(q_target)
    expected_source_labels = (
        (0, 1, 2, 3, 8, 9, 10, 11, 12),
        (3, 4, 5, 6, 7, 8, 9, 10, 11),
    )
    expected_target_labels = (
        (0, 1, 2, 3, 8, 9, 10, 11, 12),
        (3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
    )
    common_s = tuple(sorted(set(source_labels[0]) & set(target_labels[0])))
    common_t = tuple(sorted(set(source_labels[1]) & set(target_labels[1])))
    require(
        source_labels == expected_source_labels
        and target_labels == expected_target_labels
        and common_s == expected_source_labels[0]
        and common_t == expected_source_labels[1],
        "common full target labels changed",
    )
    for point, expected in (
        (q_source, expected_source_labels),
        (q_target, expected_target_labels),
    ):
        require(
            relative.full_target_labels(point - relative.Q_RADIUS)
            == expected
            and relative.full_target_labels(point + relative.Q_RADIUS)
            == expected,
            "full target labels change on the common open cylinder",
        )

    require(
        rails[8][:3] == (1, 4, 12)
        and relative.half.shallow(q_source)
        == relative.half.shallow(q_target) == 1
        and relative.half.owner(q_source)
        == relative.half.owner(q_target) == 4,
        "rail or reverse-clock witness changed",
    )
    broad_segment = (
        141869054887560, 142120818960120, 27581135604
    )
    require(
        containing_weighted_piece(q_source * T, rails[8][3])
        == (broad_segment,)
        and containing_weighted_piece(q_target * T, rails[8][3])
        == (broad_segment,),
        "semantic pair left its common rail segment",
    )
    require(
        not relative.half.strict_interval_member(
            q_source * T, present[1, 7], starts[1, 7]
        )
        and not relative.half.strict_interval_member(
            q_target * T, present[1, 7], starts[1, 7]
        ),
        "semantic pair left the relative present complement",
    )
    active_clocks = tuple(
        ell for ell in range(7)
        if relative.half.strict_interval_member(
            delayed,
            relative.half.prefix_intervals(
                pair_prefixes[0][ell][6][1]
            ),
        )
    )
    require(active_clocks == (1, 2, 3, 4, 5, 6),
            "delayed clock mask changed")

    natural_equal_rails = []
    equal_weight_unit_rails = []
    equal_weight_amplitudes = []
    rail8_data = None
    for rail_index in relative.SOURCE_ONE_RAILS:
        source_vector, target_vector, rail_pairs, details = overlap_vectors(
            module, pair_prefixes, rails, present, rail_index
        )
        if source_vector == target_vector:
            natural_equal_rails.append(rail_index)

        equal_source, equal_target, equal_pairs, equal_details = (
            overlap_vectors(
                module, pair_prefixes, rails, present, rail_index,
                equal_weight_only=True,
            )
        )
        require(
            equal_source == equal_target
            and equal_source[0] == 0
            and len(set(equal_source[1:])) == 1
            and equal_source[1] > 0
            and relative.private.is_unit(
                equal_source, SOURCE_STATE[3], CONTENT
            )
            and relative.private.is_unit(
                equal_target, TARGET_STATE[3], CONTENT
            ),
            f"equal-weight overlap clutch failed on rail {rail_index}",
        )
        equal_weight_unit_rails.append(rail_index)
        equal_weight_amplitudes.append(equal_source[1])
        if rail_index == 8:
            rail8_data = (
                source_vector, target_vector, rail_pairs, details,
                equal_source, equal_target, equal_pairs, equal_details,
            )

    require(
        tuple(natural_equal_rails)
        == (0, 4, 5, 6, 7, 8, 9, 10, 11, 12)
        and tuple(equal_weight_unit_rails) == tuple(range(14)),
        "uniform source-one rail overlap census changed",
    )
    require(rail8_data is not None, "rail-eight overlap data vanished")
    (rail8_source, rail8_target, rail8_pairs, rail8_details,
     rail8_equal_source, rail8_equal_target, _equal_pairs,
     _equal_details) = rail8_data
    amplitude = 5359949020444386606638400
    expected_vector = (0,) + (amplitude,) * 6
    require(
        rail8_source == rail8_target == expected_vector
        and rail8_equal_source == rail8_equal_target == expected_vector
        and all(source_weight == target_weight
                for _left, _right, source_weight, target_weight
                in rail8_pairs)
        and tuple(sorted({row[2] for row in rail8_pairs}))
        == (27580222516, 27581135604, 27582102210, 27582558770),
        "uncut rail-eight coefficient equality changed",
    )
    require(
        all(value % CONTENT == 0 for value in expected_vector)
        and (amplitude // CONTENT) % P == 11
        and relative.private.is_unit(expected_vector, 12, CONTENT)
        and relative.private.is_unit(expected_vector, 1, CONTENT),
        "overlap vector lost primitive-unit typing",
    )
    source_profile = normalized_profile(expected_vector, 12)
    target_profile = normalized_profile(expected_vector, 1)
    require(
        source_profile[1] == (11, 0, 0, 0, 0, 0)
        and target_profile[1] == (2, 0, 0, 0, 0, 0),
        "root-normalized overlap profiles changed",
    )
    clock_remainder = tuple(
        expected_vector[index] - expected_vector[-1]
        for index in range(6)
    )
    require(
        clock_remainder == (-amplitude, 0, 0, 0, 0, 0),
        "Phi_7 clock-character remainder changed",
    )

    source_carrier, target_carrier = rail8_details[1]
    expected_source_piece = (
        142004992589460, 142005019034340, 27581135604
    )
    expected_target_piece = (
        142005424522500, 142005450967380, 27581135604
    )
    require(
        containing_weighted_piece(q_source * T, source_carrier)
        == (expected_source_piece,)
        and containing_weighted_piece(q_target * T, target_carrier)
        == (expected_target_piece,)
        and weighted_slack(q_source * T, source_carrier)
        == weighted_slack(q_target * T, target_carrier)
        == overlap_q_margin,
        "semantic pair left the exact coefficient-overlap carrier",
    )
    local_rows = (
        {
            "q": q_source, "rail": 8, "shallow": 1, "owner": 4,
            "residue": 7, "edge": 1, "carry": 12, "root": 12,
        },
        {
            "q": q_target, "rail": 8, "shallow": 1, "owner": 4,
            "residue": 8, "edge": 0, "carry": 6, "root": 1,
        },
    )
    local_radii = tuple(
        relative.local_relative_radius(
            module, pair_prefixes, rails, present, starts, row
        ) for row in local_rows
    )
    require(local_radii == (relative.Q_RADIUS, relative.Q_RADIUS),
            "semantic pair lost its inherited open cylinder")

    # Hostile controls.  The old full half-tooth coefficient rows do not
    # intertwine; the repair depends on restricting to the chart overlap.
    full_bank = relative.coefficient_rows(
        module, pair_prefixes, rails, present, starts
    )
    full_source = next(
        row["relative"] for row in full_bank
        if row["rail"] == 8 and row["residue"] == 7
    )
    full_target = next(
        row["relative"] for row in full_bank
        if row["rail"] == 8 and row["residue"] == 8
    )
    require(
        full_source
        == (
            18516840472140598620032000,
            27674381326591285594468160,
            25848100717515875123124800,
            25848100717515875123124800,
            23249320129850353269130400,
            21941496594836135191274400,
            24109770976706672296239200,
        )
        and full_target == (0,) + (5377073458209384646915200,) * 6
        and full_source != full_target,
        "full-half hostile control unexpectedly intertwined",
    )

    right_only = Fraction(1, 2)
    left_only = Fraction(27, 2)
    require(
        right_only > 0 and right_only < 13
        and not (1 < right_only < 14)
        and 1 < left_only < 14
        and not (0 < left_only < 13),
        "outside-overlap membership control changed",
    )
    old_flank = (Fraction(181), Fraction(181) + Fraction(12, 371293))
    require(old_flank[0] >= source_overlap[0][1],
            "THM-2672 one-edge flank entered the new overlap")
    require(("right", 0) != ("left", 1),
            "chart relettering became an edge-label identity")

    print("LRC14 ROOT-ZERO OVERLAP CLUTCH AUDIT")
    print("status=FINITE-EXACT partial physical chart isomorphism with recomputed unit")
    print(f"R={R} T={T} tau_grid={SHIFT} deep_shift=14")
    print("chart_convention=script edge1:right epsilon0; edge0:left epsilon1")
    print("adjacent_chart_law=H^R_r intersect H^L_(r+1)=(14r+1,14r+13)/182 for all 13 roots")
    print("forward_reverse_overlap=H^R_12 intersect T^-1(H^L_1)=(169,181)/182 <-> (1,13)/182=H^R_0 intersect H^L_1")
    print(f"delayed={delayed} eta={eta} eta_minus_1/14={eta - Fraction(1,14)}")
    print(f"witness=n{n_source}->n{n_target} q={q_source} qprime={q_target} tau=7/{R}")
    print(f"deep_coordinates=({deep_source},{deep_target}) overlap_grid_margin={overlap_grid_margin}")
    print(f"overlap_q_margin={overlap_q_margin} inherited_radius={relative.Q_RADIUS} ratio=56447")
    print("witness_semantics=E3->D6 Q_(3,{1,2}) at both endpoints")
    print(f"rail8_metadata={rails[8][:3]} common_segment={broad_segment}")
    print(f"reverse_clock=(1,4) delayed_active_clocks={active_clocks}")
    print(f"source_full_target_labels=s{source_labels[0]} t{source_labels[1]}")
    print(f"target_full_target_labels=s{target_labels[0]} t{target_labels[1]}")
    print(f"common_full_target_labels=s{common_s} t{common_t} chosen=(0,3)")
    print(f"local_relative_radii={local_radii} coefficient_carrier_slack={overlap_q_margin}")
    print(f"natural_uncut_equal_rails={tuple(natural_equal_rails)} count={len(natural_equal_rails)}/14")
    print(f"equal_weight_restricted_unit_rails={tuple(equal_weight_unit_rails)} count=14/14")
    print(f"equal_weight_amplitude_range=({min(equal_weight_amplitudes)},{max(equal_weight_amplitudes)})")
    print(f"rail8_overlap_vector_source={expected_vector}")
    print(f"rail8_overlap_vector_target={expected_vector}")
    print(f"rail8_profiles_mod13=source_root12:{source_profile[1]} target_root1:{target_profile[1]}")
    print(f"clock_character_law=P mod Phi_7={clock_remainder}; P(zeta_7^k)=-{amplitude} for k=1..6")
    print(f"full_half_control=UNEQUAL source_support={sum(bool(v) for v in full_source)} target_support={sum(bool(v) for v in full_target)}")
    print("outside_overlap_control=R0\\L1=(0,1) and L1\\R0=(13,14) have no opposite-edge rechart")
    print(f"old_flank_control=THM2672 strip={old_flank} lies beyond overlap endpoint181")
    print("edge_label_control=physical point fixed but (right,root0)->(left,root1); not edge-preserving")
    print("mechanism=translate by 7/R; restrict to both rail/present-complement supports and the open tooth overlap; rechart; recompute before determinant")
    print("next_test=retain one common U_(s,t), e.g. (0,3), and D^-6 Q_(3,{1,2}) inside the overlap integral, then recompute equality and both units")
    print("SCOPE: all 13 geometric adjacent overlaps; natural coefficient equality on 10/14 source-one rails; equal-weight restricted units on 14/14; one explicit semantic rail8 pair; no global action, endpoint amplitude, row exclusion, or LRC14")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
