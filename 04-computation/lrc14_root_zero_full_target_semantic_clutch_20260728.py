#!/usr/bin/env python3
"""Exact two-clock full-target refinement of the rail-8 root-zero clutch.

This exact audit inserts the common lawful full-target label (s,t)=(0,3)
and the genuine E3 -> D^6 -> Q_(3,{1,2}) semantic cospan *inside* the
THM-2744 overlap integral.  The physical cut is made symmetrically: in the
source coordinate it retains both F(x) and F(x+7/13^6).  The delayed fork is
inserted by intersecting the y={13^6 x} prefix with Q, avoiding a materialized
13^6-fold pullback.

The script also exhausts all 81 full-target labels stable on both displayed
source/target open cylinders and all 49 pairs of physical-present and delayed
coefficient clocks.  It proves an object-level rank-one separation of those
two clocks.  It is a canonical-fibre computation, not a global clutch action,
endpoint current, row exclusion, or LRC(14) proof.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import gcd
from pathlib import Path
import hashlib
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

DEPENDENCIES = {
    "lrc14_root_zero_overlap_clutch_20260728.py":
        "e10fa7c9a5a238461ef422ea314dc334f7e65ec1787cf65d4e4bea12b96aefb8",
    "lrc14_two_target_present_semantic_attachment_probe_20260728.py":
        "062b352f4db12a5f01822b293cdbb10629632dacc5fa27b406d8dd321e550709",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path):
    """Return the repository's declared LF-normalized evidence image."""
    return path.read_bytes().replace(bytes((13, 10)), bytes((10,)))


for dependency, expected_hash in DEPENDENCIES.items():
    actual_hash = hashlib.sha256(lf_bytes(COMP / dependency)).hexdigest()
    require(actual_hash == expected_hash,
            f"audited dependency changed: {dependency}")

import lrc14_root_zero_overlap_clutch_20260728 as overlap
import lrc14_relative_present_semantic_lift_probe_20260728 as relative
import lrc14_two_target_present_semantic_attachment_probe_20260728 as target


P = relative.P
R = relative.R
T = relative.lift.m.T
SHIFT = overlap.SHIFT
CONTENT = relative.CANONICAL_CONTENT


def intervals_from_prefix(prefix):
    starts, lengths, _cumulative = prefix
    return tuple((left, left + length)
                 for left, length in zip(starts, lengths))


def q_restricted_pair_prefixes(module, pair_prefixes, fork):
    result = []
    for ell in range(7):
        pair = []
        for prefix in pair_prefixes[0][ell][6]:
            intervals = intervals_from_prefix(prefix)
            restricted = target.intersect_sorted(intervals, fork)
            pair.append(module.make_prefix(restricted))
        result.append(tuple(pair))
    return tuple(result)


def common_physical_cut(pieces, full_section, direction):
    """Impose F at this endpoint and at its adjacent translated endpoint.

    direction=-1 is the source coordinate, where the target condition is
    T_tau^(-1)F.  direction=+1 is the target coordinate, where the source
    condition is T_tau F.
    """
    cut = relative.private.old.intersect_weighted_union(
        pieces, full_section
    )
    return tuple(relative.private.old.intersect_weighted_union(
        cut, overlap.shift_union(full_section, direction * SHIFT)
    ))


def merge_weighted(pieces):
    """Canonicalize a disjoint weighted union, merging equal-weight joins."""
    result = []
    for left, right, weight in sorted(pieces):
        if left >= right:
            continue
        if result and left < result[-1][1]:
            require(
                weight == result[-1][2],
                "overlapping pieces with different weights",
            )
            result[-1] = (
                result[-1][0], max(result[-1][1], right), weight
            )
        elif (result and left == result[-1][1]
              and weight == result[-1][2]):
            result[-1] = (result[-1][0], right, weight)
        else:
            result.append((left, right, weight))
    return tuple(result)


def subtract_weighted(pieces, removed):
    """Remove the support of a weighted subprofile from another profile."""
    removed_support = overlap.merge_intervals(
        (left, right) for left, right, weight in removed if weight
    )
    return merge_weighted(relative.private.old.intersect_weighted_union(
        pieces, overlap.complement(removed_support)
    ))


def carry_cell_restriction(pieces, carry):
    """Independent explicit predecessor-carry cut on the canonical grid.

    With U=T/R and S=R/13, floor(13{Sx})=c exactly on grid cells
    [nU,(n+1)U] with n=c mod13.
    """
    require(T % R == 0, "carry cells do not resolve on the canonical grid")
    unit = T // R
    out = []
    for left, right, weight in pieces:
        require(isinstance(left, int) and isinstance(right, int),
                "explicit carry replay received a nonintegral endpoint")
        index = left // unit
        while index * unit < right:
            a = max(left, index * unit)
            b = min(right, (index + 1) * unit)
            if a < b and index % P == carry:
                out.append((a, b, weight))
            index += 1
    return tuple(out)


def valuation(value, prime):
    result = 0
    while value and value % prime == 0:
        value //= prime
        result += 1
    return result


def multiply_polynomials(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            result[i + j] += first * second
    return tuple(result)


def cyclic_convolution(left, right, modulus):
    result = [0] * modulus
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            result[(i + j) % modulus] += first * second
    return tuple(result)


def multiply_phi_seven(left, right):
    """Multiply six power-basis coordinates modulo (13, Phi_7)."""
    require(len(left) == len(right) == 6, "invalid Phi_7 product")
    product = [0] * 11
    for first_index, first in enumerate(left):
        for second_index, second in enumerate(right):
            product[first_index + second_index] += first * second
    for degree in range(10, 5, -1):
        coefficient = product[degree] % P
        for offset in range(6):
            product[degree - 6 + offset] -= coefficient
        product[degree] = 0
    return tuple(value % P for value in product[:6])


def safe_factor_slack(value, speed, centre):
    """Physical-q distance to the nearest shifted width-1/14 boundary."""
    phase = relative.semantic.centered(speed * value - centre)
    radius = Fraction(1, 14)
    return min(abs(phase - radius), abs(phase + radius)) / speed


def main():
    require((P, R, T, SHIFT, CONTENT) == (
        13, 4826809, 297836897838480, 431933040, 26
    ), "canonical scales changed")

    module, _prefixes, _, _, rails, present, _starts = (
        relative.lift.m.core.build_carrier_data()
    )
    pair_prefixes = relative.private.build_pair_prefixes(module)
    raw_source, raw_target, rail_pairs, overlap_details = (
        overlap.overlap_vectors(
            module, pair_prefixes, rails, present, rail_index=8
        )
    )
    require(raw_source == raw_target
            and all(a == b for _l, _r, a, b in rail_pairs),
            "rail-8 raw overlap stopped being translation symmetric")

    full_module = target.load_present_module()
    require(full_module.T == module.T, "target and overlap grids differ")
    e3 = target.exclusive_source(full_module, 3)
    fork = target.deepest_fork(full_module)
    clock_comb = tuple(
        full_module.make_comb(
            full_module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        for ell in range(7)
    )
    require(
        all(not target.intersect_sorted(e3, present[ell, 7])
            for ell in range(7)),
        "E3 stopped making every inherited present complement automatic",
    )
    fork_pair_prefixes = q_restricted_pair_prefixes(
        full_module, pair_prefixes, fork
    )
    require(
        all(prefix[2][-1] == 0 for prefix in fork_pair_prefixes[0])
        and all(fork_pair_prefixes[ell] == fork_pair_prefixes[1]
                for ell in range(1, 7))
        and
        tuple(prefix[2][-1]
              for pair in fork_pair_prefixes[1:]
              for prefix in pair)
        == (206986279500,) * 12,
        "terminal-fork delayed-clock object law changed",
    )
    delayed_window = (0, 1, 1, 1, 1, 1, 1)

    q_source = Fraction(47850889647341, 100360982066072)
    q_target = Fraction(47851035194197, 100360982066072)
    source_labels = relative.full_target_labels(q_source)
    target_labels = relative.full_target_labels(q_target)
    common_s = tuple(sorted(set(source_labels[0]) & set(target_labels[0])))
    common_t = tuple(sorted(set(source_labels[1]) & set(target_labels[1])))
    require(
        common_s == (0, 1, 2, 3, 8, 9, 10, 11, 12)
        and common_t == (3, 4, 5, 6, 7, 8, 9, 10, 11)
        and (0 in common_s and 3 in common_t),
        "common full-target bank changed",
    )
    label_factor_slacks = []
    for point in (q_source, q_target):
        for s in common_s:
            label_factor_slacks.extend((
                safe_factor_slack(
                    point, relative.semantic.UNITS[0], -Fraction(s, P)
                ),
                safe_factor_slack(
                    point, relative.semantic.BLOCKERS[1], Fraction(s, P)
                ),
            ))
        for t in common_t:
            label_factor_slacks.extend((
                safe_factor_slack(
                    point, relative.semantic.UNITS[1], -Fraction(t, P)
                ),
                safe_factor_slack(
                    point, relative.semantic.BLOCKERS[2], Fraction(t, P)
                ),
            ))
    common_label_margin = min(label_factor_slacks)
    require(
        common_label_margin
        == Fraction(1541619, 100360982066072)
        == 1541619 * relative.Q_RADIUS
        and common_label_margin > relative.Q_RADIUS,
        "common full-target bank lost whole-cylinder stability",
    )

    refined_source = []
    refined_target = []
    no_q_source = []
    no_q_target = []
    one_sided_source = []
    one_sided_target = []
    e3_clock_source = []
    e3_clock_target = []
    q_only_source = []
    q_only_target = []
    direct_source = []
    direct_target = []
    piece_counts = []
    chosen_physical_sections = []
    chosen_common_source_carriers = []
    chosen_common_target_carriers = []
    chosen_one_sided_source_carriers = []
    chosen_one_sided_target_carriers = []

    for ell, (source_pieces, target_pieces) in enumerate(overlap_details):
        full_section = target.source_present_section(
            full_module, e3, ell, 0, 3, clock_comb
        )
        chosen_physical_sections.append(tuple(full_section))
        source_cut = common_physical_cut(
            source_pieces, full_section, -1
        )
        target_cut = common_physical_cut(
            target_pieces, full_section, +1
        )
        chosen_common_source_carriers.append(source_cut)
        chosen_common_target_carriers.append(target_cut)
        require(
            overlap.shift_weighted(source_cut, SHIFT) == target_cut,
            f"symmetrized carrier failed translation at ell={ell}",
        )
        piece_counts.append((len(source_cut), len(target_cut)))

        q_pair = fork_pair_prefixes[ell]
        q_only_source.append(relative.private.delayed_carry_pair(
            source_pieces, q_pair, {}
        )[12][1])
        q_only_target.append(relative.private.delayed_carry_pair(
            target_pieces, q_pair, {}
        )[6][1])
        source_values = relative.private.delayed_carry_pair(
            source_cut, q_pair, {}
        )
        target_values = relative.private.delayed_carry_pair(
            target_cut, q_pair, {}
        )
        refined_source.append(source_values[12][1])
        refined_target.append(target_values[6][1])

        # D^-6 Q is genuinely inserted in q_pair.  On this carrier it is
        # redundant at the coefficient level; record that sharp boundary.
        no_q_source.append(relative.private.delayed_carry_pair(
            source_cut, pair_prefixes[0][ell][6], {}
        )[12][1])
        no_q_target.append(relative.private.delayed_carry_pair(
            target_cut, pair_prefixes[0][ell][6], {}
        )[6][1])

        # Hostile: omitting the pulled-target/source F factor destroys raw
        # equality even though both one-sided endpoint sections are lawful.
        source_one = tuple(relative.private.old.intersect_weighted_union(
            source_pieces, full_section
        ))
        target_one = tuple(relative.private.old.intersect_weighted_union(
            target_pieces, full_section
        ))
        chosen_one_sided_source_carriers.append(source_one)
        chosen_one_sided_target_carriers.append(target_one)
        one_sided_source.append(relative.private.delayed_carry_pair(
            source_one, q_pair, {}
        )[12][1])
        one_sided_target.append(relative.private.delayed_carry_pair(
            target_one, q_pair, {}
        )[6][1])

        # Reconcile the old unmarked overlap with the new one-sided hostile.
        # On this chosen fibre the q1/c2,q2/c3 target-safe cuts do not further
        # alter the coefficient after E3 and the physical c1 clock are added.
        e3_clock = target.intersect_sorted(e3, clock_comb[ell])
        source_clock = relative.private.old.intersect_weighted_union(
            source_pieces, e3_clock
        )
        target_clock = relative.private.old.intersect_weighted_union(
            target_pieces, e3_clock
        )
        e3_clock_source.append(relative.private.delayed_carry_pair(
            source_clock, q_pair, {}
        )[12][1])
        e3_clock_target.append(relative.private.delayed_carry_pair(
            target_clock, q_pair, {}
        )[6][1])

        # Independent route: cut the predecessor carry cell explicitly and
        # use the older weighted prefix numerator, not delayed_carry_pair.
        direct_source.append(
            relative.private.old.delayed_weighted_numerator(
                carry_cell_restriction(source_cut, 12), q_pair[1]
            )
        )
        direct_target.append(
            relative.private.old.delayed_weighted_numerator(
                carry_cell_restriction(target_cut, 6), q_pair[1]
            )
        )

    refined_source = tuple(refined_source)
    refined_target = tuple(refined_target)
    no_q_source = tuple(no_q_source)
    no_q_target = tuple(no_q_target)
    one_sided_source = tuple(one_sided_source)
    one_sided_target = tuple(one_sided_target)
    e3_clock_source = tuple(e3_clock_source)
    e3_clock_target = tuple(e3_clock_target)
    q_only_source = tuple(q_only_source)
    q_only_target = tuple(q_only_target)
    direct_source = tuple(direct_source)
    direct_target = tuple(direct_target)

    expected_vector = (
        0,
        339633525654239542165440,
        750593782703678965571520,
        719200126392878704654080,
        0, 0, 0,
    )
    expected_one_sided_source = (
        0,
        339633525654239542165440,
        750593782703678965571520,
        722054095148406001101120,
        0, 0, 0,
    )
    expected_one_sided_target = (
        0,
        345341652135823400016960,
        756301720214733558465600,
        724908063903933297548160,
        0, 0, 0,
    )
    require(
        tuple(piece_counts)
        == ((0, 0), (239, 239), (526, 526), (504, 504),
            (0, 0), (0, 0), (0, 0))
        and refined_source == refined_target == expected_vector
        and direct_source == direct_target == expected_vector,
        "chosen full-target semantic clutch vector changed",
    )
    require(
        no_q_source == no_q_target == expected_vector,
        "terminal-fork redundancy boundary changed",
    )
    require(
        q_only_source == q_only_target == raw_source == raw_target,
        "terminal fork unexpectedly changed the old overlap vector",
    )
    require(
        one_sided_source == expected_one_sided_source
        and one_sided_target == expected_one_sided_target
        and e3_clock_source == expected_one_sided_source
        and e3_clock_target == expected_one_sided_target
        and one_sided_source != one_sided_target,
        "one-sided full-target hostile changed",
    )

    # First typed repair of MISTAKE-313: if the physical clock is fixed at
    # e=1 on both the one-sided and common carriers, the source wing is
    # coefficient-null and the right wing carries the one-sided defect.
    fixed_e1_common = refined_source[1]
    fixed_e1_source = one_sided_source[1]
    fixed_e1_target = one_sided_target[1]
    fixed_e1_left_wing = fixed_e1_source - fixed_e1_common
    fixed_e1_right_wing = fixed_e1_target - fixed_e1_common
    fixed_e1_source_profile = overlap.normalized_profile(
        (0,) + (fixed_e1_source,) * 6, 12
    )[1]
    fixed_e1_target_profile = overlap.normalized_profile(
        (0,) + (fixed_e1_target,) * 6, 1
    )[1]
    fixed_e1_right_profile = overlap.normalized_profile(
        (0,) + (fixed_e1_right_wing,) * 6, 1
    )[1]
    fixed_e1_folded_gain = (
        fixed_e1_target_profile[0]
        * pow(fixed_e1_source_profile[0], -1, P) % P
    )
    require(
        fixed_e1_common == fixed_e1_source
        == 339633525654239542165440
        and fixed_e1_target == 345341652135823400016960
        and fixed_e1_left_wing == 0
        and fixed_e1_right_wing == 5708126481583857851520
        and fixed_e1_source_profile == (9, 0, 0, 0, 0, 0)
        and fixed_e1_target_profile == (8, 0, 0, 0, 0, 0)
        and fixed_e1_right_profile == (4, 0, 0, 0, 0, 0)
        and fixed_e1_folded_gain == 11,
        "fixed-e=1 MISTAKE-313 repair changed",
    )

    # Mayer--Vietoris hostile: the physical clock sections are a disjoint
    # partition of the unclocked sheet on this carrier.  Compare the union of
    # all same-clock two-sided pieces with the literal unclocked intersection,
    # including every possible cross-clock pair, before subtracting wings.
    chosen_unclocked_section = overlap.merge_intervals(
        interval
        for section in chosen_physical_sections
        for interval in section
    )
    physical_sections_disjoint = all(
        not target.intersect_sorted(
            chosen_physical_sections[first],
            chosen_physical_sections[second],
        )
        for first in range(7)
        for second in range(first + 1, 7)
    )
    common_source_disjoint = all(
        not target.intersect_sorted(
            tuple((left, right) for left, right, _weight
                  in chosen_common_source_carriers[first]),
            tuple((left, right) for left, right, _weight
                  in chosen_common_source_carriers[second]),
        )
        for first in range(7)
        for second in range(first + 1, 7)
    )
    common_target_disjoint = all(
        not target.intersect_sorted(
            tuple((left, right) for left, right, _weight
                  in chosen_common_target_carriers[first]),
            tuple((left, right) for left, right, _weight
                  in chosen_common_target_carriers[second]),
        )
        for first in range(7)
        for second in range(first + 1, 7)
    )
    reference_source_pieces, reference_target_pieces = overlap_details[1]
    full_common_source = common_physical_cut(
        reference_source_pieces, chosen_unclocked_section, -1
    )
    full_common_target = common_physical_cut(
        reference_target_pieces, chosen_unclocked_section, +1
    )
    full_one_sided_source = tuple(
        relative.private.old.intersect_weighted_union(
            reference_source_pieces, chosen_unclocked_section
        )
    )
    full_one_sided_target = tuple(
        relative.private.old.intersect_weighted_union(
            reference_target_pieces, chosen_unclocked_section
        )
    )
    same_clock_common_source = merge_weighted(
        piece for carrier in chosen_common_source_carriers for piece in carrier
    )
    same_clock_common_target = merge_weighted(
        piece for carrier in chosen_common_target_carriers for piece in carrier
    )
    partitioned_one_sided_source = merge_weighted(
        piece
        for carrier in chosen_one_sided_source_carriers
        for piece in carrier
    )
    partitioned_one_sided_target = merge_weighted(
        piece
        for carrier in chosen_one_sided_target_carriers
        for piece in carrier
    )
    require(
        physical_sections_disjoint
        and common_source_disjoint and common_target_disjoint
        and same_clock_common_source == merge_weighted(full_common_source)
        and same_clock_common_target == merge_weighted(full_common_target)
        and partitioned_one_sided_source
        == merge_weighted(full_one_sided_source)
        and partitioned_one_sided_target
        == merge_weighted(full_one_sided_target),
        "physical-clock object union stopped being the full unclocked sheet",
    )
    common_amplitude_source = relative.private.delayed_carry_pair(
        full_common_source, fork_pair_prefixes[1], {}
    )[12][1]
    common_amplitude_target = relative.private.delayed_carry_pair(
        full_common_target, fork_pair_prefixes[1], {}
    )[6][1]
    natural_amplitude_source = relative.private.delayed_carry_pair(
        full_one_sided_source, fork_pair_prefixes[1], {}
    )[12][1]
    natural_amplitude_target = relative.private.delayed_carry_pair(
        full_one_sided_target, fork_pair_prefixes[1], {}
    )[6][1]
    literal_source_wing = subtract_weighted(
        full_one_sided_source, full_common_source
    )
    literal_target_wing = subtract_weighted(
        full_one_sided_target, full_common_target
    )
    literal_source_wing_value = relative.private.delayed_carry_pair(
        literal_source_wing, fork_pair_prefixes[1], {}
    )[12][1]
    literal_target_wing_value = relative.private.delayed_carry_pair(
        literal_target_wing, fork_pair_prefixes[1], {}
    )[6][1]
    source_wing = natural_amplitude_source - common_amplitude_source
    target_wing = natural_amplitude_target - common_amplitude_target
    source_wing_components = tuple(
        natural - common
        for natural, common in zip(one_sided_source, refined_source)
    )
    target_wing_components = tuple(
        natural - common
        for natural, common in zip(one_sided_target, refined_target)
    )
    source_wing_component_residues = tuple(
        value // CONTENT * pow(12, -1, P) % P
        for value in source_wing_components
    )
    target_wing_component_residues = tuple(
        value // CONTENT % P for value in target_wing_components
    )
    source_wing_physical_profile = tuple(
        (source_wing_component_residues[index]
         - source_wing_component_residues[-1]) % P
        for index in range(6)
    )
    target_wing_physical_profile = tuple(
        (target_wing_component_residues[index]
         - target_wing_component_residues[-1]) % P
        for index in range(6)
    )
    source_wing_physical_det = (
        relative.private.old.sat.multiplication_determinant_7(
            source_wing_physical_profile
        )
    )
    target_wing_physical_det = (
        relative.private.old.sat.multiplication_determinant_7(
            target_wing_physical_profile
        )
    )
    physical_wing_ratio = (0, 2, 2, 2, 2, 6)
    physical_wing_ratio_det = (
        relative.private.old.sat.multiplication_determinant_7(
            physical_wing_ratio
        )
    )
    target_wing_factorization = multiply_phi_seven(
        multiply_phi_seven(
            (0, 2, 0, 0, 0, 0),
            (12, 1, 0, 0, 0, 0),
        ),
        (2, 1, 0, 0, 0, 0),
    )
    source_wing_vector = (0,) + (source_wing,) * 6
    target_wing_vector = (0,) + (target_wing,) * 6
    source_wing_profile = overlap.normalized_profile(
        source_wing_vector, 12
    )[1]
    target_wing_profile = overlap.normalized_profile(
        target_wing_vector, 1
    )[1]
    source_wing_det = (
        relative.private.old.sat.multiplication_determinant_7(
            source_wing_profile
        )
    )
    target_wing_det = (
        relative.private.old.sat.multiplication_determinant_7(
            target_wing_profile
        )
    )
    source_wing_root_residue = source_wing // CONTENT * pow(12, -1, P) % P
    target_wing_root_residue = target_wing // CONTENT % P
    require(
        common_amplitude_source == common_amplitude_target
        == sum(expected_vector)
        == 1809427434750797212391040
        and natural_amplitude_source == sum(expected_one_sided_source)
        == 1812281403506324508838080
        and natural_amplitude_target == sum(expected_one_sided_target)
        == 1826551436254490256030720
        and source_wing == 2853968755527296447040
        and target_wing == 17124001503693043639680
        and literal_source_wing_value == source_wing
        and literal_target_wing_value == target_wing
        and source_wing_components == (
            0, 0, 0, 2853968755527296447040, 0, 0, 0
        )
        and target_wing_components == (
            0,
            5708126481583857851520,
            5707937511054592894080,
            5707937511054592894080,
            0, 0, 0,
        )
        and source_wing_component_residues
        == (0, 0, 0, 12, 0, 0, 0)
        and target_wing_component_residues
        == (0, 9, 2, 2, 0, 0, 0)
        and source_wing_physical_profile == (0, 0, 0, 12, 0, 0)
        and target_wing_physical_profile == (0, 9, 2, 2, 0, 0)
        and source_wing_physical_det == 1
        and target_wing_physical_det == 11
        and physical_wing_ratio_det == 11
        and multiply_phi_seven(
            physical_wing_ratio, source_wing_physical_profile
        ) == target_wing_physical_profile
        and target_wing_factorization == target_wing_physical_profile
        and sum(source_wing_component_residues) % P == 12
        and sum(target_wing_component_residues) % P == 0
        and valuation(source_wing, P) == 1
        and valuation(target_wing, P) == 2
        and source_wing_root_residue == 12
        and target_wing_root_residue == 0
        and source_wing_profile == (1, 0, 0, 0, 0, 0)
        and target_wing_profile == (0, 0, 0, 0, 0, 0)
        and source_wing_det == 1 and target_wing_det == 0
        and relative.private.is_unit(source_wing_vector, 12, CONTENT)
        and not relative.private.is_unit(target_wing_vector, 1, CONTENT),
        "full-intersection wing degeneration changed",
    )

    chosen_gcd = gcd(*refined_source)
    require(
        chosen_gcd == 41337303276709440
        and valuation(chosen_gcd, P) == 1
        and all(value % CONTENT == 0 for value in refined_source),
        "chosen vector lost canonical content-26 typing",
    )
    source_profile = overlap.normalized_profile(refined_source, 12)[1]
    target_profile = overlap.normalized_profile(refined_target, 1)[1]
    require(
        source_profile == (0, 4, 10, 8, 0, 0)
        and target_profile == (0, 9, 3, 5, 0, 0)
        and target_profile
        == tuple((-value) % P for value in source_profile),
        "endpoint root-normalization sign law changed",
    )
    source_det = relative.private.old.sat.multiplication_determinant_7(
        source_profile
    )
    target_det = relative.private.old.sat.multiplication_determinant_7(
        target_profile
    )
    require(
        source_det == target_det == 1
        and relative.private.is_unit(refined_source, 12, CONTENT)
        and relative.private.is_unit(refined_target, 1, CONTENT),
        "chosen refined clutch lost a private unit",
    )

    # Exhaust the full two-clock object over every label stable on both open
    # endpoint cylinders.  Here e is the physical c1-present clock and j is
    # the delayed coefficient clock.  The old diagonal calculation used e=j.
    # We now construct all 81*7*7 cells before taking that diagonal.
    common_rows = []
    two_clock_universe = len(common_s) * len(common_t) * 7 * 7
    carrier_object_checks = 0
    carrier_translation_checks = 0
    raw_source_target_checks = 0
    outer_product_checks = 0
    rank_one_labels = 0
    diagonal_amplitude_labels = 0
    physical_clock_zero_edge_labels = 0
    fixed_thm2749_bridge_labels = 0
    chosen_two_clock_matrix = None
    thm2749_coindexed_bridge = None
    for s in common_s:
        for t in common_t:
            source_matrix = []
            target_matrix = []
            carrier_translation = True
            carrier_object_law = True
            for physical_clock in range(7):
                full_section = target.source_present_section(
                    full_module, e3, physical_clock, s, t, clock_comb
                )
                source_cuts = tuple(
                    common_physical_cut(source_pieces, full_section, -1)
                    for source_pieces, _target_pieces in overlap_details
                )
                target_cuts = tuple(
                    common_physical_cut(target_pieces, full_section, +1)
                    for _source_pieces, target_pieces in overlap_details
                )
                reference_source = source_cuts[1]
                reference_target = target_cuts[1]
                source_row = []
                target_row = []
                for delayed_clock in range(7):
                    source_cut = source_cuts[delayed_clock]
                    target_cut = target_cuts[delayed_clock]
                    same_object = (
                        source_cut == reference_source
                        and target_cut == reference_target
                    )
                    translated_object = (
                        overlap.shift_weighted(source_cut, SHIFT)
                        == target_cut
                    )
                    carrier_object_checks += int(same_object)
                    carrier_translation_checks += int(translated_object)
                    carrier_object_law = carrier_object_law and same_object
                    carrier_translation = (
                        carrier_translation and translated_object
                    )
                    q_pair = fork_pair_prefixes[delayed_clock]
                    source_value = relative.private.delayed_carry_pair(
                        source_cut, q_pair, {}
                    )[12][1]
                    target_value = relative.private.delayed_carry_pair(
                        target_cut, q_pair, {}
                    )[6][1]
                    source_row.append(source_value)
                    target_row.append(target_value)
                    raw_source_target_checks += int(
                        source_value == target_value
                    )
                source_matrix.append(tuple(source_row))
                target_matrix.append(tuple(target_row))

            source_matrix = tuple(source_matrix)
            target_matrix = tuple(target_matrix)
            source_amplitude = tuple(row[1] for row in source_matrix)
            target_amplitude = tuple(row[1] for row in target_matrix)
            expected_source_matrix = tuple(
                tuple(amplitude * value for value in delayed_window)
                for amplitude in source_amplitude
            )
            expected_target_matrix = tuple(
                tuple(amplitude * value for value in delayed_window)
                for amplitude in target_amplitude
            )
            for physical_clock in range(7):
                for delayed_clock in range(7):
                    outer_product_checks += int(
                        source_matrix[physical_clock][delayed_clock]
                        == expected_source_matrix[physical_clock][delayed_clock]
                        and target_matrix[physical_clock][delayed_clock]
                        == expected_target_matrix[physical_clock][delayed_clock]
                    )

            source_vector = tuple(
                source_matrix[clock][clock] for clock in range(7)
            )
            target_vector = tuple(
                target_matrix[clock][clock] for clock in range(7)
            )
            all_minors_zero = all(
                source_matrix[e][j] * source_matrix[f][k]
                == source_matrix[e][k] * source_matrix[f][j]
                for e in range(7) for f in range(7)
                for j in range(7) for k in range(7)
            )
            exact_rank_one = (
                source_matrix == target_matrix == expected_source_matrix
                and source_amplitude == target_amplitude
                and any(source_amplitude)
                and all_minors_zero
            )
            diagonal_is_amplitude = (
                source_vector == target_vector
                == source_amplitude == target_amplitude
            )
            rank_one_labels += int(exact_rank_one)
            diagonal_amplitude_labels += int(diagonal_is_amplitude)
            physical_clock_zero_edge_labels += int(
                source_amplitude[0] == target_amplitude[0] == 0
            )
            if s == 0:
                fixed_thm2749_bridge_labels += int(
                    source_matrix[1]
                    == (0,) + (expected_vector[1],) * 6
                )
            if (s, t) == (0, 3):
                chosen_two_clock_matrix = source_matrix
            if (s, t) == (0, 4):
                thm2749_coindexed_bridge = source_vector
            common_rows.append((
                s, t, carrier_translation and carrier_object_law,
                source_vector, target_vector, source_amplitude,
            ))

    require(len(common_rows) == 81, "common-label census size changed")
    require(
        carrier_object_checks == two_clock_universe
        and carrier_translation_checks == two_clock_universe
        and raw_source_target_checks == two_clock_universe
        and outer_product_checks == two_clock_universe
        and rank_one_labels == len(common_rows)
        and diagonal_amplitude_labels == len(common_rows)
        and physical_clock_zero_edge_labels == len(common_rows)
        and fixed_thm2749_bridge_labels == len(common_t),
        "two-clock rank-one separability census changed",
    )
    expected_chosen_matrix = tuple(
        tuple(amplitude * value for value in delayed_window)
        for amplitude in expected_vector
    )
    require(
        chosen_two_clock_matrix == expected_chosen_matrix,
        "chosen two-clock matrix changed",
    )
    require(
        thm2749_coindexed_bridge == expected_vector,
        "rail-8 (0,4) coindexed bridge to THM-2749 changed",
    )
    support_counts = Counter(
        tuple(index for index, value in enumerate(source_vector) if value)
        for _s, _t, _carrier, source_vector, _target_vector, _amplitude
        in common_rows
    )
    distinct_vectors = {
        source_vector
        for _s, _t, _carrier, source_vector, _target_vector, _amplitude
        in common_rows
    }
    gcds = tuple(
        gcd(*source_vector)
        for _s, _t, _carrier, source_vector, _target_vector, _amplitude
        in common_rows
    )
    require(
        all(carrier and source_vector == target_vector
            for _s, _t, carrier, source_vector, target_vector, _amplitude
            in common_rows)
        and all(any(source_vector)
                for _s, _t, _carrier, source_vector, _target, _amplitude
                in common_rows)
        and all(all(value % CONTENT == 0 for value in source_vector)
                for _s, _t, _carrier, source_vector, _target, _amplitude
                in common_rows),
        "common-label raw equality/content census changed",
    )
    require(
        all(relative.private.is_unit(source_vector, 12, CONTENT)
            and relative.private.is_unit(target_vector, 1, CONTENT)
            for _s, _t, _carrier, source_vector, target_vector, _amplitude
            in common_rows),
        "a common full-target label lost an endpoint unit",
    )
    require(
        support_counts == Counter({
            (1, 2, 3): 35,
            (1, 2): 21,
            (1, 3): 21,
            (1,): 4,
        })
        and len(distinct_vectors) == 15
        and min(gcds) == 5905329039529920
        and max(gcds) == 302530703523944466130560
        and all(valuation(value, P) == 1 for value in gcds),
        "common-label support/content profile changed",
    )

    # Fully mark the endpoint-dipole target coordinate.  Outside the nine
    # labels common to both strict cylinders the marked profile is zero by
    # definition; this is not a physical-emptiness assertion for unrestricted
    # full-target sections.  Clock ell=1 is positive on every common label,
    # and supplies an exact witness for all primitive target characters.
    row_lookup = {
        (s, t): source_vector
        for s, t, _carrier, source_vector, _target_vector, _amplitude
        in common_rows
    }
    zero_vector = (0,) * 7
    per_s_character_count = 0
    character_witness_clocks = Counter()
    for s in common_s:
        vector_profile = tuple(
            row_lookup.get((s, t), zero_vector)
            for t in range(P)
        )
        require(
            tuple(t for t, vector in enumerate(vector_profile) if any(vector))
            == common_t,
            f"marked target support changed at s={s}",
        )
        for frequency in range(1, P):
            coordinates_by_clock = tuple(
                target.primitive_fourier_coordinates(
                    tuple(vector[ell] for vector in vector_profile),
                    frequency,
                )
                for ell in range(7)
            )
            witnesses = tuple(
                ell for ell, coordinates in enumerate(coordinates_by_clock)
                if any(coordinates)
            )
            require(witnesses and witnesses[0] == 1,
                    "vector-valued primitive target character vanished")
            character_witness_clocks[witnesses[0]] += 1
            per_s_character_count += 1

    aggregate_target_profile = tuple(
        tuple(
            sum((row_lookup.get((s, t), zero_vector)[ell]
                 for s in common_s), 0)
            for ell in range(7)
        )
        for t in range(P)
    )
    aggregate_support = tuple(
        t for t, vector in enumerate(aggregate_target_profile) if any(vector)
    )
    aggregate_character_count = 0
    for frequency in range(1, P):
        clock_one_coordinates = target.primitive_fourier_coordinates(
            tuple(vector[1] for vector in aggregate_target_profile),
            frequency,
        )
        require(any(clock_one_coordinates),
                "aggregated primitive target character vanished")
        aggregate_character_count += 1
    aggregate_clock_one_profile = tuple(
        vector[1] for vector in aggregate_target_profile
    )
    aggregate_amplitude = 2554386600508776388555200
    require(
        per_s_character_count == len(common_s) * (P - 1)
        and character_witness_clocks == Counter({1: 108})
        and aggregate_support == common_t
        and aggregate_character_count == P - 1,
        "fully marked target-character census changed",
    )
    require(
        aggregate_clock_one_profile
        == (0, 0, 0) + (aggregate_amplitude,) * 9 + (0,),
        "aggregated clock-one target profile changed",
    )

    # The nine-label window has a positive integral inverse modulo Phi_13.
    window = (0, 0, 0) + (1,) * 9
    inverse_window = tuple(
        int(exponent in (2, 6, 10)) for exponent in range(11)
    )
    bezout_factor = (-1, 1, 0, 0, 0, 1, 0, 0, 0, 1)
    phi_thirteen = (1,) * 13
    window_times_inverse_minus_one = list(
        multiply_polynomials(window, inverse_window)
    )
    window_times_inverse_minus_one[0] -= 1
    require(
        tuple(window_times_inverse_minus_one)
        == multiply_polynomials(bezout_factor, phi_thirteen)
        and cyclic_convolution(
            window + (0,), inverse_window + (0, 0), P
        ) == (3,) + (2,) * 12,
        "nine-label target-window Bezout identity changed",
    )

    print("LRC14 ROOT-ZERO TWO-CLOCK FULL-TARGET CLUTCH AUDIT")
    print("status=FINITE-EXACT canonical rail-8 two-clock separability theorem")
    print(f"p={P} R={R} T={T} tau_grid={SHIFT} content={CONTENT}")
    print(f"dependency_hashes={tuple(sorted(DEPENDENCIES.items()))}")
    print(f"common_full_target_bank=s{common_s} t{common_t} count=81")
    print(f"common_label_factor_margin={common_label_margin} ratio_to_inherited_radius=1541619 all_81_stable_on_both_open_cylinders")
    print("chosen_label=(s,t)=(0,3)")
    print(f"chosen_piece_counts_by_ell={tuple(piece_counts)}")
    print(f"chosen_source_vector={refined_source}")
    print(f"chosen_target_vector={refined_target}")
    print(f"chosen_joint_gcd={chosen_gcd} v13={valuation(chosen_gcd,P)}")
    print(f"chosen_profiles=source_root12:{source_profile} target_root1:{target_profile}")
    print(f"chosen_determinants=source:{source_det} target:{target_det}")
    print("normalization_mechanism=root12^-1=-1 and root1^-1=1, so endpoint classes are negatives; degree6 norm preserves the determinant")
    print("translation_mechanism=y={R*x} is invariant and predecessor carry shifts 12->6 under x->x+7/R")
    print("terminal_fork_cut=inserted by pair_prefix intersect Q_(3,{1,2}); coefficient unchanged on this overlap carrier")
    print("delayed_prefix_object_law=j0 pair empty; j1=...=j6 as exact rebuilt prefix objects")
    print("carrier_independence_mechanism=E3 intersect Present_(j,7) is empty for every j (c3-danger versus inherited t0 c3-safe)")
    print("old_probe_reconciliation=Q-refined prefixes leave the unmarked old overlap vector unchanged")
    print(f"one_sided_hostile_source={one_sided_source}")
    print(f"one_sided_hostile_target={one_sided_target}")
    print("one_sided_hostile=UNEQUAL; both F(x) and F(x+tau) are load-bearing")
    print("one_sided_hostile_mechanism=E3 plus the physical c1 clock already gives the displayed unequal rows; shifted target-safe factors do not further change them on this fibre")
    print(f"fixed_e1_MISTAKE313_repair=source:{fixed_e1_source}=common:{fixed_e1_common} left_wing:{fixed_e1_left_wing}; target:{fixed_e1_target} right_wing:{fixed_e1_right_wing}; profiles=({fixed_e1_source_profile},{fixed_e1_target_profile},{fixed_e1_right_profile}) folded_gain:{fixed_e1_folded_gain}")
    print("forgotten_e_clock_blind_partition=seven F_e pairwise disjoint; same-clock two-sided union equals full unclocked A_intersect_B object; all e!=e' cross pieces empty")
    print(f"forgotten_e_full_unclocked_amplitudes=A:{natural_amplitude_source} M_full=A_intersect_B:{common_amplitude_source} B:{natural_amplitude_target}")
    print(f"forgotten_e_wings=L=A-M_full:{source_wing} R=B-M_full:{target_wing}")
    print(f"true_wing_physical_rows=source:{source_wing_components} target:{target_wing_components} literal_object_subtraction_matches")
    print(f"true_wing_physical_residues_after_content26_and_roots=source:{source_wing_component_residues} target:{target_wing_component_residues} augmentations=({sum(source_wing_component_residues)%P},{sum(target_wing_component_residues)%P}) valuations=({valuation(source_wing,P)},{valuation(target_wing,P)})")
    print(f"physical_wing_units=source_profile:{source_wing_physical_profile} det:{source_wing_physical_det} target_profile:{target_wing_physical_profile}=2z(z-1)(z+2) det:{target_wing_physical_det}")
    print(f"physical_wing_ratio=g:{physical_wing_ratio}=2z+2z^2+2z^3+2z^4+6z^5 det:{physical_wing_ratio_det}; g*source=target")
    print(f"augmented_delayed_wing_profiles=source_root12:{source_wing_root_residue}:{source_wing_profile} target_root1:{target_wing_root_residue}:{target_wing_profile} determinants=({source_wing_det},{target_wing_det})")
    print("wing_boundary=fixed-e repair has null left wing; forgotten-e physical rows are units with coefficient-algebra ratio g, but target augmentation is zero so no scalar wing gain; no physical L-to-R carrier map or cyclic operation; THM2751 gain2 frame is refuted")
    print("direct_control=explicit carry-cell cut plus delayed_weighted_numerator matches all 14 chosen coefficients")
    print(f"two_clock_universe=81*7*7={two_clock_universe} physical_clock=e delayed_clock=j")
    print(f"two_clock_object_census=carrier_independent_of_j:{carrier_object_checks}/{two_clock_universe} translated:{carrier_translation_checks}/{two_clock_universe} raw_source_target_equal:{raw_source_target_checks}/{two_clock_universe}")
    print(f"two_clock_rank_one_census=outer_product_cells:{outer_product_checks}/{two_clock_universe} exact_rank_one_labels:{rank_one_labels}/81 diagonal_equals_physical_amplitude:{diagonal_amplitude_labels}/81")
    print(f"two_clock_factorization=B_ej(s,t)=a_e(s,t)*w_j with w={delayed_window}=coefficients(z+...+z^6)=Phi_7-1")
    print(f"two_clock_edge_cases=delayed_j0_zero_all; physical_e0_amplitude_zero:{physical_clock_zero_edge_labels}/81; fixed_e1_s0_bridge_to_THM2749:{fixed_thm2749_bridge_labels}/9; coindexed_(0,4)_rail8_bridge:1/1")
    print(f"chosen_two_clock_matrix={chosen_two_clock_matrix}")
    print("common_bank_census=carrier_translation:81/81 raw_equal:81/81 content26:81/81 units_at_both_roots:81/81")
    print(f"common_bank_support_profiles={tuple(sorted(support_counts.items()))}")
    print(f"common_bank_distinct_vectors={len(distinct_vectors)} gcd_range=({min(gcds)},{max(gcds)}) all_v13=1")
    print(f"marked_target_support={aggregate_support} zero_outside_common_bank_by_definition")
    print(f"aggregate_clock1_target_profile={aggregate_clock_one_profile}")
    print(f"primitive_target_characters=per_s:{per_s_character_count}/108 aggregate:{aggregate_character_count}/12 witness_clock=1")
    print("scalar_factored_target_support_window_unit=W=u^3+...+u^11 inverse=u^2+u^6+u^10; W*inverse-1=(u^9+u^5+u-1)Phi_13; cyclic_product=(3,2,...,2)")
    print("target_window_norm=1 because W(zeta)=zeta^3(1-zeta^9)/(1-zeta), gcd(9,13)=1")
    print("target_character_type=endpoint-dipole lambda coordinate; not a physical deck character, paired left-relation current, or endpoint amplitude")
    print("SCOPE: one canonical rail-8 overlap and its 81 common full-target labels; no global clutch action, physical endpoint current, row exclusion, or LRC14 conclusion")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
