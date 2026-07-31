#!/usr/bin/env python3
"""Exact fully marked root-zero overlap clutch for THM-2749.

This inserts the exclusive E3 source, the ordinary D^6 terminal fork
Q_(3,{1,2}), and one lawful full two-target section F_(ell,0,4) inside the
THM-2744 overlap integral.  Translation by 7/13^6 identifies the resulting
source and target carriers, fixes the delayed phase, and maps carry 12 to
carry 6.  Exact enumeration then proves positivity and a primitive unit on
all fourteen equal-weight rail loci.

The seven entries are a clock-indexed family of counterfactual fibres.  They
are not a seven-step orbit, a THM-2334 relation-address current, an endpoint
amplitude, a row exclusion, or an LRC(14) proof.  The target-character
profile below is separately defined by attaching this family to the one
explicit transported open cylinder; this distinction removes the extra raw
``t=2`` coefficient, which does not contain that cylinder.
"""

from __future__ import annotations

from bisect import bisect_right
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


for dependency, expected_hash in DEPENDENCIES.items():
    payload = (COMP / dependency).read_bytes().replace(b"\r\n", b"\n")
    actual_hash = hashlib.sha256(payload).hexdigest()
    require(actual_hash == expected_hash,
            f"audited dependency changed: {dependency}")


import lrc14_root_zero_overlap_clutch_20260728 as clutch
import lrc14_two_target_present_semantic_attachment_probe_20260728 as semantic


P = 13
R = P**6
S = R // P
T = clutch.T
SHIFT = clutch.SHIFT
H = 6
KAPPA = 1
LABEL = (0, 4)
G04 = 413_915_261_760

EXPECTED_SUPPORTS = (
    (1,), (6,), (2, 3), (5,), (2, 3), (4, 6), (2, 3),
    (5, 6), (1, 2, 3), (5,), (1, 3), (5,), (2,), (5, 6),
)
EXPECTED_DETERMINANTS = (1, 1, 3, 1, 12, 8, 12, 3, 12, 12, 1, 12, 1, 1)
EXPECTED_ATTACHED_TARGET_PROFILE = (
    0, 0, 0,
    4_371_492_433_154, 4_371_492_433_154, 4_371_492_433_154,
    4_371_492_433_154, 4_371_492_433_154,
    2_633_938_414_646, 2_633_938_414_646,
    2_558_092_802_727, 2_558_092_802_727,
    0,
)
EXPECTED_ATTACHED_CROSS_PROFILE = (
    0, 0, 0,
    6_980_796_083_273_674_034_188_354,
    6_980_796_083_273_674_034_188_354,
    6_980_796_083_273_674_034_188_354,
    6_980_796_083_273_674_034_188_354,
    6_980_796_083_273_674_034_188_354,
    3_961_702_116_040_374_827_642_290,
    3_961_702_116_040_374_827_642_290,
    3_692_377_863_640_893_849_986_025,
    3_692_377_863_640_893_849_986_025,
    0,
)
EXPECTED_VECTORS = (
    (0, 587951430950255675022720, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 998922134277175124016000),
    (0, 0, 505152469728331471126080, 433803250840149059950080, 0, 0, 0),
    (0, 0, 0, 0, 0, 753472695569069236400640, 0),
    (0, 0, 753447751459206262018560, 753447751459206262018560, 0, 0, 0),
    (0, 0, 0, 0, 219755594175601826422080, 0, 753472695569069236400640),
    (0, 0, 750593782703678965571520, 750593782703678965571520, 0, 0, 0),
    (0, 0, 0, 0, 0, 750593782703678965571520, 342487588895031471091200),
    (0, 339633525654239542165440, 750593782703678965571520,
     719200126392878704654080, 0, 0, 0),
    (0, 0, 0, 0, 0, 750593782703678965571520, 0),
    (0, 750618632328277307474880, 0, 216901625420074529975040, 0, 0, 0),
    (0, 0, 0, 0, 0, 750593782703678965571520, 0),
    (0, 0, 753472695569069236400640, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 502298500972804174679040, 639310165937392079370240),
)


def build_semantic_prefixes(module, terminal_fork):
    """[sector][clock][kappa] prefixes on y={13^6 x}, at h=6."""
    result = []
    for word in clutch.relative.private.prior.sector_words(module):
        by_clock = []
        for ell in range(7):
            qell = module.subtract_comb(
                word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
            )
            qell = semantic.intersect_sorted(qell, terminal_fork)
            pair = []
            for kappa in range(2):
                digit = clutch.relative.private.old.sat.intersect_interval(
                    qell,
                    (2 * H + kappa) * T // (2 * P),
                    (2 * H + kappa + 1) * T // (2 * P),
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


def prefix_intervals(prefix):
    return tuple(
        (left, left + length)
        for left, length in zip(prefix[0], prefix[1])
    )


def strict_slack(value, intervals, denominator):
    """Exact radius of the open interval containing ``value``."""
    coordinate = value * denominator
    starts = [left for left, _right in intervals]
    index = bisect_right(starts, coordinate) - 1
    require(
        index >= 0
        and intervals[index][0] < coordinate < intervals[index][1],
        "point has no strict open-interval slack",
    )
    left, right = intervals[index]
    return min(coordinate - left, right - coordinate) / denominator


def supports_open_cylinder(value, intervals, denominator, radius):
    """Whether the open radius-``radius`` cylinder stays in one interval."""
    if not semantic.strict_member(value, intervals, denominator):
        return False
    return strict_slack(value, intervals, denominator) >= radius


def primitive_target_coordinates(values, frequency):
    """Coordinates of sum_t values[t] zeta_13^(frequency*t).

    The basis is 1,zeta,...,zeta^11 and zeta^12 is replaced using Phi_13.
    Thus a nonzero returned tuple is an exact primitive-character witness.
    """
    require(len(values) == P and 0 < frequency < P,
            "invalid primitive target-character request")
    coordinates = [0] * (P - 1)
    for target, value in enumerate(values):
        exponent = frequency * target % P
        if exponent == P - 1:
            for index in range(P - 1):
                coordinates[index] -= value
        else:
            coordinates[exponent] += value
    return tuple(coordinates)


def direct_lifted_numerator(weighted, target):
    """Enumerate R[a,b] intersect (kT+[u,v]) without a prefix formula."""
    total = 0
    hits = 0
    for left, right, weight in weighted:
        rleft = R * left
        rright = R * right
        for low, high in target:
            first_turn = (rleft - high) // T + 1
            last_turn = (rright - low - 1) // T
            for turn in range(first_turn, last_turn + 1):
                a = max(rleft, turn * T + low)
                b = min(rright, turn * T + high)
                if a < b:
                    total += weight * (b - a)
                    hits += 1
    return total, hits


def multiplication_profile(values, root, content):
    require(root and content > 0
            and all(value % content == 0 for value in values),
            "selected content does not divide a private row")
    inverse = pow(root, -1, P)
    normalized = tuple(
        (value // content) * inverse % P for value in values
    )
    reduced = tuple(
        (normalized[index] - normalized[-1]) % P
        for index in range(6)
    )
    determinant = (
        clutch.relative.private.old.sat.multiplication_determinant_7(
            reduced
        )
    )
    return normalized, reduced, determinant


def fully_marked_vectors(module, rails, present, prefixes, sections,
                         rail_index, *, common_sections=True,
                         equal_weight_only=True):
    """Recompute source/root12 and target/root1 rows on one rail.

    With common_sections=True, both source and translated-target full target
    sections are retained.  The resulting physical carriers must be literal
    translates before either delayed-carry functional is evaluated.
    """
    target_pullback = clutch.shift_weighted(
        rails[rail_index][3], -SHIFT
    )
    rail_pairs = clutch.intersect_weighted_profiles(
        rails[rail_index][3], target_pullback
    )
    if equal_weight_only:
        rail_pairs = tuple(
            row for row in rail_pairs if row[2] == row[3]
        )
    require(rail_pairs, f"rail {rail_index} has no weighted overlap")

    source_vector = []
    target_vector = []
    details = []
    for ell in range(7):
        source, target = clutch.restrict_to_relative_overlap(
            module, present, rail_pairs, ell
        )
        if common_sections:
            source_section = clutch.intersect_sorted(
                sections[ell],
                clutch.shift_union(sections[ell], -SHIFT),
            )
            target_section = clutch.intersect_sorted(
                sections[ell],
                clutch.shift_union(sections[ell], SHIFT),
            )
        else:
            source_section = sections[ell]
            target_section = sections[ell]
        source = tuple(
            clutch.relative.private.old.intersect_weighted_union(
                source, source_section
            )
        )
        target = tuple(
            clutch.relative.private.old.intersect_weighted_union(
                target, target_section
            )
        )

        source_values = clutch.relative.private.delayed_carry_pair(
            source, prefixes[0][ell], {}
        )
        target_values = clutch.relative.private.delayed_carry_pair(
            target, prefixes[0][ell], {}
        )
        source_value = source_values[12][KAPPA]
        target_value = target_values[6][KAPPA]

        source_carry = tuple(
            clutch.relative.private.old.intersect_weighted_comb(
                source, S, P, 12, 13
            )
        )
        target_carry = tuple(
            clutch.relative.private.old.intersect_weighted_comb(
                target, S, P, 6, 7
            )
        )
        direct_source = (
            clutch.relative.private.old.delayed_weighted_numerator(
                source_carry, prefixes[0][ell][KAPPA]
            )
        )
        direct_target = (
            clutch.relative.private.old.delayed_weighted_numerator(
                target_carry, prefixes[0][ell][KAPPA]
            )
        )
        require(
            (source_value, target_value)
            == (direct_source, direct_target),
            f"carry descent/direct carry-cell disagreement on rail "
            f"{rail_index}, clock {ell}",
        )
        if common_sections:
            require(
                clutch.shift_weighted(source, SHIFT) == target
                and clutch.shift_weighted(source_carry, SHIFT)
                == target_carry,
                f"marked carrier transport failed on rail {rail_index}, "
                f"clock {ell}",
            )
        source_vector.append(source_value)
        target_vector.append(target_value)
        details.append((source, target, source_carry, target_carry))
    return tuple(source_vector), tuple(target_vector), tuple(details)


def carry_at(value):
    z = clutch.relative.semantic.frac(S * value)
    return clutch.relative.half.floor_fraction(P * z)


def main():
    require((P, R, S, T, SHIFT) == (
        13, 4826809, 371293, 297836897838480, 431933040
    ), "canonical scales changed")
    require(
        R * SHIFT == 7 * T
        and S * SHIFT == 7 * T // P,
        "slope-seven phase/carry translation changed",
    )

    module, _pair_prefixes, _, _, rails, present, _starts = (
        clutch.relative.lift.m.core.build_carrier_data()
    )
    require(module.C3 == 742586 and len(rails) >= 14,
            "canonical deep blocker or source-one rail bank changed")
    source_e3 = semantic.exclusive_source(module, 3)
    terminal_fork = semantic.deepest_fork(module)
    prefixes = build_semantic_prefixes(module, terminal_fork)
    sections = semantic_sections(module, source_e3, *LABEL)

    vectors = []
    details = []
    content = 0
    for rail_index in range(14):
        source, target, rail_details = fully_marked_vectors(
            module, rails, present, prefixes, sections, rail_index
        )
        require(source == target,
                f"raw marked vectors differ on rail {rail_index}")
        vectors.append(source)
        details.append(rail_details)
        for value in source + target:
            content = gcd(content, value)
    vectors = tuple(vectors)
    details = tuple(details)
    require(vectors == EXPECTED_VECTORS,
            "uniform fully marked rail bank changed")
    require(content == G04 and content % P == 0
            and (content // P) % P != 0,
            "fully marked global content or its 13-valuation changed")

    rail_table = []
    for rail_index, vector in enumerate(vectors):
        source_profile = multiplication_profile(vector, 12, content)
        target_profile = multiplication_profile(vector, 1, content)
        support = tuple(
            index for index, value in enumerate(vector) if value
        )
        require(
            support == EXPECTED_SUPPORTS[rail_index]
            and source_profile[2] == target_profile[2]
            == EXPECTED_DETERMINANTS[rail_index]
            and source_profile[0]
            == tuple((-value) % P for value in target_profile[0])
            and source_profile[1]
            == tuple((-value) % P for value in target_profile[1]),
            f"unit/sign profile changed on rail {rail_index}",
        )
        rail_table.append((
            rail_index, support, source_profile[1], target_profile[1],
            source_profile[2],
        ))

    # A second exact integration route expands the lifted R*x intersections
    # directly on rail 8; it uses neither the delayed prefix antiderivative
    # nor the one-step carry descent.
    direct_vector = []
    direct_hits = []
    for ell in range(7):
        target = prefix_intervals(prefixes[0][ell][KAPPA])
        value, hits = direct_lifted_numerator(
            details[8][ell][2], target
        )
        direct_vector.append(value)
        direct_hits.append(hits)
    require(tuple(direct_vector) == vectors[8],
            "lifted-interval referee path disagrees on rail 8")

    # The original adjacent semantic pair is a strict whole-cylinder witness
    # for ell=1 and the chosen uniform label (0,4).
    q_source = clutch.relative.semantic.frac(
        clutch.relative.Z + Fraction(7 * 6715, R)
    )
    q_target = clutch.relative.semantic.frac(
        clutch.relative.Z + Fraction(7 * 6716, R)
    )
    require(
        q_source == Fraction(47850889647341, 100360982066072)
        and q_target == Fraction(47851035194197, 100360982066072)
        and q_target == q_source + Fraction(7, R),
        "adjacent semantic witness changed",
    )
    delayed_center = clutch.relative.semantic.frac(R * q_source)
    require(
        strict_slack(q_source, sections[1], T)
        == 56_447 * clutch.relative.Q_RADIUS
        and strict_slack(q_target, sections[1], T)
        == 8_854_585 * clutch.relative.Q_RADIUS
        and strict_slack(
            delayed_center, terminal_fork, T
        ) / R == clutch.relative.Q_RADIUS
        and strict_slack(
            delayed_center,
            prefix_intervals(prefixes[0][1][KAPPA]), T,
        ) / R == clutch.relative.Q_RADIUS
        and clutch.weighted_slack(
            q_source * T, details[8][1][0]
        ) == 56_447 * clutch.relative.Q_RADIUS
        and clutch.weighted_slack(
            q_target * T, details[8][1][1]
        ) == 56_447 * clutch.relative.Q_RADIUS,
        "strengthened whole-cylinder slack changed",
    )
    # The cylinder is open: its endpoints lie on the delayed-fork boundary.
    # Probe strict interior points after certifying the exact boundary slack.
    for delta in (-clutch.relative.Q_RADIUS / 2,
                  clutch.relative.Q_RADIUS / 2):
        source_point = q_source + delta
        target_point = q_target + delta
        delayed_source = clutch.relative.semantic.frac(R * source_point)
        delayed_target = clutch.relative.semantic.frac(R * target_point)
        require(
            semantic.strict_member(source_point, sections[1], T)
            and semantic.strict_member(target_point, sections[1], T)
            and semantic.strict_member(
                delayed_source,
                prefix_intervals(prefixes[0][1][KAPPA]), T,
            )
            and semantic.strict_member(delayed_source, terminal_fork, T)
            and delayed_source == delayed_target
            and carry_at(source_point) == 12
            and carry_at(target_point) == 6
            and clutch.containing_weighted_piece(
                source_point * T, details[8][1][0]
            )
            and clutch.containing_weighted_piece(
                target_point * T, details[8][1][1]
            )
            and clutch.relative.semantic.semantic_record(source_point)
            == (3, (1, 2))
            and clutch.relative.semantic.semantic_record(target_point)
            == (3, (1, 2))
            and LABEL[0]
            in clutch.relative.full_target_labels(source_point)[0]
            and LABEL[0]
            in clutch.relative.full_target_labels(target_point)[0]
            and LABEL[1]
            in clutch.relative.full_target_labels(source_point)[1]
            and LABEL[1]
            in clutch.relative.full_target_labels(target_point)[1],
            "strengthened whole-cylinder witness left a marked factor",
        )

    # The unqualified nine-label claim in the reservation has two distinct
    # readings.  The raw rail-8 integral is positive also at t=2.  The
    # transported-cylinder-attached pushforward is instead
    #
    #   B_t = 1_{the q_source/q_target open cylinder lies in F_(1,0,t)}
    #         * G04^{-1} * sum_ell A_(ell,t).
    #
    # The bilinear comparison is separately typed on the common clock index:
    #
    #   C_t = 1_{same cylinder condition}
    #         * G04^{-2} * sum_ell A^-_(ell,t) A^+_(ell,t).
    #
    # Both profiles have support t=3,...,11.  Since a rational length-13
    # profile has a vanishing primitive Fourier coefficient iff it is
    # constant (its polynomial must be a multiple of Phi_13), nonconstancy
    # proves all twelve primitive characters survive for both B and C.
    raw_target_rows = []
    attached_target_profile = []
    attached_cross_profile = []
    for target_label in range(P):
        target_sections = semantic_sections(
            module, source_e3, LABEL[0], target_label
        )
        target_source, target_target, target_details = fully_marked_vectors(
            module, rails, present, prefixes, target_sections, 8,
            equal_weight_only=True,
        )
        require(target_source == target_target,
                f"target scan transport failed at t={target_label}")
        cylinder_attached = (
            supports_open_cylinder(
                q_source, target_sections[1], T,
                clutch.relative.Q_RADIUS,
            )
            and supports_open_cylinder(
                q_target, target_sections[1], T,
                clutch.relative.Q_RADIUS,
            )
        )
        if cylinder_attached:
            require(
                clutch.containing_weighted_piece(
                    q_source * T, target_details[1][0]
                )
                and clutch.containing_weighted_piece(
                    q_target * T, target_details[1][1]
                ),
                f"attached target label {target_label} lost the marked pair",
            )
        raw_amplitude = sum(target_source)
        require(raw_amplitude % G04 == 0,
                f"target amplitude left selected content at t={target_label}")
        attached_amplitude = (
            raw_amplitude // G04 if cylinder_attached else 0
        )
        cross_amplitude = (
            sum(
                (source_value // G04) * (target_value // G04)
                for source_value, target_value
                in zip(target_source, target_target)
            ) if cylinder_attached else 0
        )
        raw_target_rows.append((
            target_label,
            tuple(index for index, value in enumerate(target_source) if value),
            cylinder_attached,
            raw_amplitude // G04,
        ))
        attached_target_profile.append(attached_amplitude)
        attached_cross_profile.append(cross_amplitude)
    raw_target_rows = tuple(raw_target_rows)
    attached_target_profile = tuple(attached_target_profile)
    attached_cross_profile = tuple(attached_cross_profile)
    require(
        tuple(row[0] for row in raw_target_rows if row[1])
        == tuple(range(2, 12))
        and tuple(row[0] for row in raw_target_rows if row[2])
        == tuple(range(3, 12))
        and attached_target_profile == EXPECTED_ATTACHED_TARGET_PROFILE
        and attached_cross_profile == EXPECTED_ATTACHED_CROSS_PROFILE,
        "raw or cylinder-attached target support changed",
    )
    target_character_coordinates = tuple(
        primitive_target_coordinates(attached_target_profile, frequency)
        for frequency in range(1, P)
    )
    require(
        all(any(coordinate) for coordinate in target_character_coordinates),
        "a primitive target character vanished",
    )
    cross_character_coordinates = tuple(
        primitive_target_coordinates(attached_cross_profile, frequency)
        for frequency in range(1, P)
    )
    require(
        all(any(coordinate) for coordinate in cross_character_coordinates),
        "a primitive cross character vanished",
    )
    target_character_certificates = tuple(
        (
            frequency,
            next(index for index, coordinate in enumerate(coordinates)
                 if coordinate),
            next(coordinate for coordinate in coordinates if coordinate),
        )
        for frequency, coordinates in enumerate(
            target_character_coordinates, start=1
        )
    )
    cross_character_certificates = tuple(
        (
            frequency,
            next(index for index, coordinate in enumerate(coordinates)
                 if coordinate),
            next(coordinate for coordinate in coordinates if coordinate),
        )
        for frequency, coordinates in enumerate(
            cross_character_coordinates, start=1
        )
    )

    # Sharp controls.  The t=0 full section is definitionally disjoint from
    # E3.  Conversely, retaining only each endpoint's own F_(ell,0,4), rather
    # than the common source/translated-target section, keeps positive mass
    # but destroys raw coefficient transport.
    zero_sections = semantic_sections(module, source_e3, 0, 0)
    zero_source, zero_target, _zero_details = fully_marked_vectors(
        module, rails, present, prefixes, zero_sections, 8
    )
    require(zero_source == zero_target == (0,) * 7,
            "t=0 E3 hostile unexpectedly survived")
    own_source, own_target, _own_details = fully_marked_vectors(
        module, rails, present, prefixes, sections, 8,
        common_sections=False, equal_weight_only=True,
    )
    expected_own_source = (
        0, 339633525654239542165440, 750593782703678965571520,
        722054095148406001101120, 0, 0, 0,
    )
    expected_own_target = (
        0, 345341652135823400016960, 756301720214733558465600,
        724908063903933297548160, 0, 0, 0,
    )
    require(
        own_source == expected_own_source
        and own_target == expected_own_target
        and own_source != own_target,
        "own-section/no-common hostile changed",
    )

    print("LRC14 FULLY MARKED ROOT-ZERO CLUTCH AUDIT")
    print("status=PROVED+VERIFIED-EXACT partial clock-indexed semantic clutch")
    print(f"R={R} S={S} T={T} tau=7/{R} tau_grid={SHIFT} "
          "delayed_shift=7 carry_shift=7/13 deep_shift=14")
    print("marked_carrier=both relative-present complements + both "
          "E3 intersect F_(ell,0,4) sections + R12/L1 overlap + "
          "sector0,h6,kappa1,D^-6 Q_(3,{1,2})")
    print(f"uniform_label={LABEL} exact_subbank_content={content} "
          "v13=1 equal_weight_raw_equality=14/14 units=14/14")
    print(f"rail_table={tuple(rail_table)}")
    print(f"rail_vectors={vectors}")
    print(f"rail8_direct_lifted_hits={tuple(direct_hits)} "
          f"rail8_vector={vectors[8]}")
    print(f"normalized_sign=source_root12:-target_root1 "
          f"rail8_det={EXPECTED_DETERMINANTS[8]}")
    print(f"whole_cylinder_witness=ell1 {q_source}->{q_target} "
          f"radius={clutch.relative.Q_RADIUS} label={LABEL}")
    print(f"raw_target_rows={raw_target_rows}")
    print(f"cylinder_attached_target_profile_over_G04="
          f"{attached_target_profile}")
    print(f"cylinder_attached_cross_profile_over_G04_squared="
          f"{attached_cross_profile}")
    print(f"primitive_target_character_certificates="
          f"{target_character_certificates}")
    print(f"primitive_cross_character_certificates="
          f"{cross_character_certificates}")
    print("target_profile_scope=raw integral support is t=2..11; "
          "the explicit transported-cylinder attachment removes t=2 and "
          "has support t=3..11")
    print(f"hostile_t0={zero_source}")
    print(f"hostile_own_section_source={own_source}")
    print(f"hostile_own_section_target={own_target}")
    print("transport_lemma=on any common translated weighted carrier, "
          "R*tau is integral and carry12 maps to carry6, so raw equality "
          "holds for every retained delayed target predicate")
    print("clock_scope=the seven entries use the family "
          "E3 intersect F_(ell,0,4); the explicit point is in ell1 only; "
          "the other entries are counterfactual fibres, not an orbit")
    print("NEXT: insert this carrier idempotent into one R=13^6 THM-2334 "
          "marked triangle and retain exact relation addresses before the "
          "THM-2625 Radon quotient")
    print("SCOPE: partial derived equal-weight Cech clutch; no global action, "
          "THM-2334 endpoint current, row exclusion, or LRC14")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
