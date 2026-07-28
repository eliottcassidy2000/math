#!/usr/bin/env python3
"""Exact natural single-sheet target-label bank at the root-zero clutch.

The computation retains the canonical rail-8 adjacent-chart overlap, the
complete E3 source mask, the delayed Q_(3,{1,2}) prefix, and the lawful
two-target sheet U_(s,t).  Source carry 12/root 12 and target carry 6/root 1
are evaluated independently.  All arithmetic is integral/rational.

The proved THM-2749 uses a two-sided common section and is a different
carrier.  This one-sided naturality sidecar is reflection-grade pending audit.
In particular, a nonzero coefficient here is not the global scalar-cover
endpoint current.
"""

from __future__ import annotations

from collections import Counter
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


DEPENDENCIES = {
    "lrc14_root_zero_overlap_clutch_20260728.py":
        "e10fa7c9a5a238461ef422ea314dc334f7e65ec1787cf65d4e4bea12b96aefb8",
    "lrc14_two_target_present_semantic_attachment_probe_20260728.py":
        "062b352f4db12a5f01822b293cdbb10629632dacc5fa27b406d8dd321e550709",
}
for dependency, expected in DEPENDENCIES.items():
    actual = hashlib.sha256(
        (COMP / dependency).read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    require(actual == expected,
            f"fully marked dependency changed: {dependency}")

import lrc14_root_zero_overlap_clutch_20260728 as clutch
import lrc14_two_target_present_semantic_attachment_probe_20260728 as two


relative = clutch.relative
old = relative.private.old
P = 13
CONTENT = clutch.CONTENT
SOURCE_ROOT = clutch.SOURCE_STATE[3]
TARGET_ROOT = clutch.TARGET_STATE[3]
SOURCE_CARRY = clutch.SOURCE_STATE[2]
TARGET_CARRY = clutch.TARGET_STATE[2]
ENDPOINT_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
ENDPOINT_T = tuple(range(3, 12))
CLOCKED_S = (0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12)


def build_q3_pair_prefixes(module):
    """Delayed sector-zero prefixes with the full Q_(3,{1,2}) word."""
    word = relative.private.prior.sector_words(module)[0]
    by_clock = []
    for ell in range(7):
        qell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        qell = module.intersect_comb(qell, module.C1, 182, -13, 13)
        by_h = []
        for h in range(P):
            pair = []
            for kappa in range(2):
                digit = old.sat.intersect_interval(
                    qell,
                    (2 * h + kappa) * relative.lift.m.T // (2 * P),
                    (2 * h + kappa + 1) * relative.lift.m.T // (2 * P),
                )
                pair.append(module.make_prefix(digit))
            by_h.append(tuple(pair))
        by_clock.append(tuple(by_h))
    return tuple(by_clock)


def cyclotomic_coordinates(values, frequency):
    """Power-basis coordinates of sum_t values[t] zeta_13^(frequency*t)."""
    require(len(values) == P and 0 < frequency < P,
            "invalid primitive target character")
    coordinates = [0] * (P - 1)
    for t, value in enumerate(values):
        exponent = frequency * t % P
        if exponent == P - 1:
            for slot in range(P - 1):
                coordinates[slot] -= value
        else:
            coordinates[exponent] += value
    return tuple(coordinates)


def scalar_amplitude(vector):
    """Return a when vector=(0,a,...,a), otherwise None."""
    if vector[0] == 0 and len(set(vector[1:])) == 1:
        return vector[1]
    return None


def normalized_scalar(vector, root):
    """Root-normalized scalar profile when the clock polynomial is constant."""
    require(all(value % CONTENT == 0 for value in vector),
            "coefficient left the content-26 lattice")
    profile = clutch.normalized_profile(vector, root)[1]
    require(profile[1:] == (0,) * 5,
            "clock profile stopped being scalar")
    return profile[0]


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

    source_base = {}
    target_base = {}
    require(
        all(prefix[2][-1] == 0 for prefix in pair_prefixes[0][6])
        and all(
            pair_prefixes[ell][6] == pair_prefixes[1][6]
            for ell in range(1, 7)
        ),
        "sector-zero Q_(3,{1,2}) clock law changed",
    )
    phi_cache = {}
    for ell in range(7):
        source_base[ell], target_base[ell] = clutch.restrict_to_relative_overlap(
            module, present, rail_pairs, ell
        )
    exclusive_e3 = two.exclusive_source(module, 3)
    source_clock_combs = tuple(
        module.make_comb(module.C1, 182, 26 * ell - 13, 26 * ell + 13)
        for ell in range(7)
    )

    rows = {}
    for s in range(P):
        for t in range(P):
            # This exact constructor is the MISTAKE-313 repair: unlike the
            # legacy helper, it includes the source-one clock d_(1,1).
            section = tuple(two.source_present_section(
                module, exclusive_e3, 1, s, t, source_clock_combs
            ))
            section_starts = tuple(left for left, _right in section)
            source_pieces = tuple(
                old.intersect_weighted_union(
                    source_base[ell], section, section_starts
                )
                for ell in range(7)
            )
            target_pieces = tuple(
                old.intersect_weighted_union(
                    target_base[ell], section, section_starts
                )
                for ell in range(7)
            )
            require(
                all(source_pieces[ell] == source_pieces[1]
                    for ell in range(7))
                and all(target_pieces[ell] == target_pieces[1]
                        for ell in range(7)),
                f"clocked natural carrier depends on coefficient clock at {(s,t)}",
            )
            source = source_pieces[1]
            target = target_pieces[1]
            source_mass_value = sum(
                (right - left) * weight
                for left, right, weight in source
            )
            target_mass_value = sum(
                (right - left) * weight
                for left, right, weight in target
            )
            source_values = relative.private.delayed_carry_pair(
                source, pair_prefixes[1][6], phi_cache
            )
            target_values = relative.private.delayed_carry_pair(
                target, pair_prefixes[1][6], phi_cache
            )
            source_vector = (0,) + (source_values[SOURCE_CARRY][1],) * 6
            target_vector = (0,) + (target_values[TARGET_CARRY][1],) * 6
            source_mass = (source_mass_value,) * 7
            target_mass = (target_mass_value,) * 7
            source_amplitude = scalar_amplitude(source_vector)
            target_amplitude = scalar_amplitude(target_vector)
            require(source_amplitude is not None and target_amplitude is not None,
                    f"non-scalar clock vector at {(s, t)}")
            source_scalar = normalized_scalar(source_vector, SOURCE_ROOT)
            target_scalar = normalized_scalar(target_vector, TARGET_ROOT)
            rows[s, t] = {
                "source_vector": source_vector,
                "target_vector": target_vector,
                "source_mass": source_mass,
                "target_mass": target_mass,
                "source_amplitude": source_amplitude,
                "target_amplitude": target_amplitude,
                "source_scalar": source_scalar,
                "target_scalar": target_scalar,
            }

    source_support = tuple(
        key for key in sorted(rows) if rows[key]["source_amplitude"]
    )
    target_support = tuple(
        key for key in sorted(rows) if rows[key]["target_amplitude"]
    )
    common_support = tuple(
        key for key in sorted(rows)
        if rows[key]["source_amplitude"] and rows[key]["target_amplitude"]
    )
    require(
        source_support == tuple(
            (s, t) for s in CLOCKED_S for t in ENDPOINT_T
        )
        and target_support == tuple(
            (s, t) for s in CLOCKED_S for t in range(3, 13)
        )
        and common_support == source_support,
        "source-clocked natural support laws changed",
    )
    endpoint_support = tuple(
        (s, t) for s in ENDPOINT_S for t in ENDPOINT_T
        if rows[s, t]["source_amplitude"]
        and rows[s, t]["target_amplitude"]
    )
    require(len(endpoint_support) == 81,
            "canonical endpoint window lost a clocked natural sheet")

    # The source and target witnesses have the same exact lawful support.
    # Check every primitive t-character on each fixed-s row and on their raw
    # cross product.  The latter is the finite coefficient proposed after the
    # one-sheet shear was found.
    source_characters = 0
    target_characters = 0
    cross_characters = 0
    per_s_character_counts = []
    ratios = set()
    common_ratio_histogram = Counter()
    endpoint_ratio_histogram = Counter()
    for s in range(P):
        source_profile = tuple(rows[s, t]["source_amplitude"] for t in range(P))
        target_profile = tuple(rows[s, t]["target_amplitude"] for t in range(P))
        cross_profile = tuple(
            source_profile[t] * target_profile[t] for t in range(P)
        )
        counts = [0, 0, 0]
        for frequency in range(1, P):
            for slot, profile in enumerate(
                (source_profile, target_profile, cross_profile)
            ):
                if any(cyclotomic_coordinates(profile, frequency)):
                    counts[slot] += 1
        source_characters += counts[0]
        target_characters += counts[1]
        cross_characters += counts[2]
        per_s_character_counts.append((s, *counts))
        for t in range(P):
            row = rows[s, t]
            if row["source_scalar"]:
                ratio = (row["target_scalar"]
                         * pow(row["source_scalar"], -1, P) % P)
                ratios.add(ratio)
                common_ratio_histogram[ratio] += 1
        if s in ENDPOINT_S:
            for t in ENDPOINT_T:
                source_scalar = rows[s, t]["source_scalar"]
                target_scalar = rows[s, t]["target_scalar"]
                endpoint_ratio_histogram[
                    target_scalar * pow(source_scalar, -1, P) % P
                ] += 1

    endpoint_character_counts = [0, 0, 0]
    endpoint_cross_bank = {
        (s, t): (
            rows[s, t]["source_amplitude"]
            * rows[s, t]["target_amplitude"]
            if s in ENDPOINT_S and t in ENDPOINT_T else 0
        )
        for s in range(P) for t in range(P)
    }
    for s in ENDPOINT_S:
        profiles = (
            tuple(
                rows[s, t]["source_amplitude"] if t in ENDPOINT_T else 0
                for t in range(P)
            ),
            tuple(
                rows[s, t]["target_amplitude"] if t in ENDPOINT_T else 0
                for t in range(P)
            ),
            tuple(endpoint_cross_bank[s, t] for t in range(P)),
        )
        for profile_index, profile in enumerate(profiles):
            endpoint_character_counts[profile_index] += sum(
                any(cyclotomic_coordinates(profile, frequency))
                for frequency in range(1, P)
            )

    # A complete 2D target-character check records whether the bank supplies
    # genuinely mixed (s,t) phase data, rather than only separate rowwise t
    # profiles.  A character has the SINGLE cyclotomic value
    # zeta^(fs*s+ft*t); it is not a tensor of two independent cyclotomic
    # fields.  Bucket by that exponent and reduce once modulo Phi_13.
    mixed_nonzero = 0
    mixed_total = (P - 1) ** 2
    cross_bank = {
        (s, t): rows[s, t]["source_amplitude"]
        * rows[s, t]["target_amplitude"]
        for s in range(P) for t in range(P)
    }
    for fs in range(1, P):
        for ft in range(1, P):
            buckets = [0] * P
            for s in range(P):
                for t in range(P):
                    buckets[(fs * s + ft * t) % P] += cross_bank[s, t]
            if len(set(buckets)) > 1:
                mixed_nonzero += 1

    # Retain the originally tempting tensor-field count only as a rejected
    # diagnostic.  It is not a C13xC13 character calculation because both
    # target coordinates use the same primitive thirteenth root.
    rejected_tensor_count = 0
    endpoint_mixed_nonzero = 0
    endpoint_rejected_tensor_count = 0
    for fs in range(1, P):
        for ft in range(1, P):
            t_coordinates_by_s = tuple(
                cyclotomic_coordinates(
                    tuple(cross_bank[s, t] for t in range(P)), ft
                )
                for s in range(P)
            )
            coordinate_rows = tuple(
                cyclotomic_coordinates(
                    tuple(t_coordinates_by_s[s][slot] for s in range(P)), fs
                )
                for slot in range(P - 1)
            )
            rejected_tensor_count += any(
                value for row in coordinate_rows for value in row
            )
            endpoint_buckets = [0] * P
            for s in range(P):
                for t in range(P):
                    endpoint_buckets[(fs * s + ft * t) % P] += (
                        endpoint_cross_bank[s, t]
                    )
            endpoint_mixed_nonzero += len(set(endpoint_buckets)) > 1
            endpoint_t_coordinates_by_s = tuple(
                cyclotomic_coordinates(
                    tuple(endpoint_cross_bank[s, t] for t in range(P)), ft
                )
                for s in range(P)
            )
            endpoint_coordinate_rows = tuple(
                cyclotomic_coordinates(
                    tuple(
                        endpoint_t_coordinates_by_s[s][slot]
                        for s in range(P)
                    ),
                    fs,
                )
                for slot in range(P - 1)
            )
            endpoint_rejected_tensor_count += any(
                value for row in endpoint_coordinate_rows for value in row
            )

    require(
        rows[0, 3]["source_vector"]
        == (0,) + (339633525654239542165440,) * 6
        and rows[0, 3]["target_vector"]
        == (0,) + (345341652135823400016960,) * 6
        and (rows[0, 3]["source_scalar"], rows[0, 3]["target_scalar"])
        == (9, 8),
        "one-sheet hostile control changed",
    )
    require(
        tuple(sorted(ratios)) == (6, 11)
        and tuple(sorted(common_ratio_histogram.items()))
        == ((6, 18), (11, 81))
        and tuple(sorted(endpoint_ratio_histogram.items())) == ((11, 81),),
        "clocked natural gain census changed",
    )
    require(
        source_characters == target_characters == cross_characters == 132
        and tuple(endpoint_character_counts) == (108, 108, 108)
        and mixed_nonzero == endpoint_mixed_nonzero == mixed_total == 144
        and rejected_tensor_count == endpoint_rejected_tensor_count == 144,
        "clocked natural character census changed",
    )

    print("LRC14 SOURCE-CLOCKED NATURAL ROOT-ZERO TARGET BANK")
    print(f"dependency_sha256={tuple(DEPENDENCIES.items())}")
    print("typed_state=rail8,E3,D^6,Q_(3,{1,2}),roots12->1,carries12->6")
    source_support_by_s = tuple(
        (s, tuple(t for t in range(P) if (s, t) in source_support))
        for s in range(P)
    )
    target_support_by_s = tuple(
        (s, tuple(t for t in range(P) if (s, t) in target_support))
        for s in range(P)
    )
    common_support_by_s = tuple(
        (s, tuple(t for t in range(P) if (s, t) in common_support))
        for s in range(P)
    )
    print(f"source_support_count={len(source_support)} by_s={source_support_by_s}")
    print(f"target_support_count={len(target_support)} by_s={target_support_by_s}")
    print(f"common_support_count={len(common_support)} by_s={common_support_by_s}")
    print(f"canonical_endpoint_window=S{ENDPOINT_S}xT{ENDPOINT_T} count={len(endpoint_support)}")
    print(f"normalized_target_source_ratios={tuple(sorted(ratios))}")
    print(f"common_ratio_histogram={tuple(sorted(common_ratio_histogram.items()))}")
    print(f"canonical_endpoint_ratio_histogram={tuple(sorted(endpoint_ratio_histogram.items()))}")
    print(f"per_s_(source,target,cross)_primitive_t_characters={tuple(per_s_character_counts)}")
    print(f"primitive_t_character_totals=source:{source_characters}/156,target:{target_characters}/156,cross:{cross_characters}/156")
    print(f"mixed_cross_(s,t)_characters={mixed_nonzero}/{mixed_total}")
    print(f"rejected_two_field_tensor_count={rejected_tensor_count}/{mixed_total}")
    print(f"canonical_endpoint_(source,target,cross)_primitive_t_characters={tuple(endpoint_character_counts)}/(108,108,108)")
    print(f"canonical_endpoint_mixed_cross_(s,t)_characters={endpoint_mixed_nonzero}/{mixed_total}")
    print(f"canonical_endpoint_rejected_two_field_tensor_count={endpoint_rejected_tensor_count}/{mixed_total}")
    print(f"sheet_(0,3)_source_vector={rows[0,3]['source_vector']}")
    print(f"sheet_(0,3)_target_vector={rows[0,3]['target_vector']}")
    print(f"sheet_(0,3)_normalized_scalars=({rows[0,3]['source_scalar']},{rows[0,3]['target_scalar']})")
    print("SCOPE: exact source-clock-one physical-chart coefficient bank; no two-sided common-section identification, global endpoint current, row exclusion, or LRC14 conclusion")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
