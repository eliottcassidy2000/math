#!/usr/bin/env python3
"""Exact all-81 physical-wing common-atom and factor-origin audit.

For each of the 81 lawful (s,t) labels and seven physical clocks, construct
the one-sided source/target carriers, their common carrier, and the two wing
cofibres.  Pull the target wing back to the source chart and test literal
physical co-support.  The two wings are opposite exclusive pieces, so their
intersection is empty in every cell.

The companion also disjointizes each wing by the first failed factor of the
translated full-target section.  This distinguishes the mixed A/gain-two-B
faces from the factor-pure gain-two-D corner used by the extended determinant
probe.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "lrc14_root_zero_full_target_semantic_clutch_20260728.py":
        "208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8",
    "lrc14_natural_single_sheet_root_zero_target_bank_20260728.py":
        "591241e912c452b1985c7fc700884183a6e50440a579e1d853123d276cd416a8",
}
for name, expected in PINNED.items():
    actual = sha256(
        (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    require(actual == expected, f"pinned dependency changed: {name}")

import lrc14_root_zero_full_target_semantic_clutch_20260728 as physical


P = 13
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)
ORIGIN_NAMES = (
    "shifted_E3",
    "shifted_clock",
    "shifted_q1",
    "shifted_q2",
    "shifted_c2",
    "shifted_c3",
)


def mass(pieces):
    return sum((right - left) * weight for left, right, weight in pieces)


def support(pieces):
    return physical.overlap.merge_intervals(
        (left, right) for left, right, weight in pieces if weight
    )


def weighted_intersection(pieces, intervals):
    return tuple(
        physical.relative.private.old.intersect_weighted_union(
            pieces, intervals
        )
    )


def main():
    module, _prefixes, _, _, rails, present, _starts = (
        physical.relative.lift.m.core.build_carrier_data()
    )
    pair_prefixes = physical.relative.private.build_pair_prefixes(module)
    _raw_source, _raw_target, _rail_pairs, overlap_details = (
        physical.overlap.overlap_vectors(
            module, pair_prefixes, rails, present, rail_index=8
        )
    )
    full_module = physical.target.load_present_module()
    e3 = physical.target.exclusive_source(full_module, 3)
    fork = physical.target.deepest_fork(full_module)
    clock_combs = tuple(
        full_module.make_comb(
            full_module.C1, 182, 26 * clock - 13, 26 * clock + 13
        )
        for clock in range(7)
    )
    q_pairs = physical.q_restricted_pair_prefixes(
        full_module, pair_prefixes, fork
    )
    universe = ((0, full_module.T),)

    literal_decompositions = 0
    empty_both_corners = 0
    nonempty_source_wings = 0
    nonempty_target_wings = 0
    selected = {}

    for s in COMMON_S:
        for t in COMMON_T:
            for clock, (source_base, target_base) in enumerate(
                overlap_details
            ):
                section = physical.target.source_present_section(
                    full_module, e3, clock, s, t, clock_combs
                )
                one_source = weighted_intersection(source_base, section)
                one_target = weighted_intersection(target_base, section)
                common_source = physical.common_physical_cut(
                    source_base, section, -1
                )
                common_target = physical.common_physical_cut(
                    target_base, section, +1
                )
                source_wing = physical.subtract_weighted(
                    one_source, common_source
                )
                target_wing = physical.subtract_weighted(
                    one_target, common_target
                )
                require(
                    physical.merge_weighted(
                        tuple(common_source) + tuple(source_wing)
                    ) == one_source
                    and physical.merge_weighted(
                        tuple(common_target) + tuple(target_wing)
                    ) == one_target,
                    "literal one-sided/common/wing decomposition changed",
                )
                literal_decompositions += 1

                target_pullback = physical.overlap.shift_weighted(
                    target_wing, -physical.SHIFT
                )
                common_pairs = physical.overlap.intersect_weighted_profiles(
                    source_wing, target_pullback
                )
                require(
                    not common_pairs,
                    "a physical both-wing atom unexpectedly appeared",
                )
                empty_both_corners += 1
                nonempty_source_wings += bool(support(source_wing))
                nonempty_target_wings += bool(support(target_wing))

                key = (s, t, clock)
                if key in ((8, 3, 3), (2, 3, 3), (0, 10, 3)):
                    factors = (
                        e3,
                        clock_combs[clock],
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
                    for factor in factors:
                        reconstructed = physical.target.intersect_sorted(
                            reconstructed, factor
                        )
                    require(
                        tuple(reconstructed) == tuple(section),
                        "selected full-target factorization changed",
                    )

                    def first_failures(one, direction):
                        remaining = one
                        origins = []
                        for factor in factors:
                            shifted = physical.overlap.shift_union(
                                factor, direction * physical.SHIFT
                            )
                            origins.append(
                                weighted_intersection(
                                    remaining,
                                    physical.overlap.complement(shifted),
                                )
                            )
                            remaining = weighted_intersection(
                                remaining, shifted
                            )
                        return tuple(origins), remaining

                    source_origins, source_remainder = first_failures(
                        one_source, -1
                    )
                    target_origins, target_remainder = first_failures(
                        one_target, +1
                    )
                    require(
                        source_remainder == common_source
                        and target_remainder == common_target,
                        "selected first-failure remainder changed",
                    )
                    pair = q_pairs[clock]

                    def coefficient(pieces, carry):
                        return physical.relative.private.delayed_carry_pair(
                            pieces, pair, {}
                        )[carry][1]

                    source_value = coefficient(source_wing, 12)
                    target_value = coefficient(target_wing, 6)
                    source_origin_values = tuple(
                        coefficient(origin, 12)
                        for origin in source_origins
                    )
                    target_origin_values = tuple(
                        coefficient(origin, 6)
                        for origin in target_origins
                    )
                    require(
                        sum(source_origin_values) == source_value
                        and sum(target_origin_values) == target_value,
                        "selected origin coefficients did not sum",
                    )
                    selected[key] = {
                        "source_origin_values": source_origin_values,
                        "target_origin_values": target_origin_values,
                        "source_origin_masses": tuple(
                            mass(origin) for origin in source_origins
                        ),
                        "target_origin_masses": tuple(
                            mass(origin) for origin in target_origins
                        ),
                        "source_value": source_value,
                        "target_value": target_value,
                        "source_pieces": len(source_wing),
                        "target_pieces": len(target_wing),
                        "source_mass": mass(source_wing),
                        "target_mass": mass(target_wing),
                        "chart_translate": (
                            physical.overlap.shift_weighted(
                                source_wing, physical.SHIFT
                            ) == target_wing
                        ),
                    }

    require(
        literal_decompositions == empty_both_corners == 81 * 7,
        "all-81x7 census changed",
    )
    require(
        nonempty_source_wings == nonempty_target_wings == 193,
        "physical wing-support census changed",
    )
    a = selected[(8, 3, 3)]
    b = selected[(2, 3, 3)]
    d = selected[(0, 10, 3)]
    require(
        a["source_origin_values"][0] != 0
        and all(not value for value in a["source_origin_values"][1:])
        and a["target_origin_values"][0] != 0
        and a["target_origin_values"][4] != 0
        and all(
            not value
            for index, value in enumerate(a["target_origin_values"])
            if index not in (0, 4)
        ),
        "A-face mixed-origin hostile changed",
    )
    require(
        b["source_origin_values"][0] != 0
        and all(not value for value in b["source_origin_values"][1:])
        and b["target_origin_values"][0] != 0
        and b["target_origin_values"][4] != 0,
        "gain-two B-face mixed-origin hostile changed",
    )
    require(
        d["source_origin_values"][0] != 0
        and d["target_origin_values"][0] == 2 * d["source_origin_values"][0]
        and all(not value for value in d["source_origin_values"][1:])
        and all(not value for value in d["target_origin_values"][1:])
        and not d["chart_translate"]
        and d["source_mass"] != d["target_mass"],
        "gain-two D-face pure-origin hostile changed",
    )
    require(
        a["source_origin_values"]
        == (2_853_968_755_527_296_447_040, 0, 0, 0, 0, 0)
        and a["target_origin_values"]
        == (
            5_707_937_511_054_592_894_080, 0, 0, 0,
            5_707_937_511_054_592_894_080, 0,
        )
        and d["source_origin_masses"]
        == (2_188_067_024_426_754_240, 0, 0, 0, 0, 0)
        and d["target_origin_masses"]
        == (1_458_711_349_617_836_160, 0, 0, 0, 0, 0),
        "selected exact origin values changed",
    )

    print("all-81 physical-wing common-atom/factor-origin audit")
    print(f"dependency_sha256={tuple(PINNED.items())}")
    print("source_clocked_natural_bank=fixed-e slice of the one-sided family; "
          "wing=one_sided-common")
    print(f"literal_one_sided_equals_common_plus_wing="
          f"{literal_decompositions}/567")
    print(f"physical_nonempty_wings=(source={nonempty_source_wings}/567,"
          f"target={nonempty_target_wings}/567)")
    print(f"source_wing_intersect_pullback_target_wing="
          f"empty:{empty_both_corners}/567")
    print("allocation_K4_both_wing_corner=ABSENT-PHYSICALLY")
    for name, row in (("A_8_3_3", a), ("gain2_B_2_3_3", b),
                      ("gain2_D_0_10_3", d)):
        print(f"{name}={row}")
    print(
        "CONCLUSION: independent endpoint Fourier products may form a "
        "virtual mixed face, but the physical source and target wing "
        "cofibres never occupy one common ancestry atom.  A is additionally "
        "mixed in its first-failure factor; D removes that factor ambiguity "
        "but remains disjoint, nontranslate, and mass-unequal."
    )
    print("scope=no common-wing K4 atom; one-sided/common carriers and "
          "atom-dependent allocation lifts remain unexcluded")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
