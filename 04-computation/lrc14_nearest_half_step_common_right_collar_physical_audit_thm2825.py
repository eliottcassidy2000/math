#!/usr/bin/env python3
"""Full physical sidecar companion for THM-2825.

This script independently retains native-factor, 13-twist carrier, and all
169 endpoint-address masks for the +h and +2h common/right collars.  The
primary companion handles the faster geometry and semantic-parity census.
No Python ``assert`` is used.
"""

from collections import Counter
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "lrc14_right_cofiber_positive_copy_stratification_thm2818.py":
        "85edac9bb03f1fef198343268f4fb1226bec122d45ded79a049f8fa9a73882a8",
}
for name, expected in PINNED.items():
    actual = sha256(
        (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    require(actual == expected, f"pinned dependency changed: {name}")


import lrc14_right_cofiber_positive_copy_stratification_thm2818 as copies


P = 13
H = copies.T // (2 * P**5)
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)


def cell_objects(details, full_module, e3, clocks, clock, sigma, target):
    section = copies.physical.target.source_present_section(
        full_module, e3, clock, sigma, target, clocks
    )
    source_base, target_base = details[clock]
    source = copies.weighted_intersection(source_base, section)
    target = copies.weighted_intersection(target_base, section)
    pulled = copies.physical.overlap.shift_weighted(target, -copies.SHIFT)
    aligned = copies.physical.overlap.intersect_weighted_profiles(source, pulled)
    require(all(a == b for _l, _r, a, b in aligned), "unequal common weight")
    common = tuple((left, right, a) for left, right, a, _b in aligned)
    right = copies.physical.subtract_weighted(pulled, common)
    return source, target, common, right


def section_factors(full_module, e3, clocks, clock, sigma, target):
    universe = ((0, full_module.T),)
    return (
        e3,
        clocks[clock],
        full_module.subtract_comb(
            universe, full_module.W[1], 182,
            -14 * sigma - 13, -14 * sigma + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.W[2], 182,
            -14 * target - 13, -14 * target + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.C2, 182,
            14 * sigma - 13, 14 * sigma + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.C3, 182,
            14 * target - 13, 14 * target + 13,
        ),
    )


def endpoint_row(interval, present_sets, cache):
    key = interval[:2]
    if key not in cache:
        source = copies.endpoint_mask(key, present_sets)
        target_interval = tuple(endpoint + copies.SHIFT for endpoint in key)
        target = copies.endpoint_mask(target_interval, present_sets)
        cache[key] = (source, target)
    return cache[key]


def endpoint_summary(row):
    return (
        row[0].count(0), row[0].count(1),
        row[1].count(0), row[1].count(1),
    )


def pair_endpoint_summary(left, right):
    source_distance = sum(a != b for a, b in zip(left[0], right[0]))
    target_distance = sum(a != b for a, b in zip(left[1], right[1]))
    source_translations = copies.translated_masks(left[0], right[0])
    target_translations = copies.translated_masks(left[1], right[1])
    return (
        source_distance, target_distance,
        source_translations, target_translations,
    )


def carrier_supports(source, target):
    unit = copies.T // P
    source_supports = tuple(
        copies.support_of(copies.physical.overlap.shift_weighted(
            source, -twist * unit
        ))
        for twist in range(P)
    )
    target_supports = tuple(
        copies.support_of(copies.physical.overlap.shift_weighted(
            target, twist * unit
        ))
        for twist in range(P)
    )
    return source_supports, target_supports


def carrier_mask(interval, supports):
    target_interval = tuple(endpoint + copies.SHIFT for endpoint in interval[:2])
    return (
        tuple(bool(copies.intersect_sorted((interval[:2],), support))
              for support in supports[0]),
        tuple(bool(copies.intersect_sorted((target_interval,), support))
              for support in supports[1]),
    )


def main():
    (
        _module, _rails, _present, details, full_module, e3, clocks,
        _q_pairs, _delayed, _sw, _tw, _rail,
    ) = copies.physical_setup()
    present_sets = copies.endpoint_base.present_cache()
    endpoint_cache = {}

    factor_census = {name: Counter() for name in ("R", "M1", "M2")}
    carrier_census = {name: Counter() for name in ("R", "M1", "M2")}
    endpoint_census = {name: Counter() for name in ("R", "M1", "M2")}
    endpoint_pair_census = {name: Counter() for name in ("h", "2h")}
    pairs = 0

    for clock in range(7):
        for sigma in COMMON_S:
            for target_label in COMMON_T:
                source, target, common, right = cell_objects(
                    details, full_module, e3, clocks,
                    clock, sigma, target_label,
                )
                require(bool(common) == bool(right), "support shadows split")
                if not right:
                    continue
                common_by_left = {row[0]: row for row in common}
                factors = section_factors(
                    full_module, e3, clocks, clock, sigma, target_label
                )
                supports = carrier_supports(source, target)
                for rpiece in right:
                    m1 = common_by_left.get(rpiece[0] + H)
                    m2 = common_by_left.get(rpiece[0] + 2 * H)
                    require(m1 is not None and m2 is not None, "collar piece missing")
                    rows = {"R": rpiece, "M1": m1, "M2": m2}
                    endpoint_rows = {}
                    for name, piece in rows.items():
                        factor_census[name][copies.factor_masks(
                            piece[:2], factors
                        )] += 1
                        carrier_census[name][carrier_mask(piece, supports)] += 1
                        endpoint_rows[name] = endpoint_row(
                            piece, present_sets, endpoint_cache
                        )
                        endpoint_census[name][endpoint_summary(
                            endpoint_rows[name]
                        )] += 1
                    endpoint_pair_census["h"][pair_endpoint_summary(
                        endpoint_rows["R"], endpoint_rows["M1"]
                    )] += 1
                    endpoint_pair_census["2h"][pair_endpoint_summary(
                        endpoint_rows["R"], endpoint_rows["M2"]
                    )] += 1
                    pairs += 1

    print("LRC HALF-STEP COLLAR PHYSICAL SIDECAR SCOUT")
    print(f"half_step={H};pairs={pairs};endpoint_cache={len(endpoint_cache)}")
    for name in ("R", "M1", "M2"):
        print(f"factor_census_{name}={tuple(sorted(factor_census[name].items(), key=repr))}")
        print(f"carrier_census_{name}={tuple(sorted(carrier_census[name].items(), key=repr))}")
        print(f"endpoint_census_{name}={tuple(sorted(endpoint_census[name].items()))}")
    print(f"endpoint_pair_census_h={tuple(sorted(endpoint_pair_census['h'].items(), key=repr))}")
    print(f"endpoint_pair_census_2h={tuple(sorted(endpoint_pair_census['2h'].items(), key=repr))}")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
