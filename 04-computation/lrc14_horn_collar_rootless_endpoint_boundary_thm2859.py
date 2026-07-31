#!/usr/bin/env python3
"""Exact structural explanation for the common-only Z-orbit exclusion.

This tiny companion imports the frozen full-forest audit, reconstructs only
the 98 common-only components, and de-duplicates their 8,554 labelled atoms
to 286 physical intervals.  Each source/target endpoint mask is excluded
from the vertical orbit of the q0 rectangle at its first failed invariant:

1. mass is not 81;
2. mass is 81 but the horizontal factor is wrong;
3. both agree, but the vertical 9-subset is not a translate of the q0 one.

No executable Python ``assert`` statement is used.
"""

from collections import Counter
from hashlib import sha256
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MAIN_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_horn_collar_endpoint_orbit_action_thm2859.py"
)
MAIN_SHA256 = (
    "a4c145892ec3caf1f199a13bac07472e726f5dada0b523e86274ad8c14ce2846"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


require(
    sha256(MAIN_PATH.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    == MAIN_SHA256,
    "frozen full-forest orbit audit changed",
)
spec = importlib.util.spec_from_file_location("z8_full_audit", MAIN_PATH)
full = importlib.util.module_from_spec(spec)
spec.loader.exec_module(full)


EXPECTED_VERTICAL_ROWS = (
    ((0, 1, 4, 5, 8, 9, 10, 11, 12), 11, 10),
    ((0, 1, 4, 7, 8, 9, 10, 11, 12), 23, 23),
    ((0, 1, 6, 7, 8, 9, 10, 11, 12), 23, 23),
)


def frame_audit(intervals, shift):
    base = full.endpoint.endpoint_mask(full.I)
    first_axis, second_axis, product = full.endpoint.cartesian_factor(base)
    require(product and sum(base) == 81, "q0 rectangle changed")
    orbit = frozenset(
        full.endpoint.translate_mask(base, (0, exponent))
        for exponent in range(full.P)
    )
    rows = []
    failure = Counter()
    vertical_subsets = Counter()
    descriptions = Counter()
    masses = Counter()
    for interval in intervals:
        shifted = (
            interval[0] + shift,
            interval[1] + shift,
        )
        mask = full.endpoint.endpoint_mask(shifted)
        description = full.endpoint.mask_description(mask)
        horizontal, vertical, mass, rectangle = description
        require(rectangle, "common-only endpoint mask ceased to be Cartesian")
        rows.append((interval, shifted, description))
        descriptions[description] += 1
        masses[mass] += 1
        if mass != 81:
            failure["wrong_mass"] += 1
        elif horizontal != tuple(sorted(first_axis)):
            failure["wrong_horizontal_factor"] += 1
        elif mask not in orbit:
            failure["wrong_vertical_subset"] += 1
            vertical_subsets[vertical] += 1
        else:
            failure["orbit_hit"] += 1
    return {
        "base_first": tuple(sorted(first_axis)),
        "base_second": tuple(sorted(second_axis)),
        "rows": tuple(rows),
        "unique_masks": len(descriptions),
        "masses": masses,
        "failure": failure,
        "vertical_subsets": vertical_subsets,
        "description_digest": sha256(
            repr(tuple(sorted(descriptions.items(), key=repr))).encode("ascii")
        ).hexdigest(),
        "row_digest": sha256(repr(tuple(rows)).encode("ascii")).hexdigest(),
    }


def main():
    paths, _cells, _module, _e3, _clocks, kinds = full.build_paths()
    common_only = tuple(
        path for path in paths if path["kind"] == "common_only"
    )
    intervals = tuple(sorted({
        piece[:2]
        for path in common_only
        for piece in path["vertices"]
    }))
    require(
        kinds == Counter({"cofiber_rooted": 587, "common_only": 98})
        and len(common_only) == 98
        and sum(len(path["vertices"]) for path in common_only) == 8_554
        and len(intervals) == 286,
        "common-only universe changed",
    )

    source = frame_audit(intervals, 0)
    target = frame_audit(intervals, full.SHIFT)
    expected_source_vertical = Counter({
        row[0]: row[1] for row in EXPECTED_VERTICAL_ROWS
    })
    expected_target_vertical = Counter({
        row[0]: row[2] for row in EXPECTED_VERTICAL_ROWS
    })
    require(
        source["base_first"] == target["base_first"]
        == (0, 1, 2, 3, 4, 5, 6, 7, 12)
        and source["base_second"] == target["base_second"]
        == (0, 1, 2, 3, 4, 5, 8, 9, 10)
        and source["unique_masks"] == target["unique_masks"] == 17,
        "base or common-only mask-type census changed",
    )
    require(
        source["masses"]
        == Counter({81: 68, 90: 114, 99: 29, 100: 46, 110: 29})
        and target["masses"]
        == Counter({81: 66, 90: 115, 99: 29, 100: 47, 110: 29}),
        "common-only mass census changed",
    )
    require(
        source["failure"]
        == Counter({
            "wrong_mass": 218,
            "wrong_horizontal_factor": 11,
            "wrong_vertical_subset": 57,
        })
        and target["failure"]
        == Counter({
            "wrong_mass": 220,
            "wrong_horizontal_factor": 10,
            "wrong_vertical_subset": 56,
        })
        and source["vertical_subsets"] == expected_source_vertical
        and target["vertical_subsets"] == expected_target_vertical,
        "common-only first-failure stratification changed",
    )

    print("THM-2859 ROOTLESS ENDPOINT-ORBIT BOUNDARY AUDIT")
    print(
        "universe=components:98,labelled_atoms:8554,"
        f"physical_intervals:{len(intervals)};"
        f"base_A={source['base_first']};base_B={source['base_second']}"
    )
    for name, data in (("source", source), ("target", target)):
        print(
            f"{name}_boundary="
            f"unique_masks={data['unique_masks']};"
            f"masses={tuple(sorted(data['masses'].items()))};"
            f"first_failure={tuple(sorted(data['failure'].items()))};"
            f"surviving_vertical_subsets="
            f"{tuple(sorted(data['vertical_subsets'].items()))};"
            f"description_sha256={data['description_digest']};"
            f"typed_rows_sha256={data['row_digest']}"
        )
    print(
        "MECHANISM=mass removes 218 source/220 target intervals; "
        "the horizontal factor removes 11/10 more; the remaining 57/56 "
        "have one of exactly three vertical 9-subsets, none a cyclic "
        "translate of the q0 vertical factor.  Hence common-only orbit "
        "hits are exactly zero in both endpoint frames."
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
