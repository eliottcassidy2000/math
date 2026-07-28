#!/usr/bin/env python3
"""Exact definition audit for THM-2720's global present-wall no-go.

The canonical present packet F_(ell,s) has an unshifted c3-safe factor.
The exclusive source E3 has the complementary unshifted c3-danger factor.
This companion reconstructs all 7*13 canonical present packets, checks the
factorization before/after the c3-safe cut, and supplies deletion controls.
The theorem itself is the set-theoretic factor contradiction, not a finite
extrapolation from this census.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

CARRIER_PATH = (
    ROOT / "04-computation"
    / "lrc14_replica_dichotomy_typed_row_opus_20260727.py"
)
CARRIER_SHA256 = (
    "6ba64a68a9fd008d2e06949b1f1cf75012f1f4e734f75f55ce0af58ae20ad7b9"
)
if hashlib.sha256(CARRIER_PATH.read_bytes()).hexdigest() != CARRIER_SHA256:
    raise RuntimeError("audited canonical present-packet dependency changed")

import lrc14_replica_dichotomy_typed_row_opus_20260727 as carrier


P = 13
Q7 = 7


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def intersect_sorted(left, right):
    """Intersection of two sorted disjoint half-open interval unions."""
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        start = max(left[i][0], right[j][0])
        end = min(left[i][1], right[j][1])
        if start < end:
            result.append((start, end))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return result


def interval_mass(intervals):
    return sum(right - left for left, right in intervals)


def exclusive_e3():
    """Canonical E3: common-safe, c3-danger, c1/c2-safe."""
    intervals = carrier.make_comb(
        carrier.C3, 182, -13, 13
    )
    intervals = carrier.subtract_comb(
        intervals, carrier.W[carrier.GUARD], 91, -13, 13
    )
    for index in carrier.UNIT_IDX:
        intervals = carrier.subtract_comb(
            intervals, carrier.W[index], 182, -13, 13
        )
    intervals = carrier.subtract_comb(
        intervals, carrier.C1, 182, -13, 13
    )
    intervals = carrier.subtract_comb(
        intervals, carrier.C2, 182, -13, 13
    )
    return intervals


def present_before_c3_safe(ell, shift):
    """The literal build_F definition with only its final c3-safe cut omitted."""
    intervals = carrier.make_comb(
        carrier.C1, 182, 26 * ell - 13, 26 * ell + 13
    )
    intervals = carrier.subtract_comb(
        intervals, carrier.W[carrier.GUARD], 91, -13, 13
    )
    intervals = carrier.subtract_comb(
        intervals, carrier.W[1], 182,
        -14 * shift - 13, -14 * shift + 13,
    )
    for index in (2, 3, 4, 5):
        intervals = carrier.subtract_comb(
            intervals, carrier.W[index], 182, -13, 13
        )
    intervals = carrier.subtract_comb(
        intervals, carrier.C2, 182,
        14 * shift - 13, 14 * shift + 13,
    )
    return intervals


def main():
    danger_c3 = carrier.make_comb(carrier.C3, 182, -13, 13)
    safe_c3 = carrier.subtract_comb(
        [(0, carrier.T)], carrier.C3, 182, -13, 13
    )
    require(not intersect_sorted(danger_c3, safe_c3),
            "c3 danger and its safe complement overlap")
    require(interval_mass(danger_c3) + interval_mass(safe_c3) == carrier.T,
            "c3 danger/safe partition lost mass")

    e3 = exclusive_e3()
    require(e3 and not intersect_sorted(e3, safe_c3),
            "exclusive E3 left the c3-danger factor")

    deleted_positive = []
    present_masses = []
    deleted_intersection_masses = []
    for ell in range(Q7):
        for shift in range(P):
            precursor = present_before_c3_safe(ell, shift)
            expected = carrier.subtract_comb(
                precursor, carrier.C3, 182, -13, 13
            )
            present = carrier.build_F(ell, shift)
            require(present == expected,
                    f"build_F factorization changed at {(ell, shift)}")
            require(not intersect_sorted(present, danger_c3),
                    f"present packet is not c3-safe at {(ell, shift)}")
            require(not intersect_sorted(present, e3),
                    f"present packet acquired E3 at {(ell, shift)}")
            present_masses.append(interval_mass(present))

            deleted_meet = intersect_sorted(precursor, e3)
            deleted_mass = interval_mass(deleted_meet)
            deleted_intersection_masses.append(deleted_mass)
            if deleted_mass:
                deleted_positive.append((ell, shift, deleted_mass))

    require(len(present_masses) == Q7 * P
            and min(present_masses) > 0,
            "canonical present bank acquired an empty packet")
    require(deleted_positive,
            "deleting only c3-safe did not produce a sharp control")

    print("THM-2720 UNSHIFTED DEEPEST-SOURCE/PRESENT-WALL AUDIT")
    print("row=(1,14,27,40,53,66,13,2197,742586)")
    print(f"canonical_present_dependency_sha256={CARRIER_SHA256}")
    print("identity=F_(ell,s)=F_before_c3_safe_(ell,s) intersect Safe(c3)")
    print("source_identity=E3 subset Danger(c3)")
    print(f"canonical_present_packets={Q7 * P} all_E3_intersections=0")
    print(f"E3_components={len(e3)} E3_grid_mass={interval_mass(e3)}")
    print("present_grid_mass_range="
          f"({min(present_masses)},{max(present_masses)})")
    print("delete_only_c3_safe_positive_pairs="
          f"{len(deleted_positive)}/{Q7 * P}")
    print("delete_only_c3_safe_intersection_mass_range="
          f"({min(mass for mass in deleted_intersection_masses if mass)},"
          f"{max(deleted_intersection_masses)})")
    print(f"first_deletion_control={deleted_positive[0]}")
    print("uniformity=any further rail/clock/carry/root/unit restriction "
          "remains a subset of F_(ell,s), so cannot restore E3")
    print("SCOPE: canonical typed row and unchanged build_F present grammar; "
          "no changed-present carrier, row exclusion, or LRC14 conclusion")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
