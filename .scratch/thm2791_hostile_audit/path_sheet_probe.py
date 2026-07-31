#!/usr/bin/env python3
"""Exact Boolean-path-sheet probe for the THM-2791 partial translation.

The rail-eight weight is the route-two product

    U(x) V(x-1/13),

where

    U=P_(13^5)(1_Q P_169 1_E),   V=P_(13^5)1_E.

Before marginalization, a copy of U is labelled by (a,b) in
Z/(13^5) x Z/169 and a copy of V by e' in Z/(13^5).  This companion
enumerates those literal labels from the defining interval membership
conditions.  It does not import the THM-2791 companion or use its folded
profile implementation.
"""

from bisect import bisect_right
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_base_only_bridge_opus_20260728 as base
import lrc14_b_r5_owner_clock_host_thm2581 as host
import lrc14_semantic_arm_right_wing_central_digit_thm2782 as dep


P = 13
PACK = P**2
DEPTH = P**5
ADDRESS_MODULUS = P**6
T = base.T_DEN
BASE_ADDRESS = dep.EXPECTED_BASES[1]
ADDRESS_GAP = 689364
RAIL_DISPLACEMENT = 1
SUPPLIED_PATH = (59162, 26, 56658)
INCOMING_SAME_ATOM_INTERVAL = (142004992589460, 142005019034340)

EXPECTED_DEPENDENCY_HASHES = {
    "04-computation/lrc14_base_only_bridge_opus_20260728.py":
        "99deb8891f8be1e218e7b4c5a7bad7a918ddc1da37b05b53ccbbc3ca7b8c7998",
    "04-computation/lrc14_b_r5_owner_clock_host_thm2581.py":
        "38c8eac0925d119331ef248fa7197f61bd985166c56704575cec18c361543647",
    "04-computation/lrc14_b_r5_theta_target_tensor_thm2584.py":
        "99cbb46c4cf3554fb468dda2d6bbf029f3dd7298a4f8b6009d94513484c7849d",
    "04-computation/lrc14_semantic_arm_right_wing_central_digit_thm2782.py":
        "7fbc6bb1ec303ded98eaad6e5d8205eb3d247258ada32b6f9904fc439ebb11fb",
    "04-computation/lrc14_relative_present_semantic_lift_probe_20260728.py":
        "f16754bd38ae0dfa0d7d91cc404b4447dbf359635101aa7b4223363f8064352f",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(relative_path):
    data = (ROOT / relative_path).read_bytes().replace(b"\r\n", b"\n")
    return sha256(data).hexdigest()


def scaled_intervals(intervals, denominator):
    rows = tuple(
        (denominator * left, denominator * right)
        for left, right in intervals
    )
    return rows, tuple(left for left, _right in rows)


def contains_scaled(intervals, starts, numerator):
    """Test numerator/denominator in a pre-scaled half-open interval union."""
    index = bisect_right(starts, numerator) - 1
    return index >= 0 and numerator < intervals[index][1]


def contributor_sets(
    coordinate,
    q_depth,
    q_starts,
    e_depth_pack,
    e_depth_pack_starts,
    e_depth,
    e_depth_starts,
):
    """Enumerate U's packed (a,b) labels and V's e' labels at coordinate."""
    u_labels = []
    append_u = u_labels.append
    for a in range(DEPTH):
        y_numerator = coordinate + a * T
        if not contains_scaled(q_depth, q_starts, y_numerator):
            continue
        for b in range(PACK):
            z_numerator = y_numerator + b * DEPTH * T
            if contains_scaled(
                e_depth_pack, e_depth_pack_starts, z_numerator
            ):
                append_u(a * PACK + b)

    rotated = (coordinate - RAIL_DISPLACEMENT * T // P) % T
    v_labels = []
    append_v = v_labels.append
    for e_prime in range(DEPTH):
        y_numerator = rotated + e_prime * T
        if contains_scaled(e_depth, e_depth_starts, y_numerator):
            append_v(e_prime)
    return tuple(u_labels), tuple(v_labels)


def path_digest(u_labels, v_labels):
    """Stable digest of the ordered label sets, with the two factors typed."""
    digest = sha256()
    digest.update(b"U(a*169+b):")
    for label in u_labels:
        digest.update(label.to_bytes(4, "big"))
    digest.update(b"|V(e'):")
    for label in v_labels:
        digest.update(label.to_bytes(4, "big"))
    return digest.hexdigest()


def mapped_events(intervals, multiplier, shift=0):
    """All raw contributor walls (multiplier*endpoint+shift) modulo T."""
    return tuple(sorted(set(
        (multiplier * endpoint + shift) % T
        for interval in intervals
        for endpoint in interval
    )))


def chamber(events, left, right):
    """Return consecutive raw walls strictly bracketing [left,right]."""
    index = bisect_right(events, left)
    require(index > 0 and index < len(events), "chamber wrapped unexpectedly")
    require(events[index] > right, "a contributor wall crosses the test hull")
    return events[index - 1], events[index]


def ceil_fraction(value):
    return -((-value.numerator) // value.denominator)


def main():
    actual_hashes = tuple(
        (path, lf_hash(path))
        for path in EXPECTED_DEPENDENCY_HASHES
    )
    require(
        all(
            digest == EXPECTED_DEPENDENCY_HASHES[path]
            for path, digest in actual_hashes
        ),
        "a pinned ancestry dependency changed",
    )
    require(T == 297836897838480, "base grid changed")
    require(
        dep.EXPECTED_SHIFT == 431933040,
        "physical address shift changed",
    )

    e_intervals = tuple(base.build_set(base.PAT_E, base.ZELL))
    q_intervals = tuple(base.build_set(host.PAT_QB, base.ZELL))
    q_depth, q_starts = scaled_intervals(q_intervals, DEPTH)
    e_depth_pack, e_depth_pack_starts = scaled_intervals(
        e_intervals, DEPTH * PACK
    )
    e_depth, e_depth_starts = scaled_intervals(e_intervals, DEPTH)

    half = dep.relative.Q_RADIUS * T
    cylinders = []
    for address_offset in (0, ADDRESS_GAP):
        center = dep.arm_target_center(
            BASE_ADDRESS, Fraction(address_offset, P)
        ) * T
        cylinders.append((center - half, center + half))
    cylinder_zero, cylinder_gap = cylinders

    physical_shift = Fraction(7 * ADDRESS_GAP, ADDRESS_MODULUS)
    shift_coordinate = physical_shift * T
    require(shift_coordinate.denominator == 1, "translation left the grid")
    shift_coordinate = shift_coordinate.numerator
    require(
        (
            (cylinder_zero[0] + shift_coordinate) % T,
            (cylinder_zero[1] + shift_coordinate) % T,
        )
        == cylinder_gap,
        "THM-2791 cylinders stopped translating exactly",
    )

    hull_left = min(left for left, _right in cylinders)
    hull_right = max(right for _left, right in cylinders)
    u_events = tuple(sorted(set(
        mapped_events(q_intervals, DEPTH)
        + mapped_events(e_intervals, DEPTH * PACK)
    )))
    v_events = mapped_events(
        e_intervals,
        DEPTH,
        RAIL_DISPLACEMENT * T // P,
    )
    u_chamber = chamber(u_events, hull_left, hull_right)
    v_chamber = chamber(v_events, hull_left, hull_right)
    common_chamber = (
        max(u_chamber[0], v_chamber[0]),
        min(u_chamber[1], v_chamber[1]),
    )
    require(
        common_chamber[0] < hull_left
        and hull_right < common_chamber[1],
        "selected cylinders touch an ancestry wall",
    )
    require(
        common_chamber[0] < INCOMING_SAME_ATOM_INTERVAL[0]
        and INCOMING_SAME_ATOM_INTERVAL[1] < common_chamber[1],
        "incoming same-atom interval left the common ancestry chamber",
    )

    probe_zero = ceil_fraction(cylinder_zero[0])
    require(
        cylinder_zero[0] < probe_zero < cylinder_zero[1],
        "no integral probe lies inside address-zero cylinder",
    )
    probe_gap = (probe_zero + shift_coordinate) % T
    require(
        cylinder_gap[0] < probe_gap < cylinder_gap[1],
        "translated integral probe missed the gap cylinder",
    )

    arguments = (
        q_depth,
        q_starts,
        e_depth_pack,
        e_depth_pack_starts,
        e_depth,
        e_depth_starts,
    )
    u_zero, v_zero = contributor_sets(probe_zero, *arguments)
    u_gap, v_gap = contributor_sets(probe_gap, *arguments)
    require(
        u_zero == u_gap and v_zero == v_gap,
        "literal ancestry label set changed under translation",
    )
    require(
        len(u_zero) == 966606
        and len(v_zero) == 28534
        and len(u_zero) * len(v_zero) == 27581135604,
        "ancestry census stopped factoring the rail weight",
    )

    a, b, e_prime = SUPPLIED_PATH
    require(
        a * PACK + b in u_zero and e_prime in v_zero,
        "supplied Boolean path is absent",
    )
    for coordinate in (probe_zero, probe_gap):
        arrival_root = P * coordinate // T
        deep_digit = (2 * P * coordinate // T) % P
        source_root = (arrival_root - RAIL_DISPLACEMENT) % P
        require(
            (source_root, arrival_root, deep_digit) == (5, 6, 12),
            "collision-root/deep-digit labels changed",
        )

    digest_zero = path_digest(u_zero, v_zero)
    digest_gap = path_digest(u_gap, v_gap)
    require(digest_zero == digest_gap, "ancestry digests differ")

    print("THM-2791 LITERAL BOOLEAN PATH-SHEET PROBE")
    print("status=FINITE-EXACT independent ancestry refinement")
    print("dependency_hashes=BEGIN")
    for path, digest in actual_hashes:
        print(f"{digest}  {path}")
    print("dependency_hashes=END")
    print(
        "sheet_definition="
        "U:(a,b) in Z/13^5 x Z/169; V:e' in Z/13^5; "
        "full fixed collision labels=(u,v,s,t)=(5,6,1,12)"
    )
    print(
        f"physical_shift={physical_shift}; "
        f"coordinate_shift={shift_coordinate}; exact_cylinder_map=PASS"
    )
    print(
        f"raw_wall_counts=(U:{len(u_events)},V:{len(v_events)}); "
        f"U_chamber={u_chamber}; V_chamber={v_chamber}; "
        f"common_chamber={common_chamber}"
    )
    print(
        f"probe_pair=({probe_zero},{probe_gap}); "
        f"factor_counts=({len(u_zero)},{len(v_zero)}); "
        f"product={len(u_zero) * len(v_zero)}"
    )
    print(
        f"path_set_digest={digest_zero}; identity_label_bijection=PASS"
    )
    print(
        f"supplied_path={SUPPLIED_PATH}; active_at_both_cylinders=yes; "
        "active_on_incoming_same_atom_interval=yes"
    )
    print(
        "first_lost_label=none_through_(u,v,a,b,e'); "
        "remaining_boundary=later_endpoint_origin_or_endpoint_atom_labels_"
        "are_not_constructed_by_this_probe"
    )
    print("verdict=SAME-ANCESTRY PARTIAL GERM AT THE THM-2471 RAIL SHEET")


if __name__ == "__main__":
    main()
