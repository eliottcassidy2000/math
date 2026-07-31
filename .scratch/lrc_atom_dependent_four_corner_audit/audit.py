#!/usr/bin/env python3
"""Independent hostile audit for the pending atom-dependent four-corner file.

This script does not import the witness companion or its conclusions.  It
unfolds the two Perron-transfer multiplicities into their literal inverse
branch labels at the selected physical interval, proves that no inverse-label
event occurs between that interval and its translated target copy, checks one
concrete common path and nearby hostile paths, and recomputes the displayed
two-field coefficient arithmetic.

It intentionally does not certify the absent/present allocation-to-endpoint
translation map.  That is a typing claim about the pending witness source,
not a consequence of these numerical checks.
"""

from __future__ import annotations

from bisect import bisect_right, bisect_left
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    "lrc14_b_r5_theta_target_tensor_thm2584.py":
        "99cbb46c4cf3554fb468dda2d6bbf029f3dd7298a4f8b6009d94513484c7849d",
    "lrc14_base_only_bridge_opus_20260728.py":
        "99deb8891f8be1e218e7b4c5a7bad7a918ddc1da37b05b53ccbbc3ca7b8c7998",
    "lrc14_b_r5_owner_clock_host_thm2581.py":
        "38c8eac0925d119331ef248fa7197f61bd985166c56704575cec18c361543647",
}
for name, expected in PINNED.items():
    require(lf_hash(COMP / name) == expected, f"dependency changed: {name}")


import lrc14_b_r5_theta_target_tensor_thm2584 as rail


P = 13
T = rail.GRID
RPACK = rail.RPACK
DEPTH = rail.DEPTH
SHIFT = 431933040
SOURCE_INTERVAL = (142004992589460, 142005019034340)
TARGET_INTERVAL = tuple(value + SHIFT for value in SOURCE_INTERVAL)
PATH = (59162, 26, 56658)


def inside_scaled(numerator, denominator, intervals, starts):
    """Strict membership of numerator/denominator in one interval union."""
    quotient = numerator // denominator
    index = bisect_right(starts, quotient) - 1
    return (
        index >= 0
        and intervals[index][0] * denominator < numerator
        and numerator < intervals[index][1] * denominator
    )


def labelled_fibres(y, e_set, q_set, e_starts, q_starts):
    """Literal inverse-branch labels at physical coordinate y.

    U first pulls through DEPTH with label a and then through RPACK with
    label b.  The rotated V factor has label eprime.  Packing (a,b) as
    a*RPACK+b keeps the replay reasonably small without losing information.
    """
    u_labels = []
    active_a = 0
    for a in range(DEPTH):
        q_numerator = y + a * T
        if not inside_scaled(q_numerator, DEPTH, q_set, q_starts):
            continue
        active_a += 1
        for b in range(RPACK):
            e_numerator = y + (a + b * DEPTH) * T
            if inside_scaled(
                e_numerator, DEPTH * RPACK, e_set, e_starts
            ):
                u_labels.append(a * RPACK + b)

    # Rail eight has packet displacement one, hence V(y-T/13).
    rotated_y = y - T // P
    v_labels = []
    for eprime in range(DEPTH):
        if inside_scaled(
            rotated_y + eprime * T, DEPTH, e_set, e_starts
        ):
            v_labels.append(eprime)
    return tuple(u_labels), tuple(v_labels), active_a


def label_digest(labels):
    return sha256(",".join(map(str, labels)).encode()).hexdigest()


def event_family(intervals, multiplier, offset=0):
    return tuple(sorted({
        (multiplier * endpoint + offset) % T
        for interval in intervals
        for endpoint in interval
    }))


def chamber_slacks(events, left, right):
    index = bisect_left(events, left)
    require(
        index == len(events) or events[index] >= right,
        "an inverse-label event lies inside the claimed chamber",
    )
    previous = events[index - 1] if index else events[-1] - T
    following = events[index] if index < len(events) else events[0] + T
    return left - previous, following - right


def audit_ancestry():
    e_set = tuple(rail.base.build_set(rail.base.PAT_E, rail.base.ZELL))
    q_set = tuple(rail.base.build_set(rail.host.PAT_QB, rail.base.ZELL))
    e_starts = tuple(left for left, _right in e_set)
    q_starts = tuple(left for left, _right in q_set)

    # No label can enter or leave anywhere in the whole source-to-target
    # hull.  This upgrades a midpoint comparison to a fibrewise statement.
    hull = (SOURCE_INTERVAL[0], TARGET_INTERVAL[1])
    events = {
        "U_Q": event_family(q_set, DEPTH),
        "U_E": event_family(e_set, DEPTH * RPACK),
        "V_E": event_family(e_set, DEPTH, T // P),
    }
    slacks = {
        name: chamber_slacks(values, *hull)
        for name, values in events.items()
    }
    require(all(min(pair) > 0 for pair in slacks.values()),
            "source/target intervals leave a common inverse-label chamber")

    source_y = sum(SOURCE_INTERVAL) // 2
    target_y = sum(TARGET_INTERVAL) // 2
    source_u, source_v, source_a = labelled_fibres(
        source_y, e_set, q_set, e_starts, q_starts
    )
    target_u, target_v, target_a = labelled_fibres(
        target_y, e_set, q_set, e_starts, q_starts
    )

    require(source_a == target_a == 66099, "active-a count changed")
    require(len(source_u) == len(target_u) == 966606,
            "literal U-label count changed")
    require(len(source_v) == len(target_v) == 28534,
            "literal V-label count changed")
    require(source_u == target_u and source_v == target_v,
            "source/target labelled fibres differ")
    require(len(source_u) * len(source_v) == 27581135604,
            "folded multiplicity is not the literal product fibre")

    a, b, eprime = PATH
    packed = a * RPACK + b
    require(packed in source_u and eprime in source_v,
            "selected path is not active")

    # Nearby coordinates make the path test non-vacuous.
    require(a * RPACK + (b - 1) not in source_u,
            "left-b hostile unexpectedly active")
    require(a * RPACK + (b + 1) not in source_u,
            "right-b hostile unexpectedly active")
    require((a - 1) * RPACK + b not in source_u,
            "left-a hostile unexpectedly active")
    require(eprime - 1 not in source_v,
            "left-eprime hostile unexpectedly active")

    expected_u_digest = (
        "d7b599a50d8647ac7957149d1e84588b3602fdb897dc923432102ce0e3a51cf5"
    )
    expected_v_digest = (
        "7072f55b7c54d7c32f2b8ac87ed2aa59f091b3dfdd3d1f193778e9a75c93d00f"
    )
    # These are digests of the packed integer encodings used here, not the
    # tuple-repr digests from the discovery calculation.
    require(label_digest(source_u) == expected_u_digest,
            "packed U-label digest changed")
    require(label_digest(source_v) == expected_v_digest,
            "V-label digest changed")

    return {
        "counts": (source_a, len(source_u), len(source_v),
                   len(source_u) * len(source_v)),
        "digests": (label_digest(source_u), label_digest(source_v)),
        "slacks": slacks,
    }


def hadamard(values, prime):
    v00, v10, v01, v11 = values
    return (
        (v00 + v10 + v01 + v11) % prime,
        (v00 + v10 - v01 - v11) % prime,
        (v00 - v10 + v01 - v11) % prime,
        (v00 - v10 - v01 + v11) % prime,
    )


def audit_field(prime, amplitudes, expected_products, expected_hadamard):
    p0, p1, q0, q1 = amplitudes
    products = (
        p0 * q0 % prime,
        p1 * q0 % prime,
        p0 * q1 % prime,
        p1 * q1 % prime,
    )
    coordinates = hadamard(products, prime)
    require(products == expected_products, "product square changed")
    require(coordinates == expected_hadamard, "Hadamard square changed")
    require(all(amplitudes) and all(products) and all(coordinates),
            "displayed nonzero gate failed")
    require(
        (products[0] * products[3]
         - products[1] * products[2]) % prime == 0,
        "Segre relation failed",
    )
    require(
        (coordinates[0] * coordinates[3]
         - coordinates[1] * coordinates[2]) % prime == 0,
        "Pluecker relation failed",
    )
    coefficient_face = (
        (p0 - p1) * (q0 - q1)
    ) % prime
    require(coefficient_face == coordinates[3] != 0,
            "coefficient mixed face vanished")
    return coefficient_face


def audit_arithmetic():
    face_one = audit_field(
        352341050142921841,
        (
            345283638862540358,
            180096700909697779,
            167108689435423878,
            144423867160388279,
        ),
        (
            347197813355232446,
            3438130309897238,
            320064096166488418,
            146419401152414063,
        ),
        (
            112437340698188483,
            236493496489149044,
            165063327916487722,
            170114988031260853,
        ),
    )
    face_two = audit_field(
        956354278959359281,
        (
            77115584218272234,
            620226426191367084,
            950936471790018101,
            445878683356068526,
        ),
        (
            788257626617883404,
            251268809443191686,
            161379668316310808,
            381695617330433150,
        ),
        (
            626247442748459767,
            496451150414331132,
            316672868160569376,
            757304766188814060,
        ),
    )

    source_step = (0, 1)
    target_step = (12, 0)
    determinant_face = (
        source_step[0] * target_step[1]
        - source_step[1] * target_step[0]
    ) % P
    determinant_corners = (0, 0, 0, determinant_face)
    q_corners = ((0, 0), (0, 1), (1, 0), (1, 1))
    require(determinant_face == 1, "determinant face is not normalized")
    require(determinant_corners == (0, 0, 0, 1),
            "determinant corners changed")
    require(q_corners == ((0, 0), (0, 1), (1, 0), (1, 1)),
            "target-difference corners changed")
    return face_one, face_two, determinant_face


def audit_fixed_sheet_raw_boundary():
    """Rebuild the corrected fixed-xi central and raw Mobius profiles."""
    controls = (
        (
            352341050142921841,
            161934023005863791,
            300649489018742460,
            (
                175075409527953449,
                50230041473974177,
                211754122479689528,
            ),
        ),
        (
            956354278959359281,
            155560747474790541,
            879469786637801433,
            (
                11966211344214657,
                361914377113479889,
                479268797051483842,
            ),
        ),
    )
    rows = []
    for prime, p0, q0, expected in controls:
        inverse = pow(P, -1, prime)
        p1 = p0 * inverse % prime
        q1 = q0 * inverse % prime
        base = p1 * q1 % prime
        products = (
            p0 * q0 % prime,
            p1 * q0 % prime,
            p0 * q1 % prime,
            p1 * q1 % prime,
        )
        coordinates = hadamard(products, prime)
        require(
            (p1, q1, coordinates[3]) == expected,
            "fixed-sheet field values changed",
        )
        require(
            products == tuple(
                multiplier * base % prime
                for multiplier in (P**2, P, P, 1)
            )
            and coordinates == tuple(
                multiplier * base % prime
                for multiplier in (
                    (P + 1)**2, P**2 - 1,
                    P**2 - 1, (P - 1)**2,
                )
            ),
            "fixed-sheet universal central profile failed",
        )

        raw = {}
        for source_twist in range(P):
            for target_twist in range(P):
                vector = (
                    base,
                    base if source_twist == 0 else 0,
                    base if target_twist == 0 else 0,
                    base
                    if source_twist == 0 and target_twist == 0
                    else 0,
                )
                raw[source_twist, target_twist] = (
                    vector,
                    (vector[0] - vector[1]
                     - vector[2] + vector[3]) % prime,
                )
        support_counts = tuple(
            sum(row[0][index] != 0 for row in raw.values())
            for index in range(4)
        )
        all_four = tuple(
            point for point, row in raw.items() if all(row[0])
        )
        mixed_support = tuple(
            point for point, row in raw.items() if row[1]
        )
        require(
            support_counts == (P**2, P, P, 1)
            and all_four == ((0, 0),)
            and raw[0, 0] == ((base, base, base, base), 0)
            and mixed_support
            == tuple(
                (a, b)
                for a in range(1, P)
                for b in range(1, P)
            )
            and all(raw[point][1] == base for point in mixed_support)
            and sum(row[1] for row in raw.values()) % prime
            == coordinates[3],
            "fixed-sheet raw virtual-face boundary failed",
        )
        rows.append(
            (
                prime, (p0, p1, q0, q1), products, coordinates,
                support_counts, all_four, len(mixed_support),
            )
        )
    return tuple(rows)


def main():
    ancestry = audit_ancestry()
    coefficient_faces = audit_arithmetic()
    fixed_sheet = audit_fixed_sheet_raw_boundary()
    print("independent atom-dependent four-corner hostile audit")
    print(f"ancestry_counts(active_a,U,V,UxV)={ancestry['counts']}")
    print(f"ancestry_digests(U,V)={ancestry['digests']}")
    print(f"common_chamber_slacks={ancestry['slacks']}")
    print(f"selected_path={PATH}; nearby_path_hostiles=PASS")
    print(
        "coefficient_mixed_faces(field1,field2)="
        f"{coefficient_faces[:2]}"
    )
    print(f"determinant_mixed_face={coefficient_faces[2]}")
    for row in fixed_sheet:
        print(
            "fixed_sheet_raw_boundary="
            f"prime:{row[0]};factors:{row[1]};products:{row[2]};"
            f"hadamard:{row[3]};raw_supports:{row[4]};"
            f"sole_fourfold:{row[5]};bare_only_mobius_count:{row[6]}"
        )
    print("pluecker_role=ZERO_RANK_ONE_IDENTITY_NOT_NONVANISHING")
    print(
        "unproved_gate=non_idempotent_same_atom_amplitude_then_"
        "carrier_toggles_to_derived_endpoint_translations"
    )
    print("ALL INDEPENDENT NUMERICAL AND ANCESTRY CHECKS PASSED")


if __name__ == "__main__":
    main()
