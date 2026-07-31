#!/usr/bin/env python3
"""Exact omitted-carrier constructor on the retained LRC Boolean sheet.

This scratch audit distinguishes three objects that the pending same-sheet
witness places next to each other:

1. the literal carrier-absent endpoint bank on the fixed labelled unit sheet;
2. the carrier-present zero-total-harmonic endpoint bank;
3. a freely chosen determinant-one endpoint parallelogram.

The absent bank is evaluated directly with the carrier factor omitted.  It is
then extended constantly over the primal carrier-twist variable and Fourier
transformed.  The factor 13 at total harmonic zero is consequently derived,
not used as the definition of absence.

All arithmetic decisions use ``require`` so optimized Python keeps the gate.
Only tracked computation libraries are imported; no pending witness result is
imported.
"""

from __future__ import annotations

from bisect import bisect_right
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
    "lrc14_extended_carrier_endpoint_lib.py":
        "4b3f9f195b1634e1e84a1bc8bccb878a1c8c44aec13f24d197f92547c9e36c57",
    "lrc14_root_zero_full_target_semantic_clutch_20260728.py":
        "208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8",
}
for name, expected in PINNED.items():
    actual = sha256(
        (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    require(actual == expected, f"pinned dependency changed: {name}")


import lrc14_extended_carrier_endpoint_lib as endpoint_base
import lrc14_root_zero_full_target_semantic_clutch_20260728 as physical


P = 13
SHEET = ((142004992589460, 142005019034340, 1),)
TARGET_SHEET = physical.overlap.shift_weighted(SHEET, physical.SHIFT)
ORIGIN = (0, 0)
CLAIMED_SOURCE_STEP = (0, 1)
CLAIMED_TARGET_STEP = (12, 0)


def add(left, right):
    return (
        (left[0] + right[0]) % P,
        (left[1] + right[1]) % P,
    )


def indexed_intersection(pieces, intervals, starts):
    """Intersect a tiny weighted carrier with a sorted interval union."""
    out = []
    for left, right, weight in pieces:
        index = max(0, bisect_right(starts, left) - 1)
        while index < len(intervals) and intervals[index][0] < right:
            a = max(left, intervals[index][0])
            b = min(right, intervals[index][1])
            if a < b:
                out.append((a, b, weight))
            index += 1
    return tuple(out)


def intersect_weighted(first, second):
    """Product of two sorted unit-weight profiles on the circle."""
    out = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            out.append((
                left,
                right,
                first[i][2] * second[j][2],
            ))
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def endpoint_sum(pieces, frequency):
    values = [0] * len(endpoint_base.endpoint.MODS)
    for left, right, weight in pieces:
        local = endpoint_base.endpoint.endpoint_sum(
            ((left, right),),
            frequency,
            endpoint_base.endpoint.MODS,
        )
        for index, (prime, _root) in enumerate(endpoint_base.endpoint.MODS):
            values[index] = (
                values[index] + weight * local[index]
            ) % prime
    return tuple(values)


def x_sweep(pieces, terminal, q_starts, tabs):
    values = [0] * len(endpoint_base.endpoint.MODS)
    for left, right, weight in pieces:
        local, _overlap = endpoint_base.endpoint.x_sweep(
            ((left, right),),
            terminal,
            q_starts,
            endpoint_base.X0,
            endpoint_base.endpoint.MODS,
            tabs,
        )
        for index, (prime, _root) in enumerate(endpoint_base.endpoint.MODS):
            values[index] = (
                values[index] + weight * local[index]
            ) % prime
    return tuple(values)


def hadamard(products, prime):
    v00, v10, v01, v11 = products
    return (
        (v00 + v10 + v01 + v11) % prime,
        (v00 + v10 - v01 - v11) % prime,
        (v00 - v10 + v01 - v11) % prime,
        (v00 - v10 - v01 + v11) % prime,
    )


def flag_square(factors, prime):
    p0, p1, q0, q1 = factors
    products = (
        p0 * q0 % prime,
        p1 * q0 % prime,
        p0 * q1 % prime,
        p1 * q1 % prime,
    )
    coordinates = hadamard(products, prime)
    require(
        coordinates[3] == (p0 - p1) * (q0 - q1) % prime,
        "Boolean mixed face stopped factoring",
    )
    require(
        (
            coordinates[0] * coordinates[3]
            - coordinates[1] * coordinates[2]
        ) % prime == 0,
        "Pluecker identity failed",
    )
    return products, coordinates


def selected_one_sided_carriers():
    """Rebuild the selected physical A/B carriers without pending scratch."""
    (
        module, _prefixes, _whole, _masses, rails, present, _starts
    ) = physical.relative.lift.m.core.build_carrier_data()
    pair_prefixes = physical.relative.private.build_pair_prefixes(module)
    _sv, _tv, _rail_pairs, details = physical.overlap.overlap_vectors(
        module, pair_prefixes, rails, present, rail_index=8
    )
    full_module = physical.target.load_present_module()
    e3 = physical.target.exclusive_source(full_module, 3)
    clocks = tuple(
        full_module.make_comb(
            full_module.C1, 182, 26 * clock - 13, 26 * clock + 13
        )
        for clock in range(7)
    )
    source_base, target_base = details[1]
    section = physical.target.source_present_section(
        full_module, e3, 1, 0, 4, clocks
    )
    source = tuple(
        physical.relative.private.old.intersect_weighted_union(
            source_base, section
        )
    )
    target = tuple(
        physical.relative.private.old.intersect_weighted_union(
            target_base, section
        )
    )
    source_support = physical.overlap.merge_intervals(
        (left, right)
        for left, right, weight in source
        if weight
    )
    target_support = physical.overlap.merge_intervals(
        (left, right)
        for left, right, weight in target
        if weight
    )
    return source_support, target_support


def character_covariances(bare, present, zeta, prime):
    """Find present(x)=c*zeta^(u.x)*bare(x+s) on the whole endpoint plane."""
    roots = {
        pow(zeta, exponent, prime): exponent
        for exponent in range(P)
    }
    answers = []
    failures = {}
    e1 = (1, 0)
    e2 = (0, 1)
    for shift in endpoint_base.KEYS:
        ratios = {}
        support_failure = None
        for point in endpoint_base.KEYS:
            denominator = bare[add(point, shift)]
            numerator = present[point]
            if not denominator or not numerator:
                if bool(denominator) != bool(numerator):
                    support_failure = (point, denominator, numerator)
                    break
                continue
            ratios[point] = numerator * pow(denominator, -1, prime) % prime
        if support_failure is not None:
            failures[shift] = ("support", support_failure)
            continue
        require(ORIGIN in ratios, "origin ratio unavailable")
        scalar = ratios[ORIGIN]
        lambda_one = ratios[e1] * pow(scalar, -1, prime) % prime
        lambda_two = ratios[e2] * pow(scalar, -1, prime) % prime
        if lambda_one not in roots:
            failures[shift] = (
                "first-generator-not-C13",
                (lambda_one, pow(lambda_one, P, prime) - 1),
            )
            continue
        if lambda_two not in roots:
            failures[shift] = (
                "second-generator-not-C13",
                (lambda_two, pow(lambda_two, P, prime) - 1),
            )
            continue
        u = (roots[lambda_one], roots[lambda_two])
        first = None
        for point in endpoint_base.KEYS:
            expected = (
                scalar
                * pow(
                    zeta,
                    (u[0] * point[0] + u[1] * point[1]) % P,
                    prime,
                )
            ) % prime
            if ratios[point] != expected:
                first = (point, ratios[point], expected)
                break
        if first is None:
            answers.append((shift, u, scalar))
        else:
            failures[shift] = ("first-character-failure", first)
    return tuple(answers), failures


def construct():
    terminal = tuple(
        endpoint_base.endpoint.build_set(
            endpoint_base.PAT_Q12, endpoint_base.ZERO
        )
    )
    q_starts = tuple(left for left, _right in terminal)
    tabs = endpoint_base.endpoint.make_tabs(
        terminal,
        endpoint_base.X0,
        endpoint_base.endpoint.MODS,
    )
    present_sets = endpoint_base.present_cache()
    unit = endpoint_base.T // P
    source_twists = tuple(
        physical.overlap.shift_weighted(SHEET, -twist * unit)
        for twist in range(P)
    )
    target_twists = tuple(
        physical.overlap.shift_weighted(TARGET_SHEET, twist * unit)
        for twist in range(P)
    )
    require(
        source_twists[0] == SHEET
        and target_twists[0] == TARGET_SHEET,
        "zero-twist carrier is not the retained sheet selector",
    )
    source_one_sided, target_one_sided = selected_one_sided_carriers()
    source_one_sided_twists = tuple(
        physical.overlap.shift_union(
            source_one_sided, -twist * unit
        )
        for twist in range(P)
    )
    target_one_sided_twists = tuple(
        physical.overlap.shift_union(
            target_one_sided, twist * unit
        )
        for twist in range(P)
    )
    literal_source_twists = tuple(
        indexed_intersection(
            SHEET,
            carrier,
            tuple(left for left, _right in carrier),
        )
        for carrier in source_one_sided_twists
    )
    literal_target_twists = tuple(
        indexed_intersection(
            TARGET_SHEET,
            carrier,
            tuple(left for left, _right in carrier),
        )
        for carrier in target_one_sided_twists
    )
    require(
        literal_source_twists[0] == SHEET
        and literal_target_twists[0] == TARGET_SHEET
        and all(not pieces for pieces in literal_source_twists[1:])
        and all(not pieces for pieces in literal_target_twists[1:]),
        "fixed-sheet carrier allocation is not a twist-zero delta",
    )
    source_mask = tuple(bool(pieces) for pieces in literal_source_twists)
    target_mask = tuple(bool(pieces) for pieces in literal_target_twists)
    source_mask_after_gauge = tuple(
        source_mask[(index - 1) % P] for index in range(P)
    )
    target_mask_after_gauge = tuple(
        target_mask[(index + 1) % P] for index in range(P)
    )
    require(
        source_mask == target_mask == (True,) + (False,) * (P - 1)
        and source_mask_after_gauge
        == (False, True) + (False,) * (P - 2)
        and target_mask_after_gauge
        == (False,) * (P - 1) + (True,),
        "delta mask did not transport with the representative gauge",
    )

    field_count = len(endpoint_base.endpoint.MODS)
    # These are the replacement/transported-sheet banks used by the pending
    # probe.  They are retained below as a hostile comparison.
    left_present = [dict() for _ in range(field_count)]
    right_present = [dict() for _ in range(field_count)]
    # These banks keep the labelled sheet fixed and multiply it by the
    # carrier indicator.  They are the literal same-sheet allocation.
    left_literal_present = [dict() for _ in range(field_count)]
    right_literal_present = [dict() for _ in range(field_count)]
    left_absent_raw = [dict() for _ in range(field_count)]
    right_absent_raw = [dict() for _ in range(field_count)]
    marked_gauge_left_rows = [0] * field_count
    marked_gauge_right_rows = [0] * field_count
    forgotten_sheet_left_mismatches = [0] * field_count
    forgotten_sheet_right_mismatches = [0] * field_count
    first_forgotten_sheet_left = [None] * field_count
    first_forgotten_sheet_right = [None] * field_count
    left_cache = {}
    right_cache = {}

    def left_values(ell, present, starts, pieces):
        restricted = indexed_intersection(pieces, present, starts)
        if restricted not in left_cache:
            left_cache[restricted] = x_sweep(
                restricted, terminal, q_starts, tabs
            )
        values = []
        for index, (prime, root) in enumerate(
            endpoint_base.endpoint.MODS
        ):
            zeta = pow(root, endpoint_base.NRED // P, prime)
            phase = pow(zeta, ell[endpoint_base.DEEP], prime)
            values.append(left_cache[restricted][index] * phase % prime)
        return tuple(values)

    def right_values(present, starts, pieces):
        restricted = indexed_intersection(pieces, present, starts)
        if restricted not in right_cache:
            right_cache[restricted] = endpoint_sum(
                restricted, -endpoint_base.Y0
            )
        return right_cache[restricted]

    for address, ell in endpoint_base.REPS.items():
        present = present_sets[address]
        starts = tuple(left for left, _right in present)
        shifted_ell = tuple(
            (ell[index] + endpoint_base.WMOD[index]) % P
            for index in range(9)
        )
        shifted_present = physical.overlap.shift_union(
            present, -unit
        )
        shifted_starts = tuple(
            left for left, _right in shifted_present
        )

        # Literal omission: evaluate the fixed labelled sheet with no
        # carrier-twist factor or carrier-twist variable.
        direct_left_absent = left_values(ell, present, starts, SHEET)
        direct_right_absent = right_values(
            present, starts, TARGET_SHEET
        )
        for field_index in range(field_count):
            left_absent_raw[field_index][address] = (
                direct_left_absent[field_index]
            )
            right_absent_raw[field_index][address] = (
                direct_right_absent[field_index]
            )

        # THM-2763 representative covariance is coefficientwise only when
        # the marked sheet and its twist origin move with the representative:
        # source delta_0 -> delta_1, target delta_0 -> delta_12.
        transported_left = left_values(
            shifted_ell,
            shifted_present,
            shifted_starts,
            source_twists[1],
        )
        transported_right = right_values(
            shifted_present,
            shifted_starts,
            target_twists[P - 1],
        )
        # Hostile: change representative but forget to move the marked sheet.
        forgotten_left = left_values(
            shifted_ell,
            shifted_present,
            shifted_starts,
            SHEET,
        )
        forgotten_right = right_values(
            shifted_present,
            shifted_starts,
            TARGET_SHEET,
        )
        for field_index in range(field_count):
            marked_gauge_left_rows[field_index] += (
                transported_left[field_index]
                == direct_left_absent[field_index]
            )
            marked_gauge_right_rows[field_index] += (
                transported_right[field_index]
                == direct_right_absent[field_index]
            )
            if (
                forgotten_left[field_index]
                != direct_left_absent[field_index]
            ):
                forgotten_sheet_left_mismatches[field_index] += 1
                if first_forgotten_sheet_left[field_index] is None:
                    first_forgotten_sheet_left[field_index] = (
                        address,
                        direct_left_absent[field_index],
                        forgotten_left[field_index],
                    )
            if (
                forgotten_right[field_index]
                != direct_right_absent[field_index]
            ):
                forgotten_sheet_right_mismatches[field_index] += 1
                if first_forgotten_sheet_right[field_index] is None:
                    first_forgotten_sheet_right[field_index] = (
                        address,
                        direct_right_absent[field_index],
                        forgotten_right[field_index],
                    )

        for twist in range(P):
            left_values_here = left_values(
                ell, present, starts, source_twists[twist]
            )
            right_values_here = right_values(
                present, starts, target_twists[twist]
            )
            for field_index in range(field_count):
                key = (address, twist)
                left_present[field_index][key] = (
                    left_values_here[field_index]
                )
                right_present[field_index][key] = (
                    right_values_here[field_index]
                )
                left_literal_present[field_index][key] = (
                    direct_left_absent[field_index] if twist == 0 else 0
                )
                right_literal_present[field_index][key] = (
                    direct_right_absent[field_index] if twist == 0 else 0
                )

        # The equality is a consequence of 1_sheet * 1_carrier(0)=1_sheet.
        for field_index in range(field_count):
            require(
                left_absent_raw[field_index][address]
                == left_present[field_index][address, 0],
                "direct left omitted-carrier value differs from selector identity",
            )
            require(
                right_absent_raw[field_index][address]
                == right_present[field_index][address, 0],
                "direct right omitted-carrier value differs from selector identity",
            )

    reports = []
    for field_index, (prime, root) in enumerate(
        endpoint_base.endpoint.MODS
    ):
        zeta = pow(root, endpoint_base.NRED // P, prime)
        powers = tuple(pow(zeta, exponent, prime) for exponent in range(P))
        left_absent_address = {
            address: P * left_absent_raw[field_index][address] % prime
            for address in endpoint_base.KEYS
        }
        right_absent_address = {
            address: P * right_absent_raw[field_index][address] % prime
            for address in endpoint_base.KEYS
        }
        left_literal_present_address = {
            address: sum(
                left_literal_present[field_index][address, twist]
                for twist in range(P)
            ) % prime
            for address in endpoint_base.KEYS
        }
        right_literal_present_address = {
            address: sum(
                right_literal_present[field_index][address, twist]
                for twist in range(P)
            ) % prime
            for address in endpoint_base.KEYS
        }
        left_transported_address = {
            address: sum(
                left_present[field_index][address, twist]
                for twist in range(P)
            ) % prime
            for address in endpoint_base.KEYS
        }
        right_transported_address = {
            address: sum(
                right_present[field_index][address, twist]
                for twist in range(P)
            ) % prime
            for address in endpoint_base.KEYS
        }

        def endpoint_plane(*, side, state):
            bank = {}
            sign = -1 if side == "left" else 1
            raw_absent = (
                left_absent_raw[field_index]
                if side == "left" else right_absent_raw[field_index]
            )
            raw_present = (
                left_present[field_index]
                if side == "left" else right_present[field_index]
            )
            raw_literal = (
                left_literal_present[field_index]
                if side == "left" else right_literal_present[field_index]
            )
            for point in endpoint_base.KEYS:
                total = 0
                for address in endpoint_base.KEYS:
                    pairing = (
                        address[0] * point[0]
                        + address[1] * point[1]
                    ) % P
                    phase = powers[(sign * pairing) % P]
                    if state == "absent":
                        # Direct constant extension in the omitted primal
                        # twist variable, followed by its k=0 Fourier sum.
                        for _twist in range(P):
                            total += raw_absent[address] * phase
                    elif state == "literal-present":
                        for twist in range(P):
                            total += raw_literal[address, twist] * phase
                    else:
                        require(
                            state == "transported-present",
                            "unknown endpoint-plane state",
                        )
                        for twist in range(P):
                            total += raw_present[address, twist] * phase
                bank[point] = total % prime
            return bank

        left_absent = endpoint_plane(side="left", state="absent")
        left_present_plane = endpoint_plane(
            side="left", state="literal-present"
        )
        right_absent = endpoint_plane(side="right", state="absent")
        right_present_plane = endpoint_plane(
            side="right", state="literal-present"
        )
        left_transported_plane = endpoint_plane(
            side="left", state="transported-present"
        )
        right_transported_plane = endpoint_plane(
            side="right", state="transported-present"
        )

        # Explicitly prove that the derived unnormalized central absent bank
        # is 13 times the address DFT of the directly evaluated raw bank.
        for point in endpoint_base.KEYS:
            direct_left_dft = sum(
                left_absent_raw[field_index][address]
                * powers[
                    -(
                        address[0] * point[0]
                        + address[1] * point[1]
                    ) % P
                ]
                for address in endpoint_base.KEYS
            ) % prime
            direct_right_dft = sum(
                right_absent_raw[field_index][address]
                * powers[
                    (
                        address[0] * point[0]
                        + address[1] * point[1]
                    ) % P
                ]
                for address in endpoint_base.KEYS
            ) % prime
            require(
                left_absent[point] == P * direct_left_dft % prime,
                "left factor 13 was not derived by constant extension",
            )
            require(
                right_absent[point] == P * direct_right_dft % prime,
                "right factor 13 was not derived by constant extension",
            )
            require(
                left_present_plane[point] == direct_left_dft,
                "literal left present bank is not the twist-zero delta DFT",
            )
            require(
                right_present_plane[point] == direct_right_dft,
                "literal right present bank is not the twist-zero delta DFT",
            )
            require(
                left_absent[point]
                == P * left_present_plane[point] % prime
                and right_absent[point]
                == P * right_present_plane[point] % prime,
                "literal absent/present scalar law failed",
            )

        fixed_factors = (
            left_absent[ORIGIN],
            left_present_plane[ORIGIN],
            right_absent[ORIGIN],
            right_present_plane[ORIGIN],
        )
        fixed_products, fixed_hadamard = flag_square(
            fixed_factors, prime
        )
        base_product = fixed_products[3]
        require(
            fixed_products == tuple(
                multiplier * base_product % prime
                for multiplier in (P**2, P, P, 1)
            )
            and fixed_hadamard == tuple(
                multiplier * base_product % prime
                for multiplier in (
                    (P + 1)**2,
                    P**2 - 1,
                    P**2 - 1,
                    (P - 1)**2,
                )
            ),
            "universal central allocation profile failed",
        )
        primal_common_vector = (base_product,) * 4
        primal_common_hadamard = hadamard(
            primal_common_vector, prime
        )
        require(
            primal_common_hadamard[3] == 0
            and (P - 1)**2 * base_product % prime
            == fixed_hadamard[3],
            "primal/central mixed-face boundary failed",
        )
        translated_factors = (
            left_absent[ORIGIN],
            left_present_plane[add(ORIGIN, CLAIMED_SOURCE_STEP)],
            right_absent[ORIGIN],
            right_present_plane[add(ORIGIN, CLAIMED_TARGET_STEP)],
        )
        translated_products, translated_hadamard = flag_square(
            translated_factors, prime
        )
        left_covariances, left_failures = character_covariances(
            left_absent, left_present_plane, zeta, prime
        )
        right_covariances, right_failures = character_covariances(
            right_absent, right_present_plane, zeta, prime
        )
        transported_left_covariances, _transported_left_failures = (
            character_covariances(
                left_absent, left_transported_plane, zeta, prime
            )
        )
        transported_right_covariances, _transported_right_failures = (
            character_covariances(
                right_absent, right_transported_plane, zeta, prime
            )
        )
        intrinsic_covariance = (ORIGIN, ORIGIN, pow(P, -1, prime))
        require(
            intrinsic_covariance in left_covariances
            and intrinsic_covariance in right_covariances,
            "literal allocation lost its zero-step scalar covariance",
        )
        require(
            all(
                covariance[0] != CLAIMED_SOURCE_STEP
                for covariance in left_covariances
            )
            and all(
                covariance[0] != CLAIMED_TARGET_STEP
                for covariance in right_covariances
            ),
            "a claimed determinant step became an accidental covariance",
        )
        require(
            not transported_left_covariances
            and not transported_right_covariances,
            "transported-sheet hostile unexpectedly became a covariance",
        )
        require(
            all(fixed_factors)
            and all(fixed_products)
            and all(fixed_hadamard),
            "fixed-origin allocation square lost co-support",
        )
        reports.append({
            "prime": prime,
            "raw_absent_equals_present_twist_zero": True,
            "absent_support": (
                sum(bool(value) for value in left_absent.values()),
                sum(bool(value) for value in right_absent.values()),
            ),
            "present_support": (
                sum(bool(value) for value in left_present_plane.values()),
                sum(bool(value) for value in right_present_plane.values()),
            ),
            "address_support": (
                sum(bool(value) for value in left_absent_address.values()),
                sum(
                    bool(value)
                    for value in left_literal_present_address.values()
                ),
                sum(bool(value) for value in right_absent_address.values()),
                sum(
                    bool(value)
                    for value in right_literal_present_address.values()
                ),
            ),
            "transported_address_support": (
                sum(bool(value) for value in left_transported_address.values()),
                sum(bool(value) for value in right_transported_address.values()),
            ),
            "fixed_factors": fixed_factors,
            "fixed_products": fixed_products,
            "fixed_hadamard": fixed_hadamard,
            "primal_common_vector": primal_common_vector,
            "primal_common_hadamard": primal_common_hadamard,
            "translated_factors": translated_factors,
            "translated_products": translated_products,
            "translated_hadamard": translated_hadamard,
            "left_covariances": left_covariances,
            "right_covariances": right_covariances,
            "transported_covariances": (
                transported_left_covariances,
                transported_right_covariances,
            ),
            "claimed_source_step_failure":
                left_failures[CLAIMED_SOURCE_STEP],
            "claimed_target_step_failure":
                right_failures[CLAIMED_TARGET_STEP],
            "all_shift_failure_census": (
                len(left_failures), len(right_failures)
            ),
            "marked_gauge_rows": (
                marked_gauge_left_rows[field_index],
                marked_gauge_right_rows[field_index],
            ),
            "forgotten_sheet_hostile": (
                forgotten_sheet_left_mismatches[field_index],
                forgotten_sheet_right_mismatches[field_index],
                first_forgotten_sheet_left[field_index],
                first_forgotten_sheet_right[field_index],
            ),
        })
    require(
        all(value == P**2 for value in marked_gauge_left_rows)
        and all(value == P**2 for value in marked_gauge_right_rows),
        "marked-sheet representative covariance failed",
    )
    require(
        all(value for value in forgotten_sheet_left_mismatches)
        and all(value for value in forgotten_sheet_right_mismatches),
        "forgetting the marked sheet stopped being hostile",
    )
    return (
        tuple(reports),
        (len(source_one_sided), len(target_one_sided)),
    )


def main():
    reports, carrier_piece_counts = construct()
    print("literal carrier-absent/present constructor exact audit")
    print(f"dependency_sha256={tuple(PINNED.items())}")
    print(f"retained_source_sheet={SHEET}")
    print(f"retained_target_sheet={TARGET_SHEET}")
    print(
        "selected_one_sided_carrier_piece_counts_(source,target)="
        f"{carrier_piece_counts}"
    )
    print(
        "fixed_sheet_times_full_one_sided_carrier_twists="
        "(sheet,empty,...,empty) on both endpoint sides"
    )
    print(
        "representative_transport=(ell,xi,delta0)->"
        "(ell+W,Txi,source_delta1/target_delta12);"
        "central_k0_sum_is_cyclic_reindexing_invariant"
    )
    print(
        "omitted_carrier_constructor=endpoint functional of the retained "
        "sheet with no carrier-twist factor and no twist variable"
    )
    print(
        "selector_identity=absent_raw==present_raw_at_twist_0 "
        "for all 169 addresses in both fields and both endpoint sides"
    )
    print(
        "central_absent_normalization=constant extension over 13 primal "
        "twists, then unnormalized k=0 DFT; factor 13 DERIVED"
    )
    for report in reports:
        print(f"field_prime={report['prime']}")
        print(
            "endpoint_support_(absent_left,right;present_left,right)="
            f"{(report['absent_support'], report['present_support'])}"
        )
        print(
            "inverse_address_support_(P0,P1,Q0,Q1)="
            f"{report['address_support']}"
        )
        print(
            "transported_replacement_address_support_(P1,Q1)="
            f"{report['transported_address_support']}"
        )
        print(
            "fixed_origin_factors_(P0,P1,Q0,Q1)="
            f"{report['fixed_factors']}"
        )
        print(f"fixed_origin_products={report['fixed_products']}")
        print(
            "fixed_origin_hadamard_(D0,D1,D2,D3)="
            f"{report['fixed_hadamard']}"
        )
        print(
            "primal_fourfold_cosupport_at_(a,b)=(0,0)="
            f"{report['primal_common_vector']};"
            "primal_hadamard="
            f"{report['primal_common_hadamard']}"
        )
        print(
            "primal_allocation_support_counts_(B,P,Q,H,Omega)="
            f"{(P**2, P, P, 1, (P - 1)**2)}"
        )
        print(
            "translated_candidate_factors="
            f"{report['translated_factors']}"
        )
        print(
            "translated_candidate_products="
            f"{report['translated_products']}"
        )
        print(
            "translated_candidate_hadamard="
            f"{report['translated_hadamard']}"
        )
        print(
            "literal_affine_covariances_(source,target)="
            f"{(report['left_covariances'], report['right_covariances'])}"
        )
        print(
            "transported_replacement_covariances_(source,target)="
            f"{report['transported_covariances']}"
        )
        print(
            "claimed_source_step_smallest_failure="
            f"{report['claimed_source_step_failure']}"
        )
        print(
            "claimed_target_step_smallest_failure="
            f"{report['claimed_target_step_failure']}"
        )
        print(
            "all_169_shifts_refuted_(source,target)="
            f"{report['all_shift_failure_census']}"
        )
        print(
            "marked_sheet_representative_rows_(source,target)="
            f"{report['marked_gauge_rows']}"
        )
        print(
            "forgotten_marked_sheet_hostile="
            f"{report['forgotten_sheet_hostile']}"
        )
    print(
        "true_fixed_origin_endpoint_labels=all four flags retain "
        "(L,R,q,Delta)=((0,0),(0,0),(0,0),0)"
    )
    print(
        "no_go=the literal fixed-origin Boolean allocation square has "
        "nonzero D3, but each allocation toggle has intrinsic endpoint "
        "step zero and scalar 1/13; replacing xi by translated sheets has "
        "no shift-character covariance at all.  Hence the determinant-one "
        "parallelogram is a separate carried-sample control"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
