#!/usr/bin/env python3
"""THM-2791: full arm orbit, lower-central chord, and coarse transfer.

This exact companion reconstructs the THM-2782 clock-one physical carrier,
scans the complete C_(13^6) address orbit and its C_(13^5) central coset,
classifies the positive two-cylinder chord, verifies its lower-central
digits and finite differences, and checks the transfer/decoder spectrum.

No floating point, randomized calculation, Python assertion, or scratch
dependency is used.
"""

from bisect import bisect_left
from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_semantic_arm_right_wing_central_digit_thm2782 as thm


P = 13
N = P**5
D = P**5


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def target_weighted_by_tau(context, tau):
    module, delayed, present, source, clocks, source_weight, rail_common = context
    objects = thm.carrier_objects(
        module, source, clocks, rail_common, 1, 0, tau
    )
    source_safe = thm.marked.complement(present[1, 7])
    target_safe = thm.marked.complement(
        thm.marked.shift_union(present[1, 7], thm.marked.SHIFT)
    )
    carrier = thm.marked.intersect(objects["R"], source_safe)
    carrier = thm.marked.intersect(carrier, target_safe)
    carrier = thm.marked.shift_union(carrier, -thm.marked.SHIFT)
    weighted = thm.marked.restrict_weighted(source_weight, carrier)
    weighted = thm.marked.private.old.intersect_weighted_comb(
        weighted, module.C3, 182, 1, 13
    )
    return delayed, tuple(weighted)


def candidate_indices(weighted, base, count=N, address_stride=13):
    """Return exact indices whose narrow cylinder can meet a weighted interval."""
    half = thm.relative.Q_RADIUS * thm.T
    require(half.denominator == D, "unexpected cylinder-radius denominator")
    half_num = half.numerator
    center0 = thm.arm_target_center(base, 0) * thm.T
    require(center0.denominator == D, "unexpected address-centre denominator")
    center_num = center0.numerator
    modulus = thm.T * D
    step = (address_stride * thm.EXPECTED_SHIFT * D) % modulus
    scaled = tuple((left * D, right * D) for left, right, _weight in weighted)
    rights = tuple(right for _left, right in scaled)
    hits = []
    for j in range(count):
        c = center_num % modulus
        pieces = []
        if c - half_num < 0:
            pieces.extend(((0, c + half_num), (modulus + c - half_num, modulus)))
        elif c + half_num > modulus:
            pieces.extend(((c - half_num, modulus), (0, c + half_num - modulus)))
        else:
            pieces.append((c - half_num, c + half_num))
        found = False
        for lo, hi in pieces:
            index = bisect_left(rights, lo + 1)
            if index < len(scaled) and scaled[index][0] < hi:
                found = True
                break
        if found:
            hits.append(j)
        center_num = (center_num + step) % modulus
    return tuple(hits)


def exact_cell(delayed, weighted, base, j):
    center = thm.arm_target_center(base, j)
    half = thm.relative.Q_RADIUS * thm.T
    coordinate = center * thm.T
    if coordinate - half < 0:
        cylinders = (
            (0, coordinate + half),
            (thm.T + coordinate - half, thm.T),
        )
    elif coordinate + half > thm.T:
        cylinders = (
            (coordinate - half, thm.T),
            (0, coordinate + half - thm.T),
        )
    else:
        cylinders = ((coordinate - half, coordinate + half),)
    pieces = thm.marked.restrict_weighted(weighted, cylinders)
    mass = thm.weighted_mass(pieces)
    coefficient = thm.marked.private.delayed_carry_pair(
        pieces, delayed[1], {}
    )[6][1]
    return mass, coefficient


def minimal_cyclic_period(row):
    for exponent in range(6):
        period = P**exponent
        if N % period == 0 and all(
            row[index] == row[index % period] for index in range(N)
        ):
            return period
    return N


def main():
    print("THM-2791 FULL ARM ORBIT / LOWER-CENTRAL CHORD")
    print("status=FINITE-EXACT proof-complete candidate; no LRC conclusion")
    context = thm.build_context()
    base = thm.EXPECTED_BASES[1]
    supports = {}
    coefficients = {}
    weighted_profiles = {}
    for tau in range(P):
        delayed, weighted = target_weighted_by_tau(context, tau)
        weighted_profiles[tau] = weighted
        candidates = candidate_indices(weighted, base)
        exact = tuple(
            (j, *exact_cell(delayed, weighted, base, j))
            for j in candidates
        )
        nonzero = tuple(row for row in exact if row[2])
        supports[tau] = tuple(j for j, _mass, _coefficient in nonzero)
        coefficients[tau] = {
            j: coefficient for j, _mass, coefficient in nonzero
        }
        print(
            f"tau={tau} weighted_pieces={len(weighted)} "
            f"candidates={len(candidates)} nonzero={len(nonzero)}"
        )
        print(f"support={supports[tau]}")
        print(f"values={tuple(sorted(set(coefficients[tau].values())))}")

    require(
        all(
            weighted_profiles[tau] == weighted_profiles[3]
            for tau in range(3, 12)
        ),
        "raw target columns 3,...,11 stopped sharing one physical carrier",
    )

    # Sparse period test: every proper p-power period would carry each
    # nonempty support into itself.
    for tau in range(P):
        support = set(supports[tau])
        if not support:
            continue
        periods = []
        for exponent in range(5):
            period = P**exponent
            if all(((j + period) % N) in support for j in support):
                periods.append(period)
        print(f"tau={tau} proper_support_periods={tuple(periods)}")

    union = sorted(set().union(*(set(row) for row in supports.values())))
    print(f"N={N} union_support_size={len(union)}")
    print(f"union_support={tuple(union)}")

    quotient_rows = tuple(
        tuple(
            sum(1 for j in supports[tau] if j % P == residue)
            for residue in range(P)
        )
        for tau in range(P)
    )
    print(f"quotient_pushforward_by_(tau,w)={quotient_rows}")

    # Integer target decoding commutes with the high-digit transfer.
    decoded_quotient = tuple(
        tuple(
            sum(
                quotient_rows[tau][residue]
                * thm.K_BETA[(target - tau) % P]
                for tau in range(P)
            )
            for residue in range(P)
        )
        for target in range(P)
    )
    print(f"decoded_quotient_by_(target,w)={decoded_quotient}")
    print(
        "decoded_quotient_nonconstant="
        f"{tuple(len(set(row)) > 1 for row in decoded_quotient)}"
    )

    decoded_sparse = []
    for target in range(P):
        row = {}
        for tau in range(P):
            multiplier = thm.K_BETA[(target - tau) % P]
            if multiplier:
                for j in supports[tau]:
                    row[j] = row.get(j, 0) + multiplier
        row = {j: value for j, value in row.items() if value}
        decoded_sparse.append(row)
        conductor_survival = []
        for exponent in range(1, 6):
            modulus = P**exponent
            stride = P ** (exponent - 1)
            counts = [0] * modulus
            for j, value in row.items():
                counts[j % modulus] += value
            vanishes = all(
                len({counts[a + k * stride] for k in range(P)}) == 1
                for a in range(stride)
            )
            conductor_survival.append(not vanishes)
        print(
            f"decoded_target={target} support={len(row)} "
            f"value_set={tuple(sorted(set(row.values())))} "
            f"primitive_conductors={tuple(conductor_survival)}"
        )

    # For a rational C_(13^5) profile, vanishing at one primitive
    # 13^r-character is equivalent to equality of the 13 coefficient
    # blocks modulo 13^r.  Record the exact all-or-none conductor test.
    for tau in range(P):
        support = supports[tau]
        if not support:
            continue
        conductor_survival = []
        for exponent in range(1, 6):
            modulus = P**exponent
            stride = P ** (exponent - 1)
            counts = [0] * modulus
            for j in support:
                counts[j % modulus] += 1
            vanishes = all(
                len({counts[a + k * stride] for k in range(P)}) == 1
                for a in range(stride)
            )
            conductor_survival.append(not vanishes)
        print(
            f"tau={tau} primitive_conductors_13_to_13^5="
            f"{tuple(conductor_survival)}"
        )

    delayed, weighted = target_weighted_by_tau(context, 3)
    full_candidates = candidate_indices(
        weighted, base, count=thm.ADDRESS_MODULUS, address_stride=1
    )
    full_exact = tuple(
        (
            index,
            *exact_cell(delayed, weighted, base, Fraction(index, 13)),
        )
        for index in full_candidates
    )
    # exact_cell takes an arm count and therefore adds 13*arm to the address;
    # Fraction(index,13) makes this the full +1 address orbit.
    print(f"full_address_tau3_candidates={full_candidates}")
    print(f"full_address_tau3_exact={full_exact}")
    full_metadata = tuple(
        (
            (base + index) % thm.ADDRESS_MODULUS,
            thm.omega_digits(base + index),
            thm.relative.semantic.semantic_record(
                thm.arm_target_center(base, Fraction(index, 13))
            ),
            thm.predecessor_carry(
                thm.arm_target_center(base, Fraction(index, 13))
            ),
        )
        for index in full_candidates
    )
    require(
        all(
            record == (3, (1, 2)) and carry == 6
            for _address, _omega, record, carry in full_metadata
        ),
        "full chord lost semantic/carry typing",
    )
    for index in full_candidates:
        center = thm.arm_target_center(base, Fraction(index, 13))
        sigma_labels, tau_labels = thm.relative.full_target_labels(center)
        require(
            0 in sigma_labels
            and all(tau in tau_labels for tau in range(3, 12)),
            "full chord lost one of its lawful target labels",
        )
        require(
            thm.relative.semantic_stability_radius(center)
            == thm.relative.Q_RADIUS,
            "full chord lost its sharp open semantic-radius boundary",
        )
    print(f"full_address_tau3_metadata={full_metadata}")

    require(
        full_candidates[0] == 0 and len(full_candidates) == 2,
        "full address support is not the expected two-point chord",
    )
    address_gap = full_candidates[1]
    require(address_gap % P == 0, "full address chord left the central coset")
    central_gap = address_gap // P
    correction = (central_gap - 1) // P
    require(
        central_gap == 53028
        and central_gap == 1 + P * correction
        and correction == 4079
        and correction % P == 10
        and central_gap % P == 1,
        "graded chord digits changed",
    )

    def base_p_digits(value, length):
        result = []
        for _ in range(length):
            result.append(value % P)
            value //= P
        require(value == 0, "base-p digit buffer was too short")
        return tuple(result)

    require(
        base_p_digits(central_gap, 5) == (1, 10, 1, 11, 1),
        "graded central digits changed",
    )
    print(
        "graded_chord="
        f"address_gap={address_gap}=13*{central_gap}; "
        f"central_exponent_digits={base_p_digits(central_gap, 5)}; "
        f"Z1^{central_gap}=Z1*Z2^10*Z3*Z4^11*Z5; "
        f"first_transgression=(central_gap-1)/13 mod13="
        f"{correction % P}"
    )
    print(
        "group_ring_kernel_factor="
        f"u^{central_gap}-u="
        f"u*(u^13-1)*(1+u^13+...+u^(13*{correction - 1})); "
        "divisible_by_u^13-1=yes; divisible_by_u^169-1=no"
    )

    # The two selected cylinders do carry an exact partial translation germ,
    # but its gain is the full graded chord rather than the pure Z_1 step.
    delayed3, weighted3 = target_weighted_by_tau(context, 3)

    def restricted_piece(index):
        center = thm.arm_target_center(base, Fraction(index, 13))
        half = thm.relative.Q_RADIUS * thm.T
        cylinder = ((center * thm.T - half, center * thm.T + half),)
        pieces = thm.marked.restrict_weighted(weighted3, cylinder)
        require(
            len(pieces) == 1 and pieces[0][:2] == cylinder[0],
            "selected chord cylinder stopped being a whole weighted piece",
        )
        return pieces[0]

    piece0 = restricted_piece(0)
    piece1 = restricted_piece(address_gap)
    physical_shift = Fraction(7 * address_gap, thm.ADDRESS_MODULUS)
    shift_coordinate = physical_shift * thm.T
    require(
        physical_shift == Fraction(371196, 371293)
        and physical_shift * thm.ADDRESS_MODULUS == 4825548
        and physical_shift * (thm.ADDRESS_MODULUS // P) == 371196,
        "graded chord translation lost delayed/carry integrality",
    )
    translated_left = (piece0[0] + shift_coordinate) % thm.T
    translated_right = (piece0[1] + shift_coordinate) % thm.T
    require(
        (translated_left, translated_right, piece0[2]) == piece1,
        "graded chord no longer translates the whole weighted cylinder",
    )
    require(
        thm.marked.private.delayed_carry_pair(
            (piece0,), delayed3[1], {}
        )[6]
        == thm.marked.private.delayed_carry_pair(
            (piece1,), delayed3[1], {}
        )[6]
        == (0, thm.EXPECTED_CONTENTS[1]),
        "graded chord lost its pointwise-periodic delayed coefficient",
    )
    print(
        "partial_physical_germ="
        f"shift={physical_shift}; R*shift=4825548; "
        "13^5*shift=371196; whole cylinders translate exactly with "
        f"common weight={piece0[2]} and carry-six pair="
        f"{(0, thm.EXPECTED_CONTENTS[1])}; "
        "gain is graded Z1^53028, not pure Z1"
    )

    base_residue = base % (P**2)
    other_residue = (base + address_gap) % (P**2)
    require(
        (base_residue, other_residue) == (85, 98),
        "Omega quotient support changed",
    )
    print(
        "Omega_transfer_tau3="
        f"raw support residues=({base_residue},{other_residue})="
        "((v,w)=(7,6),(7,7)); normalized central chain=1+z"
    )
    two_point = (1, 1) + (0,) * 11
    alternating = tuple(1 if index % 2 == 0 else -1 for index in range(P))
    transfer_inverse_product = tuple(
        sum(
            two_point[index]
            * alternating[(degree - index) % P]
            for index in range(P)
        )
        for degree in range(P)
    )
    require(
        transfer_inverse_product == (2,) + (0,) * 12,
        "coarse two-point transfer stopped being a unit away from char 2",
    )
    print(
        "Omega_transfer_unit="
        "(1+z)*(1-z+z^2-...+z^12)=2; "
        "inverse=(1/2)*(1-z+...+z^12) in every field of char!=2"
    )

    decoder_pairs = tuple(
        (
            sum(
                thm.K_BETA[(target - tau) % P]
                for tau in range(3, 12)
            ),
            thm.K_BETA[(target - 12) % P],
        )
        for target in range(P)
    )
    require(
        decoder_pairs
        == (
            (58, 5), (64, 5), (59, 5), (55, 7), (48, 12),
            (51, 2), (48, 8), (53, 2), (56, 9), (50, 8),
            (47, 11), (49, 0), (55, 3),
        ),
        "target decoder coefficient pairs changed",
    )
    print(
        "K_beta_full_orbit_pairs_(A_for_chord,B_for_tau12)="
        f"{decoder_pairs}; every decoded target retains primitive "
        "conductors 13,13^2,...,13^5"
    )

    # Every nonidentity lower-central translation Z_r=O^(13^r),
    # 1<=r<=5, has four-term signed finite difference on this two-point
    # current.  Z_6 is identity.
    full_support = {base % thm.ADDRESS_MODULUS,
                    (base + address_gap) % thm.ADDRESS_MODULUS}
    for depth in range(1, 7):
        shift = P**depth % thm.ADDRESS_MODULUS
        shifted = {
            (address + shift) % thm.ADDRESS_MODULUS
            for address in full_support
        }
        difference = {}
        for address in shifted:
            difference[address] = difference.get(address, 0) + 1
        for address in full_support:
            difference[address] = difference.get(address, 0) - 1
        difference = {
            address: value for address, value in difference.items() if value
        }
        expected_terms = 4 if depth <= 5 else 0
        require(
            len(difference) == expected_terms,
            f"lower-central difference Z_{depth} changed",
        )
        l1 = sum(abs(value) for value in difference.values())
        l2_squared = sum(value * value for value in difference.values())
        print(
            f"lower_central_difference_Z{depth}: "
            f"support_terms={len(difference)} l1={l1} l2sq={l2_squared}"
        )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
