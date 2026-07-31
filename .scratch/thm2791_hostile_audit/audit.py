#!/usr/bin/env python3
"""Independent hostile audit for THM-2791.

This companion does not import the THM-2791 script.  It rebuilds the physical
row from the promoted THM-2782 dependency, scans the full address orbit with
an independently written interval-overlap engine, and rechecks the group,
Fourier/conductor, quotient, decoder, and scope-critical claims.
"""

from bisect import bisect_left
from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_semantic_arm_right_wing_central_digit_thm2782 as dep


P = 13
R = P**6
N = P**5
D = P**5
T = dep.T
C = dep.EXPECTED_CONTENTS[1]
BASE = dep.EXPECTED_BASES[1]
K_BETA = dep.K_BETA


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(relative_path):
    data = (ROOT / relative_path).read_bytes().replace(b"\r\n", b"\n")
    return sha256(data).hexdigest()


EXPECTED_HASHES = {
    # Candidate.
    "04-computation/lrc14_full_arm_orbit_lower_central_chord_thm2791.py":
        "fc557d31a82ea52adf0abc7b26bfdacf1961facc99050061132d72dda231a9db",
    "05-knowledge/results/lrc14_full_arm_orbit_lower_central_chord_thm2791.out":
        "9f2b8e69b9de430f201adb7758f98fff7bf505c5bd792b03f40b3ee7c9f46edd",
    # Direct promoted physical dependency.
    "04-computation/lrc14_semantic_arm_right_wing_central_digit_thm2782.py":
        "7fbc6bb1ec303ded98eaad6e5d8205eb3d247258ada32b6f9904fc439ebb11fb",
    "05-knowledge/results/lrc14_semantic_arm_right_wing_central_digit_thm2782.out":
        "13f570d63f212808171cecdb4d8f9aa41884fbdc7ed571dbfe27122b412fadc4",
    # Direct source imports used by THM-2782.
    "04-computation/lrc14_fully_marked_root_zero_target_profile_thm2749.py":
        "d67c852c52f88feaadb2fcaa0a9a07a212f2e47018040b455855df886200595e",
    "04-computation/lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.py":
        "25cbed38026d61891173c687006250a69fe38aea56d67439406bd8bb60fa2552",
    "04-computation/lrc14_relative_present_semantic_lift_probe_20260728.py":
        "f16754bd38ae0dfa0d7d91cc404b4447dbf359635101aa7b4223363f8064352f",
    # THM-2779 coefficient/endpoint dependency.
    "04-computation/lrc14_bockstein_symplectic_heisenberg_gate_thm2779.py":
        "4c6a58c80ddd4be0fd9bdd297b310df054bbc08996eb223d519d3cce6b8ed13a",
    "05-knowledge/results/lrc14_bockstein_symplectic_heisenberg_gate_thm2779.out":
        "f7c96259777a3ab4a5e46cac8666181ae77a3be2e440cee8785997507706791a",
    "04-computation/lrc14_bockstein_symplectic_endpoint_gate_thm2779.py":
        "004e06c617f9305e2f0bc30871926e3faa7843f47dcf63af1fd8a892e63101e4",
    "05-knowledge/results/lrc14_bockstein_symplectic_endpoint_gate_thm2779.out":
        "a89c00a3830ee9ff282cc5e4557d41293af5d6f0e7feabd5d3c7e7808591e754",
    "04-computation/lrc14_heisenberg_decoder_frame_root_degree_thm2779.py":
        "ef6e9f9bcb4f11152d291342a11ae215245d1d19b96c49940a01ba9ea850cbd9",
    "05-knowledge/results/lrc14_heisenberg_decoder_frame_root_degree_thm2779.out":
        "1feb463864015035ab8d7fcfcddf9cfe8b0ec0a3ed36481f2f66d6a9149182e6",
    "04-computation/lrc14_heisenberg_decoder_frame_root_degree_hostile_audit_thm2779.py":
        "5019f87b24500a5a13825d3be01908ea983a08b360a384fd614107f476201f46",
    "05-knowledge/results/lrc14_heisenberg_decoder_frame_root_degree_hostile_audit_thm2779.out":
        "1b9ad37b35e92a14dd90d0db8c1d0cf225761e2c37ca8e2fe2120bd0f64c47d4",
    # THM-2788 convention/control.
    "04-computation/lrc14_modular_odometer_heisenberg_bockstein_thm2788.py":
        "d414bf2afb6aa3e40de9378ae20f03db1cb7bff75f59f13a60ac96e56cb95a89",
    "05-knowledge/results/lrc14_modular_odometer_heisenberg_bockstein_thm2788.out":
        "99ad33904617d45d76a285de5467b96408dc164839cb4168905c7fe678db8f66",
}


def rebuild_weighted_row(context, tau):
    """Rebuild THM-2791's target row without importing its implementation."""
    module, delayed, present, source, clocks, source_weight, rail_common = context
    objects = dep.carrier_objects(
        module, source, clocks, rail_common, 1, 0, tau
    )
    source_safe = dep.marked.complement(present[1, 7])
    target_safe = dep.marked.complement(
        dep.marked.shift_union(present[1, 7], dep.marked.SHIFT)
    )
    carrier = dep.marked.intersect(objects["R"], source_safe)
    carrier = dep.marked.intersect(carrier, target_safe)
    carrier = dep.marked.shift_union(carrier, -dep.marked.SHIFT)
    weighted = dep.marked.restrict_weighted(source_weight, carrier)
    weighted = dep.marked.private.old.intersect_weighted_comb(
        weighted, module.C3, 182, 1, 13
    )
    return delayed, tuple(weighted)


def scan_interval_hits(weighted, count, address_stride):
    """Scan an entire cyclic address orbit by scaled exact interval overlap.

    Weighted intervals are duplicated in the adjacent two periods.  A binary
    search finds the first right endpoint strictly beyond the left endpoint
    of each open cylinder.  This differs from the candidate engine in the
    audited source and checks every requested orbit point.
    """
    half = dep.relative.Q_RADIUS * T
    require(half.denominator == D, "unexpected cylinder denominator")
    half_num = half.numerator
    center = dep.arm_target_center(BASE, 0) * T
    require(center.denominator == D, "unexpected center denominator")
    center_num = center.numerator
    modulus = T * D
    step = (address_stride * dep.EXPECTED_SHIFT * D) % modulus
    base_intervals = tuple((left * D, right * D) for left, right, _ in weighted)
    require(
        all(base_intervals[i][1] <= base_intervals[i + 1][0]
            for i in range(len(base_intervals) - 1)),
        "weighted intervals are not sorted/disjoint",
    )
    extended = (
        tuple((left - modulus, right - modulus) for left, right in base_intervals)
        + base_intervals
        + tuple((left + modulus, right + modulus) for left, right in base_intervals)
    )
    rights = tuple(right for _left, right in extended)
    hits = []
    for index in range(count):
        lo = center_num - half_num
        hi = center_num + half_num
        pos = bisect_left(rights, lo + 1)
        if pos < len(extended) and extended[pos][0] < hi:
            hits.append(index)
        center_num = (center_num + step) % modulus
    return tuple(hits)


def exact_cell(delayed, weighted, address_offset):
    """Recompute mass and delayed carry-six/root-one coefficient."""
    center = dep.arm_target_center(BASE, Fraction(address_offset, P))
    half = dep.relative.Q_RADIUS * T
    coordinate = center * T
    if coordinate - half < 0:
        cylinders = (
            (0, coordinate + half),
            (T + coordinate - half, T),
        )
    elif coordinate + half > T:
        cylinders = (
            (coordinate - half, T),
            (0, coordinate + half - T),
        )
    else:
        cylinders = ((coordinate - half, coordinate + half),)
    pieces = dep.marked.restrict_weighted(weighted, cylinders)
    return (
        dep.weighted_mass(pieces),
        dep.marked.private.delayed_carry_pair(
            pieces, delayed[1], {}
        )[6][1],
        tuple(pieces),
    )


def translation_stabilizer(support, modulus):
    support = frozenset(support)
    if not support:
        return tuple(range(modulus))
    anchor = min(support)
    candidates = {(value - anchor) % modulus for value in support}
    return tuple(sorted(
        shift for shift in candidates
        if frozenset((value + shift) % modulus for value in support) == support
    ))


def conductor_survival(profile, max_exponent=5):
    """Cyclotomic block criterion for an integral C_(13^5) profile."""
    result = []
    for exponent in range(1, max_exponent + 1):
        modulus = P**exponent
        stride = P ** (exponent - 1)
        counts = [0] * modulus
        for address, value in profile.items():
            counts[address % modulus] += value
        vanishes_at_primitive = all(
            len({counts[a + digit * stride] for digit in range(P)}) == 1
            for a in range(stride)
        )
        result.append(not vanishes_at_primitive)
    return tuple(result)


def cyclic_convolution(left, right):
    return tuple(
        sum(
            left[index] * right[(degree - index) % P]
            for index in range(P)
        )
        for degree in range(P)
    )


def base_p_digits(value, length):
    digits = []
    for _ in range(length):
        digits.append(value % P)
        value //= P
    require(value == 0, "base-p digit buffer too short")
    return tuple(digits)


def main():
    print("THM-2791 INDEPENDENT HOSTILE AUDIT")

    actual_hashes = {path: lf_hash(path) for path in EXPECTED_HASHES}
    require(actual_hashes == EXPECTED_HASHES, "candidate/dependency hash drift")
    print(f"hash_pins={len(actual_hashes)} PASS")

    context = dep.build_context()
    require(
        Fraction(7 * T, R) == dep.EXPECTED_SHIFT,
        "address step was not independently recovered from 7*T/13^6",
    )
    profiles = {}
    delayed = None
    for tau in range(P):
        delayed_tau, weighted = rebuild_weighted_row(context, tau)
        if delayed is None:
            delayed = delayed_tau
        else:
            require(delayed_tau == delayed, "delayed carrier changed by target")
        profiles[tau] = weighted

    require(
        all(profiles[tau] == profiles[3] for tau in range(3, 12))
        and len(profiles[3]) == 2
        and all(not profiles[tau] for tau in range(3))
        and len(profiles[12]) == 241,
        "raw target carrier partition changed",
    )

    # Full 13^6 scan for the common main profile.
    full_hits = scan_interval_hits(profiles[3], R, 1)
    require(full_hits == (0, 689364), "full address chord changed")
    full_cells = tuple(
        (offset, *exact_cell(delayed, profiles[3], offset)[:2])
        for offset in full_hits
    )
    require(
        full_cells
        == (
            (0, Fraction(C, R), C),
            (689364, Fraction(C, R), C),
        ),
        "full address mass/coefficient changed",
    )
    print(f"full_address_scan={R} support={full_hits} PASS")

    # The exceptional target is claimed only on the central 13^5 orbit.
    exceptional_indices = scan_interval_hits(profiles[12], N, P)
    exceptional_cells = tuple(
        (index, *exact_cell(delayed, profiles[12], P * index)[:2])
        for index in exceptional_indices
    )
    require(
        len(exceptional_indices) == 121
        and all(mass == Fraction(C, R) and coefficient == C
                for _index, mass, coefficient in exceptional_cells),
        "exceptional central profile changed",
    )
    chord_indices = (0, 53028)
    require(
        tuple(offset // P for offset in full_hits) == chord_indices,
        "central chord indices changed",
    )
    print(
        f"central_scans=main:{chord_indices},exceptional:{len(exceptional_indices)}"
    )

    # Absolute address and physical typing.
    absolute = tuple((BASE + offset) % R for offset in full_hits)
    require(absolute == (3454614, 4143978), "absolute addresses changed")
    metadata = []
    for offset in full_hits:
        center = dep.arm_target_center(BASE, Fraction(offset, P))
        sigma_labels, tau_labels = dep.relative.full_target_labels(center)
        metadata.append((
            dep.omega_digits(BASE + offset),
            dep.relative.semantic.semantic_record(center),
            dep.predecessor_carry(center),
            0 in sigma_labels,
            tuple(tau in tau_labels for tau in range(3, 12)),
            dep.relative.semantic_stability_radius(center)
                == dep.relative.Q_RADIUS,
        ))
    require(
        metadata
        == [
            ((7, 6), (3, (1, 2)), 6, True, (True,) * 9, True),
            ((7, 7), (3, (1, 2)), 6, True, (True,) * 9, True),
        ],
        "semantic/label/carry/stability metadata changed",
    )
    print(f"two_point_metadata={tuple(metadata)} PASS")

    # Periods and cyclotomic conductors.
    require(
        translation_stabilizer(full_hits, R) == (0,)
        and translation_stabilizer(chord_indices, N) == (0,)
        and translation_stabilizer(exceptional_indices, N) == (0,),
        "a claimed full period has a hidden stabilizer",
    )
    chord_profile = {0: 1, 53028: 1}
    exceptional_profile = {index: 1 for index in exceptional_indices}
    require(
        R % 2 == N % 2 == 1
        and conductor_survival(chord_profile) == (True,) * 5
        and conductor_survival(exceptional_profile) == (True,) * 5,
        "period/conductor claim failed",
    )
    # Every character value 1+chi(gap)^(-1) is nonzero because its second
    # term has odd order.  Record the exact order bound and the primitive
    # conductor controls rather than using floating roots of unity.
    require(
        all(
            (R // gcd(R, character * full_hits[1])) % 2 == 1
            for character in range(R)
        ),
        "a full-address character unexpectedly has even gap order",
    )
    print("periods=13^6,13^5; all_full_characters_and_conductors=PASS")

    # Lower-central convention from THM-2788.
    gap = full_hits[1]
    a = gap // P
    correction = (a - 1) // P
    digits = base_p_digits(a, 5)
    require(
        gap == 689364
        and a == 53028
        and digits == (1, 10, 1, 11, 1)
        and a == sum(digit * P**index for index, digit in enumerate(digits))
        and correction == 4079
        and correction % P == 10
        and (a - 1) % P == 0
        and (a - 1) % (P**2) != 0,
        "lower-central chord/transgression changed",
    )
    # X O^(p^r) X^-1 O^(-p^r)=O^(p^(r+1)) for X:n->14n.
    require(
        all(
            ((1 + P) * P**depth - P**depth) % R
            == P ** (depth + 1) % R
            for depth in range(6)
        ),
        "THM-2788 commutator convention mismatch",
    )
    word_exponent = (
        P + 10 * P**2 + P**3 + 11 * P**4 + P**5
    ) % R
    require(word_exponent == gap, "little-endian Z_r word changed")
    print(
        f"lower_central=Z1*Z2^10*Z3*Z4^11*Z5; "
        f"first_transgression={correction % P} PASS"
    )

    # Exact partial physical germ, independently reconstructed.
    pieces = tuple(
        exact_cell(delayed, profiles[3], offset)[2][0]
        for offset in full_hits
    )
    require(all(len(exact_cell(delayed, profiles[3], offset)[2]) == 1
                for offset in full_hits),
            "a selected cylinder is no longer a single whole piece")
    physical_shift = Fraction(7 * gap, R)
    coordinate_shift = physical_shift * T
    half = dep.relative.Q_RADIUS * T
    expected_pieces = tuple(
        (
            dep.arm_target_center(BASE, Fraction(offset, P)) * T - half,
            dep.arm_target_center(BASE, Fraction(offset, P)) * T + half,
            27581135604,
        )
        for offset in full_hits
    )
    translated = (
        (pieces[0][0] + coordinate_shift) % T,
        (pieces[0][1] + coordinate_shift) % T,
        pieces[0][2],
    )
    require(
        physical_shift == Fraction(371196, 371293)
        and physical_shift * R == 4825548
        and physical_shift * N == 371196
        and pieces == expected_pieces
        and translated == pieces[1]
        and pieces[0][2] == pieces[1][2] == 27581135604,
        "partial graded-chord germ changed",
    )
    print(f"partial_germ_shift={physical_shift} whole_piece_transfer=PASS")

    # Quotient pushforward and the distinction from descent.
    residues = tuple(address % (P**2) for address in absolute)
    require(residues == (85, 98), "mod-169 pushforward support changed")
    require(
        0 in full_hits and 169 not in full_hits,
        "descent hostile F(0)!=F(169) disappeared",
    )
    quotient_rows = []
    supports = {
        **{tau: tuple() for tau in range(3)},
        **{tau: chord_indices for tau in range(3, 12)},
        12: exceptional_indices,
    }
    for tau in range(P):
        quotient_rows.append(tuple(
            sum(1 for index in supports[tau] if index % P == residue)
            for residue in range(P)
        ))
    require(
        quotient_rows[3] == (1, 1) + (0,) * 11
        and quotient_rows[12] == (10, 10, 10, 10) + (9,) * 9,
        "quotient pushforward normalization changed",
    )

    two_point = (1, 1) + (0,) * 11
    inverse_lift = tuple(1 if index % 2 == 0 else -1 for index in range(P))
    require(
        cyclic_convolution(two_point, inverse_lift) == (2,) + (0,) * 12,
        "1+z group-ring inverse identity failed",
    )
    value = C
    valuation = 0
    while value % P == 0:
        value //= P
        valuation += 1
    require(
        valuation == 1
        and (C // P) % P == 2
        and sum(C * entry for entry in two_point) != 1,
        "scalar/Bockstein normalization boundary changed",
    )
    print(
        "pushforward=c*z^6*(1+z); normalized_unit=PASS; "
        "raw_Z_group_ring_unit=NO; descent=NO"
    )

    decoded_quotient = tuple(
        cyclic_convolution(tuple(quotient_rows[tau][residue]
                                 for tau in range(P)), K_BETA)
        for residue in range(P)
    )
    # Transpose back to target-major rows.
    decoded_target = tuple(
        tuple(decoded_quotient[residue][target] for residue in range(P))
        for target in range(P)
    )
    require(
        all(len(set(row)) > 1 for row in decoded_target),
        "a decoded quotient target became constant",
    )

    decoder_pairs = tuple(
        (
            sum(K_BETA[(target - tau) % P] for tau in range(3, 12)),
            K_BETA[(target - 12) % P],
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
        "integral K_beta decoder convention changed",
    )

    decoded_sparse = []
    for target in range(P):
        row = {}
        for tau in range(P):
            multiplier = K_BETA[(target - tau) % P]
            for index in supports[tau]:
                row[index] = row.get(index, 0) + multiplier
        decoded_sparse.append({index: value for index, value in row.items()
                               if value})
    require(
        all(conductor_survival(row) == (True,) * 5
            for row in decoded_sparse),
        "a decoded target lost a primitive conductor",
    )
    print("exceptional_profile_and_all_13_decoded_targets=PASS")

    # Every nonidentity lower-central finite difference has four cells.
    absolute_support = frozenset(absolute)
    difference_sizes = []
    for depth in range(1, 7):
        shift = P**depth % R
        shifted = frozenset((address + shift) % R for address in absolute_support)
        signed = {}
        for address in shifted:
            signed[address] = signed.get(address, 0) + 1
        for address in absolute_support:
            signed[address] = signed.get(address, 0) - 1
        signed = {address: value for address, value in signed.items() if value}
        expected = 4 if depth <= 5 else 0
        require(
            len(signed) == expected
            and sum(abs(value) for value in signed.values()) == expected
            and sum(value * value for value in signed.values()) == expected,
            f"Z_{depth} finite difference changed",
        )
        difference_sizes.append(len(signed))
    print(f"lower_central_difference_sizes={tuple(difference_sizes)} PASS")

    print(
        "scope=TRANSFER_BY_PUSHFORWARD_PLUS_PARTIAL_GERM;"
        "NOT_DESCENT;NOT_PURE_Z1_COVARIANCE;NOT_ENDPOINT_ALLOCATION"
    )
    print("verdict=PASS;PROMOTE_AFTER_INDEPENDENT_AUDIT_METADATA")


if __name__ == "__main__":
    main()
