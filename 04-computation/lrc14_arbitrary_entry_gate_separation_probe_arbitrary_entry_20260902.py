#!/usr/bin/env python3
"""Exact small controls for the LRC(14) degree-two arbitrary-entry frontier.

The probe separates five operations that are easy to conflate: projective
fixed-pool entry, pair-adaptive Haar entry, owner-word entry, deletion, and
change of reference.  All arithmetic is exact.
"""

from fractions import Fraction
from math import gcd


POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
}
THRESHOLD = Fraction(1, 14)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_distance(speed: int, time: Fraction) -> Fraction:
    phase = speed * time
    residue = Fraction(phase.numerator % phase.denominator, phase.denominator)
    return min(residue, 1 - residue)


def safe_measure(body: tuple[int, ...]) -> Fraction:
    cuts = {Fraction(0), Fraction(1)}
    for speed in body:
        for tooth in range(speed):
            cuts.add(Fraction(14 * tooth + 1, 14 * speed))
            cuts.add(Fraction(14 * tooth + 13, 14 * speed))
    ordered = sorted(cuts)
    measure = Fraction(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if all(circle_distance(speed, midpoint) >= THRESHOLD for speed in body):
            measure += right - left
    return measure


def projective_pool_max(body: tuple[int, ...]) -> tuple[int, tuple[Fraction, ...]]:
    candidates = {Fraction(pool_label, speed) for pool_label in POOL for speed in body}
    best = -1
    scales: list[Fraction] = []
    for scale in sorted(candidates):
        image = tuple(scale * speed for speed in body)
        if any(value.denominator != 1 for value in image):
            continue
        hits = sum(int(value) in POOL for value in image)
        if hits > best:
            best = hits
            scales = [scale]
        elif hits == best:
            scales.append(scale)
    return best, tuple(scales)


def bernoulli_two(value: Fraction) -> Fraction:
    value -= value.numerator // value.denominator
    return value * value - value + Fraction(1, 6)


def cross_comb_measure(p: int, q: int) -> Fraction:
    require(0 < p < q and p % 2 == q % 2 == 1 and gcd(p, q) == 1, "tail ratio")
    return (
        Fraction(2, 49)
        + 2
        * (
            bernoulli_two(Fraction(1, 2) + Fraction(q - p, 14))
            - bernoulli_two(Fraction(1, 2) + Fraction(q + p, 14))
        )
        / (p * q)
    )


def abs_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def owner_packet(body: tuple[int, ...], clock: int) -> tuple[int, ...]:
    return tuple(
        residue
        for residue in range(clock)
        if all(14 * abs_residue(speed * residue, clock) >= clock for speed in body)
    )


def literal_universal_owner_certificate(body: tuple[int, ...], clock: int) -> bool:
    packet = owner_packet(body, clock)
    if not packet:
        return False
    odd_residues = tuple(range(1, 2 * clock, 2))
    for first_index, first in enumerate(odd_residues):
        for second in odd_residues[first_index:]:
            spoiled = True
            for residue in packet:
                for lift in (0, 1):
                    numerator = residue + lift * clock
                    first_bad = (
                        14 * abs_residue(first * numerator, 2 * clock) < 2 * clock
                    )
                    second_bad = (
                        14 * abs_residue(second * numerator, 2 * clock) < 2 * clock
                    )
                    if not first_bad and not second_bad:
                        spoiled = False
                        break
                if not spoiled:
                    break
            if spoiled:
                return False
    return True


def affine_content(velocities: tuple[int, ...]) -> int:
    base = min(velocities)
    content = 0
    for velocity in velocities:
        content = gcd(content, velocity - base)
    return content


def degree_two_row_at_reference(
    velocities: tuple[int, ...], reference: int
) -> tuple[tuple[int, ...], tuple[int, int]]:
    content = affine_content(velocities)
    row = tuple(
        sorted(abs(velocity - reference) // content for velocity in velocities if velocity != reference)
    )
    require(len(row) == len(set(row)) == 13, "collision-free reference")
    body = tuple(value // 2 for value in row if value % 2 == 0)
    tails = tuple(value for value in row if value % 2 == 1)
    require(len(body) == 11 and len(tails) == 2, "degree-two reference")
    return body, (tails[0], tails[1])


def main() -> None:
    ap_body = tuple(range(1, 12))
    separating_body = (1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14)

    ap_mass = safe_measure(ap_body)
    separating_mass = safe_measure(separating_body)
    worst_cross_mass = cross_comb_measure(1, 9)
    require(ap_mass == Fraction(10931, 194040), "AP body mass")
    require(separating_mass == Fraction(20411, 360360), "separating body mass")
    require(worst_cross_mass == Fraction(4, 63), "worst cross-comb mass")
    require(ap_mass < worst_cross_mass and separating_mass < worst_cross_mass, "mass failures")

    ap_hits, ap_scales = projective_pool_max(ap_body)
    separating_hits, separating_scales = projective_pool_max(separating_body)
    require((ap_hits, ap_scales) == (6, (Fraction(10),)), "AP projective profile")
    require(
        (separating_hits, separating_scales) == (6, (Fraction(10),)),
        "separating projective profile",
    )

    packet_23 = owner_packet(separating_body, 23)
    require(packet_23 == (4, 9, 10, 13, 14, 19), "clock-23 owner packet")
    require(
        literal_universal_owner_certificate(separating_body, 23),
        "clock-23 universal two-tail certificate",
    )
    separating_row = tuple(2 * value for value in separating_body) + (1, 9)
    require(
        min(circle_distance(speed, Fraction(2, 23)) for speed in separating_row)
        == Fraction(2, 23),
        "clock-23 concrete witness",
    )

    # Deletion reverses the useful containment: 1/11 is safe for {1,...,10}
    # but the deleted speed 11 is exactly at zero there.
    deleted_body = tuple(range(1, 11))
    delete_time = Fraction(1, 11)
    require(
        min(circle_distance(speed, delete_time) for speed in deleted_body)
        == Fraction(1, 11),
        "deleted body witness",
    )
    require(circle_distance(11, delete_time) == 0, "deleted constraint failure")

    # THM-4332's explicit outside-label point: even the entire pool does not
    # pointwise imply the singleton speed 1.
    pool_hostile = Fraction(199, 21280)
    require(
        min(circle_distance(speed, pool_hostile) for speed in POOL) >= THRESHOLD,
        "pool-safe hostile",
    )
    require(circle_distance(1, pool_hostile) < THRESHOLD, "outside singleton danger")

    # Same physical velocity set, two collision-free degree-two references.
    # Reference 22 passes the adaptive mass test; reference 0 fails it.
    velocities = tuple(range(0, 23, 2)) + (1, 9)
    body_zero, tails_zero = degree_two_row_at_reference(velocities, 0)
    body_twenty_two, tails_twenty_two = degree_two_row_at_reference(velocities, 22)
    require(body_zero == body_twenty_two == ap_body, "same AP half-body")
    require(tails_zero == (1, 9) and tails_twenty_two == (13, 21), "reference tails")
    require(cross_comb_measure(*tails_zero) == Fraction(4, 63), "reference-zero threshold")
    require(
        cross_comb_measure(*tails_twenty_two) == Fraction(2, 49),
        "reference-twenty-two threshold",
    )
    require(Fraction(2, 49) < ap_mass < Fraction(4, 63), "opposite entry verdicts")
    require(
        min(
            circle_distance(speed, Fraction(5, 24))
            for speed in tuple(2 * value for value in ap_body) + tails_zero
        )
        == Fraction(1, 12),
        "reference-zero physical witness",
    )

    last_clock_bodies = (
        (11, 14, 19, 23, 25, 26, 29, 34, 36, 37, 40),
        (12, 18, 19, 22, 25, 26, 28, 31, 34, 37, 40),
        (19, 22, 25, 26, 28, 29, 31, 34, 36, 37, 40),
    )
    expected_last_clock_masses = (
        Fraction(2331354986167, 16756505465700),
        Fraction(11062658951, 106198378650),
        Fraction(1400226987157, 10423779319800),
    )
    last_clock_profiles = []
    for body, expected_mass in zip(last_clock_bodies, expected_last_clock_masses):
        mass = safe_measure(body)
        hits, scales = projective_pool_max(body)
        require(mass == expected_mass and mass > worst_cross_mass, "last-clock mass")
        require(hits == 3, "last-clock projective hits")
        last_clock_profiles.append((body, mass, scales))

    print("LRC14 ARBITRARY-ENTRY GATE SEPARATION -- EXACT CONTROLS")
    print(
        "combined_gate_survivor_body=1,2,3,4,8,9,10,11,12,13,14 "
        f"projective_pool_max_hits={separating_hits} at_scale=10 "
        f"body_mass={separating_mass} cross_mass_1_9={worst_cross_mass}"
    )
    print(
        "owner_word_repair=clock_23 packet="
        + ",".join(map(str, packet_23))
        + " universal_over_all_odd_tail_residue_pairs=1 "
        "control_tails=1,9 witness=2/23 clearance=2/23"
    )
    print(
        f"ap_gate_survivor_body=1_to_11 projective_pool_max_hits={ap_hits} "
        f"at_scale=10 body_mass={ap_mass} cross_mass_1_9={worst_cross_mass} "
        "witness=5/24 clearance=1/12"
    )
    print(
        "deletion_direction_hostile=G_1_to_11_subset_G_1_to_10 "
        "time=1/11 deleted_body_clearance=1/11 restored_speed_11_clearance=0"
    )
    print(
        "literal_pool_nonentry=point_199/21280_in_G_POOL_and_D_1 "
        "so_full_pool_does_not_imply_outside_singleton_1"
    )
    print(
        "reference_hostile=V_even_0_to_22_plus_1,9 "
        "anchor_0_tails=1,9 cross_mass=4/63 gate=FAIL "
        "anchor_22_tails=13,21 cross_mass=2/49 gate=PASS "
        "anchor_change_does_not_transfer_target"
    )
    for body, mass, scales in last_clock_profiles:
        print(
            "clock_43_boundary_is_mass_easy=body_"
            + ",".join(map(str, body))
            + f" body_mass={mass} gt_4/63=1 projective_hits=3 "
            + "projective_scales="
            + ",".join(str(scale) for scale in scales)
        )
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
