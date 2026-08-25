#!/usr/bin/env python3
"""Exact affine-window and small-denominator sieve certificate for THM-4070.

This primary path has two independent layers inside one executable:

* it reconstructs the primitive d=2 defect intervals, compares them with
  literal strict wall cells, and checks common odd-dilation pullback; and
* it proves the complete q=2,...,14 divided-pack sieve by exact rational
  masks, then exercises the mod-14 strict gain and physical controls stated
  in the theorem.

The universal mod-14 theorem is elementary; the bounded window census is a
hostile audit of the inherited affine formula, not the source of its
quantifiers.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")

RHO = F(1, 14)
CHECKS = 0

SIEVE_K = {8: 3, 9: 2, 10: 3, 11: 2, 12: 5, 13: 2, 14: 3}
EXPECTED_DANGER_CLASSES = {
    8: ((1, 15), (7, 9), (5, 11), (3, 13)),
    9: ((1, 17), (9,), (9,), (5, 13)),
    10: ((1, 19), (9, 11), (7, 13), (3, 17)),
    11: ((1, 21), (11,), (11,), (5, 17)),
    12: ((1, 23), (11, 13), (5, 19), (7, 17)),
    13: ((1, 25), (13,), (13,), (7, 19)),
    14: ((1, 27), (13, 15), (9, 19), (5, 23)),
}


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def distance(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def danger_masks(y: F, exceptions: tuple[int, int]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            label
            for label in (0, 1)
            if distance(speed * (y + label) / 2) < RHO
        )
        for speed in exceptions
    )


def fully_spoiled(y: F, exceptions: tuple[int, int]) -> bool:
    return set().union(*(set(mask) for mask in danger_masks(y, exceptions))) == {0, 1}


def wall_points(exceptions: tuple[int, int]) -> tuple[F, ...]:
    walls = {F(0), F(1)}
    for speed in exceptions:
        for label in (0, 1):
            for integer in range(-2, speed + 3):
                for sign in (-1, 1):
                    wall = F(2, speed) * (F(integer) + sign * RHO) - label
                    if 0 <= wall <= 1:
                        walls.add(wall)
    return tuple(sorted(walls))


def wall_components(exceptions: tuple[int, int]) -> tuple[tuple[F, F], ...]:
    """Strict fully-spoiled components in the linear gauge 0<=y<=1."""
    walls = wall_points(exceptions)
    components: list[tuple[F, F]] = []
    for left, right in zip(walls, walls[1:]):
        if not fully_spoiled((left + right) / 2, exceptions):
            continue
        if (
            components
            and components[-1][1] == left
            and fully_spoiled(left, exceptions)
        ):
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    return tuple(components)


def inverse_components(
    components: tuple[tuple[F, F], ...], dilation: int
) -> tuple[tuple[F, F], ...]:
    return tuple(
        sorted(
            ((left + sheet) / dilation, (right + sheet) / dilation)
            for left, right in components
            for sheet in range(dilation)
        )
    )


def primitive_defect_windows(a: int, b: int) -> tuple[tuple[int, F, F], ...]:
    """Return the lifted z-interval J_n for each attainable signed defect."""
    require(a < b and a % 2 == b % 2 == 1 and gcd(a, b) == 1, ("primitive", a, b))
    windows = []
    for defect in range(-(a + b) + 1, a + b, 2):
        if 7 * abs(defect) >= a + b:
            continue
        target = (defect - a * b) // 2
        require(2 * target == defect - a * b, ("defect parity", a, b, defect))
        if a == 1:
            centre_a_numerator = 0
        else:
            centre_a_numerator = (target * pow(b, -1, a)) % a
        centre_b_numerator = (b * centre_a_numerator - target) // a
        require(
            b * centre_a_numerator - a * centre_b_numerator == target,
            ("Bezout reconstruction", a, b, defect),
        )
        centre_a = F(centre_a_numerator, a)
        centre_b = F(centre_b_numerator, b) - F(1, 2)
        left = max(centre_a - F(1, 14 * a), centre_b - F(1, 14 * b))
        right = min(centre_a + F(1, 14 * a), centre_b + F(1, 14 * b))
        require(left < right, ("positive defect interval", a, b, defect, left, right))
        expected_z_length = min(
            F(1, 7 * b),
            F(a + b - 7 * abs(defect), 14 * a * b),
        )
        require(right - left == expected_z_length, ("defect length", a, b, defect))
        windows.append((defect, left, right))
    return tuple(windows)


def in_doubled_open_interval(y: F, left: F, right: F) -> bool:
    """Whether y lies in {2z mod 1:left<z<right}."""
    for preimage in (y / 2, (y + 1) / 2):
        for shift in (-1, 0, 1):
            if left < preimage + shift < right:
                return True
    return False


def formula_spoiled(y: F, windows: tuple[tuple[int, F, F], ...]) -> bool:
    return any(in_doubled_open_interval(y, left, right) for _, left, right in windows)


def maximum_component_width(a: int, b: int) -> F:
    if a + b <= 7:
        return F(0)
    return min(F(2, 7 * b), F(a + b - 7, 7 * a * b))


def audit_explicit_windows() -> tuple[int, int, int]:
    profiles = 0
    positive_profiles = 0
    probes = 0
    for a in range(1, 60, 2):
        for b in range(a + 2, 60, 2):
            if gcd(a, b) != 1:
                continue
            profiles += 1
            windows = primitive_defect_windows(a, b)
            expected_positive_defects = sum(
                1 for defect in range(1, a + b, 2) if 7 * defect < a + b
            )
            require(
                len(windows) == 2 * expected_positive_defects,
                ("component count", a, b, len(windows), expected_positive_defects),
            )
            require(bool(windows) == (a + b > 7), ("existence boundary", a, b))
            if windows:
                positive_profiles += 1
                image_lengths = tuple(2 * (right - left) for _, left, right in windows)
                require(
                    max(image_lengths) == maximum_component_width(a, b),
                    ("maximum primitive width", a, b),
                )
            walls = wall_points((a, b))
            test_points = set(walls)
            test_points.update((left + right) / 2 for left, right in zip(walls, walls[1:]))
            for y in test_points:
                require(
                    fully_spoiled(y, (a, b)) == formula_spoiled(y, windows),
                    ("literal/formula membership", a, b, y),
                )
                probes += 1
    return profiles, positive_profiles, probes


def audit_pullback() -> tuple[int, int]:
    bases = 0
    scaled_profiles = 0
    dilations = (1, 3, 5, 7, 9, 11)
    for a in range(1, 32, 2):
        for b in range(a + 2, 32, 2):
            if gcd(a, b) != 1:
                continue
            bases += 1
            primitive = wall_components((a, b))
            for dilation in dilations:
                scaled = wall_components((dilation * a, dilation * b))
                require(
                    scaled == inverse_components(primitive, dilation),
                    ("odd pullback", a, b, dilation),
                )
                if primitive:
                    scaled_width = max(right - left for left, right in scaled)
                    require(
                        scaled_width == maximum_component_width(a, b) / dilation,
                        ("odd pullback width", a, b, dilation),
                    )
                scaled_profiles += 1
    return bases, scaled_profiles


def clearance(speeds: tuple[int, ...], phase: F) -> F:
    return min(distance(speed * phase) for speed in speeds)


def sieve_numerators(q: int) -> tuple[int, ...]:
    if 2 <= q <= 7:
        return (1,)
    require(q in SIEVE_K, ("sieve denominator", q))
    k = SIEVE_K[q]
    return (1, 1 + q, k, k + q)


def sieve_times(q: int) -> tuple[F, ...]:
    return tuple(F(numerator, 2 * q) for numerator in sieve_numerators(q))


def safe_sieve_times(speeds: tuple[int, ...], q: int) -> tuple[F, ...]:
    return tuple(phase for phase in sieve_times(q) if clearance(speeds, phase) >= RHO)


def audit_small_denominator_sieve() -> tuple[object, ...]:
    sieve_records = []
    triple_profiles = 0
    for q in range(2, 15):
        numerators = sieve_numerators(q)
        times = sieve_times(q)
        odd_classes = tuple(range(1, 2 * q, 2))
        maximal_pack_residues = tuple(
            residue for residue in range(2 * q) if residue % q
        )
        maximal_pack_speeds = tuple(2 * residue for residue in maximal_pack_residues)
        require(
            all(clearance(maximal_pack_speeds, phase) >= RHO for phase in times),
            ("maximal pack", q),
        )

        class_sets = tuple(
            tuple(
                speed
                for speed in odd_classes
                if distance(F(speed * numerator, 2 * q)) < RHO
            )
            for numerator in numerators
        )
        if q <= 7:
            require(class_sets == ((),), ("one-time danger classes", q, class_sets))
        else:
            require(
                class_sets == EXPECTED_DANGER_CLASSES[q],
                ("four-time danger classes", q, class_sets),
            )

        safe_count_histogram: dict[int, int] = {}
        pair_records = []
        for alpha, beta in combinations_with_replacement(odd_classes, 2):
            row = maximal_pack_speeds + (alpha, beta)
            winners = safe_sieve_times(row, q)
            require(winners, ("sieve pair", q, alpha, beta))
            if q <= 7:
                require(winners == times, ("one-time universal", q, alpha, beta))
            else:
                first_bank_spoiled = times[0] not in winners and times[1] not in winners
                if first_bank_spoiled and q % 2 == 0:
                    require(
                        times[2] in winners and times[3] in winners,
                        ("even-q bank switch", q, alpha, beta, winners),
                    )
                if first_bank_spoiled and q % 2 == 1:
                    require(
                        times[3] in winners,
                        ("odd-q bank switch", q, alpha, beta, winners),
                    )
            safe_count_histogram[len(winners)] = (
                safe_count_histogram.get(len(winners), 0) + 1
            )
            pair_records.append((alpha, beta, winners))
        require(
            len(pair_records) == q * (q + 1) // 2,
            ("pair count", q, len(pair_records)),
        )

        triple_count_histogram: dict[int, int] = {}
        triple_records = []
        if q <= 7:
            require(
                all(not danger_class for danger_class in class_sets),
                ("unbounded odd-exception capacity", q),
            )
            exception_capacity: object = "unbounded"
        elif q % 2 == 0:
            exception_capacity = 3
            for exceptions in combinations_with_replacement(odd_classes, 3):
                row = maximal_pack_speeds + exceptions
                winners = safe_sieve_times(row, q)
                require(winners, ("even-q triple capacity", q, exceptions))
                triple_count_histogram[len(winners)] = (
                    triple_count_histogram.get(len(winners), 0) + 1
                )
                triple_records.append((exceptions, winners))
            expected_triples = q * (q + 1) * (q + 2) // 6
            require(
                len(triple_records) == expected_triples,
                ("triple count", q, len(triple_records)),
            )
            triple_profiles += len(triple_records)
        else:
            exception_capacity = 2
        sieve_records.append(
            (
                q,
                SIEVE_K.get(q),
                numerators,
                class_sets,
                tuple(sorted(safe_count_histogram.items())),
                tuple(pair_records),
                exception_capacity,
                tuple(sorted(triple_count_histogram.items())),
                tuple(triple_records),
            )
        )
    require(triple_profiles == 1264, ("complete even-q triple total", triple_profiles))

    q14_record = sieve_records[-1]
    class_sets = q14_record[3]
    safe_count_histogram = dict(q14_record[4])
    pair_records = q14_record[5]

    new_allowed = tuple(residue for residue in range(56) if residue % 14)
    old_forbidden = {0, 11, 14, 22, 23, 28, 33, 34, 42, 45}
    old_allowed = tuple(residue for residue in range(56) if residue not in old_forbidden)
    require(set(old_allowed) < set(new_allowed), "strict THM-4049 containment")
    gained = tuple(sorted(set(new_allowed) - set(old_allowed)))
    require(gained == (11, 22, 23, 33, 34, 45), ("gained classes", gained))
    for residue in new_allowed:
        require(distance(F(residue, 14)) >= RHO, ("pack bank one", residue))
        require(distance(F(3 * residue, 14)) >= RHO, ("pack bank two", residue))

    clean_pack = (1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13)
    clean_row = tuple(2 * speed for speed in clean_pack) + (1, 15)
    require(len(clean_row) == len(set(clean_row)) == 13, "clean row typing")
    clean_clearances = tuple(
        (phase, clearance(clean_row, phase))
        for phase in (F(1, 28), F(15, 28), F(3, 28), F(17, 28))
    )
    require(clean_clearances[2][1] == clean_clearances[3][1] == RHO, "clean equality")
    old_bank = (F(1, 28), F(15, 28), F(5, 112), F(61, 112))
    old_clean_clearances = tuple(clearance(clean_row, phase) for phase in old_bank)
    require(
        old_clean_clearances == (F(1, 28), F(1, 28), F(1, 56), F(1, 56)),
        ("strict old-bank defeat", old_clean_clearances),
    )

    ray_witnesses = []
    for dilation in range(1, 100, 2):
        row = tuple(2 * speed for speed in clean_pack) + (dilation, 15 * dilation)
        require(len(row) == len(set(row)) == 13, ("ray typing", dilation))
        winners = safe_sieve_times(row, 14)
        require(winners, ("all-shape ray", dilation))
        ray_witnesses.append((dilation, winners))
    residual_wedge_dilations = tuple(
        dilation
        for dilation in range(1, 100, 2)
        if 15 * dilation < 12 * 13
        and 15 * dilation * dilation < 6 * 13 * 9 * dilation
    )
    require(
        residual_wedge_dilations == (1, 3, 5, 7, 9),
        ("strict THM-4052 ray gain", residual_wedge_dilations),
    )

    physical_pack = (2, 3, 4, 5, 6, 7, 8, 9, 11, 2**44, 3 * 2**44)
    physical_row = tuple(2 * speed for speed in physical_pack) + (1, 15)
    require(len(physical_row) == len(set(physical_row)) == 13, "physical row typing")
    require(all(speed % 14 for speed in physical_pack), "physical mod-14 sidecar")
    physical_clearances = tuple(
        (phase, clearance(physical_row, phase))
        for phase in (F(1, 28), F(15, 28), F(3, 28), F(17, 28))
    )
    require(
        physical_clearances
        == (
            (F(1, 28), F(1, 28)),
            (F(15, 28), F(1, 28)),
            (F(3, 28), F(1, 14)),
            (F(17, 28), F(1, 14)),
        ),
        ("physical hostile closure", physical_clearances),
    )
    old_physical_clearances = tuple(clearance(physical_row, phase) for phase in old_bank)
    require(
        old_physical_clearances
        == (F(1, 28), F(1, 28), F(1, 56), F(1, 56)),
        ("physical old-bank hostile", old_physical_clearances),
    )

    multiple_fourteen_row = (28, 1, 15)
    require(
        not safe_sieve_times(multiple_fourteen_row, 14),
        "multiple-fourteen owner must defeat this fixed bank",
    )
    even_exception_row = tuple(2 * speed for speed in range(1, 12)) + (28, 56)
    require(len(even_exception_row) == len(set(even_exception_row)) == 13, "even hostile typing")
    require(not safe_sieve_times(even_exception_row, 14), "oddness must be load-bearing")

    q15_maximal_residues = tuple(residue for residue in range(30) if residue % 15)
    for numerator in range(1, 15):
        require(
            any(distance(F(residue * numerator, 15)) < RHO for residue in q15_maximal_residues),
            ("q=15 maximal-grid hostile", numerator),
        )

    return (
        tuple(sieve_records),
        class_sets,
        tuple(sorted(safe_count_histogram.items())),
        len(old_allowed),
        len(new_allowed),
        gained,
        clean_clearances,
        old_clean_clearances,
        tuple(ray_witnesses),
        residual_wedge_dilations,
        physical_clearances,
        old_physical_clearances,
        tuple(pair_records),
        triple_profiles,
    )


def main() -> None:
    window_profiles, positive_profiles, window_probes = audit_explicit_windows()
    pullback_bases, scaled_profiles = audit_pullback()
    residue_summary = audit_small_denominator_sieve()
    semantic = sha256(
        repr(
            (
                window_profiles,
                positive_profiles,
                window_probes,
                pullback_bases,
                scaled_profiles,
                residue_summary,
            )
        ).encode("utf-8")
    ).hexdigest()

    print("LRC14_D2_SMALL_DENOMINATOR_SIEVE_AFFINE_RAY_FIREWALL_THM4070")
    print("scope=LRC14_OPEN;typed_d2_residual_sieve_only")
    print("exact_criterion=lonely_iff_G(H)_minus_[h]^-1Sigma_2(a,b)_nonempty")
    print("primitive_windows", window_profiles, "positive", positive_profiles, "wall_cell_probes", window_probes)
    print("odd_pullback_bases", pullback_bases, "scaled_profiles", scaled_profiles, "dilations=1,3,5,7,9,11")
    print("small_denominator_sieve_q=2..14")
    for record in residue_summary[0]:
        print(
            "q", record[0], "k", record[1], "numerators_over_2q", record[2],
            "danger_classes", record[3], "pair_histogram", record[4],
            "odd_exception_capacity", record[6], "triple_histogram", record[7],
        )
    print("even_q_unordered_triple_profiles", residue_summary[13])
    print("q14_bank_y=(1/14,3/14);bank_x=(1/28,15/28,3/28,17/28)")
    print("danger_classes_mod28", residue_summary[1])
    print("odd_pair_safe_time_histogram_mod28", residue_summary[2])
    print("THM4049_allowed", residue_summary[3], "new_allowed", residue_summary[4], "gained", residue_summary[5])
    print("clean_ray_base_clearances", residue_summary[6])
    print("clean_old_bank_clearances", residue_summary[7])
    print("clean_ray_odd_dilations_checked", len(residue_summary[8]))
    print("clean_ray_strict_THM4052_gain_dilations", residue_summary[9])
    print("physical_THM4049_hostile_clearances", residue_summary[10])
    print("physical_old_bank_clearances", residue_summary[11])
    print("q15_maximal_residue_grid_hostile=PASS")
    print("active_requirements", CHECKS)
    print("semantic_sha256", semantic)


if __name__ == "__main__":
    main()
