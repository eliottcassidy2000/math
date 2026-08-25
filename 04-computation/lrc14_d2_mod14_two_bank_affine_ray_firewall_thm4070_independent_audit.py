#!/usr/bin/env python3
"""No-primary-import integer residue audit for THM-4070.

The primary affine-window code is intentionally not imported.  This audit
loads every allowed pack residue at once for each q=2,...,14, and then
exhausts the redundant mod-56 odd-pair universe at q=14.
"""

from hashlib import sha256
from itertools import combinations_with_replacement
import sys


sys.stdout.reconfigure(newline="\n")

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


def circular_residue(numerator: int, denominator: int) -> int:
    residue = numerator % denominator
    return min(residue, denominator - residue)


def safe(speed: int, time_numerator: int, denominator: int) -> bool:
    """Integer form of ||speed*time_numerator/denominator|| >= 1/14."""
    return 14 * circular_residue(speed * time_numerator, denominator) >= denominator


def sieve_numerators(q: int) -> tuple[int, ...]:
    if 2 <= q <= 7:
        return (1,)
    require(q in SIEVE_K, ("sieve denominator", q))
    k = SIEVE_K[q]
    return (1, 1 + q, k, k + q)


def safe_times(speeds: tuple[int, ...], q: int) -> tuple[int, ...]:
    denominator = 2 * q
    return tuple(
        numerator
        for numerator in sieve_numerators(q)
        if all(safe(speed, numerator, denominator) for speed in speeds)
    )


def main() -> None:
    sieve_records = []
    total_pairs = 0
    total_triples = 0
    for q in range(2, 15):
        denominator = 2 * q
        numerators = sieve_numerators(q)
        allowed_pack_residues = tuple(
            residue for residue in range(denominator) if residue % q
        )
        maximal_pack_speeds = tuple(2 * residue for residue in allowed_pack_residues)
        for numerator in numerators:
            require(
                all(safe(speed, numerator, denominator) for speed in maximal_pack_speeds),
                ("maximal allowed pack", q, numerator),
            )

        odd_classes = tuple(range(1, denominator, 2))
        danger_classes = tuple(
            tuple(
                speed
                for speed in odd_classes
                if not safe(speed, numerator, denominator)
            )
            for numerator in numerators
        )
        if q <= 7:
            require(danger_classes == ((),), ("one-time danger", q, danger_classes))
        else:
            require(
                danger_classes == EXPECTED_DANGER_CLASSES[q],
                ("four-time danger", q, danger_classes),
            )

        histogram: dict[int, int] = {}
        pair_records = []
        for alpha, beta in combinations_with_replacement(odd_classes, 2):
            winners = safe_times(maximal_pack_speeds + (alpha, beta), q)
            require(winners, ("odd exception pair", q, alpha, beta))
            if q <= 7:
                require(winners == numerators, ("one-time universal", q, alpha, beta))
            else:
                first_bank_spoiled = numerators[0] not in winners and numerators[1] not in winners
                if first_bank_spoiled and q % 2 == 0:
                    require(
                        numerators[2] in winners and numerators[3] in winners,
                        ("even-q bank switch", q, alpha, beta, winners),
                    )
                if first_bank_spoiled and q % 2 == 1:
                    require(
                        numerators[3] in winners,
                        ("odd-q bank switch", q, alpha, beta, winners),
                    )
            histogram[len(winners)] = histogram.get(len(winners), 0) + 1
            pair_records.append((alpha, beta, winners))
        expected_pair_count = q * (q + 1) // 2
        require(len(pair_records) == expected_pair_count, ("pair count", q))
        total_pairs += len(pair_records)

        triple_histogram: dict[int, int] = {}
        triple_records = []
        if q <= 7:
            require(
                all(not danger_class for danger_class in danger_classes),
                ("unbounded odd-exception capacity", q),
            )
            exception_capacity: object = "unbounded"
        elif q % 2 == 0:
            exception_capacity = 3
            for exceptions in combinations_with_replacement(odd_classes, 3):
                winners = safe_times(maximal_pack_speeds + exceptions, q)
                require(winners, ("even-q triple capacity", q, exceptions))
                triple_histogram[len(winners)] = triple_histogram.get(len(winners), 0) + 1
                triple_records.append((exceptions, winners))
            expected_triples = q * (q + 1) * (q + 2) // 6
            require(len(triple_records) == expected_triples, ("triple count", q))
            total_triples += len(triple_records)
        else:
            exception_capacity = 2
        sieve_records.append(
            (
                q,
                SIEVE_K.get(q),
                numerators,
                danger_classes,
                tuple(sorted(histogram.items())),
                tuple(pair_records),
                exception_capacity,
                tuple(sorted(triple_histogram.items())),
                tuple(triple_records),
            )
        )
    require(total_pairs == 559, ("complete sieve pair total", total_pairs))
    require(total_triples == 1264, ("complete even-q triple total", total_triples))

    q = 14
    denominator = 28
    bank_numerators = sieve_numerators(q)
    allowed_pack_residues = tuple(residue for residue in range(56) if residue % 14)
    maximal_pack_speeds = tuple(2 * residue for residue in allowed_pack_residues)
    for numerator in bank_numerators:
        require(
            all(safe(speed, numerator, denominator) for speed in maximal_pack_speeds),
            ("maximal allowed pack", numerator),
        )

    odd_mod56 = tuple(range(1, 56, 2))
    mod56_histogram: dict[int, int] = {}
    mod56_pair_records = []
    for alpha, beta in combinations_with_replacement(odd_mod56, 2):
        winners = safe_times(maximal_pack_speeds + (alpha, beta), q)
        require(winners, ("odd exception pair", alpha, beta))
        mod56_histogram[len(winners)] = mod56_histogram.get(len(winners), 0) + 1
        mod56_pair_records.append((alpha, beta, winners))
    require(
        len(mod56_pair_records) == 406,
        ("mod-56 pair count", len(mod56_pair_records)),
    )

    for speed in odd_mod56:
        require(
            tuple(safe(speed, numerator, denominator) for numerator in bank_numerators)
            == tuple(
                safe(speed + 28, numerator, denominator)
                for numerator in bank_numerators
            ),
            ("mod-28 factor", speed),
        )

    danger_classes = tuple(
        tuple(
            speed
            for speed in range(1, 28, 2)
            if not safe(speed, numerator, denominator)
        )
        for numerator in bank_numerators
    )
    require(
        danger_classes == ((1, 27), (13, 15), (9, 19), (5, 23)),
        ("danger class table", danger_classes),
    )

    old_forbidden = {0, 11, 14, 22, 23, 28, 33, 34, 42, 45}
    old_allowed = tuple(residue for residue in range(56) if residue not in old_forbidden)
    require(set(old_allowed) < set(allowed_pack_residues), "old/new firewall containment")
    gained = tuple(sorted(set(allowed_pack_residues) - set(old_allowed)))
    require(gained == (11, 22, 23, 33, 34, 45), ("gained residues", gained))

    boundary_pack_residues = tuple(
        residue
        for residue in allowed_pack_residues
        if any(
            14 * circular_residue(2 * residue * numerator, denominator) == denominator
            for numerator in bank_numerators
        )
    )
    require(boundary_pack_residues, "closed pack boundary must occur")

    clean_pack = (1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13)
    clean_row = tuple(2 * speed for speed in clean_pack) + (1, 15)
    clean_winners = safe_times(clean_row, q)
    require(clean_winners == (3, 17), ("clean strict gain", clean_winners))
    old_bank_numerators_over_112 = (4, 60, 5, 61)
    old_clean_distance_numerators = tuple(
        min(
            circular_residue(speed * numerator, 112)
            for speed in clean_row
        )
        for numerator in old_bank_numerators_over_112
    )
    require(
        old_clean_distance_numerators == (4, 4, 2, 2),
        ("clean old-bank defeat", old_clean_distance_numerators),
    )

    physical_pack = (2, 3, 4, 5, 6, 7, 8, 9, 11, 2**44, 3 * 2**44)
    physical_row = tuple(2 * speed for speed in physical_pack) + (1, 15)
    physical_distances = tuple(
        min(
            circular_residue(speed * numerator, denominator)
            for speed in physical_row
        )
        for numerator in bank_numerators
    )
    require(physical_distances == (1, 1, 2, 2), ("physical distance numerators", physical_distances))

    require(
        not safe_times((28, 1, 15), q),
        "a multiple-fourteen pack owner must defeat this bank",
    )
    even_hostile = tuple(2 * speed for speed in range(1, 12)) + (28, 56)
    require(not safe_times(even_hostile, q), "even-exception hostile must defeat this bank")

    q15_maximal_residues = tuple(residue for residue in range(30) if residue % 15)
    for numerator in range(1, 15):
        require(
            any(
                14 * circular_residue(residue * numerator, 15) < 15
                for residue in q15_maximal_residues
            ),
            ("q=15 maximal-grid hostile", numerator),
        )

    semantic_payload = (
        tuple(sieve_records),
        allowed_pack_residues,
        tuple(sorted(mod56_histogram.items())),
        tuple(mod56_pair_records),
        danger_classes,
        old_allowed,
        gained,
        boundary_pack_residues,
        clean_winners,
        old_clean_distance_numerators,
        physical_distances,
    )
    semantic = sha256(repr(semantic_payload).encode("utf-8")).hexdigest()

    print("LRC14_D2_SMALL_DENOMINATOR_SIEVE_THM4070_INDEPENDENT")
    print("universe=maximal_nonmultiple_pack_residues_for_each_q=2..14")
    print("complete_sieve_pair_profiles", total_pairs)
    print("complete_even_q_triple_profiles", total_triples)
    for record in sieve_records:
        print(
            "q", record[0], "k", record[1], "numerators_over_2q", record[2],
            "danger_classes", record[3], "pair_histogram", record[4],
            "odd_exception_capacity", record[6], "triple_histogram", record[7],
        )
    print("q14_redundant_universe=all_52_nonmultiple_pack_residues_mod56_simultaneously")
    print("times_numerators_over_28", bank_numerators)
    print("danger_classes_mod28", danger_classes)
    print("odd_classes_mod56", len(odd_mod56), "unordered_pairs", len(mod56_pair_records))
    print("safe_time_count_histogram", tuple(sorted(mod56_histogram.items())))
    print("old_allowed", len(old_allowed), "new_allowed", len(allowed_pack_residues), "gained", gained)
    print("clean_winners_over_28", clean_winners)
    print("clean_old_bank_distance_numerators_over_112", old_clean_distance_numerators)
    print("physical_distance_numerators_over_28", physical_distances)
    print("multiple14_and_even_exception_hostiles=PASS")
    print("q15_maximal_residue_grid_hostile=PASS")
    print("active_requirements", CHECKS)
    print("semantic_sha256", semantic)


if __name__ == "__main__":
    main()
