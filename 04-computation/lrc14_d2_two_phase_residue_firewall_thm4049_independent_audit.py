#!/usr/bin/env python3
"""Independent integer-only audit of THM-4049's four-time firewall.

This audit deliberately does not import the discovery companion.  All
clearances are tested after multiplying by the common denominator 112.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations_with_replacement


MODULUS = 112
THRESHOLD = 8
TIME_NUMERATORS = (4, 60, 5, 61)
ODD_RESIDUES = tuple(range(1, MODULUS, 2))
FORBIDDEN_PACK_RESIDUES = (0, 11, 14, 22, 23, 28, 33, 34, 42, 45)
ALLOWED_PACK_RESIDUES = tuple(
    r for r in range(56) if r not in FORBIDDEN_PACK_RESIDUES
)
H0 = tuple(range(1, 11)) + (12,)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def circular_residue(value):
    residue = value % MODULUS
    return min(residue, MODULUS - residue)


def speed_is_safe(speed, time_numerator):
    return circular_residue(speed * time_numerator) >= THRESHOLD


def safe_time_indices_for_pair(alpha, beta):
    return tuple(
        index
        for index, numerator in enumerate(TIME_NUMERATORS)
        if speed_is_safe(alpha, numerator) and speed_is_safe(beta, numerator)
    )


def all_four_times_clear_pack_residue(h):
    return all(
        speed_is_safe(2 * h, numerator) for numerator in TIME_NUMERATORS
    )


def main():
    transcript = []

    def emit(line):
        transcript.append(line)
        print(line)

    classified_allowed = tuple(
        h for h in range(56) if all_four_times_clear_pack_residue(h)
    )
    require(
        classified_allowed == ALLOWED_PACK_RESIDUES,
        "direct full-row pack classification disagrees with the theorem",
    )
    require(
        all(all_four_times_clear_pack_residue(h) for h in H0),
        "the named eleven-pack is not safe at all four times",
    )

    pairs = tuple(combinations_with_replacement(ODD_RESIDUES, 2))
    survivor_histogram = Counter()
    for alpha, beta in pairs:
        safe_indices = safe_time_indices_for_pair(alpha, beta)
        require(safe_indices, f"odd pair {(alpha, beta)} covers all four times")
        survivor_histogram[len(safe_indices)] += 1

    first_label_spoilers = tuple(
        w for w in ODD_RESIDUES if not speed_is_safe(w, TIME_NUMERATORS[0])
    )
    second_label_spoilers = tuple(
        w for w in ODD_RESIDUES if not speed_is_safe(w, TIME_NUMERATORS[1])
    )
    require(
        {w % 28 for w in first_label_spoilers} == {1, 27},
        "the first y=1/14 lift is not exactly the +/-1 mod 28 class",
    )
    require(
        {w % 28 for w in second_label_spoilers} == {13, 15},
        "the second y=1/14 lift is not exactly the +/-15 mod 28 class",
    )
    require(
        all(
            speed_is_safe(w, TIME_NUMERATORS[2])
            and speed_is_safe(w, TIME_NUMERATORS[3])
            for w in second_label_spoilers
        ),
        "a +/-15 mod 28 spoiler reaches the second bank",
    )

    require(
        circular_residue(2 * 12 * TIME_NUMERATORS[2]) == THRESHOLD,
        "the closed pack endpoint control changed",
    )
    require(
        not safe_time_indices_for_pair(22, 28),
        "the even-exception parity hostile no longer covers the bank",
    )
    require(
        not any(
            all(
                speed_is_safe(speed, numerator)
                for speed in (2 * 11, 1, 15)
            )
            for numerator in TIME_NUMERATORS
        ),
        "the h=11 pack-side hostile no longer covers the bank",
    )
    require(
        not any(
            all(
                speed_is_safe(speed, numerator)
                for speed in (2 * 14, 1, 11)
            )
            for numerator in TIME_NUMERATORS
        ),
        "the h=14 pack-side hostile no longer covers the bank",
    )

    emit("THM4049 INDEPENDENT INTEGER-RESIDUE AUDIT")
    emit("common_denominator=112;closed_threshold_numerator=8")
    emit(f"time_numerators={TIME_NUMERATORS}")
    emit(
        f"pack_residue_universe=56;allowed={len(ALLOWED_PACK_RESIDUES)};"
        f"forbidden={FORBIDDEN_PACK_RESIDUES}"
    )
    emit(
        f"odd_exception_residues={len(ODD_RESIDUES)};"
        f"unordered_pairs_with_repetition={len(pairs)}"
    )
    emit(
        "safe_time_count_histogram="
        + str(tuple(sorted(survivor_histogram.items())))
    )
    emit("first_bank_spoilers_mod28=((+/-1),(+/-15))")
    emit("plusminus15_mod28_safe_at_both_second_bank_times=True")
    emit("closed_endpoint=h12_at_x5_over112_has_clearance_1_over14=True")
    emit("pack_h11_hostile=True;pack_h14_hostile=True;even_pair_hostile=True")
    emit(
        "universal_consequence=every_allowed_pack_residue_family_and_every_"
        "odd_exception_pair_has_a_safe_displayed_time"
    )
    digest = sha256(("\n".join(transcript) + "\n").encode()).hexdigest()
    print(f"semantic_digest={digest}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
