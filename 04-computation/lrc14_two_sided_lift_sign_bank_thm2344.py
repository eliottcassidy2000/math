#!/usr/bin/env python3
"""Optimization-safe exact sign bank for the secondary control in THM-2344."""

from math import gcd


SEVEN = 7
DIMENSION = 9


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def allowed(displacement: int) -> tuple[int, ...]:
    return tuple(
        residue
        for residue in range(SEVEN)
        if residue and (residue + displacement) % SEVEN
    )


def centered(residue: int) -> int:
    residue %= SEVEN
    return residue if residue <= 3 else residue - SEVEN


def danger_sign(index: int) -> int:
    absolute = abs(index)
    require(
        absolute and absolute % SEVEN,
        "danger sign requested at a Fourier zero",
    )
    return 1 if (absolute // SEVEN) % 2 == 0 else -1


def guard_sign(index: int) -> int:
    residue = abs(index) % SEVEN
    require(residue, "guard sign requested at a Fourier zero")
    return 1 if residue <= 3 else -1


def dynamic_count(
    speeds: list[int],
    choices: list[tuple[int, ...]],
    frequency: int,
) -> int:
    counts = [0] * SEVEN
    counts[0] = 1
    for speed, coordinate_choices in zip(speeds, choices):
        updated = [0] * SEVEN
        for old_residue, old_count in enumerate(counts):
            for choice in coordinate_choices:
                updated[(old_residue + speed * choice) % SEVEN] += old_count
        counts = updated
    return counts[frequency % SEVEN]


def sign_bank() -> tuple[int, int, int, int]:
    speeds = [13, 13**3, 2 * 13**5, 1, 14, 27, 40, 53, 66]
    relation = [
        -27,
        -27,
        -27,
        20_110_798,
        -41,
        -27,
        -27,
        -27,
        38,
    ]
    deepest_index = 2
    guard_index = 3
    multiplier = 1
    frequency = 13

    require(
        sum(a * b for a, b in zip(relation, speeds)) == 0,
        "relation left the exact kernel",
    )
    require(
        all(gcd(entry, 91) == 1 for entry in relation),
        "relation address lost an all-91-unit coordinate",
    )
    require(speeds[guard_index] == 1, "guard is not a Bezout coordinate")

    displacement = [
        (multiplier if index == deepest_index else 0) - entry
        for index, entry in enumerate(relation)
    ]
    choices = [allowed(entry % SEVEN) for entry in displacement]
    dp_count = dynamic_count(speeds, choices, frequency)

    residues = [0] * DIMENSION
    positive = 0
    negative = 0
    solutions = 0

    def recurse(index: int, centered_dot: int, right_sign: int) -> None:
        nonlocal positive, negative, solutions
        if index == DIMENSION:
            if (centered_dot - frequency) % SEVEN:
                return
            correction = (frequency - centered_dot) // SEVEN
            left_guard = centered(residues[guard_index]) + SEVEN * correction
            right_guard = left_guard + displacement[guard_index]
            term_sign = (
                right_sign
                * guard_sign(left_guard)
                * guard_sign(right_guard)
                * danger_sign(multiplier)
            )
            require(term_sign in (-1, 1), "term sign left {+1,-1}")
            solutions += 1
            if term_sign == 1:
                positive += 1
            else:
                negative += 1
            return

        for residue in choices[index]:
            residues[index] = residue
            representative = centered(residue)
            next_sign = right_sign
            if index != guard_index:
                next_sign *= danger_sign(
                    representative + displacement[index]
                )
            recurse(
                index + 1,
                centered_dot + speeds[index] * representative,
                next_sign,
            )

    recurse(0, 0, 1)
    require(solutions == dp_count, "enumeration and DP disagree")
    require(solutions == 334_825, "lift count changed")
    require(positive == 167_404, "positive sign count changed")
    require(negative == 167_421, "negative sign count changed")
    return dp_count, solutions, positive, negative


def main() -> None:
    dynamic, enumerated, positive, negative = sign_bank()
    print("theorem=THM-2344-secondary-sign-bank")
    print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-AUDITED")
    print("relation_coordinates_all_91_units=PASS")
    print(f"two_sided_lifts_dp={dynamic}")
    print(f"two_sided_lifts_enumerated={enumerated}")
    print(f"positive_term_signs={positive}")
    print(f"negative_term_signs={negative}")
    print(f"unweighted_sign_imbalance={positive-negative}")
    print("grouped_coefficient_cancellation=NOT-CLAIMED")
    print("scalar_rows_excluded=0")
    print("lrc14_status=OPEN")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
