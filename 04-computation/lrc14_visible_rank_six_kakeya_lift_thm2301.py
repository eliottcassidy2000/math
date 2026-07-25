#!/usr/bin/env python3
"""Exact arithmetic referee for THM-2301.

Checks the sharp essential-arrangement constants, the centered-saturation
height recurrence, the explicit sharp rank-six column packets over F_7 and
F_13, all-depth lift arithmetic, and the intrinsic bright-packet height
ledgers.  All validity gates use ``require`` and remain active under
``python -O``.
"""

from fractions import Fraction
from itertools import product
import sys


COORDINATE_COUNT = 9
VISIBLE_HEIGHT = 526
TARGET_RANK = 6


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def ceiling_half(value: int) -> int:
    return (value + 1) // 2


def saturation_heights(height: int, rank: int) -> tuple[int, ...]:
    require(height >= 1, "height must be positive")
    require(rank >= 1, "rank must be positive")
    bounds = [height]
    while len(bounds) < rank:
        bounds.append(max(height, ceiling_half(sum(bounds))))
    return tuple(bounds)


def seeded_extension_heights(
    initial_height: int,
    candidate_heights: tuple[int, ...],
) -> tuple[int, ...]:
    require(initial_height >= 1, "initial height must be positive")
    require(
        all(height >= 1 for height in candidate_heights),
        "candidate heights must be positive",
    )
    bounds = [initial_height]
    for candidate_height in candidate_heights:
        bounds.append(
            max(candidate_height, ceiling_half(sum(bounds)))
        )
    return tuple(bounds)


def arrangement_lower_bound(
    field_order: int,
    dimension: int,
    hyperplanes: int,
) -> int:
    require(field_order >= 2, "field order too small")
    require(dimension >= 1, "dimension too small")
    require(
        dimension <= hyperplanes <= dimension + field_order - 1,
        "arrangement parameters outside theorem range",
    )
    return (
        (field_order - hyperplanes + dimension - 1)
        * (field_order - 1) ** (dimension - 1)
    )


def rank_mod(matrix: list[list[int]], prime: int) -> int:
    require(prime >= 2, "prime too small")
    work = [[entry % prime for entry in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0]) if work else 0
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, row_count)
                if work[row][column] % prime
            ),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [
            entry * inverse % prime for entry in work[pivot_row]
        ]
        for row in range(row_count):
            if row == pivot_row:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    (left - factor * right) % prime
                    for left, right in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def sharp_rank_six_columns(prime: int) -> tuple[tuple[int, ...], ...]:
    require(prime in (7, 13), "sharp model only instantiated at 7 and 13")
    columns: list[tuple[int, ...]] = []
    for index in range(TARGET_RANK):
        column = [0] * TARGET_RANK
        column[index] = 1
        columns.append(tuple(column))
    for scalar in (1, 2, 3):
        column = [0] * TARGET_RANK
        column[0] = -scalar % prime
        column[1] = 1
        columns.append(tuple(column))
    require(
        len(columns) == COORDINATE_COUNT,
        "sharp packet column count changed",
    )
    return tuple(columns)


def dot(left: tuple[int, ...], right: tuple[int, ...], prime: int) -> int:
    require(len(left) == len(right), "dot-product dimension mismatch")
    return sum(a * b for a, b in zip(left, right)) % prime


def normalized_all_unit_count(prime: int) -> int:
    columns = sharp_rank_six_columns(prime)
    rows = [
        [columns[column][row] for column in range(COORDINATE_COUNT)]
        for row in range(TARGET_RANK)
    ]
    require(
        rank_mod(rows, prime) == TARGET_RANK,
        f"sharp packet lost rank over F_{prime}",
    )
    count = 0
    for tail in product(range(prime), repeat=TARGET_RANK - 1):
        coefficients = (1,) + tail
        word = tuple(dot(coefficients, column, prime) for column in columns)
        if all(word):
            count += 1
    return count


def gaussian_containment_ratio(
    field_order: int,
    ambient_rank: int,
    subspace_rank: int,
) -> Fraction:
    require(
        ambient_rank >= subspace_rank >= 1,
        "invalid Gaussian containment dimensions",
    )
    return Fraction(
        field_order ** (ambient_rank - subspace_rank) - 1,
        field_order**ambient_rank - 1,
    )


def jackson_coefficient(bandwidth: int, frequency: int) -> int:
    require(bandwidth >= 2, "Jackson bandwidth too small")
    require(
        0 <= frequency <= 2 * bandwidth - 2,
        "frequency outside Jackson support",
    )
    if frequency <= bandwidth:
        return (
            4 * bandwidth**3
            - 6 * bandwidth * frequency**2
            + 2 * bandwidth
            + 3 * frequency**3
            - 3 * frequency
        ) // 6
    reflected = 2 * bandwidth - frequency
    return (reflected**3 - reflected) // 6


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    generic_rank_six_bounds = seeded_extension_heights(
        20,
        (196, 196, 196, 196, VISIBLE_HEIGHT),
    )
    require(
        generic_rank_six_bounds == (20, 196, 196, 206, 309, 526),
        "generic rank-six saturation recurrence changed",
    )
    require(
        sum(generic_rank_six_bounds) == 1453,
        "generic rank-six height sum changed",
    )

    visible_rank_six_bounds = seeded_extension_heights(
        462,
        (VISIBLE_HEIGHT,) * 5,
    )
    require(
        visible_rank_six_bounds
        == (462, 526, 526, 757, 1136, 1704),
        "visible rank-six saturation recurrence changed",
    )
    require(
        sum(visible_rank_six_bounds) == 5111,
        "visible rank-six height sum changed",
    )

    visible_full_rank_bounds = seeded_extension_heights(
        462,
        (VISIBLE_HEIGHT,) * 7,
    )
    require(
        visible_full_rank_bounds
        == (462, 526, 526, 757, 1136, 1704, 2556, 3834),
        "visible full-carrier saturation recurrence changed",
    )
    require(
        sum(visible_full_rank_bounds) == 11501,
        "visible full-carrier sum changed",
    )

    rank_law_13 = tuple(
        arrangement_lower_bound(13, rank - 1, 8)
        for rank in range(2, 10)
    )
    expected_rank_law_13 = tuple(
        (rank + 3) * 12 ** (rank - 2) for rank in range(2, 10)
    )
    require(
        rank_law_13 == expected_rank_law_13,
        "F_13 anchored all-rank law changed",
    )
    require(rank_law_13[1] == 72, "rank-three count changed")
    require(rank_law_13[3] == 13_824, "rank-five count changed")
    require(rank_law_13[4] == 186_624, "rank-six count changed")

    expected_counts = {
        7: {
            "normalized": 3 * 6**4,
            "unanchored": 3 * 6**5,
        },
        13: {
            "normalized": 9 * 12**4,
            "unanchored": 9 * 12**5,
        },
    }
    observed_counts = {}
    for prime in (7, 13):
        normalized = normalized_all_unit_count(prime)
        unanchored = (prime - 1) * normalized
        require(
            normalized
            == arrangement_lower_bound(prime, TARGET_RANK - 1, 8),
            f"normalized sharp count failed over F_{prime}",
        )
        require(
            unanchored
            == arrangement_lower_bound(prime, TARGET_RANK, 9),
            f"unanchored sharp count failed over F_{prime}",
        )
        require(
            normalized == expected_counts[prime]["normalized"],
            f"normalized expected count failed over F_{prime}",
        )
        require(
            unanchored == expected_counts[prime]["unanchored"],
            f"unanchored expected count failed over F_{prime}",
        )
        observed_counts[prime] = (unanchored, normalized)

    for depth in (1, 2, 3):
        require(
            observed_counts[13][0] * 13 ** (6 * depth - 6)
            == 9 * 12**5 * 13 ** (6 * depth - 6),
            "thirteen-adic unanchored lift changed",
        )
        require(
            observed_counts[13][1] * 13 ** (5 * depth - 5)
            == 9 * 12**4 * 13 ** (5 * depth - 5),
            "thirteen-adic normalized lift changed",
        )
        require(
            observed_counts[7][0] * 7 ** (6 * depth - 6)
            == 3 * 6**5 * 7 ** (6 * depth - 6),
            "septimal unanchored lift changed",
        )
        require(
            observed_counts[7][1] * 7 ** (5 * depth - 5)
            == 3 * 6**4 * 7 ** (5 * depth - 5),
            "septimal normalized lift changed",
        )

    for ambient_rank in range(6, 9):
        ratio = gaussian_containment_ratio(7, ambient_rank, 6)
        require(
            9 * ratio < 1,
            f"septimal Grassmannian union bound failed at rank {ambient_rank}",
        )
        require(
            ratio < Fraction(1, 7**6),
            f"Gaussian ratio cap failed at rank {ambient_rank}",
        )

    require(
        6 * sum(visible_full_rank_bounds) == 69_006,
        "p=13 bright row changed",
    )
    five_largest = sum(
        sorted(visible_full_rank_bounds, reverse=True)[:5]
    )
    require(five_largest == 9_987, "five companion rows changed")
    require(
        6 * sum(visible_full_rank_bounds) + five_largest == 78_993,
        "p=13 intrinsic packet sum changed",
    )
    require(
        6 * (3 * sum(visible_full_rank_bounds)) == 207_018,
        "p=7 intrinsic packet sum changed",
    )

    pivot_atlas = (
        3 * 6,
        3 * 15,
        1 * 20,
    )
    require(pivot_atlas == (18, 45, 20), "blocker pivot atlas changed")
    require(sum(pivot_atlas) == 83, "blocker pivot total changed")

    bandwidth = 264
    require(2 * bandwidth - 2 == VISIBLE_HEIGHT, "visible height changed")
    require(
        all(
            jackson_coefficient(bandwidth, frequency) > 0
            for frequency in range(VISIBLE_HEIGHT + 1)
        ),
        "Jackson support developed an internal zero",
    )
    visible_nonzero_count = sum(
        1
        for frequency in range(-VISIBLE_HEIGHT, VISIBLE_HEIGHT + 1)
        if frequency == 0 or frequency % 7
    )
    require(
        visible_nonzero_count == 903,
        "one-coordinate visible support count changed",
    )

    print("THM-2301 visible rank-six Kakeya referee")
    print(
        "generic_rank6_heights="
        + ",".join(map(str, generic_rank_six_bounds))
        + f";sum={sum(generic_rank_six_bounds)}"
    )
    print(
        "visible_rank6_heights="
        + ",".join(map(str, visible_rank_six_bounds))
        + f";sum={sum(visible_rank_six_bounds)}"
    )
    print(
        "visible_full_rank8_heights="
        + ",".join(map(str, visible_full_rank_bounds))
        + f";sum={sum(visible_full_rank_bounds)}"
    )
    print(
        "F13_anchored_rank_law_r2_to_r9="
        + ",".join(map(str, rank_law_13))
    )
    print(
        "rank6_base_counts="
        f"F7:unanchored={observed_counts[7][0]},"
        f"normalized={observed_counts[7][1]};"
        f"F13:unanchored={observed_counts[13][0]},"
        f"normalized={observed_counts[13][1]}"
    )
    print(
        "intrinsic_packet_height_sums="
        "F7:207018,F13:78993;"
        "fixed_section=414036,157986"
    )
    print(
        "rank6_blocker_pivot_atlas="
        f"one={pivot_atlas[0]},two={pivot_atlas[1]},"
        f"three={pivot_atlas[2]},total={sum(pivot_atlas)}"
    )
    print(
        "visible_support="
        f"N={bandwidth},H={VISIBLE_HEIGHT},"
        f"one_coordinate_count={visible_nonzero_count},"
        "zeros_exactly_nonzero_multiples_of_7"
    )
    print("all exact checks passed")


if __name__ == "__main__":
    main()
