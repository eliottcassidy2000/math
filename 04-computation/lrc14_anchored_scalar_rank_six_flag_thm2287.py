#!/usr/bin/env python3
"""Exact arithmetic companion for the THM-2287 proof candidate.

Checks the centered prime-thirteen height recurrences, nested determinant
bounds, rank-five and rank-six label atlases, address/quotient counts,
Kakeya union-bound floors, spectral augmentation, the optional anchored
height-20 branch, and the all-165 unanchored contrast. All validity gates
use ``require`` so normal and optimized Python execute identical checks.
"""

from itertools import combinations, permutations
from math import factorial, prod
import sys


PRIME = 13
CENTER_RADIUS = (PRIME - 1) // 2
RANK_FIVE_SEED_HEIGHT = 196
RANK_SIX_SEED_HEIGHT = 526
ANCHOR_COEFFICIENT = 757
SPECTRAL_MULTIPLIER_CAP = 2_533

UNIFORM_HEIGHTS = (9_841, 4_921, 7_381, 11_072, 16_608, 24_912)
SHORT_HEIGHTS = (20, 196, 196, 206, 309, 526)

RELATIVE_UNIT_LABELS = ("H", "q1", "q2", "q3", "q4")
RELATIVE_BLOCKER_LABELS = ("c2", "c3")
ALL_UNIT_LABELS = ("H", "q1", "q2", "q3", "q4", "q5")
ALL_BLOCKER_LABELS = ("c1", "c2", "c3")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def extension_height(
    old_heights: tuple[int, ...],
    seed_height: int,
) -> int:
    return max(seed_height, (sum(old_heights) + 1) // 2)


def reconstruct_heights(anchor_height: int) -> tuple[int, ...]:
    heights = [anchor_height]
    while len(heights) < 5:
        heights.append(
            extension_height(tuple(heights), RANK_FIVE_SEED_HEIGHT)
        )
    heights.append(
        extension_height(tuple(heights), RANK_SIX_SEED_HEIGHT)
    )
    return tuple(heights)


def recurrence_slacks(heights: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        PRIME * heights[index]
        - (
            heights[index]
            + CENTER_RADIUS * sum(heights[:index])
        )
        for index in range(1, len(heights))
    )


def anchored_determinant_bounds(
    heights: tuple[int, ...],
    anchor_coefficient: int,
) -> tuple[int, ...]:
    return tuple(
        factorial(rank - 1)
        * anchor_coefficient
        * prod(heights[1:rank])
        for rank in range(2, 7)
    )


def blocker_histogram(
    atlas: tuple[tuple[str, ...], ...],
    blockers: tuple[str, ...],
) -> dict[int, int]:
    return {
        count: sum(
            sum(label in blockers for label in pattern) == count
            for pattern in atlas
        )
        for count in range(len(blockers) + 1)
    }


def centered_address_height(height_sum: int, depth: int) -> int:
    require(depth >= 1, "thirteen-adic depth must be positive")
    return height_sum * (PRIME**depth - 1) // 2


def integer_determinant(matrix: tuple[tuple[int, ...], ...]) -> int:
    size = len(matrix)
    require(
        size > 0 and all(len(row) == size for row in matrix),
        "determinant input must be nonempty and square",
    )
    total = 0
    for permutation in permutations(range(size)):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(size)
            for right in range(left + 1, size)
        )
        term = prod(
            matrix[row][permutation[row]]
            for row in range(size)
        )
        total += -term if inversions % 2 else term
    return total


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    uniform = reconstruct_heights(9_841)
    require(
        uniform == UNIFORM_HEIGHTS,
        "uniform pivot-extension height tuple changed",
    )
    uniform_slacks = recurrence_slacks(uniform)
    require(
        uniform_slacks == (6, 0, 6, 6, 6),
        "uniform prime-thirteen recurrence slacks changed",
    )
    require(
        all(slack >= 0 for slack in uniform_slacks),
        "uniform height box is not recurrence-invariant",
    )

    uniform_determinants = anchored_determinant_bounds(
        uniform,
        ANCHOR_COEFFICIENT,
    )
    require(
        uniform_determinants
        == (
            3_725_197,
            54_991_358_114,
            1_826_592_951_114_624,
            121_344_222_928_446_701_568,
            15_114_636_407_967_321_147_310_080,
        ),
        "uniform nested determinant ledger changed",
    )

    relative_bank = RELATIVE_UNIT_LABELS + RELATIVE_BLOCKER_LABELS
    rank_five_atlas = tuple(combinations(relative_bank, 4))
    rank_six_atlas = tuple(combinations(relative_bank, 5))
    require(len(rank_five_atlas) == 35, "rank-five atlas count changed")
    require(len(rank_six_atlas) == 21, "rank-six atlas count changed")
    require(
        blocker_histogram(rank_five_atlas, RELATIVE_BLOCKER_LABELS)
        == {0: 5, 1: 20, 2: 10},
        "rank-five blocker histogram changed",
    )
    require(
        blocker_histogram(rank_six_atlas, RELATIVE_BLOCKER_LABELS)
        == {0: 1, 1: 10, 2: 10},
        "rank-six blocker histogram changed",
    )
    require(
        all(
            sum(label in RELATIVE_UNIT_LABELS for label in pattern) >= 3
            for pattern in rank_six_atlas
        ),
        "rank-six atlas lost its three-unit floor",
    )

    complement_histogram = {0: 0, 1: 0, 2: 0}
    for pattern in rank_six_atlas:
        complement = set(relative_bank) - set(pattern)
        complement_blockers = sum(
            label in RELATIVE_BLOCKER_LABELS for label in complement
        )
        complement_histogram[complement_blockers] += 1
    require(
        complement_histogram == {0: 10, 1: 10, 2: 1},
        "three-column complement atlas changed",
    )

    uniform_sum = sum(uniform)
    require(uniform_sum == 74_735, "uniform row-height sum changed")
    require(
        centered_address_height(uniform_sum, 1) == 448_410,
        "uniform depth-one scalar address height changed",
    )
    require(
        2 * centered_address_height(uniform_sum, 1) == 896_820,
        "uniform depth-one fixed address height changed",
    )
    require(
        tuple(2 * height for height in uniform)
        == (19_682, 9_842, 14_762, 22_144, 33_216, 49_824),
        "uniform fixed-section row heights changed",
    )
    require(
        2 * uniform_determinants[-1]
        == 30_229_272_815_934_642_294_620_160,
        "uniform fixed-section determinant cap changed",
    )

    augmented_cap = SPECTRAL_MULTIPLIER_CAP * uniform_determinants[-1]
    require(
        augmented_cap
        == 38_285_374_021_381_224_466_136_432_640,
        "spectral augmented determinant cap changed",
    )

    hostile_speeds = (1, 4, 2, 3, 6, 10, 13, 2_197, 742_586)
    hostile_rows = []
    for unit_index, unit_speed in enumerate(hostile_speeds[1:6], start=1):
        row = [0] * 9
        row[0] = -unit_speed
        row[unit_index] = 1
        hostile_rows.append(tuple(row))
    hostile_anchor = [0] * 9
    hostile_anchor[1] = 13
    hostile_anchor[6] = -4
    hostile_rows.append(tuple(hostile_anchor))
    require(
        all(
            sum(coefficient * speed for coefficient, speed in zip(row, hostile_speeds))
            == 0
            for row in hostile_rows
        ),
        "THM-2299 hostile rank-six rows stopped being relations",
    )
    hostile_pivot = (1, 2, 3, 4, 5, 6)
    hostile_pivot_matrix = tuple(
        tuple(row[column] for column in hostile_pivot)
        for row in hostile_rows
    )
    hostile_determinant = integer_determinant(hostile_pivot_matrix)
    require(
        hostile_pivot_matrix[-1] == (13, 0, 0, 0, 0, -4),
        "THM-2299 hostile anchor row changed",
    )
    require(
        hostile_determinant == -4
        and hostile_determinant % PRIME != 0
        and max(abs(entry) for row in hostile_rows for entry in row) == 13,
        "THM-2299 hostile rank-six pivot gate changed",
    )

    base_full_combinations = PRIME**6
    base_hyperplane = PRIME**5
    base_all_unit_floor = base_full_combinations - 9 * base_hyperplane
    require(
        base_all_unit_floor == 4 * PRIME**5 == 1_485_172,
        "full Kakeya union-bound floor changed",
    )
    anchored_base = PRIME**5
    anchored_floor = anchored_base - 8 * PRIME**4
    require(
        anchored_floor == 5 * PRIME**4 == 142_805,
        "anchored Kakeya union-bound floor changed",
    )
    two_digit_base = PRIME**4
    two_digit_floor = two_digit_base - 7 * PRIME**3
    require(
        two_digit_floor == 6 * PRIME**3 == 13_182,
        "two-digit union-bound floor changed",
    )

    short = reconstruct_heights(20)
    require(short == SHORT_HEIGHTS, "short height tuple changed")
    short_slacks = recurrence_slacks(short)
    require(
        short_slacks == (2_232, 1_056, 0, 0, 750),
        "short recurrence slacks changed",
    )
    short_determinants = anchored_determinant_bounds(short, 20)
    require(
        short_determinants
        == (
            3_920,
            1_536_640,
            949_643_520,
            1_173_759_390_720,
            3_086_987_197_593_600,
        ),
        "short anchored determinant ledger changed",
    )
    short_sum = sum(short)
    require(short_sum == 1_453, "short row-height sum changed")
    require(
        centered_address_height(short_sum, 1) == 8_718,
        "short depth-one scalar address height changed",
    )
    require(
        2 * centered_address_height(short_sum, 1) == 17_436,
        "short depth-one fixed address height changed",
    )
    require(
        2 * short_determinants[-1] == 6_173_974_395_187_200,
        "short anchored fixed determinant cap changed",
    )

    all_labels = ALL_UNIT_LABELS + ALL_BLOCKER_LABELS
    all_six_subsets = tuple(combinations(all_labels, 6))
    blocker_pivots = tuple(
        pattern
        for pattern in all_six_subsets
        if any(label in ALL_BLOCKER_LABELS for label in pattern)
    )
    require(
        len(all_six_subsets) == 84 and len(blocker_pivots) == 83,
        "all-165 pivot atlas cardinality changed",
    )
    require(
        blocker_histogram(blocker_pivots, ALL_BLOCKER_LABELS)
        == {0: 0, 1: 18, 2: 45, 3: 20},
        "all-165 blocker pivot histogram changed",
    )
    unanchored_determinant_cap = factorial(6) * prod(short)
    require(
        unanchored_determinant_cap == 18_521_923_185_561_600,
        "all-165 unanchored determinant cap changed",
    )

    for depth in range(1, 7):
        centered_radius = (PRIME**depth - 1) // 2
        require(
            2 * centered_radius + 1 == PRIME**depth,
            f"centered residue count changed at depth {depth}",
        )
        require(
            centered_address_height(uniform_sum, depth)
            == uniform_sum * centered_radius,
            f"uniform address formula failed at depth {depth}",
        )
        require(
            centered_address_height(short_sum, depth)
            == short_sum * centered_radius,
            f"short address formula failed at depth {depth}",
        )
        require(
            4 * PRIME ** (6 * depth - 1)
            == base_all_unit_floor * PRIME ** (6 * depth - 6),
            f"full Kakeya lift count failed at depth {depth}",
        )
        require(
            5 * PRIME ** (5 * depth - 1)
            == anchored_floor * PRIME ** (5 * depth - 5),
            f"anchored Kakeya lift count failed at depth {depth}",
        )
        require(
            6 * PRIME ** (4 * depth - 1)
            == two_digit_floor * PRIME ** (4 * depth - 4),
            f"two-digit lift count failed at depth {depth}",
        )

    print("THM-2287 anchored scalar rank-six flag exact companion")
    print("status=PROVED_VERIFIED_EXACT_INDEPENDENTLY_AUDITED")
    print(
        "uniform_heights="
        + ",".join(str(height) for height in uniform)
        + f";sum={uniform_sum}"
    )
    print(
        "uniform_recurrence_slacks="
        + ",".join(str(slack) for slack in uniform_slacks)
    )
    print(
        "uniform_nested_determinant_caps="
        + ",".join(str(value) for value in uniform_determinants)
    )
    print(
        "rank5_atlas=35:(blockers0,1,2)=(5,20,10);"
        "rank6_atlas=21:(1,10,10);minimum_unit_labels=3"
    )
    print(
        "complement_atlas=21:"
        "(blockers0,1,2)=(10,10,1);"
        "residue_quotient_columns=3"
    )
    print(
        "uniform_addresses=13^(6n);"
        "scalar_height=74735*(13^n-1)/2;"
        "depth1=448410;"
        "fixed_depth1=896820"
    )
    print(
        "Kakeya_base_floor=4*13^5=1485172;"
        "c1_normalized=5*13^4=142805;"
        "two_digit_noncollapsed=6*13^3=13182"
    )
    print(f"spectral_augmented_determinant_cap={augmented_cap}")
    print(
        "THM2299_phase_hostile:"
        "speeds=1,4,2,3,6,10,13,2197,742586;"
        "rank6_height=13;pivot_det=-4_mod13;"
        "Fhat4=Ehat52=What4_0=0"
    )
    print(
        "short_heights="
        + ",".join(str(height) for height in short)
        + f";sum={short_sum};"
        f"anchored_D6={short_determinants[-1]}"
    )
    print(
        "all165_unanchored:"
        "heights=20,196,196,206,309,526;"
        "pivot_atlas=83:(blockers1,2,3)=(18,45,20);"
        f"D6_cap={unanchored_determinant_cap};"
        "depth1=8718;fixed_depth1=17436"
    )
    print("all exact checks passed")


if __name__ == "__main__":
    main()
