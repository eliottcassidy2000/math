#!/usr/bin/env python3
"""Exact companion for THM-2436.

Normalize the guard on Z/91Z to G={0,...,25}.  A top blocker is one
vertical coset modulo seven, and each of the five ordinary top words is
an unoriented length-thirteen arithmetic progression with unit step.

For each guard-stabilizer orbit of blocker cosets this script:

* proves by an exact memoized set-cover search that five distinct
  unsigned ordinary steps cannot cover the punctured stalk;
* classifies every same-step pair which extends to a five-word cover;
* extracts the exact coarse and sharp parent-phase banks; and
* unions the banks over the moving blocker coset.

The searches are existence searches, so dynamic uncovered-point pivots
are safe: no candidate-index ordering is imposed.
"""

from __future__ import annotations

from functools import lru_cache
from fractions import Fraction
from math import gcd


N = 91
FULL = (1 << N) - 1
GUARD = sum(1 << x for x in range(26))
UNIT_STEPS = tuple(d for d in range(1, N // 2 + 1) if gcd(d, N) == 1)
LINE_REPRESENTATIVES = (0, 1, 2, 5)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def progression(center: int, step: int) -> int:
    mask = 0
    for offset in range(-6, 7):
        mask |= 1 << ((center + offset * step) % N)
    return mask


def vertical_line(residue: int) -> int:
    return sum(1 << x for x in range(residue, N, 7))


def points(mask: int) -> tuple[int, ...]:
    return tuple(x for x in range(N) if (mask >> x) & 1)


def reflect(mask: int) -> int:
    answer = 0
    for x in points(mask):
        answer |= 1 << ((25 - x) % N)
    return answer


ALL_APS = {
    (step, center): progression(center, step)
    for step in UNIT_STEPS
    for center in range(N)
}

require(len(UNIT_STEPS) == 36, "wrong unsigned unit-step count")
require(len(set(ALL_APS.values())) == 3276, "wrong unit-AP universe")
require(reflect(GUARD) == GUARD, "reflection does not stabilize guard")
for residue in range(7):
    require(
        reflect(vertical_line(residue)) == vertical_line((4 - residue) % 7),
        "wrong blocker-line reflection",
    )


POSITIVE_COVERS = {
    0: ((5, 56), (5, 57), (5, 58), (5, 59), (5, 60)),
    1: ((5, 56), (5, 57), (5, 58), (5, 59), (5, 60)),
    2: ((4, 51), (2, 64), (2, 38), (4, 53), (1, 84)),
    5: ((1, 84), (1, 71), (1, 33), (1, 47), (1, 58)),
}


EXPECTED_LINE_BANKS = {
    0: {
        1: (62, 74, 69),
        2: (36, 52, 45),
        3: (10, 14, 12),
        4: (12, 18, 15),
        5: (12, 20, 16),
        44: (12, 16, 14),
        45: (12, 20, 16),
    },
    1: {
        1: (64, 70, 67),
        2: (30, 48, 39),
        3: (4, 8, 6),
        4: (12, 18, 15),
        5: (14, 22, 18),
        44: (8, 12, 10),
        45: (16, 24, 20),
    },
    2: {
        1: (56, 64, 61),
        2: (30, 48, 39),
        3: (8, 20, 14),
        4: (6, 12, 9),
        5: (18, 26, 22),
        44: (10, 14, 12),
        45: (8, 16, 12),
    },
    5: {
        1: (68, 76, 75),
        2: (48, 58, 53),
        3: (10, 14, 12),
        4: (12, 18, 15),
        5: (12, 20, 16),
        44: (12, 16, 14),
        45: (16, 24, 20),
    },
}


EXPECTED_MOVING_BANKS = {
    1: (72, 76, 75),
    2: (52, 62, 57),
    3: (14, 26, 20),
    4: (12, 18, 15),
    5: (18, 26, 22),
    44: (12, 16, 14),
    45: (16, 24, 20),
}


def target_for_line(residue: int) -> int:
    return FULL & ~(GUARD | vertical_line(residue))


def no_distinct_step_cover(residue: int) -> tuple[bool, int]:
    """Return whether five distinct unsigned steps cover the target.

    The state records the uncovered target and used step classes.  The
    exact waste budget counts points lost to the guard/blocker or to
    overlaps.  Branching on any uncovered pivot is complete because
    every putative cover contains some still-available candidate
    through that pivot.
    """

    target = target_for_line(residue)
    target_size = target.bit_count()
    final_waste = 65 - target_size

    candidates: list[tuple[int, int, int, int]] = []
    for step_index, step in enumerate(UNIT_STEPS):
        seen: set[int] = set()
        for center in range(N):
            mask = ALL_APS[(step, center)] & target
            if mask in seen:
                continue
            seen.add(mask)
            if 13 - mask.bit_count() <= final_waste:
                candidates.append((mask, step_index, step, center))

    by_point: list[list[int]] = [[] for _ in range(N)]
    for index, (mask, _, _, _) in enumerate(candidates):
        remainder = mask
        while remainder:
            bit = remainder & -remainder
            by_point[bit.bit_length() - 1].append(index)
            remainder -= bit

    @lru_cache(maxsize=None)
    def exists(uncovered: int, used_steps: int) -> bool:
        if uncovered == 0:
            return True

        depth = used_steps.bit_count()
        remaining = 5 - depth
        uncovered_size = uncovered.bit_count()
        if remaining == 0 or uncovered_size > 13 * remaining:
            return False

        covered_size = target_size - uncovered_size
        waste = 13 * depth - covered_size
        slack = final_waste - waste
        if slack < 0:
            return False

        pivot_candidates: list[int] | None = None
        remainder = uncovered
        while remainder:
            bit = remainder & -remainder
            point = bit.bit_length() - 1
            remainder -= bit
            compatible: list[int] = []
            for index in by_point[point]:
                mask, step_index, _, _ = candidates[index]
                if (used_steps >> step_index) & 1:
                    continue
                gain = (mask & uncovered).bit_count()
                if 13 - gain <= slack:
                    compatible.append(index)
            if not compatible:
                return False
            if (
                pivot_candidates is None
                or len(compatible) < len(pivot_candidates)
            ):
                pivot_candidates = compatible

        require(pivot_candidates is not None, "missing distinct-step pivot")
        pivot_candidates.sort(
            key=lambda index: (
                -(candidates[index][0] & uncovered).bit_count(),
                candidates[index][2],
                candidates[index][3],
            )
        )
        for index in pivot_candidates:
            mask, step_index, _, _ = candidates[index]
            if exists(uncovered & ~mask, used_steps | (1 << step_index)):
                return True
        return False

    answer = exists(target, 0)
    return answer, exists.cache_info().misses


def extendable_same_step_banks(
    residue: int,
) -> tuple[dict[int, set[int]], int]:
    """Classify directed centre differences of extendable same-step pairs."""

    target = target_for_line(residue)

    # For completion, only the target restriction of an AP matters.
    restricted_representative: dict[int, tuple[int, int]] = {}
    for key, mask in ALL_APS.items():
        restricted_representative.setdefault(mask & target, key)
    restricted_masks = sorted(
        restricted_representative,
        key=lambda mask: (-mask.bit_count(), mask),
    )

    by_point: list[list[int]] = [[] for _ in range(N)]
    for mask in restricted_masks:
        remainder = mask
        while remainder:
            bit = remainder & -remainder
            by_point[bit.bit_length() - 1].append(mask)
            remainder -= bit

    @lru_cache(maxsize=None)
    def can_cover(uncovered: int, remaining: int) -> bool:
        if uncovered == 0:
            return True
        uncovered_size = uncovered.bit_count()
        if remaining == 0 or uncovered_size > 13 * remaining:
            return False

        minimum_gain = max(1, uncovered_size - 13 * (remaining - 1))
        pivot_candidates: list[int] | None = None
        rest = uncovered
        while rest:
            bit = rest & -rest
            point = bit.bit_length() - 1
            rest -= bit
            compatible = [
                mask
                for mask in by_point[point]
                if (mask & uncovered).bit_count() >= minimum_gain
            ]
            if not compatible:
                return False
            if (
                pivot_candidates is None
                or len(compatible) < len(pivot_candidates)
            ):
                pivot_candidates = compatible

        require(pivot_candidates is not None, "missing completion pivot")
        pivot_candidates.sort(
            key=lambda mask: (
                -(mask & uncovered).bit_count(),
                -mask.bit_count(),
                mask,
            )
        )
        for mask in pivot_candidates:
            if can_cover(uncovered & ~mask, remaining - 1):
                return True
        return False

    banks: dict[int, set[int]] = {}
    for step in UNIT_STEPS:
        differences: set[int] = set()
        # Include zero explicitly.  If two same-step blocks could coincide
        # and still extend to a cover, the physical phase argument would
        # have a stationary-pair boundary.  The exact search proves that
        # this never happens.
        for difference in range(0, N // 2 + 1):
            extendable = False
            for center in range(N):
                pair = (
                    ALL_APS[(step, center)]
                    | ALL_APS[(step, (center + difference) % N)]
                )
                if can_cover(target & ~pair, 3):
                    extendable = True
                    break
            if extendable:
                differences.add(difference)
                differences.add((-difference) % N)
        if differences:
            banks[step] = differences

    return banks, can_cover.cache_info().misses


def bank_sizes(step: int, differences: set[int]) -> tuple[int, int, int]:
    inverse = pow(step, -1, N)
    quotient = {(inverse * difference) % N for difference in differences}
    coarse = {
        (value + error) % N
        for value in quotient
        for error in (-1, 0, 1)
    }
    sharp = quotient | {(value - 1) % N for value in quotient}
    return len(differences), len(coarse), len(sharp)


line_banks: dict[int, dict[int, set[int]]] = {}
distinct_states: dict[int, int] = {}
completion_states: dict[int, int] = {}

for residue in LINE_REPRESENTATIVES:
    target = target_for_line(residue)
    expected_target_size = 56 if residue <= 4 else 55
    require(
        target.bit_count() == expected_target_size,
        "wrong punctured target size",
    )

    witness_union = GUARD | vertical_line(residue)
    for step, center in POSITIVE_COVERS[residue]:
        witness_union |= ALL_APS[(step, center)]
    require(witness_union == FULL, "positive punctured cover failed")

    distinct_exists, distinct_misses = no_distinct_step_cover(residue)
    require(not distinct_exists, "found a five-distinct-step cover")
    distinct_states[residue] = distinct_misses

    banks, completion_misses = extendable_same_step_banks(residue)
    observed = {
        step: bank_sizes(step, differences)
        for step, differences in banks.items()
    }
    require(
        observed == EXPECTED_LINE_BANKS[residue],
        "wrong fixed-line repeated-pair bank",
    )
    line_banks[residue] = banks
    completion_states[residue] = completion_misses


moving_banks: dict[int, set[int]] = {}
for banks in line_banks.values():
    for step, differences in banks.items():
        moving_banks.setdefault(step, set()).update(differences)

observed_moving = {
    step: bank_sizes(step, differences)
    for step, differences in moving_banks.items()
}
require(observed_moving == EXPECTED_MOVING_BANKS, "wrong moving-line bank")

d1_quotient = moving_banks[1]
d1_sharp = d1_quotient | {(value - 1) % N for value in d1_quotient}
require(
    d1_quotient == set(range(9, 83)) - {22, 69},
    "wrong exact d=1 quotient bank",
)
require(d1_sharp == set(range(8, 83)), "wrong exact d=1 sharp bank")

repeatable_steps = tuple(sorted(moving_banks))
require(
    repeatable_steps == (1, 2, 3, 4, 5, 44, 45),
    "wrong repeatable step classes",
)

expected_threshold_survivors = {
    52: (1, 2),
    65: (1,),
    78: (),
}
threshold_survivors = {
    threshold: tuple(
        step
        for step, (_, _, sharp_size) in sorted(observed_moving.items())
        if sharp_size >= threshold
    )
    for threshold in expected_threshold_survivors
}
require(
    threshold_survivors == expected_threshold_survivors,
    "wrong physical threshold survivors",
)

threshold_ratio_classes = {
    threshold: tuple(
        sorted(
            {
                min(pow(step, -1, N), N - pow(step, -1, N))
                for step in survivors
            }
        )
    )
    for threshold, survivors in threshold_survivors.items()
}
require(
    threshold_ratio_classes == {52: (1, 45), 65: (1,), 78: ()},
    "wrong guard-ratio survivors",
)
require(
    max(size[2] for size in observed_moving.values()) == 75,
    "wrong sharp universal pair cap",
)
require(
    Fraction(1, 2) + Fraction(1, 3) < Fraction(80, 91),
    "wrong distinct two-comb density threshold",
)
d1_tail = Fraction(8, 91) - Fraction(1, 14)
length_upper = Fraction(1, 7) / d1_tail
max_b_from_length = length_upper.numerator // length_upper.denominator
require(
    (d1_tail, max_b_from_length)
    == (Fraction(3, 182), 8),
    "wrong terminal tooth contradiction",
)
require(
    Fraction(13, 14 * max_b_from_length) > Fraction(8, 91),
    "first positive tooth reaches the wide-comb tail",
)


print("THM-2436 exact companion")
print(f"unit_ap13={len(ALL_APS)}")
print("blocker_line_orbits=0/4,1/3,2,5/6")
for residue in LINE_REPRESENTATIVES:
    encoded = ",".join(
        f"{step}:{sizes[0]}/{sizes[1]}/{sizes[2]}"
        for step, sizes in sorted(EXPECTED_LINE_BANKS[residue].items())
    )
    print(
        f"line_{residue}_target={target_for_line(residue).bit_count()} "
        f"distinct_cover=0 distinct_states={distinct_states[residue]} "
        f"completion_states={completion_states[residue]} banks={encoded}"
    )
print(
    "moving_banks="
    + ",".join(
        f"{step}:{sizes[0]}/{sizes[1]}/{sizes[2]}"
        for step, sizes in sorted(observed_moving.items())
    )
)
print("d1_quotient_bank=9..82_without_22,69")
print("d1_sharp_cells=8..82")
print("repeatable_steps=" + ",".join(map(str, repeatable_steps)))
print("parent_mass_52_over_91_surviving_steps=1,2")
print("parent_mass_65_over_91_surviving_steps=1")
print("parent_mass_78_over_91_surviving_steps=NONE")
print("parent_mass_52_over_91_signed_guard_ratio_classes=1,45")
print("parent_mass_65_over_91_signed_guard_ratio_classes=1")
print("sharp_single_pair_no_go=75/91")
print("ALL CHECKS PASSED")
