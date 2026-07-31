#!/usr/bin/env python3
"""Exact finite first-drift bank after the aligned-pair excess floor.

Scratch artifact only.  It scans every canonical six-body root and every
nonaligned integral first drift allowed by the rootwise strict cap from
``two_drift_excess_cap_audit.py``.
"""

from __future__ import annotations

import math
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import exact_carrier as X
from full_cell_tail_cap_audit import root_cap


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def nonstrict_integer_cap(bound: F) -> int:
    return bound.numerator // bound.denominator


def main() -> None:
    candidate_rows = 0
    exceptional: list[
        tuple[tuple[int, ...], int, F, F, int, int, F, int]
    ] = []
    exceptional_root_counts: list[tuple[int, tuple[int, ...]]] = []
    closest_nonexception: tuple[F, tuple[object, ...]] | None = None
    maximum_second_bound: tuple[F, tuple[object, ...]] | None = None
    closed_subcritical_rows = 0
    finite_subcritical_rows = 0

    for body in combinations(range(1, 15), 6):
        carrier = X.integer_carrier(body)
        h = F(
            sum(right - left for left, right in carrier),
            X.RULER,
        )
        components = len(carrier)
        canonical_L = 14 * math.lcm(*body)
        delta_floor = F(88, 1365) * h
        gamma = F(6 * components, 49)
        first_cap, second_floor, _ = root_cap(
            5,
            h,
            components,
            canonical_L,
        )
        require(first_cap is not None, f"{body}: refined first cap disappeared")
        root_exception_count = 0

        for label in range(15, first_cap + 1):
            # d=1 is an aligned clock and belongs to an earlier branch.
            if label % canonical_L == 0:
                continue
            candidate_rows += 1
            delta = X.integer_singleton_coverage(carrier, label) - h / 7
            if delta >= delta_floor:
                exceptional.append(
                    (
                        body,
                        label,
                        delta,
                        delta_floor,
                        canonical_L,
                        components,
                        h,
                        first_cap,
                    )
                )
                root_exception_count += 1
                continue

            gap = delta_floor - delta
            second_bound = gamma / gap
            second_cap = nonstrict_integer_cap(second_bound)
            second_lower = max(label + 1, second_floor)
            if second_lower > second_cap:
                closed_subcritical_rows += 1
            else:
                finite_subcritical_rows += 1
            row = (
                body,
                label,
                delta,
                delta_floor,
                canonical_L,
                components,
                h,
                first_cap,
                gap,
                gamma,
                second_lower,
                second_cap,
            )
            if (
                closest_nonexception is None
                or gap < closest_nonexception[0]
                or (gap == closest_nonexception[0] and row < closest_nonexception[1])
            ):
                closest_nonexception = (gap, row)
            if (
                maximum_second_bound is None
                or second_bound > maximum_second_bound[0]
                or (
                    second_bound == maximum_second_bound[0]
                    and row < maximum_second_bound[1]
                )
            ):
                maximum_second_bound = (second_bound, row)

        if root_exception_count:
            exceptional_root_counts.append((root_exception_count, body))

    require(candidate_rows == 371_622, "first-label bank size changed")
    require(len(exceptional) == 4_084, "high-excess row count changed")
    require(
        len(exceptional_root_counts) == 2_309,
        "high-excess root count changed",
    )
    require(
        max(count for count, _ in exceptional_root_counts) == 5,
        "per-root high-excess maximum changed",
    )
    require(
        max(row[1] for row in exceptional) == 66,
        "largest high-excess first label changed",
    )
    require(closest_nonexception is not None, "nonexception bank disappeared")
    require(maximum_second_bound is not None, "second-label bound disappeared")
    require(
        closed_subcritical_rows == 352_334
        and finite_subcritical_rows == 15_204,
        "subcritical full-cell split changed",
    )
    require(
        closest_nonexception
        == (
            F(20_807, 3_902_338_440),
            (
                (2, 8, 9, 11, 13, 14),
                17,
                F(421_819, 20_011_992),
                F(54_997, 2_608_515),
                1_009_008,
                30,
                F(54_997, 168_168),
                174,
                F(20_807, 3_902_338_440),
                F(180, 49),
                40_041,
                688_956,
            ),
        ),
        "closest subcritical row changed",
    )
    require(
        maximum_second_bound[0] == F(14_335_120_800, 20_807),
        "uniform second-label bound changed",
    )

    label_histogram = tuple(sorted(Counter(row[1] for row in exceptional).items()))
    root_count_histogram = tuple(
        sorted(Counter(count for count, _ in exceptional_root_counts).items())
    )
    expected_label_histogram = (
        (15, 76),
        (16, 441),
        (17, 170),
        (18, 421),
        (19, 166),
        (20, 244),
        (21, 220),
        (22, 1114),
        (23, 65),
        (24, 203),
        (25, 15),
        (26, 349),
        (27, 163),
        (28, 24),
        (30, 20),
        (31, 2),
        (32, 15),
        (33, 188),
        (36, 53),
        (39, 7),
        (40, 25),
        (44, 96),
        (46, 2),
        (52, 4),
        (66, 1),
    )
    require(
        label_histogram == expected_label_histogram,
        "high-excess label histogram changed",
    )
    require(
        root_count_histogram == ((1, 1097), (2, 759), (3, 356), (4, 84), (5, 13)),
        "high-excess per-root histogram changed",
    )

    print("LRC14 five-aligned/two-drift first-label scratch bank")
    print("engine=self-contained integer-ruler endpoint/primitive sweep")
    print(
        f"candidate_rows={candidate_rows};"
        f"high_excess_rows={len(exceptional)};"
        f"high_excess_roots={len(exceptional_root_counts)};"
        f"closed_subcritical_rows={closed_subcritical_rows};"
        f"finite_subcritical_rows={finite_subcritical_rows}"
    )
    print(f"high_excess_label_histogram={label_histogram}")
    print(f"high_excess_per_root_histogram={root_count_histogram}")
    print(
        "closest_subcritical="
        f"gap={closest_nonexception[0]};row={closest_nonexception[1]}"
    )
    print(
        "uniform_second_bound_on_subcritical_bank="
        f"{maximum_second_bound[0]};"
        f"nonstrict_integer_cap={nonstrict_integer_cap(maximum_second_bound[0])};"
        f"row={maximum_second_bound[1]}"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
