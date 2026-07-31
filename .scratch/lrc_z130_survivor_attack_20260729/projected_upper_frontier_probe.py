#!/usr/bin/env python3
"""Cheap exact-first / analytic-second probe after the projected k=5 wall."""

from fractions import Fraction as F
from itertools import combinations
from math import lcm

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENGINE = HERE.parent / "lrc_j7_two_drift_transition_child_20260729"
sys.path.insert(0, str(ENGINE))
import exact_carrier as X  # noqa: E402


ETA = F(88, 1365)
ALPHA = F(2275, 18627)


def main() -> None:
    rows = []
    for body in combinations(range(1, 15), 6):
        carrier = X.carrier(body)
        h = X.mass(carrier)
        components = len(carrier)
        canonical_l = 14 * lcm(*body)
        lower = h * ETA
        wall_q = ALPHA * canonical_l
        wall = max(15, wall_q.numerator // wall_q.denominator + 1)
        for first in range(15, 131):
            if first % canonical_l == 0:
                continue
            delta_first = X.singleton_coverage(carrier, first) - h / 7
            second_min = max(first + 1, wall)
            upper = delta_first + F(6 * components, 49 * second_min)
            if upper >= lower:
                rows.append(
                    (
                        first,
                        upper - lower,
                        body,
                        second_min,
                        delta_first,
                        lower,
                        canonical_l,
                        components,
                        h,
                    )
                )

    maximum = max(row[0] for row in rows)
    frontiers = sorted(row for row in rows if row[0] == maximum)
    histogram = tuple(
        (first, sum(row[0] == first for row in rows))
        for first in sorted({row[0] for row in rows}, reverse=True)
    )
    print("k5 projected-wall analytic-second upper probe")
    print(f"rows={len(rows)};max_first={maximum};frontier_rows={len(frontiers)}")
    for row in frontiers:
        print(f"FRONTIER={row}")
    print(f"top20_histogram={histogram[:20]}")
    print("scope=exact delta(first), analytic delta(second) upper only")


if __name__ == "__main__":
    main()
