#!/usr/bin/env python3
"""Exact all-root audit of the LRC full-cell lower wall plus BV tail cap.

Scratch artifact only; it does not decide any residual finite tree.
"""

from __future__ import annotations

import math
from fractions import Fraction as F
from itertools import combinations

import exact_carrier as X


ETA = {
    2: F(1, 91),
    3: F(3, 91),
    4: F(51, 1183),
    5: F(88, 1365),
    6: F(22, 273),
}
EXPECTED_TOP = {
    2: (6_649, (1, 8, 10, 12, 13, 14), F(811, 2548), 38, 152_880, 1_821),
    3: (1_576, (2, 8, 9, 10, 12, 14), F(1993, 5880), 36, 35_280, 688),
    4: (860, (2, 8, 9, 10, 12, 14), F(1993, 5880), 36, 35_280, 1_009),
    5: (305, (2, 6, 8, 10, 12, 14), F(363, 980), 36, 11_760, 467),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def strict_floor(value: F) -> int:
    """Smallest integer strictly greater than ``value``."""
    return value.numerator // value.denominator + 1


def maximum_reciprocal_sum(first: int, count: int, largest_floor: int) -> F:
    """Maximum sum for distinct first=z1<...<zd and zd>=largest_floor."""
    require(count >= 1 and first >= 15, "invalid drift packet")
    largest_floor = max(largest_floor, first + count - 1)
    if count == 1:
        require(first >= largest_floor, "one-drift label misses lower wall")
        return F(1, first)
    if largest_floor == first + count - 1:
        labels = tuple(range(first, first + count))
    else:
        labels = tuple(range(first, first + count - 1)) + (largest_floor,)
    require(
        len(labels) == len(set(labels)) == count
        and labels[0] == first
        and labels[-1] >= largest_floor,
        "extremal reciprocal packet is mistyped",
    )
    return sum((F(1, label) for label in labels), F(0))


def root_cap(
    k: int,
    h: F,
    components: int,
    denominator: int,
) -> tuple[int | None, int, int]:
    """Return refined cap, full-cell largest-label floor, and crude cap."""
    drift_count = 7 - k
    gamma = F(6 * components, 49)
    target = ETA[k] * h
    full_cell_wall = F(k * denominator, 14 * (14 - k))
    largest_floor = max(15, strict_floor(full_cell_wall))
    crude_bound = F(drift_count) * gamma / target
    crude_cap = crude_bound.numerator // crude_bound.denominator

    if drift_count == 1:
        if largest_floor > crude_cap:
            return None, largest_floor, crude_cap
        return crude_cap, largest_floor, crude_cap

    def feasible(first: int) -> bool:
        return (
            gamma
            * maximum_reciprocal_sum(first, drift_count, largest_floor)
            >= target
        )

    low = 15
    high = crude_cap
    best: int | None = None
    while low <= high:
        middle = (low + high) // 2
        if feasible(middle):
            best = middle
            low = middle + 1
        else:
            high = middle - 1
    if best is not None:
        require(feasible(best), f"k={k}: accepted cap is infeasible")
        require(
            not feasible(best + 1),
            f"k={k}: cap is not the last feasible first label",
        )
    return best, largest_floor, crude_cap


def main() -> None:
    roots: list[tuple[tuple[int, ...], F, int, int]] = []
    for body in combinations(range(1, 15), 6):
        carrier = X.carrier(body)
        roots.append(
            (
                body,
                X.mass(carrier),
                len(carrier),
                14 * math.lcm(*body),
            )
        )
    require(len(roots) == 3_003, "six-body universe changed")

    print("LRC14 full-cell/BV-tail smallest-drift scratch audit")
    print(
        "full_cell_wall=M>kL/[14(14-k)];"
        "tail=delta(z)<=6r/(49z)"
    )
    for k in range(2, 6):
        rows: list[tuple[int, tuple[int, ...], F, int, int, int, int]] = []
        for body, h, components, denominator in roots:
            cap, largest_floor, crude_cap = root_cap(
                k,
                h,
                components,
                denominator,
            )
            require(cap is not None, f"k={k}: root unexpectedly scalar-closed")
            rows.append(
                (
                    cap,
                    body,
                    h,
                    components,
                    denominator,
                    largest_floor,
                    crude_cap,
                )
            )
        maximum_cap = max(row[0] for row in rows)
        top = tuple(row for row in rows if row[0] == maximum_cap)
        require(len(top) == 1, f"k={k}: maximum cap is not unique")
        require(top[0][:6] == EXPECTED_TOP[k], f"k={k}: top row changed")
        print(
            f"k={k};d={7-k};maximum_refined_cap={maximum_cap};"
            f"unique_top={top[0][1]};h={top[0][2]};r={top[0][3]};"
            f"L={top[0][4]};largest_drift_floor={top[0][5]};"
            f"root_crude_cap={top[0][6]}"
        )

    k = 6
    survivor_rows: list[tuple[int, tuple[int, ...], F, int, int, int, int]] = []
    for body, h, components, denominator in roots:
        cap, largest_floor, crude_cap = root_cap(
            k,
            h,
            components,
            denominator,
        )
        if cap is not None:
            survivor_rows.append(
                (
                    cap,
                    body,
                    h,
                    components,
                    denominator,
                    largest_floor,
                    crude_cap,
                )
            )
    require(len(survivor_rows) == 62, "k=6 scalar survivor count changed")
    maximum_cap = max(row[0] for row in survivor_rows)
    top = tuple(row for row in survivor_rows if row[0] == maximum_cap)
    require(
        maximum_cap == 91
        and tuple(row[1] for row in top)
        == (
            (1, 2, 3, 8, 10, 12),
            (1, 3, 4, 8, 10, 12),
        ),
        "k=6 top scalar rows changed",
    )
    print(
        "k=6;d=1;"
        f"scalar_closed_roots={3_003-len(survivor_rows)};"
        f"scalar_survivor_roots={len(survivor_rows)};"
        f"maximum_survivor_cap={maximum_cap};top_roots="
        f"{tuple(row[1] for row in top)}"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
