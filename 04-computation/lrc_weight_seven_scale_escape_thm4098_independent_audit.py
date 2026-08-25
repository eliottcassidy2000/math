#!/usr/bin/env python3
"""Independent wall-cell audit for THM-4098.

Unlike the primary interval-merger, this script classifies the cells cut out
by every literal danger wall and evaluates safety on walls and midpoints.  It
also checks the odd/even dilation pullback pointwise on exact rational grids.
"""

from __future__ import annotations

from fractions import Fraction as Q
from itertools import combinations


DELTA = Q(1, 14)
ROWS = {
    7: (Q(4, 35), Q(13, 98), 4, 3, 55),
    6: (Q(4, 35), Q(1, 7), 5, 2, 35),
    5: (Q(4, 35), Q(1, 7), 6, 1, 35),
}


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def norm(x: Q) -> Q:
    r = x % 1
    return min(r, 1 - r)


def weight(v: int) -> int:
    return 2 if v & 1 else 1


def safe2(v: int, t: Q) -> bool:
    return norm(v * t) >= DELTA and norm(v * (t + Q(1, 2))) >= DELTA


def safe1(v: int, t: Q) -> bool:
    return norm(v * t) >= DELTA


def wall_points(speed: int, left: Q, right: Q, two_phase: bool) -> set[Q]:
    points: set[Q] = set()
    root_count = (weight(speed) if two_phase else 1) * speed
    radius = DELTA / speed
    low = root_count * (left - radius)
    high = root_count * (right + radius)
    lo = low.numerator // low.denominator - 2
    hi = high.numerator // high.denominator + 3
    for index in range(lo, hi + 1):
        centre = Q(index, root_count)
        for point in (centre - radius, centre + radius):
            if left <= point <= right:
                points.add(point)
    return points


def cell_survivor(
    bank: tuple[int, ...], left: Q, right: Q, *, two_phase: bool
) -> tuple[Q, Q] | None:
    """Return a safe point and the length of one safe wall cell."""
    walls = {left, right}
    for speed in bank:
        walls |= wall_points(speed, left, right, two_phase)
    ordered = sorted(walls)
    predicate = safe2 if two_phase else safe1
    for a, b in zip(ordered, ordered[1:]):
        mid = (a + b) / 2
        if all(predicate(speed, mid) for speed in bank):
            return mid, b - a
    for point in ordered:
        if all(predicate(speed, point) for speed in bank):
            # A safe wall might be isolated, so length zero is honest here.
            return point, Q(0)
    return None


def positive_safe_cell(bank: tuple[int, ...], *, two_phase: bool) -> Q:
    walls = {Q(0), Q(1)}
    for speed in bank:
        walls |= wall_points(speed, Q(0), Q(1), two_phase)
    ordered = sorted(walls)
    predicate = safe2 if two_phase else safe1
    best = Q(0)
    for a, b in zip(ordered, ordered[1:]):
        mid = (a + b) / 2
        if all(predicate(speed, mid) for speed in bank):
            best = max(best, b - a)
    return best


def main() -> None:
    pullback_gates = 0
    census_banks = 0
    cell_gates = 0
    minimum_cell: Q | None = None

    # An independent census over banks in [1,9], with three scales per row.
    for core, (left, right, count, odd_count, first_q) in ROWS.items():
        length = right - left
        check(first_q * length >= 1 > (first_q - 1) * length, "threshold arithmetic")
        for bank in combinations(range(1, 10), count):
            if sum(v & 1 for v in bank) != odd_count:
                continue
            check(sum(weight(v) for v in bank) == 7, "weight-seven typing")
            odd_cell = positive_safe_cell(bank, two_phase=True)
            even_cell = positive_safe_cell(bank, two_phase=False)
            check(odd_cell > 0, "odd base has no positive cell")
            check(even_cell > 0, "even base has no positive cell")

            for q in (first_q, first_q + 1, first_q + 2):
                base_cell = odd_cell if q & 1 else even_cell
                check((1 - base_cell) / q < length, "preimage mesh gap")
                scaled = tuple(q * v for v in bank)
                survivor = cell_survivor(scaled, left, right, two_phase=True)
                check(survivor is not None, "literal wall classifier found a cover")
                point, cell_length = survivor
                check(all(safe2(v, point) for v in tuple(range(1, core + 1)) + scaled), "body safety")
                if cell_length > 0:
                    minimum_cell = cell_length if minimum_cell is None else min(minimum_cell, cell_length)
                cell_gates += 1
            census_banks += 1

            # Exact pointwise pullback, including walls and generic grid points.
            for q in (3, 4, first_q, first_q + 1):
                for denominator in (17, 29):
                    for numerator in range(denominator):
                        theta = Q(numerator, denominator)
                        phi = (q * theta) % 1
                        for speed in bank:
                            target = safe2(q * speed, theta)
                            source = safe2(speed, phi) if q & 1 else safe1(speed, phi)
                            check(target == source, "dilation pullback mismatch")
                            pullback_gates += 1

    # The old AP7 bank is a true q=1 hostile, but both parity branches escape.
    hostile = (8, 9, 11, 13)
    left, right, *_ = ROWS[7]
    check(cell_survivor(hostile, left, right, two_phase=True) is None, "hostile no longer covers")
    for q in (55, 56):
        check(
            cell_survivor(tuple(q * v for v in hostile), left, right, two_phase=True) is not None,
            "hostile scale did not escape",
        )

    check(minimum_cell is not None and minimum_cell > 0, "no positive target cells recorded")
    print("THM-4098 independent exact wall-cell audit: PASS")
    print(f"declared bank census={census_banks}; target cell gates={cell_gates}")
    print(f"pointwise pullback gates={pullback_gates}")
    print(f"smallest positive target wall cell={minimum_cell}")
    print("q=1 AP7 hostile=cover; q=55 odd and q=56 even=survive")


if __name__ == "__main__":
    main()
