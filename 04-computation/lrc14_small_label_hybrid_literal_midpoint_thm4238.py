#!/usr/bin/env python3
"""Independent common-grid midpoint ledger for the hybrid exceptions."""

from fractions import Fraction
from math import gcd


BODY = (170, 176, 190, 193, 240, 252, 264, 286, 290)
RANGES = ((6, 590, 613), (25, 590, 597))


def lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def geometry(speeds: tuple[int, ...]) -> tuple[int, int, int]:
    common = 1
    for speed in speeds:
        common = lcm(common, 14 * speed)
    walls = {0, common}
    for speed in speeds:
        unit = common // (14 * speed)
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    ordered = sorted(walls)
    safe_cells = []
    mass = 0
    period = 2 * common
    for left, right in zip(ordered, ordered[1:]):
        midpoint_twice = left + right
        safe = all(
            common <= 7 * ((speed * midpoint_twice) % period) <= 13 * common
            for speed in speeds
        )
        safe_cells.append(safe)
        if safe:
            mass += right - left
    components = sum(
        safe_cells[index] and not safe_cells[index - 1]
        for index in range(len(safe_cells))
    )
    return common, mass, components


def main() -> None:
    print("LRC14_SMALL_LABEL_HYBRID_LITERAL_LEDGER")
    print("BODY " + ",".join(map(str, BODY)))
    best = None
    rows = 0
    for q, low, high in RANGES:
        for r in range(low, high + 1):
            denominator, mass, components = geometry((q,) + BODY + (r,))
            margin = 63 * mass - 4 * denominator
            assert margin > 0
            scaled = Fraction(margin, denominator)
            candidate = (scaled, q, r)
            if best is None or candidate < best:
                best = candidate
            rows += 1
            print(
                f"ROW Q {q} R {r} DENOMINATOR {denominator} "
                f"MASS_TICKS {mass} COMPONENTS {components} MARGIN63 {margin}"
            )
    assert best is not None and rows == 32
    print(
        f"SUMMARY ROWS {rows} NONPOSITIVE 0 MIN_SCALED_MARGIN "
        f"{best[0].numerator}/{best[0].denominator} AT_Q {best[1]} AT_R {best[2]}"
    )


if __name__ == "__main__":
    main()

