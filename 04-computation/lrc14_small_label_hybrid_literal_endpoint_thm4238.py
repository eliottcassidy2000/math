#!/usr/bin/env python3
"""Exact endpoint-toggle ledger for the two hybrid small-q exceptions."""

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
    delta: dict[int, int] = {}
    for speed in speeds:
        unit = common // (14 * speed)
        for tooth in range(speed):
            leave = (14 * tooth + 1) * unit
            enter = (14 * tooth + 13) * unit
            delta[leave] = delta.get(leave, 0) - 1
            delta[enter] = delta.get(enter, 0) + 1
    failed = len(speeds)
    previous = 0
    previous_safe = False
    mass = 0
    components = 0
    for tick in sorted(delta):
        safe = failed == 0
        if safe:
            mass += tick - previous
            if not previous_safe:
                components += 1
        previous_safe = safe
        failed += delta[tick]
        previous = tick
    safe = failed == 0
    if safe:
        mass += common - previous
        if not previous_safe:
            components += 1
    assert failed == len(speeds)
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

