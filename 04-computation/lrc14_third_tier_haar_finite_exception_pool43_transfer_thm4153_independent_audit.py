#!/usr/bin/env python3
"""No-import interval audit for THM-4153.

This referee does not import the primary certificate.  It constructs danger
arcs, doubles their literal cross intersection to the quotient, intersects
safe teeth one speed at a time, and recomputes the clock and family ledgers.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, gcd
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

D = F(1, 14)
TIER = F(4, 91)
LENGTH = F(1, 399)
P40 = (
    1, 2, 4, 5, 8, 10, 16, 17, 19, 20, 23,
    25, 29, 31, 32, 34, 38, 40, 41, 43, 47,
    50, 51, 53, 58, 62, 64, 67, 69, 71, 73,
    75, 76, 80, 82, 86, 89, 93, 95, 141,
)
P43 = tuple(sorted(P40 + (111, 159, 285)))
Y = F(57, 742)
X = F(799, 1484)
EXPECTED_SEMANTIC = "fec0ab982776064029be96988978f783e365914e2615e03cc5bd0e1a17f3a998"

EXCEPTIONS = {
    (1, 9): (F(4, 63), F(2, 63), (1, 3, 5, 7, 9, 11)),
    (1, 11): (F(4, 77), F(2, 77), (1, 3, 5, 7, 9)),
    (3, 11): (F(4, 77), F(2, 77), (1, 3, 5, 7, 9)),
    (1, 23): (F(8, 161), F(2, 161), (1, 3)),
    (5, 11): (F(18, 385), F(9, 385), (1, 3, 5, 7, 9)),
    (1, 37): (F(12, 259), F(2, 259), (1, 3)),
    (1, 25): (F(8, 175), F(2, 175), (1, 3)),
    (3, 25): (F(8, 175), F(2, 175), (1, 3)),
    (3, 23): (F(22, 483), F(2, 161), (1, 3)),
    (1, 51): (F(16, 357), F(2, 357), (1,)),
    (5, 9): (F(2, 45), F(1, 45), (1, 3, 5, 7)),
}
ABOVE = tuple(sorted(EXCEPTIONS))
EQUAL = ((1, 13), (1, 39), (1, 65), (3, 13), (5, 13))

TESTS = 0


def verify(predicate: bool, label: object) -> None:
    global TESTS
    TESTS += 1
    if not predicate:
        raise ArithmeticError(f"independent audit failed: {label}")


def pair(value: F) -> tuple[int, int]:
    return value.numerator, value.denominator


def merge(intervals: list[tuple[F, F]]) -> tuple[tuple[F, F], ...]:
    result: list[tuple[F, F]] = []
    for left, right in sorted(
        (left, right) for left, right in intervals if left < right
    ):
        if result and left <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], right))
        else:
            result.append((left, right))
    return tuple(result)


def intersect(
    first: tuple[tuple[F, F], ...],
    second: tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    cells: list[tuple[F, F]] = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            cells.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        elif second[j][1] < first[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return merge(cells)


def measure(intervals: tuple[tuple[F, F], ...]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def circle_arc(center: F, radius: F) -> list[tuple[F, F]]:
    center %= 1
    left, right = center - radius, center + radius
    if left < 0:
        return [(F(0), right), (left + 1, F(1))]
    if right > 1:
        return [(F(0), right - 1), (left, F(1))]
    return [(left, right)]


def danger(speed: int, shift: F = F(0)) -> tuple[tuple[F, F], ...]:
    arcs: list[tuple[F, F]] = []
    for integer in range(speed):
        arcs.extend(circle_arc(F(integer, speed) - shift, D / speed))
    result = merge(arcs)
    verify(measure(result) == F(1, 7), ("danger measure", speed, shift))
    return result


def double_image(intervals: tuple[tuple[F, F], ...]) -> tuple[tuple[F, F], ...]:
    pieces: list[tuple[F, F]] = []
    for left, right in intervals:
        verify(right - left < F(1, 2), ("short double-image interval", left, right))
        doubled_left = (2 * left) % 1
        length = 2 * (right - left)
        doubled_right = doubled_left + length
        if doubled_right <= 1:
            pieces.append((doubled_left, doubled_right))
        else:
            pieces.append((doubled_left, F(1)))
            pieces.append((F(0), doubled_right - 1))
    return merge(pieces)


def cross_geometry(p: int, q: int) -> tuple[F, int, F]:
    verify(0 < p < q and p % 2 == q % 2 == 1 and gcd(p, q) == 1,
           ("primitive ratio", p, q))
    physical_cross = intersect(danger(p), danger(q, F(1, 2)))
    quotient_cross = double_image(physical_cross)
    return (
        measure(quotient_cross),
        len(quotient_cross),
        max((right - left for left, right in quotient_cross), default=F(0)),
    )


def b2(value: F) -> F:
    value %= 1
    return value * value - value + F(1, 6)


def formula(p: int, q: int) -> F:
    minus = (F(1, 2) + F(q - p, 14)) % 1
    plus = (F(1, 2) + F(q + p, 14)) % 1
    return 2 * (
        F(1, 49) + (b2(minus) - b2(plus)) / (p * q)
    )


def safe_teeth(speed: int) -> tuple[tuple[F, F], ...]:
    return tuple(
        ((F(integer) + D) / speed, (F(integer + 1) - D) / speed)
        for integer in range(speed)
    )


def safe_set(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    current = ((F(0), F(1)),)
    for speed in speeds:
        current = intersect(current, safe_teeth(speed))
    return current


def gap(speed: int, phase: F) -> F:
    residue = (speed * phase) % 1
    return min(residue, 1 - residue)


def main() -> None:
    verify(TIER - F(2, 49) == F(2, 637), "analytic threshold surplus")
    verify(F(1, 322) < F(2, 637), "large-product strict bound")

    small = tuple(
        (p, q)
        for q in range(3, 160, 2)
        for p in range(1, q, 2)
        if gcd(p, q) == 1 and p * q < 160
    )
    direct_values = {ratio: cross_geometry(*ratio)[0] for ratio in small}
    for ratio, direct in direct_values.items():
        verify(direct == formula(*ratio), ("direct/formula", ratio))
    above = tuple(sorted(ratio for ratio, value in direct_values.items()
                         if value > TIER))
    equal = tuple(sorted(ratio for ratio, value in direct_values.items()
                         if value == TIER))
    verify(len(small) == 117, "finite product universe")
    verify(above == ABOVE, "above-threshold exception list")
    verify(equal == EQUAL, "threshold equality list")

    exception_rows = []
    for ratio in ABOVE:
        direct_measure, components, width = cross_geometry(*ratio)
        expected_measure, expected_width, expected_scales = EXCEPTIONS[ratio]
        scales = tuple(t for t in range(1, 100, 2) if t * LENGTH < width)
        verify(direct_measure == expected_measure, ("exception measure", ratio))
        verify(width == expected_width, ("exception width", ratio))
        verify(scales == expected_scales, ("exception scale bank", ratio))
        exception_rows.append(
            (ratio, pair(direct_measure), components, pair(width), scales)
        )

    safe = safe_set(P43)
    safe_measure = measure(safe)
    largest = max(right - left for left, right in safe)
    maximizers = tuple(cell for cell in safe if cell[1] - cell[0] == largest)
    verify(len(P43) == len(set(P43)) == 43, "pool cardinality")
    verify(set(P40) < set(P43), "strict P40 inheritance")
    verify(tuple(value for value in P43 if value not in P40) == (111, 159, 285),
           "new labels")
    verify(len(safe) == 46, "safe component count")
    verify(safe_measure ==
           F(10080921463555906580413, 211196778145191767531400),
           "safe measure")
    verify(safe_measure > TIER, "strict measure threshold")
    verify(largest == LENGTH, "largest safe component")
    verify(maximizers ==
           ((F(911, 3990), F(307, 1330)),
            (F(1023, 1330), F(3079, 3990))),
           "largest reflection pair")

    verify(X == (Y + 1) / 2, "physical lift")
    body_gap = min(gap(value, Y) for value in P43)
    owners = tuple(value for value in P43 if gap(value, Y) == body_gap)
    verify(body_gap == D and owners == (53,), "body clock and owner")
    clock_rows = []
    tail_minimum = F(1, 2)
    for ratio in ABOVE:
        p, q = ratio
        rows = []
        for scale in EXCEPTIONS[ratio][2]:
            gaps = gap(p * scale, X), gap(q * scale, X)
            verify(min(gaps) > D, ("strict tail clock", ratio, scale))
            tail_minimum = min(tail_minimum, *gaps)
            rows.append((scale, pair(gaps[0]), pair(gaps[1])))
        clock_rows.append((ratio, tuple(rows)))
    verify(tail_minimum == F(113, 1484), "finite-bank minimum")
    verify(tail_minimum - D == F(1, 212), "finite-bank surplus")

    hostile_measure, _, hostile_width = cross_geometry(5, 13)
    hostile_scales = tuple(
        scale for scale in range(1, 100, 2)
        if scale * LENGTH < hostile_width
    )
    hostile_gaps = tuple(
        (scale, pair(gap(5 * scale, X)), pair(gap(13 * scale, X)))
        for scale in hostile_scales
    )
    verify(hostile_measure == TIER, "next-tier hostile measure")
    verify(hostile_width == F(2, 91), "next-tier hostile width")
    verify(hostile_scales == (1, 3, 5, 7), "next-tier hostile scales")
    verify(all(
        min(F(*first), F(*second)) < D
        for _, first, second in hostile_gaps
    ), "current clock fails every hostile scale")

    total, inherited = comb(43, 11), comb(40, 11)
    width = 0
    for i, low in enumerate(P43):
        for j in range(i + 10, len(P43)):
            high = P43[j]
            if 27 * (13 * low - high) >= 4 * low * high:
                width += comb(j - i - 1, 9)
    verify(total == 5_752_004_349, "family count")
    verify(inherited == 2_311_801_440, "P40 family count")
    verify(total - inherited == 3_440_202_909, "P40 increment")
    verify(width == 253_100, "width overlap")
    verify(total - width == 5_751_751_249, "beyond-width count")

    semantic = {
        "threshold": (pair(TIER), len(small), above, equal),
        "exceptions": tuple(exception_rows),
        "pool": (
            P43, len(safe), pair(safe_measure), pair(largest),
            tuple((pair(left), pair(right)) for left, right in maximizers),
        ),
        "clock": (
            pair(Y), pair(X), pair(body_gap), owners,
            pair(tail_minimum), tuple(clock_rows),
        ),
        "next_hostile": (
            (5, 13), pair(hostile_measure), pair(hostile_width),
            hostile_scales, hostile_gaps,
        ),
        "families": (
            total, inherited, total - inherited, width, total - width,
        ),
    }
    semantic_hash = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC != "TO_BE_FILLED":
        verify(semantic_hash == EXPECTED_SEMANTIC, "semantic digest")

    print("LRC14_THIRD_TIER_POOL43_INDEPENDENT_INTERVAL_AUDIT_20260825")
    print("status=PASS;implementation=direct interval intersections;LRC14=OPEN")
    print(f"threshold=4/91;finite_pairs={len(small)};above={above};equal={equal}")
    print(f"exception_geometry={tuple(exception_rows)}")
    print(f"pool=(size={len(P43)},components={len(safe)},measure={safe_measure},max_width={largest})")
    print(f"max_components={maximizers}")
    print(f"clock=({Y},{X},{body_gap},{owners},{tail_minimum});rows={tuple(clock_rows)}")
    print(f"next_hostile={semantic['next_hostile']}")
    print(f"families={total};P40={inherited};increment={total-inherited}")
    print(f"width_covered={width};beyond_width={total-width}")
    print(f"checks={TESTS}")
    print(f"semantic_sha256={semantic_hash}")


if __name__ == "__main__":
    main()
