#!/usr/bin/env python3
"""No-import interval-union audit for the second-tier P40 transfer."""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, gcd
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

D = F(1, 14)
TIER = F(4, 77)
CAP = F(2, 63)
P33 = (
    1, 2, 4, 5, 8, 10, 16, 17, 19, 20, 23,
    25, 29, 31, 32, 34, 38, 40, 41, 43, 47,
    50, 51, 53, 58, 62, 64, 69, 71, 73, 75, 76, 80,
)
P40 = (
    1, 2, 4, 5, 8, 10, 16, 17, 19, 20, 23,
    25, 29, 31, 32, 34, 38, 40, 41, 43, 47,
    50, 51, 53, 58, 62, 64, 67, 69, 71, 73,
    75, 76, 80, 82, 86, 89, 93, 95, 141,
)
Y = F(57, 742)
X = F(799, 1484)
EXPECTED_SEMANTIC = "ce2f02310458b21dbc247590ff969065f64008ef4d93241dab052118b09553ee"
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
    for left, right in sorted((left, right) for left, right in intervals if left < right):
        if result and left <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], right))
        else:
            result.append((left, right))
    return tuple(result)


def intersect(
    first: tuple[tuple[F, F], ...], second: tuple[tuple[F, F], ...]
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


def cross_measure(p: int, q: int) -> F:
    verify(0 < p < q and p % 2 == q % 2 == 1 and gcd(p, q) == 1,
           ("primitive ratio", p, q))
    return 2 * measure(intersect(danger(p), danger(q, F(1, 2))))


def b2(value: F) -> F:
    value %= 1
    return value * value - value + F(1, 6)


def formula(p: int, q: int) -> F:
    minus = (F(1, 2) + F(q - p, 14)) % 1
    plus = (F(1, 2) + F(q + p, 14)) % 1
    return 2 * (F(1, 49) + (b2(minus) - b2(plus)) / (p * q))


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
    verify(b2(0) - b2(F(1, 2)) == F(1, 4), "B2 oscillation")
    verify(F(1, 90) < TIER - F(2, 49), "product threshold")

    ratios = 0
    maximum = F(0)
    maximizers: list[tuple[int, int]] = []
    above: list[tuple[int, int]] = []
    for q in range(3, 152, 2):
        for p in range(1, q, 2):
            if gcd(p, q) != 1:
                continue
            ratios += 1
            direct = cross_measure(p, q)
            verify(direct == formula(p, q), ("direct/formula", p, q))
            if direct > TIER:
                above.append((p, q))
                continue
            if direct > maximum:
                maximum, maximizers = direct, [(p, q)]
            elif direct == maximum:
                maximizers.append((p, q))
    verify(ratios == 2350, "ratio census")
    verify(above == [(1, 9)], "only second-tier exception")
    verify(maximum == TIER and maximizers == [(1, 11), (3, 11)],
           "sharp second-tier equality")

    safe = safe_set(P40)
    safe_measure = measure(safe)
    largest = max(right - left for left, right in safe)
    maximizer_intervals = tuple(interval for interval in safe if interval[1] - interval[0] == largest)
    verify(len(safe) == 44, "safe component count")
    verify(safe_measure == F(23518182747542658511201, 441420293060220631236800),
           "safe measure")
    verify(safe_measure > TIER and safe_measure < F(4, 63), "tier separation")
    verify(largest == F(137, 32900), "component width")
    verify(maximizer_intervals ==
           ((F(299, 658), F(321, 700)), (F(379, 700), F(359, 658))),
           "component locations")
    scales = tuple(t for t in range(1, 100, 2) if t * largest < CAP)
    verify(scales == (1, 3, 5, 7), "finite scale bank")
    verify(9 * largest > CAP, "scale-nine escape")

    verify(X == (Y + 1) / 2, "upper physical lift")
    body_gap = min(gap(2 * speed, X) for speed in P40)
    verify(body_gap == D, "body clock gap")
    tail_ledger = []
    for scale in scales:
        first, second = gap(scale, X), gap(9 * scale, X)
        verify(min(first, second) > D, ("tail clock", scale))
        tail_ledger.append((scale, pair(first), pair(second)))

    verify(set(P33) < set(P40), "P33 strict subset")
    total, old = comb(40, 11), comb(33, 11)
    width = 0
    for i, low in enumerate(P40):
        for j in range(i + 10, len(P40)):
            high = P40[j]
            if 27 * (13 * low - high) >= 4 * low * high:
                width += comb(j - i - 1, 9)
    verify(total == 2_311_801_440, "family total")
    verify(old == 193_536_720 and total - old == 2_118_264_720, "P33 increment")
    verify(width == 253_100 and total - width == 2_311_548_340, "width split")

    semantic = {
        "ratios": (ratios, pair(maximum), tuple(maximizers), tuple(above)),
        "pool": (len(safe), pair(safe_measure), pair(largest),
                 tuple((pair(left), pair(right)) for left, right in maximizer_intervals)),
        "clock": (pair(Y), pair(X), pair(body_gap), tuple(tail_ledger)),
        "families": (total, old, total - old, width, total - width),
    }
    semantic_hash = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC != "TO_BE_FILLED":
        verify(semantic_hash == EXPECTED_SEMANTIC, "semantic digest")

    print("LRC14_SECOND_TIER_POOL40_INDEPENDENT_INTERVAL_AUDIT_20260825")
    print("status=PASS;implementation=direct periodic-tooth intersections;LRC14=OPEN")
    print(f"ratios={ratios};above_4/77={tuple(above)};equal={tuple(maximizers)}")
    print(f"pool=(components={len(safe)},measure={safe_measure},max_width={largest})")
    print(f"max_components={maximizer_intervals}")
    print(f"finite_scales={scales};clock=({Y},{X},{body_gap},{tuple(tail_ledger)})")
    print(f"families={total};P33={old};increment={total-old}")
    print(f"width_covered={width};beyond_width={total-width}")
    print(f"checks={TESTS}")
    print(f"semantic_sha256={semantic_hash}")


if __name__ == "__main__":
    main()
