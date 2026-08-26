#!/usr/bin/env python3
"""Clean-room interval-union audit for THM-4150.

Unlike the primary wall-cell implementation, this audit constructs every
danger comb as a union of clipped circular teeth and constructs every body
safe set by repeated intersection of its closed safe teeth.
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
BOUND = F(4, 63)
WIDTH = F(2, 189)
POOL = (
    1, 2, 4, 5, 8, 10, 16, 17, 19, 20, 23,
    25, 29, 31, 32, 34, 38, 40, 41, 43, 47,
    50, 51, 53, 58, 62, 64, 69, 71, 73, 75, 76, 80,
)
BODY = (1, 17, 31, 32, 41, 47, 50, 51, 58, 62, 71)
EXPECTED_SEMANTIC = "1686df22845483de1378e8c91e9fb1260bfd826fb56e3652e607e8eb67ccf81e"
TESTS = 0


def check(predicate: bool, label: object) -> None:
    global TESTS
    TESTS += 1
    if not predicate:
        raise ArithmeticError(f"failed independent check: {label}")


def frac_pair(value: F) -> tuple[int, int]:
    return value.numerator, value.denominator


def merge(intervals: list[tuple[F, F]]) -> tuple[tuple[F, F], ...]:
    """Merge positive-length intervals; isolated points do not affect Haar measure."""
    ordered = sorted((left, right) for left, right in intervals if left < right)
    result: list[tuple[F, F]] = []
    for left, right in ordered:
        if result and left <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], right))
        else:
            result.append((left, right))
    return tuple(result)


def intersect(
    first: tuple[tuple[F, F], ...], second: tuple[tuple[F, F], ...]
) -> tuple[tuple[F, F], ...]:
    result: list[tuple[F, F]] = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            result.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        elif second[j][1] < first[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return merge(result)


def measure(intervals: tuple[tuple[F, F], ...]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def circular_arc(center: F, radius: F) -> list[tuple[F, F]]:
    center %= 1
    left, right = center - radius, center + radius
    if left < 0:
        return [(F(0), right), (left + 1, F(1))]
    if right > 1:
        return [(F(0), right - 1), (left, F(1))]
    return [(left, right)]


def danger_teeth(speed: int, shift: F = F(0)) -> tuple[tuple[F, F], ...]:
    arcs: list[tuple[F, F]] = []
    for integer in range(speed):
        arcs.extend(circular_arc(F(integer, speed) - shift, D / speed))
    result = merge(arcs)
    check(measure(result) == F(1, 7), ("danger measure", speed, shift))
    return result


def safe_teeth(speed: int) -> tuple[tuple[F, F], ...]:
    return tuple(
        ((F(integer) + D) / speed, (F(integer + 1) - D) / speed)
        for integer in range(speed)
    )


def full_safe_set(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    current = ((F(0), F(1)),)
    for speed in speeds:
        current = intersect(current, safe_teeth(speed))
    return current


def cross_measure(p: int, q: int) -> F:
    check(0 < p < q and p % 2 == q % 2 == 1 and gcd(p, q) == 1,
          ("primitive pair", p, q))
    # B=D_p intersect (D_q-1/2); the other cross term is B+1/2 and
    # doubling identifies the pair, so mu(C)=2mu(B).
    packet = intersect(danger_teeth(p), danger_teeth(q, F(1, 2)))
    return 2 * measure(packet)


def b2(value: F) -> F:
    value %= 1
    return value * value - value + F(1, 6)


def bernoulli_cross_measure(p: int, q: int) -> F:
    u = (F(1, 2) + F(q - p, 14)) % 1
    v = (F(1, 2) + F(q + p, 14)) % 1
    return 2 * (F(1, 49) + (b2(u) - b2(v)) / (p * q))


def main() -> None:
    check(b2(0) - b2(F(1, 2)) == F(1, 4), "Bernoulli oscillation")
    check(F(1, 2 * 23) < BOUND - F(2, 49), "product-23 strict gate")

    maximum = F(0)
    maximizers: list[tuple[int, int]] = []
    ratios = 0
    for q in range(3, 152, 2):
        for p in range(1, q, 2):
            if gcd(p, q) != 1:
                continue
            ratios += 1
            tooth_value = cross_measure(p, q)
            formula_value = bernoulli_cross_measure(p, q)
            check(tooth_value == formula_value, ("tooth/Bernoulli", p, q))
            check(tooth_value <= BOUND, ("measure bound", p, q))
            if tooth_value > maximum:
                maximum, maximizers = tooth_value, [(p, q)]
            elif tooth_value == maximum:
                maximizers.append((p, q))
    check(ratios == 2350, "ratio census")
    check(maximum == BOUND and maximizers == [(1, 9)], "sharp unique ratio")

    pool_set = full_safe_set(POOL)
    pool_measure = measure(pool_set)
    pool_width = max(right - left for left, right in pool_set)
    check(len(pool_set) == 46, "pool components")
    check(pool_measure == F(110551382435042260737, 1702610555154297252800),
          "pool measure")
    check(pool_measure > BOUND, "pool threshold")
    check(pool_width == F(3, 700) < WIDTH, "pool no-wide-component hostile")
    check(tuple((1 - right, 1 - left) for left, right in reversed(pool_set)) == pool_set,
          "pool reflection")

    body_set = full_safe_set(BODY)
    body_measure = measure(body_set)
    body_width = max(right - left for left, right in body_set)
    check(len(body_set) == 100, "body components")
    check(body_measure == F(1142176622583, 5854727790800), "body measure")
    check(body_width == F(21, 4100) < WIDTH, "body fragmented width")

    total = comb(33, 11)
    width_count = 0
    for i, smallest in enumerate(POOL):
        for j in range(i + 10, len(POOL)):
            largest = POOL[j]
            if smallest >= 3 and 27 * (13 * smallest - largest) >= 4 * smallest * largest:
                width_count += comb(j - i - 1, 9)
    check(total == 193_536_720, "family total")
    check(width_count == 208_000, "width subfamily")
    check(total - width_count == 193_328_720, "beyond-width subfamily")

    semantic = {
        "ratios": ratios,
        "maximum": frac_pair(maximum),
        "maximizers": maximizers,
        "pool": (len(pool_set), frac_pair(pool_measure), frac_pair(pool_width)),
        "body": (len(body_set), frac_pair(body_measure), frac_pair(body_width)),
        "families": (total, width_count, total - width_count),
    }
    semantic_hash = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC != "TO_BE_FILLED":
        check(semantic_hash == EXPECTED_SEMANTIC, "semantic digest")

    print("THM4150_INDEPENDENT_INTERVAL_UNION_AUDIT_20260825")
    print("status=PASS;implementation=direct danger-tooth and safe-tooth intersections")
    print(f"ratios={ratios};cross_bound={maximum};maximizers={tuple(maximizers)}")
    print(f"pool=(components={len(pool_set)},measure={pool_measure},max_width={pool_width})")
    print(f"fragmented_body=(components={len(body_set)},measure={body_measure},max_width={body_width})")
    print(f"families={total};width_covered={width_count};beyond_width={total-width_count}")
    print(f"checks={TESTS}")
    print(f"semantic_sha256={semantic_hash}")


if __name__ == "__main__":
    main()
