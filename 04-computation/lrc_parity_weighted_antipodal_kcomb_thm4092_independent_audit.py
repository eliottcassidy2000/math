#!/usr/bin/env python3
"""Independent original-theta/cell-sweep audit for THM-4092.

This implementation never merges tooth intervals.  It builds the exact wall
arrangement, classifies open cells and wall points directly with the original
two-phase inequalities, and independently replays finite threshold banks.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations


Q = Fraction
RADIUS = Q(1, 14)
ROWS = {
    8: (Q(11, 49), Q(13, 56), 98),
    7: (Q(4, 35), Q(13, 98), 98),
    6: (Q(4, 35), Q(1, 7), 70),
    5: (Q(4, 35), Q(1, 7), 70),
}
FLOORS = {
    (8, 0): 85, (8, 1): 106, (8, 2): 150, (8, 3): 281,
    (7, 0): 63, (7, 1): 90, (7, 2): 172,
    (6, 0): 76, (6, 1): 146,
    (5, 0): 181,
}
NEGATIVE_ROWS = {
    (8, 0): (10, 22, 26), (8, 1): (9, 10, 26),
    (8, 2): (9, 10, 11), (8, 3): (9, 11, 13),
    (7, 0): (8, 10, 12, 26), (7, 1): (8, 9, 10, 12),
    (7, 2): (8, 9, 10, 11),
    (6, 0): (8, 10, 14, 22, 26), (6, 1): (7, 8, 10, 12, 26),
    (5, 0): (6, 8, 10, 14, 22, 26),
}


def check(flag: bool, text: str) -> None:
    if not flag:
        raise ArithmeticError(text)


def norm(x: Q) -> Q:
    x %= 1
    return min(x, 1 - x)


def bad(v: int, t: Q) -> bool:
    return min(norm(v * t), norm(v * (t + Q(1, 2)))) < RADIUS


def weight(v: int) -> int:
    return 2 if v & 1 else 1


def walls(v: int, a: Q, b: Q) -> set[Q]:
    number = weight(v) * v
    half_width = Q(1, 14 * v)
    lo = (number * (a - half_width)).numerator // (number * (a - half_width)).denominator - 2
    hi = (number * (b + half_width)).numerator // (number * (b + half_width)).denominator + 3
    answer: set[Q] = set()
    for j in range(lo, hi + 1):
        centre = Q(j, number)
        for x in (centre - half_width, centre + half_width):
            if a <= x <= b:
                answer.add(x)
    return answer


def arrangement(speeds: tuple[int, ...], a: Q, b: Q) -> list[Q]:
    points = {a, b}
    for v in speeds:
        points.update(walls(v, a, b))
    return sorted(points)


def occupied(speeds: tuple[int, ...], a: Q, b: Q) -> Q:
    cuts = arrangement(speeds, a, b)
    total = Q(0)
    for left, right in zip(cuts, cuts[1:]):
        if any(bad(v, (left + right) / 2) for v in speeds):
            total += right - left
    return total


def survivor(speeds: tuple[int, ...], a: Q, b: Q) -> Q | None:
    for point in arrangement(speeds, a, b):
        if all(not bad(v, point) for v in speeds):
            return point
    return None


def is_open_cover(speeds: tuple[int, ...], a: Q, b: Q) -> bool:
    cuts = arrangement(speeds, a, b)
    probes = cuts + [(x + y) / 2 for x, y in zip(cuts, cuts[1:])]
    return all(any(bad(v, point) for v in speeds) for point in probes)


def even_denominator(x: Q) -> int:
    return x.denominator if x.denominator % 2 == 0 else 2 * x.denominator


def mod_distance(x: int, n: int) -> int:
    x %= n
    return min(x, n - x)


def main() -> None:
    # Original two-phase tests on every core arrangement cell and wall.
    core_tests = 0
    for m, (a, b, _) in ROWS.items():
        speeds = tuple(range(1, m + 1))
        cuts = arrangement(speeds, a, b)
        probes = cuts + [(x + y) / 2 for x, y in zip(cuts, cuts[1:])]
        for point in probes:
            for v in speeds:
                check(not bad(v, point), f"core interval failed AP{m},v={v}")
                core_tests += 1

    # Independent cell-sweep replay of the sharp discrepancy inequality.
    window_tests = 0
    equality_tests = 0
    for v in range(1, 61):
        w = weight(v)
        period = Q(1, w * v)
        density = Q(w, 7)
        exact_length = density * period
        exact_occupancy = occupied((v,), -exact_length / 2, exact_length / 2)
        exact_bound = density * exact_length + Q(7 - w, 49 * v)
        check(exact_occupancy == exact_bound, f"sharpness failed at v={v}")
        equality_tests += 1
        for denominator in range(7, 24):
            for numerator in range(denominator):
                a = Q(numerator, denominator)
                length = Q((5 * numerator + 3) % denominator + 1, 4 * denominator)
                b = a + length
                upper = density * length + Q(7 - w, 49 * v)
                check(occupied((v,), a, b) <= upper, f"cell discrepancy failed at v={v}")
                window_tests += 1

    # Exhaustive small banks above every displayed threshold.
    bank_rows = 0
    bank_combinations = 0
    owner_tests = 0
    witnesses: list[tuple[int, int, tuple[int, ...], Q, int]] = []
    for (m, odd_count), floor in sorted(FLOORS.items()):
        a, b, q0 = ROWS[m]
        length = b - a
        count = 11 - m
        candidates = range(floor, floor + 2 * count + 4)
        row_witness: tuple[int, ...] | None = None
        for speeds in combinations(candidates, count):
            if sum(v & 1 for v in speeds) != odd_count:
                continue
            total_weight = sum(weight(v) for v in speeds)
            boundary = sum((Q(7 - weight(v), v) for v in speeds), Q(0))
            check(total_weight < 7, "admitted bank crossed the weighted ceiling")
            check(boundary < 7 * length * (7 - total_weight), "uniform threshold failed")
            check(occupied(speeds, a, b) < length, "bank row covers despite strict budget")
            point = survivor(speeds, a, b)
            check(point is not None, "bank row has no arrangement survivor")
            body = tuple(range(1, m + 1)) + speeds
            check(all(not bad(v, point) for v in body), "survivor fails original theta test")
            n = even_denominator(point)
            check(n <= max(q0, 14 * max(speeds)), "even clock bound failed")
            if row_witness is None:
                row_witness = speeds
                label = (point * n).numerator % n
                other = (label + n // 2) % n
                for z in range(1, 2 * n, 2):
                    check(
                        not (7 * mod_distance(z * label, n) < n and 7 * mod_distance(z * other, n) < n),
                        "complementary-label owner test failed",
                    )
                    owner_tests += 1
                witnesses.append((m, odd_count, speeds, point, n))
            bank_combinations += 1
        check(row_witness is not None, f"empty threshold replay AP{m}/o{odd_count}")
        bank_rows += 1

    # Literal wall-point tests distinguish open coverage from a.e. coverage.
    hostile_tests = 0
    for (m, odd_count), speeds in NEGATIVE_ROWS.items():
        a, b, _ = ROWS[m]
        check(sum(v & 1 for v in speeds) == odd_count, "negative parity changed")
        check(is_open_cover(speeds, a, b), f"negative row no longer literally covers AP{m}")
        hostile_tests += 1
    ceiling = (8, 9, 11, 13)
    check(sum(weight(v) for v in ceiling) == 7, "ceiling weight changed")
    check(is_open_cover(ceiling, ROWS[7][0], ROWS[7][1]), "weighted-ceiling control failed")

    print("THM-4092 independent original-theta cell audit: PASS")
    print(f"core tests={core_tests}; window tests={window_tests}; sharp equalities={equality_tests}")
    print(f"threshold banks={bank_rows}; combinations={bank_combinations}; owner tests={owner_tests}")
    print(f"literal low-height covers={hostile_tests}; weighted-ceiling cover={ceiling}")
    print("independent witnesses:")
    for m, odd_count, speeds, point, n in witnesses:
        print(f"  AP{m}/o{odd_count} speeds={speeds} theta={point} N={n}")


if __name__ == "__main__":
    main()
