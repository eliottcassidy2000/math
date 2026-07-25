#!/usr/bin/env python3
"""Independent interval-geometry referee for THM-2255.

This replay does not call the folded overlap formula.  It constructs every
danger comb as an exact union of rational intervals and directly intersects
the 70 reduced pairs in the finite bank ab<=345.
"""

from fractions import Fraction as Q
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def merge(intervals: list[tuple[Q, Q]]) -> list[tuple[Q, Q]]:
    ordered = sorted(intervals)
    result: list[tuple[Q, Q]] = []
    for left, right in ordered:
        if not result or left > result[-1][1]:
            result.append((left, right))
        elif right > result[-1][1]:
            result[-1] = (result[-1][0], right)
    return result


def danger_intervals(speed: int) -> list[tuple[Q, Q]]:
    radius = Q(1, 14 * speed)
    pieces: list[tuple[Q, Q]] = []
    for index in range(speed):
        center = Q(index, speed)
        left = center - radius
        right = center + radius
        if left < 0:
            pieces.append((Q(0), right))
            pieces.append((Q(1) + left, Q(1)))
        elif right > 1:
            pieces.append((left, Q(1)))
            pieces.append((Q(0), right - 1))
        else:
            pieces.append((left, right))
    result = merge(pieces)
    require(sum(right - left for left, right in result) == Q(1, 7),
            f"danger mass failed at speed {speed}")
    return result


def intersection_mass(a: int, b: int) -> Q:
    left_intervals = danger_intervals(a)
    right_intervals = danger_intervals(b)
    i = 0
    j = 0
    total = Q(0)
    while i < len(left_intervals) and j < len(right_intervals):
        left = max(left_intervals[i][0], right_intervals[j][0])
        right = min(left_intervals[i][1], right_intervals[j][1])
        if right > left:
            total += right - left
        if left_intervals[i][1] < right_intervals[j][1]:
            i += 1
        elif right_intervals[j][1] < left_intervals[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return total


CAP = Q(25, 1183)
DELTA = Q(961, 6930)
bank: list[tuple[int, int, Q]] = []
for a in range(1, 346):
    for b in range(a + 1, 346):
        if a * b > 345:
            break
        if gcd(a, b) != 1:
            continue
        if (a % 13 == 0) == (b % 13 == 0):
            continue
        bank.append((a, b, intersection_mass(a, b)))

require(len(bank) == 70, "independent bank-size mismatch")
maximum = max(value for _, _, value in bank)
equalities = [(a, b) for a, b, value in bank if value == maximum]
require(maximum == CAP, "independent cap mismatch")
require(equalities == [(1, 169)], "independent equality mismatch")

# Independent common-dilate checks use interval geometry at the unreduced speeds.
for scale in (1, 2, 3, 5, 8):
    require(intersection_mass(scale, 169 * scale) == CAP,
            f"independent scaled equality failed at scale {scale}")

strict_floor = DELTA - 3 * CAP
double_floor = DELTA - Q(1, 14) - 2 * CAP
require(strict_floor == Q(88159, 1171170), "strict floor replay mismatch")
require(double_floor == Q(14627, 585585), "double floor replay mismatch")

print("THM-2255 independent exact interval referee")
print(f"directly intersected finite pairs: {len(bank)}")
print(f"direct interval maximum: {maximum}")
print(f"direct reduced equality classes: {equalities}")
print(f"strict-depth exclusive-owner floor: {strict_floor}")
print(f"double-depth exclusive-owner floor: {double_floor}")
print("all independent checks passed")
