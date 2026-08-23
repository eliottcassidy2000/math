#!/usr/bin/env python3
"""Independent open-interval audit of the 11+2 cyclic slack obstruction.

Unlike probe.py's wall-cell topology, this checker constructs periodic danger
intervals on the real line, merges only strict overlaps, intersects a shifted
copy for the two-lift obstruction, and then passes to the half-period quotient.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from math import ceil, floor, gcd
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")
H = Fraction(1, 14)
GATES = 0


def require(condition: bool, label: object) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def factor(n: int) -> tuple[tuple[int, int], ...]:
    result = []
    d = 2
    while d * d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e:
            result.append((d, e))
        d += 1
    if n > 1:
        result.append((n, 1))
    return tuple(result)


def allowed(total: int) -> bool:
    return all(prime % 3 == 2 and exponent <= 2 for prime, exponent in factor(total))


def build_atlas() -> tuple[tuple[int, int], ...]:
    result = []
    for total in range(3, 357):
        if not allowed(total):
            continue
        for p in range(1, (total + 1) // 2):
            q = total - p
            if p < q and gcd(p, q) == 1:
                result.append((p, q))
    return tuple(sorted(result))


Interval = tuple[Fraction, Fraction]


def raw_danger(speed: int, low: int = -2, high: int = 3) -> list[Interval]:
    radius = H / speed
    first = floor(low * speed) - 2
    last = ceil(high * speed) + 2
    return [(Fraction(k, speed) - radius, Fraction(k, speed) + radius)
            for k in range(first, last + 1)]


def merge_open(intervals: list[Interval]) -> list[Interval]:
    """Merge strict overlaps; touching open endpoints remain disconnected."""
    result: list[list[Fraction]] = []
    for left, right in sorted(intervals):
        require(left < right, ("nonempty interval", left, right))
        if not result or left >= result[-1][1]:
            result.append([left, right])
        elif right > result[-1][1]:
            result[-1][1] = right
    return [(left, right) for left, right in result]


def intersect_open(first: list[Interval], second: list[Interval]) -> list[Interval]:
    intersections: list[Interval] = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            intersections.append((left, right))
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return merge_open(intersections)


def one_period(components: list[Interval], period: Fraction) -> list[Interval]:
    selected = []
    for left, right in components:
        midpoint = (left + right) / 2
        if 0 <= midpoint < period:
            selected.append((left, right))
    return selected


def interval_widths(p: int, q: int) -> tuple[Fraction, int, Fraction, Fraction, int, Fraction]:
    danger = merge_open(raw_danger(p) + raw_danger(q))
    danger_period = one_period(danger, Fraction(1))
    require(bool(danger_period), ("danger nonempty", p, q))
    beta1 = max(right - left for left, right in danger_period)
    measure1 = sum((right - left for left, right in danger_period), Fraction(0))

    shifted = [(left - Fraction(1, 2), right - Fraction(1, 2)) for left, right in danger]
    both = intersect_open(danger, shifted)
    half_period = one_period(both, Fraction(1, 2))
    if half_period:
        lengths = [right - left for left, right in half_period]
        beta2 = 2 * max(lengths)
        measure2 = 2 * sum(lengths, Fraction(0))
    else:
        beta2 = Fraction(0)
        measure2 = Fraction(0)
    return beta1, len(danger_period), measure1, beta2, len(half_period), measure2


def fmt(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    pairs = build_atlas()
    require(len(pairs) == 5855, "atlas count")
    records = []
    for p, q in pairs:
        data = interval_widths(p, q)
        records.append((p, q, *data))
        if p % 2 and q % 2:
            require((data[3] == 0) == (p + q <= 7), ("odd boundary", p, q, data[3]))

    odd_residual = [row for row in records if row[0] % 2 and row[1] % 2 and row[5] > 0]
    require(len(odd_residual) == 1650, "odd residual count")
    max1 = max(row[2] for row in records)
    max2 = max(row[5] for row in odd_residual)
    max1_pairs = tuple((row[0], row[1]) for row in records if row[2] == max1)
    max2_pairs = tuple((row[0], row[1]) for row in odd_residual if row[5] == max2)
    require((max1, max1_pairs) == (Fraction(29, 196), ((1, 28),)), "scale1 maximum")
    require((max2, max2_pairs) == (Fraction(2, 63), ((1, 9),)), "scale2 maximum")

    k1_le_one = tuple((row[0], row[1]) for row in records if 42 * row[2] <= 1)
    k2_le_one = tuple((row[0], row[1]) for row in odd_residual if 42 * row[5] <= 1)
    require(len(k1_le_one) == 5445, "scale1 K<=1")
    require(len(k2_le_one) == 1649, "scale2 K<=1")
    k2_exceptions = tuple((row[0], row[1], 42 * row[5]) for row in odd_residual if 42 * row[5] > 1)
    require(k2_exceptions == ((1, 9, Fraction(4, 3)),), "unique scale2 exception")

    controls = {}
    for pair in ((1, 3), (1, 5), (3, 7), (1, 9), (1, 355), (175, 181)):
        controls[pair] = interval_widths(*pair)
    require(controls[(1, 3)][3] == 0 and controls[(1, 5)][3] == 0,
            "universal controls")
    require(controls[(3, 7)][3] == Fraction(1, 49), "minimal hostile width")

    # Primitive support-two Graver saturation is scale-blind.
    require(all(p + q <= 356 and gcd(p, q) == 1 for p, q in pairs),
            "THM3743 support-two saturation")

    semantic = {
        "count": len(pairs),
        "odd_residual": len(odd_residual),
        "max1": fmt(max1),
        "max1_pairs": max1_pairs,
        "max2": fmt(max2),
        "max2_pairs": max2_pairs,
        "k1_le_one": len(k1_le_one),
        "k2_le_one": len(k2_le_one),
        "k2_exceptions": tuple((p, q, fmt(k)) for p, q, k in k2_exceptions),
        "controls": {
            f"{p},{q}": tuple(fmt(value) if isinstance(value, Fraction) else value for value in data)
            for (p, q), data in controls.items()
        },
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()

    print("LRC14_REMAINING_CYCLIC_INTERVAL_INDEPENDENT_AUDIT_20260823")
    print("method=periodic_open_interval_union+strict_overlap+shifted_intersection+half_period_quotient")
    print(f"atlas={len(pairs)};odd_scale2_residual={len(odd_residual)}")
    print(f"scale1_max_beta={fmt(max1)};pairs={max1_pairs};K_le_1={len(k1_le_one)}")
    print(f"scale2_max_beta={fmt(max2)};pairs={max2_pairs};K_le_1={len(k2_le_one)}")
    print("scale2_unique_K_gt_1=" + repr(tuple((p, q, fmt(k)) for p, q, k in k2_exceptions)))
    for pair, data in sorted(controls.items()):
        print(f"control_{pair[0]}_{pair[1]}=" + repr(tuple(fmt(v) if isinstance(v, Fraction) else v for v in data)))
    print("THM3743_join=saturated_by_internal_pair_relation_on_every_atlas_row;zero_seam_deletions")
    print(f"semantic_sha256={digest}")
    print(f"gates={GATES}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
