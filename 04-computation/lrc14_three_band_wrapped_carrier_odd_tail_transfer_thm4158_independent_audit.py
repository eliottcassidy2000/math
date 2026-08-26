#!/usr/bin/env python3
"""Independent exact controls for THM-4158's wrapped carrier.

This implementation imports no primary code. It reconstructs the complete
common-safe alphabet, the full m=7 Haar-safe set, the counted families, and
the asymptotic density from normalized rational arithmetic.
"""

from fractions import Fraction as Q
from math import comb


DELTA = Q(1, 14)
HAAR_GATE = Q(4, 63)


def require(condition, label):
    if not condition:
        raise RuntimeError("FAIL: " + str(label))


def gap(speed, phase):
    residue = (speed * phase) % 1
    return min(residue, 1 - residue)


def endpoints(m):
    return Q(1, 14 * m), Q(8, 7 * (12 * m + 1))


def formula_pool(m):
    rows = []
    for k in range(3 if m >= 2 else 4):
        lo = m * (14 * k + 1)
        hi = ((12 * m + 1) * (14 * k + 13)) // 16
        rows.extend(range(lo, hi + 1))
    return tuple(rows)


def interval_safe(m, speed):
    left, right = endpoints(m)
    lower = speed * left
    upper = speed * right
    k = lower.numerator // lower.denominator
    return (
        k == upper.numerator // upper.denominator
        and lower - k >= DELTA
        and upper - k <= 1 - DELTA
    )


def exhaustive_pool(m):
    # If h*|J|>6/7, the connected image cannot fit in one closed safe band.
    width = endpoints(m)[1] - endpoints(m)[0]
    cap_fraction = Q(6, 7) / width
    cap = cap_fraction.numerator // cap_fraction.denominator
    return tuple(h for h in range(1, cap + 1) if interval_safe(m, h))


def safe_components(speeds):
    walls = {Q(0), Q(1)}
    for speed in speeds:
        for integer in range(speed):
            walls.add(Q(14 * integer + 1, 14 * speed))
            walls.add(Q(14 * integer + 13, 14 * speed))
    walls = sorted(walls)
    cells = []
    for left, right in zip(walls, walls[1:]):
        middle = (left + right) / 2
        if all(gap(speed, middle) >= DELTA for speed in speeds):
            require(
                all(gap(speed, left) >= DELTA for speed in speeds),
                ("left endpoint", left),
            )
            require(
                all(gap(speed, right) >= DELTA for speed in speeds),
                ("right endpoint", right),
            )
            cells.append((left, right))
    merged = []
    for left, right in cells:
        if merged and merged[-1][1] == left:
            merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    return tuple(merged)


def count_at_bound(m, bound):
    return sum(1 for speed in formula_pool(m) if speed <= bound)


def family_count(bound):
    return sum(
        comb(count_at_bound(m, bound) - 1, 10)
        for m in range(2, bound + 1)
        if count_at_bound(m, bound) >= 11
    )


def integrate_power(constant, slope, left, right):
    return (
        (constant + slope * right) ** 11
        - (constant + slope * left) ** 11
    ) / (11 * slope)


for first in range(1, 251):
    left, right = endpoints(first)
    require(
        right - left == Q(4 * first - 1, 14 * first * (12 * first + 1)),
        ("width", first),
    )
    require(exhaustive_pool(first) == formula_pool(first), ("complete pool", first))

P7 = tuple(range(7, 70)) + tuple(range(105, 144)) + tuple(range(203, 218))
require(formula_pool(7) == P7 and len(P7) == 117, "P7")
components = safe_components(P7)
require(
    components
    == (
        (Q(1, 98), Q(13, 966)),
        (Q(953, 966), Q(97, 98)),
    ),
    "P7 full safe set",
)
haar = sum((right - left for left, right in components), Q(0))
require(haar == Q(22, 3381) < HAAR_GATE, "P7 Haar hostile")

expected_counts = {
    20: 75_582,
    40: 812_850_987,
    80: 3_595_550_244_611,
    120: 397_529_462_747_261,
    160: 10_616_582_432_233_990,
    200: 132_777_517_674_540_845,
}
actual_counts = {bound: family_count(bound) for bound in expected_counts}
require(actual_counts == expected_counts, "finite family counts")

pieces = (
    (Q(0), Q(4, 123), Q(0), Q(63, 4)),
    (Q(4, 123), Q(1, 29), Q(1), Q(-15)),
    (Q(1, 29), Q(4, 81), Q(0), Q(14)),
    (Q(4, 81), Q(1, 15), Q(1), Q(-25, 4)),
    (Q(1, 15), Q(4, 39), Q(0), Q(35, 4)),
    (Q(4, 39), Q(1), Q(1), Q(-1)),
)
density = 11 * sum(
    (
        integrate_power(constant, slope, left, right)
        for left, right, constant, slope in pieces
    ),
    Q(0),
)
expected_density = Q(
    848953086913769850118498851618778832628468542103282298683365532079,
    2481088067163593416217816176836483026480276818419826456353950662656,
)
require(density == expected_density, "density")

anchors = {7, 120, 126, 143}
require(anchors <= set(P7), "anchors in P7")
require(comb(len(P7) - len(anchors), 7) == 38_620_298_376, "m7 family")
for modulus in range(2, 15):
    require(
        any(anchor % modulus == 0 for anchor in anchors),
        ("divisor owner", modulus),
    )

print("THM4158 INDEPENDENT EXACT AUDIT")
print("result=ACCEPT")
print("complete_pool_reconstruction=m1..250")
print("m1_pool=[1,10]|[15,21]|[29,33]|[43,44]")
print("m7_pool=[7,69]|[105,143]|[203,217];size=117")
print("m7_safe_components=[1/98,13/966]|[953/966,97/98]")
print("m7_haar=22/3381<4/63")
print("m7_anchored_families=38620298376")
print("finite_counts=" + ",".join(f"{n}:{value}" for n, value in actual_counts.items()))
print("density=" + str(density))
