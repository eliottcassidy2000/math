#!/usr/bin/env python3
"""Exact audit of the cell estimate used in THM-1195.

The audited family is the 13-speed tight family
    {1,2,...,11,13,24}.
On the consecutive-zero cell [1/24,1/13], the lower envelope has only
three active affine pieces.  All arithmetic below is rational.
"""

from fractions import Fraction as Q
from itertools import combinations


V = tuple(range(1, 12)) + (13, 24)
LEFT, RIGHT = Q(1, 24), Q(1, 13)


def circle_distance(v: int, t: Q) -> Q:
    residue = (v * t) % 1
    return min(residue, 1 - residue)


def lower_envelope(t: Q) -> Q:
    return min(circle_distance(v, t) for v in V)


def arrangement_points() -> list[Q]:
    points = {LEFT, RIGHT}
    for v in V:
        points.update(
            Q(k, v) for k in range(v + 1) if LEFT <= Q(k, v) <= RIGHT
        )
        points.update(
            Q(2 * k + 1, 2 * v)
            for k in range(v)
            if LEFT <= Q(2 * k + 1, 2 * v) <= RIGHT
        )
    for a, b in combinations(V, 2):
        for denominator in (a + b, abs(a - b)):
            if denominator:
                points.update(
                    Q(k, denominator)
                    for k in range(denominator + 1)
                    if LEFT <= Q(k, denominator) <= RIGHT
                )
    return sorted(points)


points = arrangement_points()
raw_segments = []
area = Q(0)
height = Q(0)
for x, y in zip(points, points[1:]):
    gx, gy = lower_envelope(x), lower_envelope(y)
    midpoint = (x + y) / 2
    values = {v: circle_distance(v, midpoint) for v in V}
    active = min(v for v, value in values.items() if value == min(values.values()))
    slope = (gy - gx) / (y - x)
    raw_segments.append((x, y, active, slope, gx, gy))
    area += (gx + gy) * (y - x) / 2
    height = max(height, gx, gy)

# Merge adjacent arrangement pieces carrying the same affine function.
segments = []
for x, y, active, slope, gx, gy in raw_segments:
    if segments and segments[-1][2] == active and segments[-1][3] == slope:
        old = segments[-1]
        segments[-1] = (old[0], y, active, slope, old[4], gy)
    else:
        segments.append((x, y, active, slope, gx, gy))

expected = [
    (Q(1, 24), Q(1, 23), 24, Q(24), Q(0), Q(1, 23)),
    (Q(1, 23), Q(1, 14), 1, Q(1), Q(1, 23), Q(1, 14)),
    (Q(1, 14), Q(1, 13), 13, Q(-13), Q(1, 14), Q(0)),
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


require(segments == expected, "unexpected lower-envelope decomposition")
require(
    all(a[3] >= b[3] for a, b in zip(segments, segments[1:])),
    "cell slopes are not nonincreasing",
)

cell_length = RIGHT - LEFT
claimed_upper_bound = height * cell_length / 2
require(area == Q(185, 100464), "unexpected exact cell area")
require(height == Q(1, 14), "unexpected cell height")
require(claimed_upper_bound == Q(11, 8736), "unexpected former tent bound")
require(
    area - claimed_upper_bound == Q(3, 5152) > 0,
    "counterexample surplus did not verify",
)
require(
    area / claimed_upper_bound == Q(370, 253),
    "unexpected counterexample ratio",
)

print("THM-1195 CELL-DIRECTION AUDIT (exact rational arithmetic)")
print(f"speeds={V}")
print(f"consecutive-zero cell=[{LEFT},{RIGHT}], length={cell_length}")
print("lower-envelope pieces: left right active_speed slope g(left) g(right)")
for segment in segments:
    print(*segment)
print(f"slopes_nonincreasing={all(a[3] >= b[3] for a, b in zip(segments, segments[1:]))}")
print(f"cell_area={area}")
print(f"cell_height={height}")
print(f"height_times_length_over_2={claimed_upper_bound}")
print(f"area_minus_claimed_upper_bound={area - claimed_upper_bound}")
print(f"area_over_claimed_upper_bound={area / claimed_upper_bound}")
print("VERIFIED: the claimed upper bound has the wrong direction on this 13-speed cell")
