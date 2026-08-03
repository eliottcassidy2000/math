#!/usr/bin/env python3
"""Exact controls for THM-3211's signed one-over-g correction.

The uniform proof is analytic in the theorem.  This certificate evaluates the
endpoint barycenter exactly, checks every ray through Q=100, and independently
matches the common residue-linear coefficient through Q=30.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from math import gcd

from lrc_uniform_channel_limit_engine_thm3211 import (
    CELL,
    EDGES,
    H,
    L,
    bernoulli_limit,
    cleared_num,
    exact_period_scout,
    interpolate_g_poly,
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


RADIUS = F(1, 14)


def clipped_comb_intervals(P, phase):
    intervals = []
    # The periodic comb makes this deliberately oversized integer bank safe.
    for k in range(-2, P + 3):
        lo = F(k, P) + phase / P - RADIUS / P
        hi = F(k, P) + phase / P + RADIUS / P
        lo, hi = max(lo, F(0)), min(hi, F(1))
        if lo < hi:
            intervals.append((lo, hi))
    intervals.sort()
    return intervals


def overlap_barycenter(P, alpha, Q, beta):
    left = clipped_comb_intervals(P, alpha)
    right = clipped_comb_intervals(Q, beta)
    i = j = 0
    centered = F(0)
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            centered += (hi * hi - lo * lo) / 2 - (hi - lo) / 2
        if left[i][1] < right[j][1]:
            i += 1
        elif right[j][1] < left[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return centered


def correction_barycenter(P, Q, e, f):
    R, S = e * CELL % L, f * CELL % L
    initial = overlap_barycenter(P, F(R, L), Q, F(S, L))
    terminal = overlap_barycenter(P, F(R + e, L), Q, F(S + f, L))
    return terminal - initial


def stabilized_poly_one_residue(P, Q, e, f):
    period = abs(Q * e - P * f) or 1
    points = []
    for j in range(4):
        g = 1 + period * (10000 + j)
        points.append((g, cleared_num(e, g * P, f, g * Q)))
    polynomial = interpolate_g_poly(points[:3])
    x, y = points[3]
    require(
        polynomial[0] * x * x + polynomial[1] * x + polynomial[2] == y,
        ("unstable fourth point", P, Q, e, f, period, points, polynomial),
    )
    return polynomial, period


lane_order = []
for i, j in EDGES:
    lane_order.extend(((H[i], H[j]), (H[j], H[i])))

records = []
counts = {"channels": 0, "rays": 0, "positive": 0, "negative": 0, "zero": 0}
witness = {lane: {"positive": None, "negative": None} for lane in lane_order}
largest = None
for Q in range(2, 101):
    for P in range((Q + 1) // 2, Q):
        if gcd(P, Q) > 1 or P + Q < 8:
            continue
        counts["channels"] += 1
        for e, f in lane_order:
            limit = bernoulli_limit(e, P, f, Q)
            correction = correction_barycenter(P, Q, e, f)
            polynomial, period = stabilized_poly_one_residue(P, Q, e, f)
            d2 = L * L * P * Q
            d1 = -L * (P * f + Q * e)
            require(polynomial[0] == d2 * limit, ("leading mismatch", P, Q, e, f))
            require(
                polynomial[1] == d2 * correction + d1 * limit,
                ("linear mismatch", P, Q, e, f, polynomial[1], correction),
            )
            sign = "positive" if correction > 0 else "negative" if correction < 0 else "zero"
            counts[sign] += 1
            counts["rays"] += 1
            if sign != "zero" and witness[(e, f)][sign] is None:
                witness[(e, f)][sign] = (P, Q, correction)
            row = (abs(correction), P, Q, e, f, correction)
            if largest is None or row[0] > largest[0]:
                largest = row
            records.append((P, Q, e, f, limit, correction, polynomial[0], polynomial[1], period))

expected_counts = {
    "channels": 1519,
    "rays": 27342,
    "positive": 15213,
    "negative": 12129,
    "zero": 0,
}
require(counts == expected_counts, ("count mismatch", counts))
require(
    all(witness[lane]["positive"] is not None and witness[lane]["negative"] is not None for lane in lane_order),
    ("missing mixed-sign lane", witness),
)
require(
    largest == (F(8213, 1411200), 3, 5, 6, 1, F(-8213, 1411200)),
    ("largest correction mismatch", largest),
)

# Exhaust every residue state in a smaller exact bank.  This checks the
# theorem's common linear coefficient against the independent floor-sum path.
residue_classes = 0
for Q in range(2, 31):
    for P in range((Q + 1) // 2, Q):
        if gcd(P, Q) > 1 or P + Q < 8:
            continue
        for e, f in lane_order:
            period, _, polynomials = exact_period_scout(e, P, f, Q, h0=20)
            residue_classes += period
            linear = {row[1] for row in polynomials}
            require(len(linear) == 1, ("residue linear variation", P, Q, e, f, period))
            correction = correction_barycenter(P, Q, e, f)
            limit = bernoulli_limit(e, P, f, Q)
            expected = L * L * P * Q * correction - L * (P * f + Q * e) * limit
            require(next(iter(linear)) == expected, ("residue correction mismatch", P, Q, e, f))

equality = tuple((e, f, correction_barycenter(3, 5, e, f)) for e, f in ((1, 2), (3, 1), (2, 3)))
require(
    equality == ((1, 2, F(71, 264600)), (3, 1, F(23, 12600)), (2, 3, F(1, 5400))),
    ("equality correction mismatch", equality),
)
hostile = correction_barycenter(3, 5, 6, 1)
require(hostile == F(-8213, 1411200), ("hostile mismatch", hostile))

tree = ast.parse(open(__file__, encoding="utf-8").read(), filename=__file__)
require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), ("assert node", __file__))

digest = sha256("\n".join(map(str, records)).encode()).hexdigest()
print("LRC SIGNED ONE-OVER-G CORRECTION CERTIFICATE")
print(f"counts={counts};record_sha256={digest}")
print(f"exact_residue_classes_Q_le_30={residue_classes}")
print(f"equality_corrections={equality}")
print(f"hostile_correction={hostile};largest_abs_Q_le_100={largest}")
print(f"mixed_sign_lanes={len(lane_order)}")
print("FAILED_CHECKS=NONE")
