#!/usr/bin/env python3
"""No-import exact H99 audit for the reserved THM-4431 candidate.

This is finite corroboration only.  The all-height colored-basis and
three-direction classification are proved separately in the audit report.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import gcd
import sys


CHECKS = 0
TARGET = Q(6, 77)


def need(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(payload)


def strict_bound(numerator):
    return (numerator - 1) // 14


def carriers_middle(w):
    """Range the first/third coordinates and solve exactly for the middle."""
    a, b, c = w
    bx, by, bz = (strict_bound(3 * (b + c)),
                  strict_bound(3 * (a + c)),
                  strict_bound(3 * (a + b)))
    result = set()
    for x in range(-bx, bx + 1):
        if x % 3 == 0:
            continue
        for z in range(-bz, bz + 1):
            if z % 3 == 0:
                continue
            numerator = -a * x - c * z
            if numerator % b:
                continue
            y = numerator // b
            if abs(y) <= by and y % 3:
                result.add((x, y, z))
    return result


def primitive_ray(point):
    divisor = gcd(gcd(abs(point[0]), abs(point[1])), abs(point[2]))
    ray = tuple(value // divisor for value in point)
    return ray if next(value for value in ray if value) > 0 else tuple(-value for value in ray)


def determinant(left, right):
    return left[0] * right[1] - left[1] * right[0]


def owner_color(w, ray):
    residues = tuple(w[index] * ray[index] % 3 for index in range(3))
    need(len(set(residues)) == 1 and residues[0], ("owner color", w, ray, residues))
    return residues[0]


def oriented_color_one(w, ray):
    return ray if owner_color(w, ray) == 1 else tuple(-value for value in ray)


def classify_three(w, directions):
    native = 3 * w[2]
    basis = next((pair for pair in combinations(directions, 2)
                  if abs(determinant(*pair)) == native), None)
    need(basis is not None, ("native Gamma basis", w, directions))
    u, v = (oriented_color_one(w, ray) for ray in basis)
    third = next(ray for ray in directions
                 if ray not in (u, tuple(-x for x in u), v, tuple(-x for x in v)))
    denominator = determinant(u, v)
    m = Q(determinant(third, v), denominator)
    n = Q(determinant(u, third), denominator)
    need(m.denominator == n.denominator == 1, ("integral basis coordinates", w, m, n))
    m, n = int(m), int(n)
    need(gcd(abs(m), abs(n)) == 1 and (m + n) % 3,
         ("primitive live third coordinate", w, m, n))
    need(abs(m) == 1 or abs(n) == 1, ("unit-coordinate classification", w, m, n))
    if abs(n) != 1:
        m, n = n, m
    if n < 0:
        m, n = -m, -n
    need(n == 1 and m in (1, -2), ("two normal forms", w, m, n))
    return m


def projections(w, points):
    a, b, c = w
    denominator = 14 * a * b * c
    cap = 6 * a * b
    totals = []
    for index, speed in enumerate(w):
        roof_sum = 0
        for point in points:
            roof = 3 * (sum(w) - speed) - 14 * abs(point[index])
            need(roof > 0, ("positive roof", w, point, index))
            roof_sum += min(cap, speed * roof)
        totals.append(Q(roof_sum, denominator))
    return tuple(totals)


def analyze_three(w, points, directions):
    shape = classify_three(w, directions)
    maxima = tuple(max(abs(value) for value in ray) for ray in directions)
    for ray, maximum in zip(directions, maxima):
        norm = sum(abs(value) for value in ray)
        need(norm % 2 == 0 and norm >= 16 and maximum >= 7,
             ("short-relation floor", w, ray, norm, maximum))
    threshold = Q(3 * w[2], 2)
    need(all(left * right >= threshold for left, right in combinations(maxima, 2)),
         ("pairwise determinant products", w, maxima, threshold))
    reciprocal = sum((Q(1, maximum) for maximum in maxima), Q(0))
    need(reciprocal <= Q(1, 7) + Q(14, threshold) or reciprocal * reciprocal <= Q(9, threshold),
         ("reciprocal dichotomy", w, maxima, reciprocal))
    values = projections(w, points)
    need(max(values) < Q(12, 49) * reciprocal + Q(12, 7 * w[2]),
         ("three-ray count envelope", w, values, reciprocal))
    need(max(values) < TARGET, ("finite three-ray target", w, values))
    return shape, values, maxima


def main():
    sys.stdout.reconfigure(newline="\n")
    values = tuple(value for value in range(1, 99, 2) if value % 3)
    universe = tuple(w for w in combinations(values, 3)
                     if gcd(gcd(w[0], w[1]), w[2]) == 1)
    need(len(universe) == 5409, ("H99 universe", len(universe)))
    multi = selected = 0
    types = Counter()
    leader = (Q(0), None, None)
    digest = sha256()
    for w in universe:
        points = carriers_middle(w)
        directions = tuple(sorted({primitive_ray(point) for point in points}))
        if len(directions) < 2:
            continue
        multi += 1
        need(any(abs(determinant(left, right)) == 3 * w[2]
                 for left, right in combinations(directions, 2)),
             ("finite native basis", w, directions))
        if len(directions) != 3:
            continue
        selected += 1
        shape, values_row, maxima = analyze_three(w, points, directions)
        types[shape] += 1
        candidate = (max(values_row), w, values_row)
        if candidate > leader:
            leader = candidate
        digest.update(repr((w, tuple(sorted(points)), directions, shape, values_row, maxima)).encode("ascii"))

    need(multi == 3500, ("H99 multi-direction count", multi))
    need(selected == 1791 and types == Counter({1: 1107, -2: 684}),
         ("H99 exactly-three count/types", selected, types))
    need(leader == (Q(18, 301), (5, 37, 43),
                    (Q(240, 11137), Q(18, 301), Q(2822, 55685))),
         ("H99 leader", leader))

    linear_99 = Q(12, 343) + Q(4, 99)
    need(linear_99 == Q(2560, 33957) < TARGET,
         ("linear tail branch", linear_99))
    need(Q(96, 26411) < Q(4, 1089), "radical tail branch after positive squaring")

    wide = []
    for w in ((19, 23, 29), (5, 191, 199), (7, 611, 613)):
        points = carriers_middle(w)
        directions = tuple(sorted({primitive_ray(point) for point in points}))
        need(any(abs(determinant(left, right)) == 3 * w[2]
                 for left, right in combinations(directions, 2)),
             ("wide native basis", w, directions))
        result = analyze_three(w, points, directions) if len(directions) == 3 else None
        wide.append((w, len(points), len(directions), result))
    need(wide[1][1] == 16 and wide[1][2] == 3, ("unbounded-count three-ray hostile", wide[1]))
    need(wide[2][2] == 13, ("four-plus-direction frontier", wide[2]))

    print("RESERVED THM-4431 H99 -- NO-IMPORT MIDDLE-COORDINATE AUDIT")
    print("status=FINITE_EXACT_PASS; all_height_requires_separate_colored_basis_and_classification_proof")
    print("universe_rows=%d; multi_direction_rows=%d; exactly_three_direction_rows=%d" %
          (len(universe), multi, selected))
    print("normal_forms=%s" % dict(sorted(types.items())))
    print("leader=%s" % (leader,))
    print("linear_tail_at_99=%s; radical_squared=%s<%s; target=%s" %
          (linear_99, Q(96, 26411), Q(4, 1089), TARGET))
    for row in wide:
        print("wide_control=%s" % (row,))
    print("semantic_sha256=" + digest.hexdigest())
    print("checks=%d" % CHECKS)
    print("verdict=PASS_FINITE_ONLY; four_or_more_directions=OPEN; LRC14=OPEN")


if __name__ == "__main__":
    main()
