#!/usr/bin/env python3
"""Exact standalone extension of the one-zero-mod-3 relation atlas to l1=20.

No repository mathematics is imported.  All proof gates use explicit
RuntimeError checks and therefore remain live under ``python -O``.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd
import sys


Q = Fraction
R = Q(3, 14)
CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError("CHECK FAILED: " + label)


def gcd_all(values):
    answer = 0
    for value in values:
        answer = gcd(answer, abs(value))
    return answer


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def normalize_ray(v):
    g = gcd_all(v)
    check(g > 0, "nonzero ray")
    v = tuple(x // g for x in v)
    if next(x for x in v if x) < 0:
        v = tuple(-x for x in v)
    return v


def signed_rays(pattern):
    answer = set()
    for magnitudes in set(permutations(pattern)):
        for signs in product((-1, 1), repeat=3):
            c = tuple(a * s for a, s in zip(magnitudes, signs))
            if min(c) < 0 < max(c):
                answer.add(normalize_ray(c))
    return tuple(sorted(answer))


def patterns_through(limit, zero_count):
    answer = []
    for norm in range(4, limit + 1, 2):
        for a in range(1, norm):
            for b in range(a, norm):
                c = norm - a - b
                if c < b or c <= 0:
                    continue
                if gcd_all((a, b, c)) != 1:
                    continue
                if sum(x % 3 == 0 for x in (a, b, c)) != zero_count:
                    continue
                answer.append((a, b, c))
    return tuple(answer)


def floor_q(x):
    return x.numerator // x.denominator


def ceil_q(x):
    return -floor_q(-x)


def egcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x - (a // b) * y


def bezout(v):
    a, b, c = v
    g12, x, y = egcd(abs(a), abs(b))
    g, s, t = egcd(g12, abs(c))
    answer = (
        x * s * (1 if a >= 0 else -1),
        y * s * (1 if b >= 0 else -1),
        t * (1 if c >= 0 else -1),
    )
    check(dot(answer, v) == g, "Bezout identity")
    return g, answer


def relation_defects(c):
    norm = sum(abs(x) for x in c)
    bound = R * norm
    unit = all(x % 3 for x in c)
    answer = []
    for delta in range(-ceil_q(bound), ceil_q(bound) + 1):
        if abs(delta) >= bound:
            continue
        if unit and delta % 3 == 0:
            answer.append(delta)
        if not unit and delta % 3 != 0:
            answer.append(delta)
    return tuple(answer)


def roof(w, C):
    a, b, c = w
    A, B, Cc = C
    return max(Q(0), min(
        2 * R / a,
        2 * R / b,
        2 * R / c,
        R / a + R / b - Q(abs(Cc), a * b),
        R / a + R / c - Q(abs(B), a * c),
        R / b + R / c - Q(abs(A), b * c),
    ))


def integer_open_interval(lo, hi):
    if lo >= hi:
        return range(0)
    return range(floor_q(lo) + 1, ceil_q(hi))


def quotient_components(w, c):
    """Enumerate the exact affine carrier fibres C_delta+k*c."""
    check(dot(w, c) == 0, "relation annihilates speed")
    g, v = bezout(c)
    check(g == 1, "primitive relation")
    expected_live = 2 if all(x % 3 for x in c) else 1
    components = {}
    masses = {}
    for delta in relation_defects(c):
        C0 = cross(w, tuple(delta * x for x in v))
        check(dot(C0, w) == 0, "base carrier in raw lattice")
        residue_rows = tuple(
            tuple((C0[i] + k * c[i]) % 3 for i in range(3))
            for k in range(3)
        )
        check(sum(all(x for x in row) for row in residue_rows) == expected_live,
              "owner-live affine classes " + repr((w, c, delta)))

        lower = None
        upper = None
        for i in range(3):
            j, k = tuple(index for index in range(3) if index != i)
            bound = R * (w[j] + w[k])
            lo = (-bound - C0[i]) / c[i]
            hi = (bound - C0[i]) / c[i]
            if lo > hi:
                lo, hi = hi, lo
            lower = lo if lower is None else max(lower, lo)
            upper = hi if upper is None else min(upper, hi)

        mass = Q(0)
        for shift in integer_open_interval(lower, upper):
            C = tuple(C0[i] + shift * c[i] for i in range(3))
            if any(x % 3 == 0 for x in C):
                continue
            value = roof(w, C)
            if value <= 0:
                continue
            check(C not in components, "disjoint defect fibres")
            components[C] = value
            mass += value
        masses[delta] = mass
    return components, masses


def presentation_universe(pattern, height):
    """Solve every labeled relation for the final sorted coordinate."""
    speeds = tuple(x for x in range(1, height + 1, 2) if x % 3)
    speed_set = set(speeds)
    answer = defaultdict(set)
    for c in signed_rays(pattern):
        for w0 in speeds:
            for w1 in speeds:
                if w1 <= w0:
                    continue
                numerator = -(c[0] * w0 + c[1] * w1)
                if numerator % c[2]:
                    continue
                w2 = numerator // c[2]
                if w2 <= w1 or w2 not in speed_set:
                    continue
                w = (w0, w1, w2)
                if gcd_all(w) != 1:
                    continue
                check(dot(c, w) == 0, "presentation solver")
                answer[w].add(c)
    return {w: tuple(sorted(cs)) for w, cs in sorted(answer.items())}


def brute_universe(pattern, height):
    speeds = tuple(x for x in range(1, height + 1, 2) if x % 3)
    rays = signed_rays(pattern)
    answer = {}
    for w in combinations(speeds, 3):
        if gcd_all(w) != 1:
            continue
        hits = tuple(c for c in rays if dot(c, w) == 0)
        if hits:
            answer[w] = hits
    return answer


@lru_cache(maxsize=None)
def danger_list(a):
    radius = R / a
    return tuple(
        (max(Q(0), Q(n, a) - radius), min(Q(1), Q(n, a) + radius), n)
        for n in range(a + 1)
    )


def intersect_annotated(lists):
    indices = [0, 0, 0]
    while all(indices[i] < len(lists[i]) for i in range(3)):
        current = tuple(lists[i][indices[i]] for i in range(3))
        lo = max(row[0] for row in current)
        hi = min(row[1] for row in current)
        if lo < hi:
            yield lo, hi, tuple(row[2] for row in current)
        first_end = min(row[1] for row in current)
        for i, row in enumerate(current):
            if row[1] == first_end:
                indices[i] += 1


@lru_cache(maxsize=None)
def physical_components(w):
    """Construct the failure comb directly on the physical circle."""
    answer = defaultdict(Fraction)
    pieces = 0
    for lo, hi, n in intersect_annotated(tuple(danger_list(a) for a in w)):
        owners = tuple((-w[i] * n[i]) % 3 for i in range(3))
        if set(owners) != {0, 1, 2}:
            continue
        C = cross(w, n)
        answer[C] += hi - lo
        pieces += 1
    return tuple(sorted(answer.items())), pieces


def physical_dict(w):
    items, pieces = physical_components(tuple(w))
    return dict(items), pieces


def clip(poly, A, B, D):
    if not poly:
        return []
    answer = []
    previous = poly[-1]
    previous_value = A * previous[0] + B * previous[1] - D
    for current in poly:
        current_value = A * current[0] + B * current[1] - D
        previous_inside = previous_value <= 0
        current_inside = current_value <= 0
        if previous_inside != current_inside:
            t = previous_value / (previous_value - current_value)
            answer.append((
                previous[0] + t * (current[0] - previous[0]),
                previous[1] + t * (current[1] - previous[1]),
            ))
        if current_inside:
            answer.append(current)
        previous = current
        previous_value = current_value
    return answer


def area(poly):
    if len(poly) < 3:
        return Q(0)
    twice = sum(
        poly[i][0] * poly[(i + 1) % len(poly)][1]
        - poly[(i + 1) % len(poly)][0] * poly[i][1]
        for i in range(len(poly))
    )
    return abs(twice) / 2


def bulk_polygon(pattern, delta):
    a, b, c = pattern
    square = [(-R, -R), (R, -R), (R, R), (-R, R)]
    square = clip(square, a, b, c * R - delta)
    square = clip(square, -a, -b, c * R + delta)
    return area(square) / c


def bulk_convolution(pattern, delta):
    answer = Q(0)
    for mask in product((0, 1), repeat=3):
        shift = sum(pattern[i] for i in range(3) if mask[i])
        x = abs(delta) + R * sum(pattern) - 2 * R * shift
        answer += (-1 if sum(mask) % 2 else 1) * max(Q(0), x) ** 2
    return answer / (2 * pattern[0] * pattern[1] * pattern[2])


OLD_PATTERNS = (
    (1, 2, 3),
    (1, 1, 6), (1, 3, 4),
    (2, 3, 5),
    (1, 2, 9), (1, 3, 8), (1, 5, 6), (2, 3, 7), (3, 4, 5),
    (1, 1, 12), (1, 3, 10), (1, 4, 9), (1, 6, 7), (3, 4, 7),
)

PATTERNS = (
    (2, 3, 11), (2, 5, 9), (3, 5, 8), (5, 5, 6),
    (1, 2, 15), (1, 3, 14), (1, 5, 12), (1, 6, 11),
    (1, 8, 9), (2, 3, 13), (2, 7, 9), (3, 4, 11),
    (3, 5, 10), (3, 7, 8), (4, 5, 9), (5, 6, 7),
    (1, 1, 18), (1, 3, 16), (1, 4, 15), (1, 6, 13),
    (1, 7, 12), (1, 9, 10), (3, 4, 13), (3, 7, 10),
    (4, 7, 9), (6, 7, 7),
)

EXPECTED = {
    (2, 3, 11): dict(height=31, triples=19, rays=19, positive=15, maximum=Q(58, 833), winner=(5, 7, 17), bulks={1: Q(9, 539), 2: Q(83, 6468)}, baseline=Q(191, 9702), threshold=Q(282744, 8237), masses={-2: Q(8, 833), -1: Q(3, 119), 1: Q(3, 119), 2: Q(8, 833)}),
    (2, 5, 9): dict(height=29, triples=17, rays=19, positive=15, maximum=Q(12, 161), winner=(1, 11, 23), bulks={1: Q(41, 2205), 2: Q(1, 105)}, baseline=Q(124, 6615), threshold=Q(65205, 2122), masses={-2: Q(3, 161), -1: Q(3, 161), 1: Q(3, 161), 2: Q(3, 161)}),
    (3, 5, 8): dict(height=25, triples=11, rays=11, positive=8, maximum=Q(6, 77), winner=(1, 5, 11), bulks={1: Q(221, 11760), 2: Q(33, 3920)}, baseline=Q(8, 441), threshold=Q(4158, 145), masses={-2: Q(0), -1: Q(3, 77), 1: Q(3, 77), 2: Q(0)}),
    (5, 5, 6): dict(height=31, triples=13, rays=13, positive=12, maximum=Q(8, 119), winner=(5, 11, 17), bulks={1: Q(281, 14700), 2: Q(1, 147)}, baseline=Q(127, 7350), threshold=Q(214200, 6241), masses={-2: Q(1, 119), -1: Q(3, 119), 1: Q(3, 119), 2: Q(1, 119)}),
    (1, 2, 15): dict(height=37, triples=15, rays=15, positive=14, maximum=Q(60, 1001), winner=(1, 11, 13), bulks={1: Q(3, 245), 2: Q(3, 245)}, baseline=Q(4, 245), threshold=Q(15015, 382), masses={-2: Q(8, 1001), -1: Q(2, 91), 1: Q(2, 91), 2: Q(8, 1001)}),
    (1, 3, 14): dict(height=31, triples=13, rays=13, positive=9, maximum=Q(12, 175), winner=(1, 13, 25), bulks={1: Q(9, 686), 2: Q(9, 686)}, baseline=Q(6, 343), threshold=Q(2450, 73), masses={-2: Q(3, 175), -1: Q(3, 175), 1: Q(3, 175), 2: Q(3, 175)}),
    (1, 5, 12): dict(height=31, triples=15, rays=16, positive=12, maximum=Q(64, 931), winner=(7, 13, 19), bulks={1: Q(3, 196), 2: Q(23, 1960)}, baseline=Q(53, 2940), threshold=Q(95760, 2833), masses={-2: Q(11, 931), -1: Q(3, 133), 1: Q(3, 133), 2: Q(11, 931)}),
    (1, 6, 11): dict(height=35, triples=19, rays=20, positive=17, maximum=Q(6, 91), winner=(1, 7, 13), bulks={1: Q(107, 6468), 2: Q(23, 2156)}, baseline=Q(8, 441), threshold=Q(4914, 137), masses={-2: Q(0), -1: Q(3, 91), 1: Q(3, 91), 2: Q(0)}),
    (1, 8, 9): dict(height=37, triples=25, rays=25, positive=19, maximum=Q(118, 1925), winner=(7, 11, 25), bulks={1: Q(37, 2352), 2: Q(23, 2352)}, baseline=Q(5, 294), threshold=Q(138600, 3581), masses={-2: Q(26, 1925), -1: Q(3, 175), 1: Q(3, 175), 2: Q(26, 1925)}),
    (2, 3, 13): dict(height=29, triples=8, rays=8, positive=7, maximum=Q(12, 161), winner=(1, 11, 23), bulks={1: Q(9, 637), 2: Q(2, 147)}, baseline=Q(106, 5733), threshold=Q(113022, 3695), masses={-2: Q(3, 161), -1: Q(3, 161), 1: Q(3, 161), 2: Q(3, 161)}),
    (2, 7, 9): dict(height=37, triples=25, rays=25, positive=23, maximum=Q(10, 161), winner=(5, 13, 23), bulks={1: Q(17, 1029), 2: Q(10, 1029)}, baseline=Q(6, 343), threshold=Q(3381, 88), masses={-2: Q(2, 161), -1: Q(3, 161), 1: Q(3, 161), 2: Q(2, 161)}),
    (3, 4, 11): dict(height=31, triples=17, rays=18, positive=14, maximum=Q(46, 665), winner=(5, 7, 19), bulks={1: Q(215, 12936), 2: Q(19, 1617)}, baseline=Q(367, 19404), threshold=Q(3160080, 92647), masses={-2: Q(8, 665), -1: Q(3, 133), 1: Q(3, 133), 2: Q(8, 665)}),
    (3, 5, 10): dict(height=31, triples=16, rays=16, positive=12, maximum=Q(58, 833), winner=(5, 7, 17), bulks={1: Q(127, 7350), 2: Q(51, 4900)}, baseline=Q(407, 22050), threshold=Q(642600, 19181), masses={-2: Q(8, 833), -1: Q(3, 119), 1: Q(3, 119), 2: Q(8, 833)}),
    (3, 7, 8): dict(height=25, triples=14, rays=15, positive=13, maximum=Q(6, 77), winner=(1, 5, 11), bulks={1: Q(93, 5488), 2: Q(51, 5488)}, baseline=Q(6, 343), threshold=Q(539, 19), masses={-2: Q(3, 77), -1: Q(0), 1: Q(0), 2: Q(3, 77)}),
    (4, 5, 9): dict(height=37, triples=25, rays=25, positive=22, maximum=Q(10, 161), winner=(5, 13, 23), bulks={1: Q(311, 17640), 2: Q(1, 105)}, baseline=Q(479, 26460), threshold=Q(1043280, 26783), masses={-2: Q(2, 161), -1: Q(3, 161), 1: Q(3, 161), 2: Q(2, 161)}),
    (5, 6, 7): dict(height=31, triples=23, rays=25, positive=19, maximum=Q(64, 931), winner=(7, 13, 19), bulks={1: Q(53, 2940), 2: Q(169, 20580)}, baseline=Q(6, 343), threshold=Q(5586, 167), masses={-2: Q(11, 931), -1: Q(3, 133), 1: Q(3, 133), 2: Q(11, 931)}),
    (1, 1, 18): dict(height=65, triples=23, rays=23, positive=20, maximum=Q(2, 37), winner=(1, 19, 37), bulks={1: Q(1, 98), 2: Q(1, 98), 4: Q(1, 441)}, baseline=Q(20, 1323), threshold=Q(62937, 953), masses={-4: Q(1, 259), -2: Q(3, 259), -1: Q(3, 259), 1: Q(3, 259), 2: Q(3, 259), 4: Q(1, 259)}),
    (1, 3, 16): dict(height=59, triples=41, rays=42, positive=37, maximum=Q(12, 203), winner=(5, 17, 29), bulks={1: Q(9, 784), 2: Q(9, 784), 4: Q(1, 1176)}, baseline=Q(1, 63), threshold=Q(4698, 79), masses={-4: Q(0), -2: Q(3, 203), -1: Q(3, 203), 1: Q(3, 203), 2: Q(3, 203), 4: Q(0)}),
    (1, 4, 15): dict(height=47, triples=29, rays=29, positive=25, maximum=Q(58, 833), winner=(5, 7, 17), bulks={1: Q(3, 245), 2: Q(3, 245), 4: Q(1, 1470)}, baseline=Q(37, 2205), threshold=Q(13770, 283), masses={-4: Q(0), -2: Q(8, 833), -1: Q(3, 119), 1: Q(3, 119), 2: Q(8, 833), 4: Q(0)}),
    (1, 6, 13): dict(height=41, triples=24, rays=25, positive=21, maximum=Q(6, 77), winner=(1, 5, 11), bulks={1: Q(9, 637), 2: Q(29, 2548), 4: Q(1, 1911)}, baseline=Q(199, 11466), threshold=Q(324324, 7639), masses={-4: Q(0), -2: Q(0), -1: Q(3, 77), 1: Q(3, 77), 2: Q(0), 4: Q(0)}),
    (1, 7, 12): dict(height=49, triples=42, rays=42, positive=40, maximum=Q(8, 119), winner=(5, 11, 17), bulks={1: Q(125, 8232), 2: Q(29, 2744), 4: Q(1, 2058)}, baseline=Q(6, 343), threshold=Q(7497, 145), masses={-4: Q(0), -2: Q(1, 119), -1: Q(3, 119), 1: Q(3, 119), 2: Q(1, 119), 4: Q(0)}),
    (1, 9, 10): dict(height=59, triples=59, rays=59, positive=57, maximum=Q(12, 203), winner=(7, 11, 29), bulks={1: Q(43, 2940), 2: Q(29, 2940), 4: Q(1, 2205)}, baseline=Q(22, 1323), threshold=Q(49329, 815), masses={-4: Q(0), -2: Q(3, 203), -1: Q(3, 203), 1: Q(3, 203), 2: Q(3, 203), 4: Q(0)}),
    (3, 4, 13): dict(height=41, triples=27, rays=27, positive=21, maximum=Q(6, 77), winner=(1, 5, 11), bulks={1: Q(9, 637), 2: Q(191, 15288), 4: Q(1, 3822)}, baseline=Q(137, 7644), threshold=Q(216216, 5045), masses={-4: Q(0), -2: Q(3, 77), -1: Q(0), 1: Q(0), 2: Q(3, 77), 4: Q(0)}),
    (3, 7, 10): dict(height=53, triples=47, rays=47, positive=45, maximum=Q(6, 91), winner=(1, 7, 13), bulks={1: Q(47, 2940), 2: Q(69, 6860), 4: Q(1, 5145)}, baseline=Q(6, 343), threshold=Q(637, 12), masses={-4: Q(0), -2: Q(0), -1: Q(3, 91), 1: Q(3, 91), 2: Q(0), 4: Q(0)}),
    (4, 7, 9): dict(height=41, triples=37, rays=37, positive=31, maximum=Q(6, 77), winner=(1, 5, 11), bulks={1: Q(101, 6174), 2: Q(10, 1029), 4: Q(1, 6174)}, baseline=Q(6, 343), threshold=Q(1617, 38), masses={-4: Q(0), -2: Q(0), -1: Q(3, 77), 1: Q(3, 77), 2: Q(0), 4: Q(0)}),
    (6, 7, 7): dict(height=49, triples=32, rays=32, positive=31, maximum=Q(64, 931), winner=(7, 13, 19), bulks={1: Q(124, 7203), 2: Q(64, 7203), 4: Q(1, 7203)}, baseline=Q(6, 343), threshold=Q(8379, 167), masses={-4: Q(0), -2: Q(11, 931), -1: Q(3, 133), 1: Q(3, 133), 2: Q(11, 931), 4: Q(0)}),
}


EXPECTED_UNIT112_INTERSECTIONS = {
    (2, 3, 11): ((1, 5, 7), (1, 13, 25), (5, 7, 17), (7, 13, 19)),
    (2, 5, 9): ((1, 7, 13), (1, 11, 23), (5, 7, 19)),
    (3, 5, 8): ((1, 5, 11), (7, 13, 19)),
    (5, 5, 6): ((1, 5, 7), (5, 11, 17)),
    (1, 2, 15): ((1, 11, 13), (1, 17, 19), (5, 13, 31), (5, 17, 29)),
    (1, 3, 14): ((1, 13, 25), (5, 13, 31), (5, 17, 29), (7, 11, 29)),
    (1, 5, 12): ((7, 11, 25), (7, 11, 29), (7, 13, 19), (11, 17, 23)),
    (1, 6, 11): ((1, 7, 13), (5, 13, 23), (11, 17, 23)),
    (1, 8, 9): ((1, 17, 19), (7, 11, 25)),
    (2, 3, 13): ((1, 11, 23), (1, 17, 19), (5, 7, 17), (7, 11, 29)),
    (2, 7, 9): ((5, 13, 23), (7, 11, 25)),
    (3, 4, 11): ((1, 7, 13), (1, 17, 19), (5, 7, 19), (7, 11, 25)),
    (3, 5, 10): ((1, 13, 25), (5, 7, 17), (5, 13, 23), (7, 11, 25)),
    (3, 7, 8): ((1, 5, 11), (1, 11, 13), (1, 11, 23), (1, 17, 19), (5, 13, 23)),
    (4, 5, 9): ((1, 17, 19), (5, 13, 23)),
    (5, 6, 7): ((1, 17, 19), (5, 11, 17), (7, 13, 19)),
    (1, 1, 18): ((1, 17, 35), (1, 19, 37)),
    (1, 3, 16): ((1, 5, 7), (1, 17, 35), (5, 13, 31), (5, 17, 29), (7, 19, 31)),
    (1, 4, 15): ((1, 7, 13), (5, 7, 17), (5, 13, 23), (7, 11, 29), (7, 19, 31)),
    (1, 6, 13): ((1, 5, 11), (7, 11, 25), (13, 19, 25)),
    (1, 7, 12): ((5, 11, 17), (5, 13, 23), (5, 13, 31), (13, 19, 25)),
    (1, 9, 10): ((1, 17, 19), (7, 11, 29)),
    (3, 4, 13): ((1, 5, 7), (1, 5, 11), (5, 7, 19), (5, 17, 29), (11, 17, 23)),
    (3, 7, 10): ((1, 7, 13), (11, 17, 23)),
    (4, 7, 9): ((1, 5, 7), (1, 5, 11), (1, 11, 23), (1, 13, 25), (5, 11, 17)),
    (6, 7, 7): ((1, 5, 7), (7, 13, 19)),
}


EXPECTED_WINNER_GROUPS = (
    ((1, 5, 11), Q(6, 77), ((3, 5, 8), (3, 7, 8), (1, 6, 13), (3, 4, 13), (4, 7, 9))),
    ((1, 7, 13), Q(6, 91), ((1, 6, 11), (3, 7, 10))),
    ((1, 11, 13), Q(60, 1001), ((1, 2, 15),)),
    ((1, 11, 23), Q(12, 161), ((2, 5, 9), (2, 3, 13))),
    ((1, 13, 25), Q(12, 175), ((1, 3, 14),)),
    ((1, 19, 37), Q(2, 37), ((1, 1, 18),)),
    ((5, 7, 17), Q(58, 833), ((2, 3, 11), (3, 5, 10), (1, 4, 15))),
    ((5, 7, 19), Q(46, 665), ((3, 4, 11),)),
    ((5, 11, 17), Q(8, 119), ((5, 5, 6), (1, 7, 12))),
    ((5, 13, 23), Q(10, 161), ((2, 7, 9), (4, 5, 9))),
    ((5, 17, 29), Q(12, 203), ((1, 3, 16),)),
    ((7, 11, 25), Q(118, 1925), ((1, 8, 9),)),
    ((7, 11, 29), Q(12, 203), ((1, 9, 10),)),
    ((7, 13, 19), Q(64, 931), ((1, 5, 12), (5, 6, 7), (6, 7, 7))),
)

EXPECTED_PROOF_TOTALS = (636, 646, 1638, 370, 340)

EXPECTED_NN_HIST = (
    ((1, 2), 96), ((2, 2), 2388), ((2, 3), 11), ((3, 3), 303),
    ((3, 4), 4), ((4, 4), 40), ((5, 5), 6), ((5, 6), 1),
    ((6, 6), 2), ((6, 7), 1), ((7, 8), 1), ((8, 9), 1),
)

EXPECTED_NU_HIST = (
    ((1, 1, 1, 1), 2411), ((1, 1, 1, 2), 2), ((1, 1, 2, 2), 27),
    ((1, 1, 3, 3), 9), ((1, 1, 4, 4), 1), ((1, 2, 1, 1), 21),
    ((1, 2, 3, 3), 1), ((2, 2, 1, 1), 980), ((2, 2, 1, 2), 4),
    ((2, 2, 2, 2), 34), ((2, 2, 3, 3), 10), ((2, 2, 4, 4), 1),
    ((2, 3, 1, 1), 3), ((2, 3, 2, 2), 1), ((3, 3, 1, 1), 158),
    ((3, 3, 1, 2), 1), ((3, 3, 2, 2), 21), ((3, 3, 3, 3), 6),
    ((3, 4, 1, 1), 4), ((4, 4, 1, 1), 27), ((4, 4, 1, 2), 1),
    ((4, 4, 2, 2), 6), ((4, 4, 2, 3), 1), ((4, 4, 3, 3), 4),
    ((4, 4, 4, 4), 1), ((5, 5, 1, 1), 3), ((5, 5, 3, 3), 2),
    ((5, 5, 4, 4), 1), ((5, 6, 1, 1), 1), ((6, 6, 3, 3), 1),
    ((6, 6, 8, 8), 1), ((6, 7, 1, 1), 1), ((7, 8, 1, 1), 1),
    ((8, 9, 1, 1), 1),
)

EXPECTED_NN_DIGEST = "c82f9fe5aef6c9154ed204ba940a81333654e9b250e9d0a6d0f85b3e857ecf62"
EXPECTED_NU_DIGEST = "a45070cfa80b9fcbb62d3aa20d77fd7a1e75b9dd5d0dec86837a02ebfde62c75"
EXPECTED_ALL_ONE_ZERO_DIGEST = "ef58ff98a9d69afc70f31c3c12a0bb0890f73b32102bd451e6a85756ca572e11"
EXPECTED_EXTENSION_ONE_ZERO_DIGEST = "ab30482ff3b66194f41c513215cf2e93f7abaf00dec44f0b8f7473fd5849c537"
EXPECTED_ALL_ONE_ZERO_UNIT_DIGEST = "b98fd07fa93a140d6cb2566a1d7e2bf2558c64095e5ba56d8333be86f8a06231"
EXPECTED_EXTENSION_ONE_ZERO_UNIT_DIGEST = "44a7694c5ff4cc678ccc5c09583e3ae8d1f702685f78509d5bece1a31ac59cab"


def canonical_speed(c, d):
    w = cross(c, d)
    if not all(w) or not (all(x > 0 for x in w) or all(x < 0 for x in w)):
        return None
    if all(x < 0 for x in w):
        w = tuple(-x for x in w)
    g = gcd_all(w)
    w = tuple(sorted(x // g for x in w))
    if len(set(w)) != 3 or any(x % 2 == 0 or x % 3 == 0 for x in w):
        return None
    return w


def intersections_between(left_patterns, right_patterns, same_pool=False):
    left = tuple((p, c) for p in left_patterns for c in signed_rays(p))
    right = left if same_pool else tuple((p, c) for p in right_patterns for c in signed_rays(p))
    candidates = set()
    pairs = combinations(left, 2) if same_pool else product(left, right)
    for (_, c), (_, d) in pairs:
        w = canonical_speed(c, d)
        if w is not None:
            candidates.add(w)
    return tuple(sorted(candidates)), left, right


def incidence_records(candidates, left_patterns, right_patterns=None):
    left = tuple((p, signed_rays(p)) for p in left_patterns)
    right = tuple((p, signed_rays(p)) for p in right_patterns) if right_patterns else ()
    records = []
    for w in candidates:
        left_hits = tuple((p, tuple(c for c in cs if dot(c, w) == 0)) for p, cs in left)
        left_hits = tuple((p, cs) for p, cs in left_hits if cs)
        right_hits = tuple((p, tuple(c for c in cs if dot(c, w) == 0)) for p, cs in right)
        right_hits = tuple((p, cs) for p, cs in right_hits if cs)
        records.append((w, left_hits, right_hits))
    return tuple(records)


def digest(records):
    return sha256((repr(records) + "\n").encode("ascii")).hexdigest()


def main():
    if "--tripwire" in sys.argv:
        check(False, "optimization-live tripwire")

    all_one_zero = patterns_through(20, 1)
    unit_patterns = patterns_through(20, 0)
    check(all_one_zero == OLD_PATTERNS + PATTERNS, "complete one-zero atlas through 20")
    check(tuple(p for p in PATTERNS if sum(p) == 16) == PATTERNS[:4], "norm-16 partition")
    check(tuple(p for p in PATTERNS if sum(p) == 18) == PATTERNS[4:16], "norm-18 partition")
    check(tuple(p for p in PATTERNS if sum(p) == 20) == PATTERNS[16:], "norm-20 partition")
    check(len(unit_patterns) == 33, "33 unit magnitude patterns through 20")

    result_rows = []
    proof_union = set()
    proof_total_components = 0
    unit_winner_rows = []
    winner_groups = defaultdict(list)
    winner_profiles = defaultdict(list)

    for pattern in PATTERNS:
        expected = EXPECTED[pattern]
        height = expected["height"]
        solved = presentation_universe(pattern, height)
        brute = brute_universe(pattern, height)
        check(solved == brute, "independent universe routes " + repr(pattern))
        check(len(solved) == expected["triples"], "triple count " + repr(pattern))
        ray_count = sum(len(cs) for cs in solved.values())
        check(ray_count == expected["rays"], "relation-chart count " + repr(pattern))

        defects = relation_defects(signed_rays(pattern)[0])
        expected_defects = (-2, -1, 1, 2) if sum(pattern) < 20 else (-4, -2, -1, 1, 2, 4)
        check(defects == expected_defects, "strict defect states " + repr(pattern))
        bulks = {}
        for delta in sorted(set(abs(x) for x in defects)):
            convolution = bulk_convolution(pattern, delta)
            polygon = bulk_polygon(pattern, delta)
            check(convolution == polygon, "two-route cube slice " + repr((pattern, delta)))
            check(convolution == expected["bulks"][delta], "cube slice value " + repr((pattern, delta)))
            bulks[delta] = convolution
        baseline = sum(bulk_convolution(pattern, delta) for delta in defects) / 3
        check(baseline == expected["baseline"], "bulk baseline " + repr(pattern))
        error_numerator = Q(3 * len(defects), 7)
        threshold = error_numerator / (expected["maximum"] - baseline)
        check(threshold == expected["threshold"], "tail threshold " + repr(pattern))
        admissible_below = tuple(x for x in range(1, floor_q(threshold) + 1) if x % 2 and x % 3)
        check(height == max(admissible_below), "smallest admissible proof cutoff " + repr(pattern))
        next_admissible = next(x for x in range(height + 1, height + 7) if x % 2 and x % 3)
        check(next_admissible > threshold, "next admissible height crosses threshold " + repr(pattern))
        check(baseline + error_numerator / height >= expected["maximum"], "cutoff minimality lower side " + repr(pattern))
        check(baseline + error_numerator / next_admissible < expected["maximum"], "strict admissible tail exclusion " + repr(pattern))

        measures = {}
        positive = 0
        components_here = 0
        for w, cs in solved.items():
            proof_union.add(w)
            direct, pieces = physical_dict(w)
            check(pieces == len(direct), "physical interval/carrier bijection " + repr(w))
            value = sum(direct.values(), Q(0))
            measures[w] = value
            positive += value > 0
            components_here += pieces
            for c in cs:
                quotient, masses = quotient_components(w, c)
                check(quotient == direct, "quotient equals physical comb " + repr((pattern, w, c)))
                check(sum(masses.values(), Q(0)) == value, "defect total " + repr((pattern, w, c)))

        maximum = max(measures.values())
        winners = tuple(w for w in sorted(measures) if measures[w] == maximum)
        check(positive == expected["positive"], "positive count " + repr(pattern))
        check(maximum == expected["maximum"], "sharp mass " + repr(pattern))
        check(winners == (expected["winner"],), "unique unordered winner " + repr(pattern))
        check(any(candidate[-1] == height for candidate in solved), "cutoff height actually occupied " + repr(pattern))
        w = winners[0]
        profiles = tuple((c, quotient_components(w, c)[1]) for c in solved[w])
        check(all(masses == expected["masses"] for _, masses in profiles), "winning defect profiles " + repr(pattern))
        winner_profiles[w].append((pattern, profiles))
        winner_groups[(w, maximum)].append(pattern)

        unit_hits = tuple(
            (p, c, quotient_components(w, c)[1])
            for p in unit_patterns for c in signed_rays(p) if dot(c, w) == 0
        )
        check(unit_hits and unit_hits[0][0] == (1, 1, 2), "winner lies in 112 unit sector " + repr(pattern))
        for p, c, masses in unit_hits:
            quotient, again = quotient_components(w, c)
            direct, _ = physical_dict(w)
            check(quotient == direct and masses == again, "unit chart equals physical comb " + repr((p, w, c)))
        unit_winner_rows.append((pattern, w, maximum, unit_hits))

        proof_total_components += components_here
        result_rows.append((pattern, sum(pattern), defects, height, len(solved), ray_count,
                            positive, components_here, bulks, baseline, threshold,
                            maximum, winners, profiles))

    groups = tuple((w, mass, tuple(patterns)) for (w, mass), patterns in sorted(winner_groups.items()))
    check(groups == EXPECTED_WINNER_GROUPS, "all physical equality-comb groups")
    overall = max(row[11] for row in result_rows)
    overall_rows = tuple(row[0] for row in result_rows if row[11] == overall)
    check(overall == Q(6, 77), "global new-row sharp value")
    check(overall_rows == ((3, 5, 8), (3, 7, 8), (1, 6, 13), (3, 4, 13), (4, 7, 9)), "five global equality presentations")
    check({EXPECTED[p]["winner"] for p in overall_rows} == {(1, 5, 11)}, "unique physical global equality comb")

    new_reallocations = []
    for w, rows in sorted(winner_profiles.items()):
        normalized = {
            tuple((delta, masses.get(delta, Q(0))) for delta in (-4, -2, -1, 1, 2, 4))
            for _, profiles in rows for _, masses in profiles
        }
        if len(normalized) > 1:
            new_reallocations.append((w, tuple(sorted(normalized))))
    check(tuple(w for w, _ in new_reallocations) == ((1, 5, 11),),
          "only repeated new winner with chart-dependent defect split")

    # Height-free cross-product compression against the smallest unit sector.
    unit112 = ((1, 1, 2),)
    compression_rows = []
    for pattern in PATTERNS:
        candidates, _, _ = intersections_between((pattern,), unit112)
        check(candidates == EXPECTED_UNIT112_INTERSECTIONS[pattern], "complete 112 intersection list " + repr(pattern))
        check(EXPECTED[pattern]["winner"] in candidates, "winner on 112 intersection list " + repr(pattern))
        for w in candidates:
            direct, _ = physical_dict(w)
            for p in (pattern, (1, 1, 2)):
                hits = tuple(c for c in signed_rays(p) if dot(c, w) == 0)
                check(hits, "cross-product relation survives canonical sorting")
                for c in hits:
                    quotient, _ = quotient_components(w, c)
                    check(quotient == direct, "112 compression chart/physical equality")
        compression_rows.append((pattern, candidates, EXPECTED[pattern]["winner"]))

    # Reconstruct every one-zero chart of the global equality comb, including
    # the old <=14 rows, to expose the chart-dependent defect reallocation.
    equality_w = (1, 5, 11)
    equality_direct, _ = physical_dict(equality_w)
    equality_charts = []
    for pattern in all_one_zero:
        hits = tuple(c for c in signed_rays(pattern) if dot(c, equality_w) == 0)
        if not hits:
            continue
        rows = []
        for c in hits:
            quotient, masses = quotient_components(equality_w, c)
            check(quotient == equality_direct, "global equality one-zero chart/physical equality")
            rows.append((c, masses))
        equality_charts.append((pattern, tuple(rows)))
    check(tuple(p for p, _ in equality_charts) == (
        (1, 1, 6), (1, 3, 4), (2, 3, 5), (2, 3, 7), (1, 4, 9),
        (3, 5, 8), (3, 7, 8), (1, 6, 13), (3, 4, 13), (4, 7, 9),
    ), "ten one-zero equality presentation shapes")
    normalized_profiles = {
        tuple((delta, masses.get(delta, Q(0))) for delta in (-4, -2, -1, 1, 2, 4))
        for _, rows in equality_charts for _, masses in rows
    }
    check(normalized_profiles == {
        ((-4, Q(0)), (-2, Q(0)), (-1, Q(3, 77)), (1, Q(3, 77)), (2, Q(0)), (4, Q(0))),
        ((-4, Q(0)), (-2, Q(3, 77)), (-1, Q(0)), (1, Q(0)), (2, Q(3, 77)), (4, Q(0))),
    }, "exact two-profile nonunit defect reallocation")

    equality_unit_charts = []
    for pattern in unit_patterns:
        for c in signed_rays(pattern):
            if dot(c, equality_w):
                continue
            quotient, masses = quotient_components(equality_w, c)
            check(quotient == equality_direct, "global equality unit chart/physical equality")
            equality_unit_charts.append((pattern, c, masses))
    check(tuple((p, c) for p, c, _ in equality_unit_charts) == (
        ((1, 1, 2), (1, 2, -1)),
        ((1, 1, 16), (16, -1, -1)),
        ((1, 2, 17), (17, 1, -2)),
        ((1, 5, 14), (14, -5, 1)),
    ), "all unit charts through 20 at global equality")

    # Complete, height-free presentation-overlap atlases.  Completeness is by
    # cross product: two independent labeled relations determine the primitive
    # positive kernel generator, and all rays in the finite shape lists occur.
    nn_candidates, new_rays, _ = intersections_between(PATTERNS, PATTERNS, same_pool=True)
    nn_records = incidence_records(nn_candidates, PATTERNS)
    nn_records = tuple(row for row in nn_records if sum(len(cs) for _, cs in row[1]) >= 2)
    check(len(new_rays) == 441, "441 new signed relation rays")
    check(len(nn_records) == 2854, "complete new/new physical overlap count")
    nn_cross = sum(len(hits) >= 2 for _, hits, _ in nn_records)
    nn_within_only = sum(len(hits) == 1 for _, hits, _ in nn_records)
    check((nn_cross, nn_within_only) == (2758, 96), "cross-pattern/within-pattern overlap split")

    nu_candidates, _, unit_rays = intersections_between(PATTERNS, unit_patterns)
    nu_records = incidence_records(nu_candidates, PATTERNS, unit_patterns)
    check(len(unit_rays) == 495, "495 unit signed relation rays through 20")
    check(len(nu_records) == 3747, "complete new/unit physical overlap count")

    nn_hist = tuple(sorted(Counter((len(left), sum(len(cs) for _, cs in left)) for _, left, _ in nn_records).items()))
    nu_hist = tuple(sorted(Counter((len(left), sum(len(cs) for _, cs in left), len(right), sum(len(cs) for _, cs in right)) for _, left, right in nu_records).items()))
    check(nn_hist == EXPECTED_NN_HIST, "new/new overlap multiplicity histogram")
    check(nu_hist == EXPECTED_NU_HIST, "new/unit overlap multiplicity histogram")
    check(digest(nn_records) == EXPECTED_NN_DIGEST, "new/new overlap atlas digest")
    check(digest(nu_records) == EXPECTED_NU_DIGEST, "new/unit overlap atlas digest")
    overlap_summary = (
        len(nn_records), nn_cross, nn_within_only, max(w[-1] for w, _, _ in nn_records),
        digest(nn_records), nn_hist,
        len(nu_records), max(w[-1] for w, _, _ in nu_records), digest(nu_records), nu_hist,
    )

    # Full <=20 presentation atlas, so old/new overlaps are not silently lost.
    all_candidates, all_rays, _ = intersections_between(all_one_zero, all_one_zero, same_pool=True)
    all_records = incidence_records(all_candidates, all_one_zero)
    all_records = tuple(row for row in all_records if sum(len(cs) for _, cs in row[1]) >= 2)
    new_pattern_set = set(PATTERNS)
    extension_records = tuple(
        row for row in all_records if any(pattern in new_pattern_set for pattern, _ in row[1])
    )
    check(len(all_rays) == 675, "675 one-zero signed rays through 20")
    check((len(all_records), len(extension_records)) == (4374, 4218), "full/extension one-zero overlaps")
    check((sum(len(hits) >= 2 for _, hits, _ in all_records),
           sum(len(hits) == 1 for _, hits, _ in all_records)) == (4263, 111),
          "full cross-shape/within-shape overlap split")
    check(digest(all_records) == EXPECTED_ALL_ONE_ZERO_DIGEST, "full one-zero overlap digest")
    check(digest(extension_records) == EXPECTED_EXTENSION_ONE_ZERO_DIGEST, "extension one-zero overlap digest")

    all_unit_candidates, _, _ = intersections_between(all_one_zero, unit_patterns)
    all_unit_records = incidence_records(all_unit_candidates, all_one_zero, unit_patterns)
    extension_unit_records = tuple(
        row for row in all_unit_records if any(pattern in new_pattern_set for pattern, _ in row[1])
    )
    check((len(all_unit_records), len(extension_unit_records)) == (4247, 3747),
          "full/extension one-zero-unit overlaps")
    check(digest(all_unit_records) == EXPECTED_ALL_ONE_ZERO_UNIT_DIGEST, "full one-zero/unit overlap digest")
    check(digest(extension_unit_records) == EXPECTED_EXTENSION_ONE_ZERO_UNIT_DIGEST,
          "extension one-zero/unit overlap digest")
    all_hist = tuple(sorted(Counter((len(left), sum(len(cs) for _, cs in left)) for _, left, _ in all_records).items()))
    all_unit_hist = tuple(sorted(Counter((len(left), sum(len(cs) for _, cs in left), len(right), sum(len(cs) for _, cs in right)) for _, left, right in all_unit_records).items()))
    through20_overlap_summary = (
        len(all_records), len(extension_records), 4263, 111,
        max(w[-1] for w, _, _ in all_records), digest(all_records), digest(extension_records), all_hist,
        len(all_unit_records), len(extension_unit_records),
        max(w[-1] for w, _, _ in all_unit_records),
        digest(all_unit_records), digest(extension_unit_records), all_unit_hist,
    )

    proof_positive_union = sum(sum(physical_dict(w)[0].values(), Q(0)) > 0 for w in proof_union)
    proof_totals = (
        sum(row[4] for row in result_rows), sum(row[5] for row in result_rows),
        proof_total_components, len(proof_union), proof_positive_union,
    )
    check(proof_totals == EXPECTED_PROOF_TOTALS, "exact proof-universe totals")

    print("status=PASS")
    print("scope=26_additional_one_zero_mod3_presentation_shapes_norms_16_18_20; no_minimality_filter; LRC(14)=OPEN")
    print("patterns=" + repr(PATTERNS))
    print("unit_patterns_through20=" + repr(unit_patterns))
    print("result_rows=" + repr(tuple(result_rows)))
    print("proof_totals=" + repr(proof_totals))
    print("physical_equality_groups=" + repr(groups))
    print("global_sharp=" + repr((overall, overall_rows, equality_w)))
    print("new_winner_defect_reallocations=" + repr(tuple(new_reallocations)))
    print("winner_unit_charts=" + repr(tuple(unit_winner_rows)))
    print("unit112_compression=" + repr(tuple(compression_rows)))
    print("global_equality_one_zero_charts=" + repr(tuple(equality_charts)))
    print("global_equality_unit_charts=" + repr(tuple(equality_unit_charts)))
    print("complete_overlap_summary=" + repr(overlap_summary))
    print("complete_through20_overlap_summary=" + repr(through20_overlap_summary))
    print("physical_cache=" + repr(physical_components.cache_info()))
    print("optimization_safe_checks=yes")
    print("explicit_checks=" + str(CHECKS))


if __name__ == "__main__":
    main()
