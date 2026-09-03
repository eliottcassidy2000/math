#!/usr/bin/env python3
"""Clean-room exact verifier for the minimal ternary-unit l1=20 LRC shell.

The program imports no repository or producer implementation.  All theorem
checks use integers/Fractions and explicit exceptions that remain live under
python -O.  There is no floating-point theorem evidence.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd
import json


R = Fraction(3, 14)
CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


def gcd_many(values):
    g = 0
    for value in values:
        g = gcd(g, abs(value))
    return g


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def primitive_direction(v):
    g = gcd_many(v)
    check(g > 0, f"nonzero primitive direction {v}")
    q = tuple(x // g for x in v)
    first = next(x for x in q if x)
    return tuple(-x for x in q) if first < 0 else q


def extended_gcd(a, b):
    old_r, r = abs(a), abs(b)
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    if a < 0:
        old_s = -old_s
    if b < 0:
        old_t = -old_t
    return old_r, old_s, old_t


def bezout_vector(c):
    g12, x, y = extended_gcd(c[0], c[1])
    g, s, t = extended_gcd(g12, c[2])
    check(g == 1, f"primitive coefficient vector {c}")
    v = (s * x, s * y, t)
    check(dot(c, v) == 1, f"Bezout identity {c}, {v}")
    return v


def shell_patterns(limit, exact=False):
    """Sorted primitive ternary-unit magnitude patterns of allowable parity."""
    rows = []
    for a in range(1, limit + 1):
        for b in range(a, limit + 1):
            for c in range(b, limit + 1):
                total = a + b + c
                if total > limit:
                    break
                if exact and total != limit:
                    continue
                if total % 2:
                    continue
                if any(x % 3 == 0 for x in (a, b, c)):
                    continue
                if gcd_many((a, b, c)) != 1:
                    continue
                rows.append((a, b, c))
    return tuple(rows)


@lru_cache(maxsize=None)
def relation_vectors(pattern):
    """All signed coordinate placements, modulo simultaneous negation."""
    rows = set()
    for order in set(permutations(pattern)):
        for signs in product((-1, 1), repeat=3):
            c = tuple(signs[i] * order[i] for i in range(3))
            if all(x > 0 for x in c) or all(x < 0 for x in c):
                continue
            rows.add(primitive_direction(c))
    return tuple(sorted(rows))


def admissible_speeds(height):
    return tuple(x for x in range(1, height + 1, 2) if x % 3)


def primitive_triples(height):
    for w in combinations(admissible_speeds(height), 3):
        if gcd_many(w) == 1:
            yield w


@lru_cache(maxsize=None)
def generated_relation_map(pattern, height):
    """All sorted primitive speed triples on a signed coefficient ray."""
    values = admissible_speeds(height)
    value_set = set(values)
    rows = {}
    for c in relation_vectors(pattern):
        for w1 in values:
            for w2 in values:
                if w2 <= w1:
                    continue
                numerator = -(c[0] * w1 + c[1] * w2)
                if numerator % c[2]:
                    continue
                w3 = numerator // c[2]
                if w3 <= w2 or w3 > height or w3 not in value_set:
                    continue
                w = (w1, w2, w3)
                if gcd_many(w) != 1:
                    continue
                check(dot(c, w) == 0, f"generated relation {pattern}, {w}, {c}")
                rows.setdefault(w, set()).add(c)
    return rows


def direct_relations(w, patterns):
    rows = {}
    for pattern in patterns:
        relations = tuple(c for c in relation_vectors(pattern) if dot(c, w) == 0)
        if relations:
            rows[pattern] = relations
    return rows


def component_length(w, C):
    return max(Fraction(0), min(
        2 * R / w[0],
        2 * R / w[1],
        2 * R / w[2],
        R / w[0] + R / w[1] - Fraction(abs(C[2]), w[0] * w[1]),
        R / w[0] + R / w[2] - Fraction(abs(C[1]), w[0] * w[2]),
        R / w[1] + R / w[2] - Fraction(abs(C[0]), w[1] * w[2]),
    ))


def floor_fraction(value):
    return value.numerator // value.denominator


def ceil_fraction(value):
    return -floor_fraction(-value)


def real_line_support(w, c, C0):
    lower = []
    upper = []
    for i in range(3):
        j, k = tuple(index for index in range(3) if index != i)
        bound = R * (w[j] + w[k])
        left = Fraction(-bound - C0[i], c[i])
        right = Fraction(bound - C0[i], c[i])
        lower.append(min(left, right))
        upper.append(max(left, right))
    return max(lower), min(upper)


def integer_line_support(w, c, C0):
    left, right = real_line_support(w, c, C0)
    return floor_fraction(left) + 1, ceil_fraction(right) - 1


def raw_relation_components(w, c):
    """Exact affine defect-fibre carriers with (delta,k) chart metadata."""
    v = bezout_vector(c)
    components = {}
    metadata = {}
    for delta in (-3, 0, 3):
        n0 = tuple(delta * x for x in v)
        C0 = cross(w, n0)
        check(dot(C0, w) == 0, f"base carrier lies in Lambda {w}, {c}, {delta}")
        check(all(x % 3 == 0 for x in C0), f"base carrier is zero mod three {w}, {c}, {delta}")
        lo, hi = integer_line_support(w, c, C0)
        for k in range(lo, hi + 1):
            C = tuple(C0[i] + k * c[i] for i in range(3))
            check(dot(C, w) == 0, f"affine carrier relation {w}, {c}, {C}")
            gate = all(x % 3 for x in C)
            check(gate == (k % 3 != 0), f"two residue classes on fibre {w}, {c}, {delta}, {k}")
            if not gate:
                continue
            length = component_length(w, C)
            if not length:
                continue
            check(C not in components, f"unique defect address {w}, {c}, {C}")
            components[C] = length
            metadata[C] = (delta, k)
    return components, metadata


def literal_physical_components(w):
    """Definition-level nearest-integer interval sweep on y in [0,1]."""
    interval_lists = []
    for speed in w:
        inverse = pow(speed, -1, 3)
        intervals = []
        for nearest in range(speed + 1):
            left = max(Fraction(0), (Fraction(nearest) - R) / speed)
            right = min(Fraction(1), (Fraction(nearest) + R) / speed)
            if left < right:
                intervals.append((left, right, nearest, (-inverse * nearest) % 3))
        interval_lists.append(intervals)

    indices = [0, 0, 0]
    components = {}
    representatives = {}
    while all(indices[i] < len(interval_lists[i]) for i in range(3)):
        current = tuple(interval_lists[i][indices[i]] for i in range(3))
        left = max(row[0] for row in current)
        right = min(row[1] for row in current)
        if left < right and len({row[3] for row in current}) == 3:
            n = tuple(row[2] for row in current)
            owners = tuple(row[3] for row in current)
            C = cross(w, n)
            for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
                expected = w[j] * w[k] * (owners[j] - owners[k])
                check((C[i] - expected) % 3 == 0,
                      f"literal owner/carrier gate {w}, {C}, {i}")
            check(all(x % 3 for x in C), f"literal full owner gate {w}, {C}")
            if C in representatives:
                difference = tuple(n[i] - representatives[C][i] for i in range(3))
                check(cross(w, difference) == (0, 0, 0),
                      f"literal repeated lift differs by w {w}, {C}")
            else:
                representatives[C] = n
            components[C] = components.get(C, Fraction(0)) + right - left

        earliest_right = min(row[1] for row in current)
        for i in range(3):
            if current[i][1] == earliest_right:
                indices[i] += 1
    return components, representatives


def slice_integral(pattern, delta):
    """Integral of the affine roof from the three-box convolution formula."""
    target = Fraction(abs(delta)) + R * sum(pattern)
    alternating_sum = Fraction(0)
    for mask in range(8):
        removed = 2 * R * sum(pattern[i] for i in range(3) if mask & (1 << i))
        positive = max(Fraction(0), target - removed)
        alternating_sum += (-1) ** mask.bit_count() * positive * positive
    return alternating_sum / (2 * pattern[0] * pattern[1] * pattern[2])


def affine_roof_integral(w, c, delta):
    """Independent exact integration of the piecewise-affine lower envelope."""
    v = bezout_vector(c)
    C0 = cross(w, tuple(delta * x for x in v))
    support_left, support_right = real_line_support(w, c, C0)
    check(support_left < support_right, f"nonempty real slice {w}, {c}, {delta}")
    coarse = {support_left, support_right}
    for a, b in zip(C0, c):
        kink = Fraction(-a, b)
        if support_left < kink < support_right:
            coarse.add(kink)
    coarse = sorted(coarse)
    pair_data = (
        (0, R / w[1] + R / w[2], w[1] * w[2]),
        (1, R / w[0] + R / w[2], w[0] * w[2]),
        (2, R / w[0] + R / w[1], w[0] * w[1]),
    )
    total = Fraction(0)
    for coarse_left, coarse_right in zip(coarse, coarse[1:]):
        midpoint = (coarse_left + coarse_right) / 2
        lines = [(2 * R / speed, Fraction(0)) for speed in w]
        for index, base, denominator in pair_data:
            sign = 1 if C0[index] + midpoint * c[index] >= 0 else -1
            alpha = base - Fraction(sign * C0[index], denominator)
            beta = -Fraction(sign * c[index], denominator)
            lines.append((alpha, beta))
        refined = {coarse_left, coarse_right}
        for (a1, b1), (a2, b2) in combinations(lines, 2):
            if b1 == b2:
                continue
            crossing = Fraction(a2 - a1, b1 - b2)
            if coarse_left < crossing < coarse_right:
                refined.add(crossing)
        refined = sorted(refined)
        for left, right in zip(refined, refined[1:]):
            mid = (left + right) / 2
            alpha, beta = min(lines, key=lambda line: line[0] + line[1] * mid)
            check(alpha + beta * mid > 0,
                  f"positive lower envelope {w}, {c}, {delta}, {left}, {right}")
            total += alpha * (right - left) + beta * (right * right - left * left) / 2
    return total


SHORT_EXPECTED = (
    (1, 1, 2), (1, 1, 4), (1, 1, 8), (1, 1, 10), (1, 1, 14), (1, 1, 16),
    (1, 2, 5), (1, 2, 7), (1, 2, 11), (1, 2, 13),
    (1, 4, 5), (1, 4, 7), (1, 4, 11), (1, 4, 13),
    (1, 5, 8), (1, 5, 10), (1, 7, 8), (1, 7, 10),
    (2, 5, 5), (2, 5, 7), (2, 5, 11), (2, 7, 7),
    (4, 5, 5), (4, 5, 7), (4, 7, 7), (5, 5, 8),
)

SHELL20 = (
    (1, 2, 17), (1, 5, 14), (1, 8, 11), (2, 5, 13),
    (2, 7, 11), (4, 5, 11), (5, 7, 8),
)

CLAIMS = {
    (1, 2, 17): {
        "integrals": (Fraction(9, 833), Fraction(9, 833)),
        "bulk": Fraction(18, 833), "threshold": Fraction(22015, 53), "height": 415,
        "all": 1029, "minimal": 1012, "rays": 1013, "positive": 1012,
        "components": 21652, "zeros": (), "maximum": Fraction(36, 1295),
        "maximizers": ((1, 11, 185),),
        "relation": (2, -17, 1), "delta_counts": (4, 4, 4),
        "masses": (Fraction(12, 1295), Fraction(12, 1295), Fraction(12, 1295)),
    },
    (1, 5, 14): {
        "integrals": (Fraction(9, 686), Fraction(9, 1372)),
        "bulk": Fraction(6, 343), "threshold": Fraction(21903, 47), "height": 466,
        "all": 1567, "minimal": 1549, "rays": 1550, "positive": 1549,
        "components": 37274, "zeros": (), "maximum": Fraction(24, 1043),
        "maximizers": ((1, 11, 149),),
        "relation": (5, -14, 1), "delta_counts": (2, 4, 2),
        "masses": (Fraction(6, 1043), Fraction(12, 1043), Fraction(6, 1043)),
    },
    (1, 8, 11): {
        "integrals": (Fraction(9, 539), Fraction(45, 8624)),
        "bulk": Fraction(39, 2156), "threshold": Fraction(1158696, 3385), "height": 342,
        "all": 1071, "minimal": 1056, "rays": 1057, "positive": 1056,
        "components": 21946, "zeros": (), "maximum": Fraction(412, 16093),
        "maximizers": ((11, 19, 121),),
        "relation": (8, -11, 1), "delta_counts": (2, 4, 2),
        "masses": (Fraction(92, 16093), Fraction(12, 847), Fraction(92, 16093)),
    },
    (2, 5, 13): {
        "integrals": (Fraction(9, 637), Fraction(18, 3185)),
        "bulk": Fraction(54, 3185), "threshold": Fraction(4009005, 11332), "height": 353,
        "all": 977, "minimal": 964, "rays": 965, "positive": 963,
        "components": 18460, "zeros": ((1, 19, 41),), "maximum": Fraction(166, 6853),
        "maximizers": ((7, 11, 89),),
        "relation": (5, 13, -2), "delta_counts": (2, 2, 2),
        "masses": (Fraction(50, 6853), Fraction(6, 623), Fraction(50, 6853)),
    },
    (2, 7, 11): {
        "integrals": (Fraction(9, 539), Fraction(18, 3773)),
        "bulk": Fraction(6, 343), "threshold": Fraction(1547469, 4369), "height": 354,
        "all": 1160, "minimal": 1147, "rays": 1147, "positive": 1147,
        "components": 24512, "zeros": (), "maximum": Fraction(608, 24563),
        "maximizers": ((11, 29, 121),),
        "relation": (7, -11, 2), "delta_counts": (2, 4, 2),
        "masses": (Fraction(130, 24563), Fraction(12, 847), Fraction(130, 24563)),
    },
    (4, 5, 11): {
        "integrals": (Fraction(9, 539), Fraction(81, 21560)),
        "bulk": Fraction(87, 5390), "threshold": Fraction(118580, 249), "height": 476,
        "all": 2092, "minimal": 2079, "rays": 2080, "positive": 2078,
        "components": 55166, "zeros": ((1, 19, 41),), "maximum": Fraction(894, 41503),
        "maximizers": ((11, 49, 121),),
        "relation": (5, -11, 4), "delta_counts": (2, 4, 2),
        "masses": (Fraction(153, 41503), Fraction(12, 847), Fraction(153, 41503)),
    },
    (5, 7, 8): {
        "integrals": (Fraction(279, 13720), Fraction(81, 27440)),
        "bulk": Fraction(6, 343), "threshold": Fraction(185955, 499), "height": 372,
        "all": 1549, "minimal": 1538, "rays": 1538, "positive": 1538,
        "components": 32424, "zeros": (), "maximum": Fraction(216, 8855),
        "maximizers": ((13, 23, 55),),
        "relation": (7, 8, -5), "delta_counts": (1, 2, 1),
        "masses": (Fraction(8, 1771), Fraction(136, 8855), Fraction(8, 1771)),
    },
}

COMMON_COUNTS = {
    (1, 2, 17): (1349, 1332, 1333),
    (1, 5, 14): (1645, 1627, 1628),
    (1, 8, 11): (2090, 2075, 2076),
    (2, 5, 13): (1769, 1756, 1757),
    (2, 7, 11): (2095, 2082, 2082),
    (4, 5, 11): (2092, 2079, 2080),
    (5, 7, 8): (2552, 2541, 2541),
}

MULTIPLE_RELATIONS = {
    (1, 19, 41): (((2, 5, 13), (13, -5, 2)), ((4, 5, 11), (4, -11, 5))),
    (5, 13, 49): (((1, 2, 17), (17, 1, -2)), ((1, 8, 11), (11, -8, 1))),
    (7, 11, 65): (((1, 2, 17), (17, 1, -2)), ((1, 8, 11), (8, -11, 1))),
    (7, 11, 97): (((1, 2, 17), (1, 17, -2)), ((1, 2, 17), (17, -2, -1))),
    (11, 29, 41): (((2, 7, 11), (11, -7, 2)), ((5, 7, 8), (5, 8, -7))),
    (11, 29, 79): (((1, 2, 17), (17, -1, -2)), ((1, 5, 14), (1, -14, 5))),
    (13, 23, 41): (((1, 5, 14), (14, 1, -5)), ((1, 8, 11), (11, -8, 1))),
    (13, 23, 47): (((4, 5, 11), (5, -11, 4)), ((4, 5, 11), (11, 4, -5))),
    (13, 23, 67): (((1, 5, 14), (1, 14, -5)), ((1, 5, 14), (14, -5, -1))),
    (13, 31, 37): (((2, 7, 11), (11, -7, 2)), ((5, 7, 8), (8, 5, -7))),
    (17, 23, 53): (((2, 5, 13), (2, -13, 5)), ((2, 5, 13), (13, -5, -2))),
    (17, 37, 53): (((1, 5, 14), (14, -5, -1)), ((1, 8, 11), (1, 11, -8))),
    (25, 29, 43): (((1, 8, 11), (1, 11, -8)), ((1, 8, 11), (11, -8, -1))),
}

OVERLAP_PROFILES = {
    (1, 19, 41): (Fraction(0), {}),
    (5, 13, 49): (Fraction(32, 4459), {(-3, 0): (1, Fraction(16, 4459)), (3, 0): (1, Fraction(16, 4459))}),
    (7, 11, 65): (Fraction(6, 455), {(-3, 0): (1, Fraction(3, 455)), (3, 0): (1, Fraction(3, 455))}),
    (7, 11, 97): (Fraction(18, 679), {
        (-3, -3): (1, Fraction(3, 679)), (-3, 0): (1, Fraction(3, 679)),
        (0, -3): (1, Fraction(3, 679)), (0, 3): (1, Fraction(3, 679)),
        (3, 0): (1, Fraction(3, 679)), (3, 3): (1, Fraction(3, 679)),
    }),
    (11, 29, 41): (Fraction(138, 8323), {
        (-3, 0): (1, Fraction(1, 203)), (0, -3): (1, Fraction(4, 1189)),
        (0, 3): (1, Fraction(4, 1189)), (3, 0): (1, Fraction(1, 203)),
    }),
    (11, 29, 79): (Fraction(360, 16037), {
        (-3, -3): (1, Fraction(50, 16037)), (-3, 0): (1, Fraction(3, 553)),
        (0, -3): (1, Fraction(43, 16037)), (0, 3): (1, Fraction(43, 16037)),
        (3, 0): (1, Fraction(3, 553)), (3, 3): (1, Fraction(50, 16037)),
    }),
    (13, 23, 41): (Fraction(38, 6601), {(-3, 0): (1, Fraction(19, 6601)), (3, 0): (1, Fraction(19, 6601))}),
    (13, 23, 47): (Fraction(102, 7567), {
        (-3, 0): (1, Fraction(4, 1081)), (0, -3): (1, Fraction(1, 329)),
        (0, 3): (1, Fraction(1, 329)), (3, 0): (1, Fraction(4, 1081)),
    }),
    (13, 23, 67): (Fraction(282, 20033), {
        (-3, 0): (1, Fraction(37, 10787)), (0, -3): (1, Fraction(22, 6097)),
        (0, 3): (1, Fraction(22, 6097)), (3, 0): (1, Fraction(37, 10787)),
    }),
    (13, 31, 37): (Fraction(142, 8029), {
        (-3, 0): (1, Fraction(46, 8029)), (0, -3): (1, Fraction(25, 8029)),
        (0, 3): (1, Fraction(25, 8029)), (3, 0): (1, Fraction(46, 8029)),
    }),
    (17, 23, 53): (Fraction(2592, 145061), {
        (-3, -3): (1, Fraction(11, 2737)), (-3, 0): (1, Fraction(1, 371)),
        (0, -3): (1, Fraction(2, 901)), (0, 3): (1, Fraction(2, 901)),
        (3, 0): (1, Fraction(1, 371)), (3, 3): (1, Fraction(11, 2737)),
    }),
    (17, 37, 53): (Fraction(90, 6307), {
        (-3, 0): (1, Fraction(4, 901)), (0, -3): (1, Fraction(1, 371)),
        (0, 3): (1, Fraction(1, 371)), (3, 0): (1, Fraction(4, 901)),
    }),
    (25, 29, 43): (Fraction(120, 8729), {
        (-3, 0): (1, Fraction(31, 8729)), (0, -3): (1, Fraction(1, 301)),
        (0, 3): (1, Fraction(1, 301)), (3, 0): (1, Fraction(31, 8729)),
    }),
}


def audit_pattern_completeness():
    short = shell_patterns(18)
    shell = shell_patterns(20, exact=True)
    check(short == SHORT_EXPECTED and len(short) == 26, f"complete <=18 patterns {short}")
    check(shell == SHELL20 and len(shell) == 7, f"complete l1=20 patterns {shell}")
    check(all(sum(p) % 2 == 0 for p in short + shell), "even relation parity")
    check(all(len(relation_vectors(p)) == 18 for p in SHELL20), "18 rays per distinct pattern")

    # A second, triple-first route validates the coefficient-first generator.
    height = 79
    universe = tuple(primitive_triples(height))
    for pattern in short + shell:
        generated = set(generated_relation_map(pattern, height))
        brute = {
            w for w in universe
            if any(dot(c, w) == 0 for c in relation_vectors(pattern))
        }
        check(generated == brute, f"generator/triple-first equality H={height}, {pattern}")
    return short, shell


def audit_multiple_relations(short_patterns, shell_patterns_exact):
    labelled = [(pattern, c) for pattern in shell_patterns_exact for c in relation_vectors(pattern)]
    candidates = set()
    for (_, c), (_, d) in combinations(labelled, 2):
        raw = cross(c, d)
        if raw == (0, 0, 0):
            continue
        if not (all(x > 0 for x in raw) or all(x < 0 for x in raw)):
            continue
        w = tuple(sorted(primitive_direction(raw)))
        if len(set(w)) != 3:
            continue
        if any(x <= 0 or x % 2 == 0 or x % 3 == 0 for x in w):
            continue
        if gcd_many(w) != 1 or direct_relations(w, short_patterns):
            continue
        relations = direct_relations(w, shell_patterns_exact)
        if sum(len(cs) for cs in relations.values()) >= 2:
            candidates.add(w)

    check(candidates == set(MULTIPLE_RELATIONS), f"complete cross-product overlap list {candidates}")
    summaries = []
    for w in sorted(candidates):
        relations = direct_relations(w, shell_patterns_exact)
        flat = tuple((pattern, c) for pattern, cs in relations.items() for c in cs)
        check(flat == MULTIPLE_RELATIONS[w], f"exact multiple relation rays {w}, {flat}")
        check(len(flat) == 2, f"exactly two l1=20 rays {w}")
        check(cross(flat[0][1], flat[1][1]) != (0, 0, 0), f"independent relations {w}")

        literal, representatives = literal_physical_components(w)
        expected_measure, expected_profile = OVERLAP_PROFILES[w]
        check(sum(literal.values(), Fraction(0)) == expected_measure,
              f"overlap physical measure {w}")
        for _, c in flat:
            raw, metadata = raw_relation_components(w, c)
            check(raw == literal, f"all charts give same carrier dictionary {w}, {c}")
            for C, n in representatives.items():
                check(metadata[C][0] == dot(c, n), f"chart defect sidecar {w}, {c}, {C}")

        profile = {}
        for C, length in literal.items():
            defects = tuple(dot(c, representatives[C]) for _, c in flat)
            check(defects != (0, 0), f"independent zero-defect physical overlap {w}, {C}")
            count, mass = profile.get(defects, (0, Fraction(0)))
            profile[defects] = (count + 1, mass + length)
        check(profile == expected_profile, f"joint defect profile {w}, {profile}")
        summaries.append((w, flat, expected_measure, tuple(sorted(profile.items()))))

    cross_count = sum(len({pattern for pattern, _ in MULTIPLE_RELATIONS[w]}) == 2
                      for w in MULTIPLE_RELATIONS)
    within_count = len(MULTIPLE_RELATIONS) - cross_count
    check((cross_count, within_count) == (8, 5), "cross/within overlap split")
    check(sum(OVERLAP_PROFILES[w][0] == 0 for w in MULTIPLE_RELATIONS) == 1,
          "unique empty multiply-presented ray")
    return tuple(summaries), cross_count, within_count


def audit_atlas(short_patterns, shell_patterns_exact):
    common_height = 476
    short_union = set()
    for pattern in short_patterns:
        short_union.update(generated_relation_map(pattern, common_height))

    all_rows = {}
    minimal_rows = {}
    for pattern in shell_patterns_exact:
        all_rows[pattern] = generated_relation_map(pattern, common_height)
        minimal_rows[pattern] = {
            w: cs for w, cs in all_rows[pattern].items() if w not in short_union
        }
        expected_all, expected_minimal, expected_rays = COMMON_COUNTS[pattern]
        check(len(all_rows[pattern]) == expected_all, f"common all count {pattern}")
        check(len(minimal_rows[pattern]) == expected_minimal, f"common minimal count {pattern}")
        check(sum(len(cs) for cs in minimal_rows[pattern].values()) == expected_rays,
              f"common ray count {pattern}")
        for w in all_rows[pattern]:
            has_short_direct = bool(direct_relations(w, short_patterns))
            check(has_short_direct == (w in short_union), f"two-route shorter filter {pattern}, {w}")

    physical_cache = {}
    row_summaries = []
    for pattern in shell_patterns_exact:
        claim = CLAIMS[pattern]
        height = claim["height"]
        all_proof = {w: cs for w, cs in all_rows[pattern].items() if max(w) <= height}
        proof = {w: cs for w, cs in minimal_rows[pattern].items() if max(w) <= height}
        check(len(all_proof) == claim["all"], f"finite all window {pattern}")
        check(len(proof) == claim["minimal"], f"finite minimal window {pattern}")
        check(sum(len(cs) for cs in proof.values()) == claim["rays"],
              f"finite relation rays {pattern}")

        I0 = slice_integral(pattern, 0)
        I3 = slice_integral(pattern, 3)
        check((I0, I3) == claim["integrals"], f"slice integrals {pattern}")
        bulk = Fraction(2, 3) * (I0 + 2 * I3)
        check(bulk == claim["bulk"], f"bulk constant {pattern}")

        values = {}
        component_count = 0
        zero_triples = []
        for w, relations in proof.items():
            if w not in physical_cache:
                physical_cache[w] = literal_physical_components(w)
            literal, representatives = physical_cache[w]
            measure = sum(literal.values(), Fraction(0))
            values[w] = measure
            component_count += len(literal)
            if not measure:
                zero_triples.append(w)
            for c in relations:
                raw, metadata = raw_relation_components(w, c)
                check(raw == literal, f"raw/literal dictionary {pattern}, {w}, {c}")
                for C, n in representatives.items():
                    delta = dot(c, n)
                    check(delta in (-3, 0, 3), f"exact defect range {pattern}, {w}, {c}, {delta}")
                    check(metadata[C][0] == delta,
                          f"raw/literal defect label {pattern}, {w}, {c}, {C}")

        maximum = max(values.values())
        maximizers = tuple(sorted(w for w, value in values.items() if value == maximum))
        positive = sum(value > 0 for value in values.values())
        check(maximum == claim["maximum"], f"sharp sector maximum {pattern}, {maximum}")
        check(maximizers == claim["maximizers"], f"all sector winners {pattern}, {maximizers}")
        check(positive == claim["positive"], f"positive count {pattern}")
        check(component_count == claim["components"], f"component count {pattern}")
        check(tuple(sorted(zero_triples)) == claim["zeros"], f"empty triples {pattern}")

        maximizing_w = maximizers[0]
        maximizing_relations = proof[maximizing_w]
        check(len(maximizing_relations) == 1, f"unique winning relation chart {pattern}")
        maximizing_c = next(iter(maximizing_relations))
        check(maximizing_c == claim["relation"], f"winning signed relation {pattern}, {maximizing_c}")
        raw, metadata = raw_relation_components(maximizing_w, maximizing_c)
        delta_counts = tuple(sum(metadata[C][0] == delta for C in raw) for delta in (-3, 0, 3))
        masses = tuple(
            sum((raw[C] for C in raw if metadata[C][0] == delta), Fraction(0))
            for delta in (-3, 0, 3)
        )
        check(delta_counts == claim["delta_counts"], f"winner defect counts {pattern}, {delta_counts}")
        check(masses == claim["masses"], f"winner defect masses {pattern}, {masses}")
        for delta, expected in ((-3, I3), (0, I0), (3, I3)):
            actual = affine_roof_integral(maximizing_w, maximizing_c, delta)
            check(actual == expected, f"two-route roof integral {pattern}, {delta}, {actual}")

        threshold = Fraction(18, 7) / (maximum - bulk)
        check(threshold == claim["threshold"], f"tail threshold {pattern}, {threshold}")
        check(height == floor_fraction(threshold), f"proof height is threshold floor {pattern}")
        check(bulk + Fraction(18, 7 * (height + 1)) < maximum,
              f"strict tail exclusion above proof height {pattern}")
        check(bulk + Fraction(18, 7 * height) >= maximum,
              f"coarse proof height cannot be lowered {pattern}")
        row_summaries.append({
            "pattern": pattern, "height": height, "all": len(all_proof),
            "minimal": len(proof), "rays": sum(len(cs) for cs in proof.values()),
            "positive": positive, "components": component_count,
            "zeros": tuple(sorted(zero_triples)), "I0": I0, "I3": I3,
            "bulk": bulk, "threshold": threshold, "maximum": maximum,
            "maximizers": maximizers, "relation": maximizing_c,
            "delta_counts": delta_counts, "masses": masses,
        })

    physical_union = set().union(*(set(row) for row in minimal_rows.values()))
    sector_incidences = sum(len(row) for row in minimal_rows.values())
    relation_rays = sum(sum(len(cs) for cs in row.values()) for row in minimal_rows.values())
    multiplicities = {}
    cross_overlaps = {}
    empty = []
    for w in physical_union:
        labels = tuple(pattern for pattern in shell_patterns_exact if w in minimal_rows[pattern])
        multiplicities[len(labels)] = multiplicities.get(len(labels), 0) + 1
        if len(labels) > 1:
            cross_overlaps[w] = labels
        first_pattern = labels[0]
        first_relation = next(iter(minimal_rows[first_pattern][w]))
        components, _ = raw_relation_components(w, first_relation)
        if not components:
            empty.append(w)
    expected_cross = {
        w: tuple(dict.fromkeys(pattern for pattern, _ in MULTIPLE_RELATIONS[w]))
        for w in MULTIPLE_RELATIONS
        if len({pattern for pattern, _ in MULTIPLE_RELATIONS[w]}) == 2
    }
    check(sector_incidences == 13492, "common sector incidences")
    check(relation_rays == 13497, "common relation rays")
    check(len(physical_union) == 13484, "common physical triples")
    check(multiplicities == {1: 13476, 2: 8}, f"common multiplicities {multiplicities}")
    check(cross_overlaps == expected_cross, f"common cross-pattern overlaps {cross_overlaps}")
    check(tuple(sorted(empty)) == ((1, 19, 41),), f"common empty physical triples {empty}")

    global_maximum = max(row["maximum"] for row in row_summaries)
    global_maximizers = tuple(
        (row["pattern"], w)
        for row in row_summaries if row["maximum"] == global_maximum
        for w in row["maximizers"]
    )
    check(global_maximum == Fraction(36, 1295), "global minimal-l20 maximum")
    check(global_maximizers == (((1, 2, 17), (1, 11, 185)),),
          f"all global winners {global_maximizers}")
    incidence = {
        "height": common_height,
        "sector_incidences": sector_incidences,
        "relation_rays": relation_rays,
        "physical_triples": len(physical_union),
        "positive_combs": len(physical_union) - len(empty),
        "empty_combs": tuple(sorted(empty)),
        "multiplicities": multiplicities,
    }
    return tuple(row_summaries), incidence, global_maximum, global_maximizers


def ff(value):
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def json_ready(value):
    if isinstance(value, Fraction):
        return ff(value)
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (tuple, list, set)):
        return [json_ready(item) for item in value]
    return value


def main():
    short, shell = audit_pattern_completeness()
    overlaps, cross_count, within_count = audit_multiple_relations(short, shell)
    rows, incidence, global_maximum, global_maximizers = audit_atlas(short, shell)
    semantic = {
        "short": short, "shell": shell, "rows": rows, "incidence": incidence,
        "overlaps": overlaps, "cross_count": cross_count, "within_count": within_count,
        "global_maximum": global_maximum, "global_maximizers": global_maximizers,
    }
    digest = sha256(json.dumps(json_ready(semantic), sort_keys=True).encode()).hexdigest()

    print("LRC14 MINIMAL TERNARY-UNIT L1=20 CLEAN-ROOM SHELL")
    print("status=PROVED_ELEMENTARY_CANDIDATE+VERIFIED_EXACT; LRC14=OPEN")
    print("minimality=no primitive full-support ternary-unit relation of l1<=18")
    print("patterns=" + ",".join("-".join(map(str, pattern)) for pattern in shell))
    print("defects=-3,0,3; affine_carrier=C_delta+k*c; allowed_k_mod3=1,2")
    for row in rows:
        print(
            "row=" + "-".join(map(str, row["pattern"]))
            + f" I0={ff(row['I0'])} I3={ff(row['I3'])} bulk={ff(row['bulk'])}"
            + f" threshold={ff(row['threshold'])} height={row['height']}"
            + f" all={row['all']} minimal={row['minimal']} rays={row['rays']}"
            + f" positive={row['positive']} components={row['components']}"
            + f" maximum={ff(row['maximum'])}"
            + " maximizers=" + ",".join("{" + ",".join(map(str, w)) + "}" for w in row["maximizers"])
            + " relation=" + ",".join(map(str, row["relation"]))
            + " delta_counts=" + ",".join(map(str, row["delta_counts"]))
            + " masses=" + ",".join(ff(value) for value in row["masses"])
            + " zeros=" + (",".join("{" + ",".join(map(str, w)) + "}" for w in row["zeros"]) or "none")
        )
    print(
        f"incidence_H={incidence['height']}"
        + f" sector_incidences={incidence['sector_incidences']}"
        + f" relation_rays={incidence['relation_rays']}"
        + f" physical_triples={incidence['physical_triples']}"
        + f" positive_combs={incidence['positive_combs']}"
        + " multiplicities=" + ",".join(f"{k}:{v}" for k, v in sorted(incidence["multiplicities"].items()))
    )
    print(f"multiple_relations={len(overlaps)} cross_pattern={cross_count} within_pattern={within_count}")
    print("multiple_empty={1,19,41}; nonempty_joint_defects_never_include=(0,0)")
    print(f"global_maximum={ff(global_maximum)} at="
          + ",".join("-".join(map(str, pattern)) + ":{" + ",".join(map(str, w)) + "}"
                     for pattern, w in global_maximizers))
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()

