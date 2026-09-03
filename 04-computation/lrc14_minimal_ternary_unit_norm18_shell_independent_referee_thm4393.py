#!/usr/bin/env python3
"""Independent exact audit of the minimal ternary-unit l1=18 LRC shell.

This program deliberately imports no repository mathematics and does not
import the producer verifier.  All proof checks use exact Fraction arithmetic
and explicit RuntimeError gates that remain active under ``python -O``.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd
import sys


Q = Fraction
R = Q(3, 14)
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("CHECK FAILED: " + label)


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


def normalize_ray(c):
    c = tuple(c)
    g = gcd_many(c)
    check(g > 0, "nonzero ray")
    c = tuple(x // g for x in c)
    first = next(x for x in c if x)
    if first < 0:
        c = tuple(-x for x in c)
    return c


def magnitude_patterns(total):
    """Primitive sorted positive ternary-unit magnitude partitions."""
    out = []
    for a in range(1, total + 1):
        for b in range(a, total + 1):
            c = total - a - b
            if c < b:
                continue
            if c <= 0:
                continue
            if any(x % 3 == 0 for x in (a, b, c)):
                continue
            if gcd_many((a, b, c)) != 1:
                continue
            out.append((a, b, c))
    return tuple(out)


def signed_rays(pattern):
    """All coordinate-labelled signed rays, modulo simultaneous negation."""
    out = set()
    for p in set(permutations(pattern)):
        for signs in product((-1, 1), repeat=3):
            c = tuple(s * x for s, x in zip(signs, p))
            if all(x > 0 for x in c) or all(x < 0 for x in c):
                continue
            out.add(normalize_ray(c))
    return tuple(sorted(out))


def relation_rays(w, rays):
    return tuple(c for c in rays if dot(c, w) == 0)


def floor_q(x):
    return x.numerator // x.denominator


def ceil_q(x):
    return -floor_q(-x)


def egcd_nonnegative(a, b):
    if b == 0:
        return (a, 1, 0)
    g, x1, y1 = egcd_nonnegative(b, a % b)
    return (g, y1, x1 - (a // b) * y1)


def bezout_vector(v):
    """Return z with z dot v = gcd(v), with signs handled exactly."""
    a, b, c = v
    g12, x, y = egcd_nonnegative(abs(a), abs(b))
    g, s, t = egcd_nonnegative(g12, abs(c))
    z = (
        x * s * (1 if a >= 0 else -1),
        y * s * (1 if b >= 0 else -1),
        t * (1 if c >= 0 else -1),
    )
    check(dot(z, v) == g, "Bezout identity")
    return g, z


def roof(w, C):
    w1, w2, w3 = w
    C1, C2, C3 = C
    values = (
        2 * R / w1,
        2 * R / w2,
        2 * R / w3,
        R / w1 + R / w2 - Q(abs(C3), w1 * w2),
        R / w1 + R / w3 - Q(abs(C2), w1 * w3),
        R / w2 + R / w3 - Q(abs(C1), w2 * w3),
    )
    return max(Q(0), min(values))


def strict_integer_range(lower, upper):
    """All integers k with lower < k < upper."""
    if lower >= upper:
        return range(0)
    return range(floor_q(lower) + 1, ceil_q(upper))


def quotient_components(w, c):
    """Exact carrier dictionary and defect masses for one relation chart."""
    check(dot(w, c) == 0, "relation annihilates speed")
    g, v = bezout_vector(c)
    check(g == 1, "primitive coefficient ray")
    components = {}
    masses = {}
    for delta in (-3, 0, 3):
        C0 = cross(w, tuple(delta * x for x in v))
        check(dot(C0, w) == 0, "base carrier lies in speed kernel")
        residue_types = []
        for k0 in range(3):
            residues = tuple((C0[i] + k0 * c[i]) % 3 for i in range(3))
            residue_types.append(residues)
        check(sum(all(x != 0 for x in row) for row in residue_types) == 2,
              "two owner-live residue classes")
        check(sum(all(x == 0 for x in row) for row in residue_types) == 1,
              "one owner-collision residue class")
        check(all((all(x != 0 for x in row) or all(x == 0 for x in row))
                  for row in residue_types), "no partial owner residue class")
        lower = None
        upper = None
        for i in range(3):
            other = [j for j in range(3) if j != i]
            bound = R * (w[other[0]] + w[other[1]])
            q = c[i]
            lo = (-bound - C0[i]) / q
            hi = (bound - C0[i]) / q
            if lo > hi:
                lo, hi = hi, lo
            lower = lo if lower is None else max(lower, lo)
            upper = hi if upper is None else min(upper, hi)
        mass = Q(0)
        for k in strict_integer_range(lower, upper):
            C = tuple(C0[i] + k * c[i] for i in range(3))
            if any(x % 3 == 0 for x in C):
                continue
            value = roof(w, C)
            if value == 0:
                continue
            check(C not in components, "defect fibres are carrier-disjoint")
            components[C] = value
            mass += value
        masses[delta] = mass
    return components, masses


@lru_cache(maxsize=None)
def danger_intervals(a, sheet):
    """Intervals in [0,1] for ||a(x+sheet/3)||<1/14, annotated by integer."""
    radius = Q(1, 14 * a)
    lower_m = floor_q(Q(a * sheet, 3) - Q(1, 14)) + 1
    upper_m = ceil_q(Q(a * (3 + sheet), 3) + Q(1, 14)) - 1
    out = []
    for m in range(lower_m, upper_m + 1):
        centre = Q(m, a) - Q(sheet, 3)
        lo = max(Q(0), centre - radius)
        hi = min(Q(1), centre + radius)
        if lo < hi:
            out.append((lo, hi, m))
    out.sort()
    return tuple(out)


def intersect_three_annotated(lists):
    pos = [0, 0, 0]
    out = []
    while all(pos[i] < len(lists[i]) for i in range(3)):
        current = [lists[i][pos[i]] for i in range(3)]
        lo = max(x[0] for x in current)
        hi = min(x[1] for x in current)
        if lo < hi:
            out.append((lo, hi, tuple(x[2] for x in current)))
        next_end = min(x[1] for x in current)
        for i in range(3):
            if current[i][1] == next_end:
                pos[i] += 1
    return out


@lru_cache(maxsize=None)
def six_sheet_components(w):
    """Literal six-sheet x-circle construction, grouped by raw carrier."""
    components = defaultdict(Fraction)
    literal_pieces = 0
    for sheets in permutations((0, 1, 2)):
        lists = tuple(danger_intervals(w[i], sheets[i]) for i in range(3))
        for lo, hi, nearest in intersect_three_annotated(lists):
            n = tuple(3 * nearest[i] - w[i] * sheets[i] for i in range(3))
            C = cross(w, n)
            components[C] += hi - lo
            literal_pieces += 1
    return tuple(sorted(components.items())), literal_pieces


@lru_cache(maxsize=None)
def body_intervals(a):
    """Intervals in the physical y-circle for ||a*y-n||<3/14."""
    radius = R / a
    out = []
    for n in range(0, a + 1):
        centre = Q(n, a)
        lo = max(Q(0), centre - radius)
        hi = min(Q(1), centre + radius)
        if lo < hi:
            out.append((lo, hi, n))
    return tuple(out)


@lru_cache(maxsize=None)
def physical_components(w):
    """Definition-level y-circle comb with the distinct-owner gate."""
    components = defaultdict(Fraction)
    literal_pieces = 0
    lists = tuple(body_intervals(a) for a in w)
    for lo, hi, nearest in intersect_three_annotated(lists):
        owners = tuple((-w[i] * nearest[i]) % 3 for i in range(3))
        if set(owners) != {0, 1, 2}:
            continue
        C = cross(w, nearest)
        components[C] += hi - lo
        literal_pieces += 1
    return tuple(sorted(components.items())), literal_pieces


def physical_dict(w):
    items, pieces = physical_components(tuple(w))
    return dict(items), pieces


def clip_halfplane(poly, A, B, D):
    if not poly:
        return []
    out = []
    for P, S in zip(poly, poly[-1:] + poly[:-1]):
        fP = A * P[0] + B * P[1] - D
        fS = A * S[0] + B * S[1] - D
        inP = fP <= 0
        inS = fS <= 0
        if inP != inS:
            t = fS / (fS - fP)
            X = (S[0] + t * (P[0] - S[0]), S[1] + t * (P[1] - S[1]))
            out.append(X)
        if inP:
            out.append(P)
    return out


def polygon_area(poly):
    if len(poly) < 3:
        return Q(0)
    twice = sum(
        poly[i][0] * poly[(i + 1) % len(poly)][1]
        - poly[(i + 1) % len(poly)][0] * poly[i][1]
        for i in range(len(poly))
    )
    return abs(twice) / 2


def slice_bulk_polygon(pattern, delta):
    """Independent strip-clipping derivation of the cube-slice density."""
    a, b, c = pattern
    poly = [(-R, -R), (R, -R), (R, R), (-R, R)]
    # |delta + a*x + b*y| <= c*r.
    poly = clip_halfplane(poly, a, b, c * R - delta)
    poly = clip_halfplane(poly, -a, -b, c * R + delta)
    return polygon_area(poly) / c


def positive_part(x):
    return max(Q(0), x)


def slice_bulk_convolution(pattern, delta):
    total = Q(0)
    for mask in product((0, 1), repeat=3):
        removed = sum(pattern[i] for i in range(3) if mask[i])
        sign = -1 if sum(mask) % 2 else 1
        total += sign * positive_part(abs(delta) + R * sum(pattern) - 2 * R * removed) ** 2
    return total / (2 * pattern[0] * pattern[1] * pattern[2])


SHORT_EXPECTED = (
    (1, 1, 2), (1, 1, 4), (1, 1, 8), (1, 1, 10), (1, 1, 14),
    (1, 2, 5), (1, 2, 7), (1, 2, 11), (1, 2, 13),
    (1, 4, 5), (1, 4, 7), (1, 4, 11),
    (1, 5, 8), (1, 5, 10), (1, 7, 8),
    (2, 5, 5), (2, 5, 7), (2, 7, 7), (4, 5, 5), (4, 5, 7),
)

PATTERNS = (
    (1, 1, 16),
    (1, 4, 13),
    (1, 7, 10),
    (2, 5, 11),
    (4, 7, 7),
    (5, 5, 8),
)

LABEL = {
    (1, 1, 16): "116",
    (1, 4, 13): "1413",
    (1, 7, 10): "1710",
    (2, 5, 11): "2511",
    (4, 7, 7): "477",
    (5, 5, 8): "558",
}

EXPECTED = {
    (1, 1, 16): {
        "height": 400, "all": 501, "minimal": 497, "rays": 497,
        "positive": 497, "pieces": 10730, "maximum": Q(36, 1225),
        "maximizers": ((1, 11, 175),),
        "masses": {-3: Q(12, 1225), 0: Q(12, 1225), 3: Q(12, 1225)},
        "bulk": (Q(9, 784), Q(9, 784)), "baseline": Q(9, 392),
        "threshold": Q(400),
    },
    (1, 4, 13): {
        "height": 335, "all": 876, "minimal": 866, "rays": 867,
        "positive": 864, "pieces": 15456, "maximum": Q(12, 497),
        "maximizers": ((5, 7, 71),),
        "masses": {-3: Q(3, 497), 0: Q(6, 497), 3: Q(3, 497)},
        "bulk": (Q(9, 637), Q(27, 5096)), "baseline": Q(3, 182),
        "threshold": Q(3692, 11),
    },
    (1, 7, 10): {
        "height": 222, "all": 495, "minimal": 485, "rays": 486,
        "positive": 485, "pieces": 6962, "maximum": Q(12, 413),
        "maximizers": ((1, 7, 59),),
        "masses": {-3: Q(3, 413), 0: Q(6, 413), 3: Q(3, 413)},
        "bulk": (Q(9, 490), Q(27, 6860)), "baseline": Q(6, 343),
        "threshold": Q(2891, 13),
    },
    (2, 5, 11): {
        "height": 395, "all": 1439, "minimal": 1429, "rays": 1430,
        "positive": 1428, "pieces": 32904, "maximum": Q(318, 14399),
        "maximizers": ((11, 17, 121),),
        "masses": {-3: Q(57, 14399), 0: Q(12, 847), 3: Q(57, 14399)},
        "bulk": (Q(9, 539), Q(9, 2695)), "baseline": Q(6, 385),
        "threshold": Q(10285, 26),
    },
    (4, 7, 7): {
        "height": 357, "all": 789, "minimal": 785, "rays": 785,
        "positive": 784, "pieces": 15756, "maximum": Q(144, 5831),
        "maximizers": ((11, 17, 49),),
        "masses": {-3: Q(22, 5831), 0: Q(100, 5831), 3: Q(22, 5831)},
        "bulk": (Q(54, 2401), Q(9, 4802)), "baseline": Q(6, 343),
        "threshold": Q(357),
    },
    (5, 5, 8): {
        "height": 376, "all": 855, "minimal": 851, "rays": 851,
        "positive": 851, "pieces": 17668, "maximum": Q(172, 7175),
        "maximizers": ((1, 25, 41),),
        "masses": {-3: Q(22, 7175), 0: Q(128, 7175), 3: Q(22, 7175)},
        "bulk": (Q(27, 1225), Q(9, 4900)), "baseline": Q(3, 175),
        "threshold": Q(18450, 49),
    },
}

H400_EXPECTED = {
    (1, 1, 16): (501, 497, 497),
    (1, 4, 13): (1237, 1227, 1228),
    (1, 7, 10): (1612, 1602, 1603),
    (2, 5, 11): (1463, 1453, 1454),
    (4, 7, 7): (981, 977, 977),
    (5, 5, 8): (972, 968, 968),
}

OVERLAP_EXPECTED = {
    (1, 19, 35): ("116", "477"),
    (1, 25, 41): ("116", "558"),
    (1, 37, 53): ("116", "1710"),
    (5, 13, 41): ("1413", "1710"),
    (7, 25, 29): ("1413", "477"),
    (11, 29, 43): ("1413", "1710"),
    (13, 19, 47): ("1413", "2511"),
}

WITHIN_EXPECTED = {
    (1, 17, 55): "1413",
    (13, 23, 31): "1710",
    (1, 17, 37): "2511",
}


def speed_universe(height):
    speeds = [x for x in range(1, height + 1, 2) if x % 3]
    return tuple(w for w in combinations(speeds, 3) if gcd_many(w) == 1)


def main():
    if "--tripwire" in sys.argv:
        check(False, "optimization-live tripwire")

    generated_short = tuple(
        p for total in range(4, 17, 2) for p in magnitude_patterns(total)
    )
    check(set(generated_short) == set(SHORT_EXPECTED), "complete l1<=16 pattern shell")
    check(magnitude_patterns(18) == PATTERNS, "six l1=18 patterns")

    short_rays = tuple(sorted({c for p in generated_short for c in signed_rays(p)}))
    pattern_rays = {p: signed_rays(p) for p in PATTERNS}
    all18_rays = tuple(sorted({c for p in PATTERNS for c in pattern_rays[p]}))

    all_records = {p: {} for p in PATTERNS}
    minimal_records = {p: {} for p in PATTERNS}
    universe = speed_universe(400)
    for w in universe:
        hits = {}
        for p in PATTERNS:
            rays = relation_rays(w, pattern_rays[p])
            if rays:
                hits[p] = rays
                all_records[p][w] = rays
        if not hits:
            continue
        shorter = any(dot(c, w) == 0 for c in short_rays)
        if not shorter:
            for p, rays in hits.items():
                minimal_records[p][w] = rays

    print("status=running")
    print("short_patterns=" + repr(generated_short))
    print("norm18_patterns=" + repr(PATTERNS))

    proof_rows = []
    for p in PATTERNS:
        expected = EXPECTED[p]
        height = expected["height"]
        all_box = {w: rays for w, rays in all_records[p].items() if w[-1] <= height}
        min_box = {w: rays for w, rays in minimal_records[p].items() if w[-1] <= height}

        ray_count = sum(len(rays) for rays in min_box.values())
        check(len(all_box) == expected["all"], LABEL[p] + " all count")
        check(len(min_box) == expected["minimal"], LABEL[p] + " minimal count")
        check(ray_count == expected["rays"], LABEL[p] + " ray count")

        I0_formula = slice_bulk_convolution(p, 0)
        I3_formula = slice_bulk_convolution(p, 3)
        I0_clip = slice_bulk_polygon(p, 0)
        I3_clip = slice_bulk_polygon(p, 3)
        check(I0_formula == I0_clip == expected["bulk"][0], LABEL[p] + " I0")
        check(I3_formula == I3_clip == expected["bulk"][1], LABEL[p] + " I3")
        baseline = Q(2, 3) * (I0_formula + 2 * I3_formula)
        check(baseline == expected["baseline"], LABEL[p] + " baseline")
        threshold = Q(18, 7) / (expected["maximum"] - baseline)
        check(threshold == expected["threshold"], LABEL[p] + " tail threshold")
        check(floor_q(threshold) == height, LABEL[p] + " proof height")
        check(baseline + Q(18, 7 * (height + 1)) < expected["maximum"], LABEL[p] + " tail strictness")

        measures = {}
        positive = 0
        literal_components = 0
        for w, rays in sorted(min_box.items()):
            direct, pieces = physical_dict(w)
            literal_components += pieces
            check(pieces == len(direct), LABEL[p] + " one y-component per carrier " + repr(w))
            measure = sum(direct.values(), Q(0))
            measures[w] = measure
            if measure > 0:
                positive += 1
            for c in rays:
                quotient, masses = quotient_components(w, c)
                check(quotient == direct, LABEL[p] + " quotient/physical " + repr((w, c)))
                check(sum(masses.values(), Q(0)) == measure, LABEL[p] + " defect total " + repr((w, c)))

        maximum = max(measures.values())
        maximizers = tuple(w for w in sorted(measures) if measures[w] == maximum)
        check(positive == expected["positive"], LABEL[p] + " positive count")
        check(literal_components == expected["pieces"], LABEL[p] + " literal component count " + repr((literal_components, expected["pieces"])))
        check(maximum == expected["maximum"], LABEL[p] + " sharp maximum")
        check(maximizers == expected["maximizers"], LABEL[p] + " unique maximizer")
        max_ray = min_box[maximizers[0]][0]
        _, max_masses = quotient_components(maximizers[0], max_ray)
        check(max_masses == expected["masses"], LABEL[p] + " maximizing defect masses")

        proof_rows.append((LABEL[p], height, len(all_box), len(min_box), ray_count,
                           positive, literal_components, maximum, maximizers,
                           I0_formula, I3_formula, baseline, threshold, max_masses))

    # Common-height incidence audit.
    union_patterns = defaultdict(set)
    union_rays = defaultdict(set)
    h400_rows = []
    for p in PATTERNS:
        all_count = len(all_records[p])
        min_count = len(minimal_records[p])
        rays_count = sum(len(rays) for rays in minimal_records[p].values())
        check((all_count, min_count, rays_count) == H400_EXPECTED[p], LABEL[p] + " H400 row")
        h400_rows.append((LABEL[p], all_count, min_count, rays_count))
        for w, rays in minimal_records[p].items():
            union_patterns[w].add(LABEL[p])
            union_rays[w].update(rays)

    sector_incidences = sum(len(x) for x in union_patterns.values())
    relation_ray_count = sum(len(x) for x in union_rays.values())
    check(sector_incidences == 6724, "H400 sector incidences")
    check(relation_ray_count == 6727, "H400 relation rays")
    check(len(union_patterns) == 6717, "H400 physical triples")

    positive_union = {}
    for w in sorted(union_patterns):
        direct, _ = physical_dict(w)
        positive_union[w] = sum(direct.values(), Q(0))
    check(sum(value > 0 for value in positive_union.values()) == 6714, "H400 positive combs")

    cross_overlaps = {
        w: tuple(sorted(labels)) for w, labels in union_patterns.items() if len(labels) > 1
    }
    check(cross_overlaps == OVERLAP_EXPECTED, "seven cross-pattern overlaps")
    within = {}
    for p in PATTERNS:
        for w, rays in minimal_records[p].items():
            if len(rays) > 1:
                within[w] = LABEL[p]
    check(within == WITHIN_EXPECTED, "three within-pattern double presentations")

    # Height-free completeness of the multiple-relation list: every such
    # speed ray is the primitive positive cross product of two coefficient rays.
    global_multi_candidates = set()
    for c, d in combinations(all18_rays, 2):
        q = cross(c, d)
        if q == (0, 0, 0):
            continue
        if all(x < 0 for x in q):
            q = tuple(-x for x in q)
        if not all(x > 0 for x in q):
            continue
        q = tuple(x // gcd_many(q) for x in q)
        order = tuple(sorted(range(3), key=lambda i: q[i]))
        w = tuple(q[i] for i in order)
        if len(set(w)) != 3 or any(x % 2 == 0 or x % 3 == 0 for x in w):
            continue
        if gcd_many(w) != 1:
            continue
        if any(dot(s, w) == 0 for s in short_rays):
            continue
        actual = tuple(r for r in all18_rays if dot(r, w) == 0)
        if len(actual) >= 2:
            global_multi_candidates.add(w)

    observed_multi = {w for w, rays in union_rays.items() if len(rays) >= 2}
    expected_multi = set(OVERLAP_EXPECTED) | set(WITHIN_EXPECTED)
    check(global_multi_candidates == expected_multi, "height-free multiple relation classification")
    check(observed_multi == expected_multi, "H400 contains all multiple relation rays")

    # Defect-sidecar audit on all ten multiply-related triples.
    multi_rows = []
    nonempty_count = 0
    for w in sorted(global_multi_candidates):
        relations = tuple(r for r in all18_rays if dot(r, w) == 0)
        check(len(relations) == 2, "exactly two l1=18 rays " + repr(w))
        direct, _ = physical_dict(w)
        measure = sum(direct.values(), Q(0))
        if measure > 0:
            nonempty_count += 1
        g, u = bezout_vector(w)
        check(g == 1, "primitive speed Bezout " + repr(w))
        defect_profile = defaultdict(lambda: [0, Q(0)])
        for C, length in direct.items():
            n = cross(C, u)
            check(cross(w, n) == C, "inverse carrier lift " + repr((w, C)))
            defects = tuple(dot(c, n) for c in relations)
            check(all(delta in (-3, 0, 3) for delta in defects), "defect range " + repr((w, defects)))
            check(not all(delta == 0 for delta in defects), "no double-zero component " + repr((w, C)))
            defect_profile[defects][0] += 1
            defect_profile[defects][1] += length
        multi_rows.append((w, tuple(LABEL[tuple(sorted(abs(x) for x in c))] for c in relations),
                           measure, tuple((key, value[0], value[1]) for key, value in sorted(defect_profile.items()))))
    check(nonempty_count == 9, "nine of ten multiply-related triples are nonempty")
    empty_multi = tuple(w for w in sorted(global_multi_candidates) if positive_union[w] == 0)
    check(empty_multi == ((7, 25, 29),), "unique empty multiple relation triple")

    # Independent normalization check against the six shifted x-sheets.  The
    # map x -> y=3x has three physical lifts; grouping by C preserves mass.
    lift_controls = tuple(sorted(set(
        w for row in proof_rows for w in row[8]
    ) | global_multi_candidates))
    lift_rows = []
    for w in lift_controls:
        y_direct, _ = physical_dict(w)
        x_items, x_pieces = six_sheet_components(w)
        x_direct = dict(x_items)
        check(x_direct == y_direct, "six-sheet/y carrier dictionary " + repr(w))
        check(x_pieces == 3 * len(y_direct), "three x-lifts per y carrier " + repr(w))
        lift_rows.append((w, len(y_direct), x_pieces, sum(y_direct.values(), Q(0))))

    overall_max = max(row[7] for row in proof_rows)
    overall_winners = tuple((row[0], row[8]) for row in proof_rows if row[7] == overall_max)
    check(overall_max == Q(36, 1225), "overall shell maximum")
    check(overall_winners == (("116", ((1, 11, 175),)),), "overall shell unique winner")

    print("proof_rows=" + repr(tuple(proof_rows)))
    print("h400_rows=" + repr(tuple(h400_rows)))
    print("h400_totals=" + repr((sector_incidences, relation_ray_count, len(union_patterns),
                                  sum(value > 0 for value in positive_union.values()))))
    print("cross_pattern_overlaps=" + repr(tuple(sorted((w, labels, positive_union[w]) for w, labels in cross_overlaps.items()))))
    print("within_pattern_duplicates=" + repr(tuple(sorted((w, label, positive_union[w]) for w, label in within.items()))))
    print("height_free_multiple_relation_rows=" + repr(tuple(multi_rows)))
    print("multiple_relation_nonempty_count=" + repr(nonempty_count))
    print("multiple_relation_empty=" + repr(empty_multi))
    print("six_sheet_lift_controls=" + repr(tuple(lift_rows)))
    print("overall_sharp=" + repr((overall_max, overall_winners)))
    print("physical_cache=" + repr(physical_components.cache_info()))
    print("six_sheet_cache=" + repr(six_sheet_components.cache_info()))
    print("optimization_safe_checks=yes")
    print("explicit_checks=" + str(CHECKS))
    print("scope=minimal_primitive_full_support_ternary_unit_l1_18_shell_only; LRC(14)=OPEN")
    print("status=PASS")


if __name__ == "__main__":
    main()

