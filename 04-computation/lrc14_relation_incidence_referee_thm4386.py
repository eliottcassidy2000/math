#!/usr/bin/env python3
"""Optimization-safe exact relation-incidence referee for THM-4386.

No repository-local module is imported.  Every load-bearing check uses the
explicit ``check`` function, so ``python -O`` leaves the audit live.
"""

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd


HEIGHT = 199
RADIUS = F(1, 14)
MASTER_PATTERNS = (
    (1, 2), (1, 4), (1, 8), (1, 10),
    (2, 5), (2, 7), (2, 11),
    (4, 5), (4, 7), (5, 8),
)
MASTER_MAGNITUDES = {tuple(sorted((1, h, m))) for h, m in MASTER_PATTERNS}
EXPECTED_SECTOR_COUNTS = {
    (1, 2): 1023, (1, 4): 515, (1, 8): 257, (1, 10): 203,
    (2, 5): 825, (2, 7): 588, (2, 11): 370,
    (4, 5): 800, (4, 7): 589, (5, 8): 516,
}
EXPECTED_CROSS_OVERLAPS = {
    ((1, 4), (5, 8)): {(1, 7, 11)},
    ((1, 8), (2, 5)): {(1, 5, 13)},
    ((1, 8), (4, 7)): {(1, 5, 13), (1, 11, 19)},
    ((1, 10), (2, 5)): {(1, 7, 17)},
    ((1, 10), (4, 7)): {(1, 13, 23)},
    ((2, 5), (4, 7)): {(1, 5, 13), (5, 7, 11)},
    ((2, 7), (2, 11)): {(1, 7, 25)},
    ((2, 7), (5, 8)): {(1, 5, 17), (7, 11, 19)},
    ((2, 11), (5, 8)): {(5, 7, 31), (7, 19, 29)},
}
EXPECTED_WITHIN_DUPLICATES = {
    (2, 7): {(1, 5, 17)},
    (2, 11): {(5, 7, 41)},
    (5, 8): {(11, 13, 23)},
}


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def mod_one(x):
    return x - x // 1


def add_circle_interval(out, center, radius):
    center = mod_one(center)
    lo, hi = center - radius, center + radius
    if lo < 0:
        out.extend(((F(0), hi), (lo + 1, F(1))))
    elif hi > 1:
        out.extend(((lo, F(1)), (F(0), hi - 1)))
    else:
        out.append((lo, hi))


def merge(intervals):
    ans = []
    for lo, hi in sorted(intervals):
        if not ans or lo > ans[-1][1]:
            ans.append([lo, hi])
        elif hi > ans[-1][1]:
            ans[-1][1] = hi
    return tuple((lo, hi) for lo, hi in ans)


@lru_cache(maxsize=None)
def danger(w, sheet):
    pieces = []
    radius = F(1, 14 * w)
    for n in range(w):
        add_circle_interval(pieces, F(n, w) - F(sheet, 3), radius)
    return merge(pieces)


def intersect(left, right):
    ans = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            ans.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(ans)


@lru_cache(maxsize=None)
def pair_intersection(a, b, b_sheet):
    return intersect(danger(a, 0), danger(b, b_sheet))


def physical_intervals(triple):
    a, b, c = triple
    positive = intersect(pair_intersection(a, b, 1), danger(c, 2))
    negative = intersect(pair_intersection(a, b, 2), danger(c, 1))
    return positive, negative


def physical_measure(triple):
    positive, negative = physical_intervals(triple)
    mu = 3 * sum((hi - lo for lo, hi in positive + negative), F(0))
    return mu, 3 * (len(positive) + len(negative))


def canonical_vector(v):
    d = gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2]))
    check(d > 0, ("zero-vector", v))
    v = tuple(x // d for x in v)
    for x in v:
        if x:
            return v if x > 0 else tuple(-y for y in v)
    raise RuntimeError(("zero-vector-after-normalization", v))


def coefficient_vectors(h, m):
    ans = set()
    for magnitudes in set(permutations((1, h, m))):
        for signs in product((-1, 1), repeat=3):
            ans.add(canonical_vector(tuple(a * b for a, b in zip(magnitudes, signs))))
    return tuple(sorted(ans))


MASTER_VECTORS = {pattern: coefficient_vectors(*pattern) for pattern in MASTER_PATTERNS}


def presentation_vectors(triple, pattern):
    return tuple(v for v in MASTER_VECTORS[pattern]
                 if sum(a * w for a, w in zip(v, triple)) == 0)


def memberships(triple):
    return tuple(pattern for pattern in MASTER_PATTERNS
                 if presentation_vectors(triple, pattern))


def positive_primitive_cross(u, v):
    w = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    if w == (0, 0, 0) or 0 in w:
        return None
    if not (all(x > 0 for x in w) or all(x < 0 for x in w)):
        return None
    if w[0] < 0:
        w = tuple(-x for x in w)
    d = gcd(gcd(w[0], w[1]), w[2])
    w = tuple(x // d for x in w)
    if len(set(w)) != 3 or any(x % 2 == 0 or x % 3 == 0 for x in w):
        return None
    return tuple(sorted(w))


def relation_overlap_atlas():
    cross = {}
    for i, left in enumerate(MASTER_PATTERNS):
        for right in MASTER_PATTERNS[i + 1:]:
            rays = {
                w for u in MASTER_VECTORS[left] for v in MASTER_VECTORS[right]
                if (w := positive_primitive_cross(u, v)) is not None
            }
            if rays:
                cross[(left, right)] = rays
    within = {}
    for pattern in MASTER_PATTERNS:
        vectors = MASTER_VECTORS[pattern]
        rays = {
            w for i, u in enumerate(vectors) for v in vectors[i + 1:]
            if (w := positive_primitive_cross(u, v)) is not None
        }
        if rays:
            within[pattern] = rays
    return cross, within


def primitive_triples(height):
    units = tuple(w for w in range(1, height + 1, 2) if w % 3)
    for triple in combinations(units, 3):
        if gcd(gcd(triple[0], triple[1]), triple[2]) == 1:
            yield triple


def nearest_integer(z):
    return (z + F(1, 2)) // 1


def primitive_carrier(w, n):
    carrier = (
        w[1] * n[2] - w[2] * n[1],
        w[2] * n[0] - w[0] * n[2],
        w[0] * n[1] - w[1] * n[0],
    )
    carrier = canonical_vector(carrier)
    check(sum(a * x for a, x in zip(carrier, w)) == 0,
          ("carrier-speed", w, n, carrier))
    check(sum(a * x for a, x in zip(carrier, n)) == 0,
          ("carrier-lift", w, n, carrier))
    return carrier


def physical_component_data(triple):
    positive, negative = physical_intervals(triple)
    ans = []
    for intervals, sheets in ((positive, (0, 1, 2)), (negative, (0, 2, 1))):
        for lo, hi in intervals:
            x = (lo + hi) / 2
            y = mod_one(3 * x)
            n = tuple(nearest_integer(w * y) for w in triple)
            errors = tuple(w * y - q for w, q in zip(triple, n))
            owners = tuple((-pow(w, -1, 3) * q) % 3 for w, q in zip(triple, n))
            check(all(abs(e) < F(3, 14) for e in errors),
                  ("component-eligibility", triple, lo, hi, y, errors))
            check(len(set(owners)) == 3, ("component-owners", triple, y, owners))
            gauge = {(owners[i] - sheets[i]) % 3 for i in range(3)}
            check(len(gauge) == 1, ("component-sheet-gauge", triple, sheets, owners))
            carrier = primitive_carrier(triple, n)
            check(all(a % 3 for a in carrier),
                  ("carrier-coefficient-unit", triple, n, owners, carrier))
            ans.append((lo, hi, y, n, owners, carrier))
    return tuple(ans)


def small_pattern_split():
    all_patterns = tuple(
        (h, m) for h in range(1, 14) for m in range(h + 1, 14 - h)
        if (h + m) % 2 == 1
    )
    zero = tuple(p for p in all_patterns if p[0] % 3 and p[1] % 3)
    defect = tuple(p for p in all_patterns
                   if (p[0] % 3 == 0) != (p[1] % 3 == 0))
    impossible = tuple(p for p in all_patterns if p[0] % 3 == p[1] % 3 == 0)
    return all_patterns, zero, defect, impossible


def coefficient_shell():
    return tuple(sorted(
        (
            (a, b, c)
            for a in range(1, 15)
            for b in range(a, 15)
            for c in range(b, 15)
            if a + b + c <= 14
            and (a + b + c) % 2 == 0
            and all(x % 3 for x in (a, b, c))
            and gcd(gcd(a, b), c) == 1
        ),
        key=lambda row: (sum(row), row),
    ))


def main():
    all_patterns, zero_patterns, defect_patterns, impossible_patterns = small_pattern_split()
    check(zero_patterns == MASTER_PATTERNS, ("zero-pattern-split", zero_patterns))
    check(defect_patterns == (
        (1, 6), (1, 12), (2, 3), (2, 9), (3, 4),
        (3, 8), (3, 10), (4, 9), (5, 6), (6, 7),
    ), ("defect-pattern-split", defect_patterns))
    check(impossible_patterns == ((3, 6),), ("impossible-pattern-split", impossible_patterns))

    shell = coefficient_shell()
    expected_shell = (
        (1, 1, 2),
        (1, 1, 4),
        (1, 2, 5),
        (1, 1, 8), (1, 2, 7), (1, 4, 5),
        (1, 1, 10), (1, 4, 7), (2, 5, 5),
        (1, 2, 11), (1, 5, 8), (2, 5, 7), (4, 5, 5),
    )
    check(shell == expected_shell, ("coefficient-shell", shell))
    missing_shell = tuple(row for row in shell if row not in MASTER_MAGNITUDES)
    check(missing_shell == ((2, 5, 5), (2, 5, 7), (4, 5, 5)),
          ("missing-carrier-shell", missing_shell))

    cross, within = relation_overlap_atlas()
    check(cross == EXPECTED_CROSS_OVERLAPS, ("cross-overlaps", cross))
    check(within == EXPECTED_WITHIN_DUPLICATES, ("within-duplicates", within))
    unique_cross = sorted(set().union(*cross.values()))
    double_relation_rays = sorted(set(unique_cross) | set().union(*within.values()))
    for triple in double_relation_rays:
        check(physical_measure(triple) == (F(0), 0),
              ("two-short-relations-not-empty", triple, physical_measure(triple)))

    records = []
    for triple in primitive_triples(HEIGHT):
        inc = memberships(triple)
        mu, cells = physical_measure(triple)
        records.append((triple, inc, mu, cells))
    check(len(records) == 47499, ("universe-size", len(records)))
    multiplicity = Counter(len(row[1]) for row in records)
    check(multiplicity == Counter({0: 41825, 1: 5663, 2: 10, 3: 1}),
          ("multiplicity", multiplicity))
    sector_counts = {p: sum(p in row[1] for row in records) for p in MASTER_PATTERNS}
    check(sector_counts == EXPECTED_SECTOR_COUNTS, ("sector-counts", sector_counts))
    covered = [row for row in records if row[1]]
    uncovered = [row for row in records if not row[1]]
    positive = [row for row in records if row[2] > 0]
    check((len(covered), len(uncovered), len(positive)) == (5674, 41825, 47463),
          ("coverage-counts", len(covered), len(uncovered), len(positive)))
    covered_positive = sum(row[2] > 0 for row in covered)
    uncovered_positive = sum(row[2] > 0 for row in uncovered)
    check((covered_positive, uncovered_positive) == (5650, 41813),
          ("positive-coverage", covered_positive, uncovered_positive))

    expected_heights = {
        11: (4, 4, 2, 2),
        23: (56, 52, 41, 37),
        43: (454, 245, 422, 221),
        79: (2910, 883, 2874, 859),
        127: (12231, 2287, 12195, 2263),
        199: (47499, 5674, 47463, 5650),
    }
    height_rows = {}
    for height in expected_heights:
        sub = [row for row in records if row[0][2] <= height]
        height_rows[height] = (
            len(sub), sum(bool(row[1]) for row in sub),
            sum(row[2] > 0 for row in sub),
            sum(bool(row[1]) and row[2] > 0 for row in sub),
        )
    check(height_rows == expected_heights, ("height-table", height_rows))

    earliest_height = min(row[0][2] for row in uncovered)
    earliest = tuple(row[:1] + row[2:] for row in uncovered if row[0][2] == earliest_height)
    check(earliest == (
        ((5, 19, 23), F(2, 665), 6),
        ((7, 19, 23), F(4, 437), 6),
        ((11, 19, 23), F(4, 437), 6),
        ((17, 19, 23), F(10, 2261), 6),
    ), ("earliest-uncovered", earliest_height, earliest))
    peak = max(row[2] for row in uncovered)
    peak_rows = tuple(row[0] for row in uncovered if row[2] == peak)
    check((peak, peak_rows) == (F(36, 1225), ((1, 11, 175),)),
          ("uncovered-peak", peak, peak_rows))

    covered_zero = tuple(row[0] for row in covered if row[2] == 0)
    check(len(covered_zero) == 24, ("covered-zero-count", covered_zero))

    hostile_expected = {
        (1, 11, 25): (((2, 5, 5), (5, -5, 2)), F(38, 1925),
                      (((F(431, 462), F(983, 1050)),),
                       ((F(67, 1050), F(31, 462)),))),
        (5, 17, 25): (((2, 5, 7), (7, -5, 2)), F(4, 425),
                      (((F(641, 1050), F(437, 714)),),
                       ((F(277, 714), F(409, 1050)),))),
        (5, 19, 23): (((4, 5, 5), (4, 5, -5)), F(2, 665),
                      (((F(69, 70), F(787, 798)),),
                       ((F(11, 798), F(1, 70)),))),
    }
    hostile_rows = {}
    for triple, (carrier_expected, mu_expected, intervals_expected) in hostile_expected.items():
        mu, cells = physical_measure(triple)
        components = physical_component_data(triple)
        carriers = {tuple(sorted(abs(x) for x in row[-1])): row[-1] for row in components}
        check((mu, cells) == (mu_expected, 6), ("hostile-measure", triple, mu, cells))
        check(physical_intervals(triple) == intervals_expected,
              ("hostile-intervals", triple, physical_intervals(triple)))
        check(carriers == {carrier_expected[0]: carrier_expected[1]},
              ("hostile-carrier", triple, carriers))
        hostile_rows[triple] = (mu, components)

    # Small coefficient-divisible-by-three relation: defects are forced away
    # from zero on the earliest raw coverage miss.
    div3_relation = (3, -2, 1)
    check(sum(a * w for a, w in zip(div3_relation, (5, 19, 23))) == 0,
          "divisible relation is not a relation")
    div3_defects = sorted({sum(a * q for a, q in zip(div3_relation, row[3]))
                           for row in hostile_rows[(5, 19, 23)][1]})
    check(div3_defects == [-1, 1], ("divisible-relation-defects", div3_defects))

    # The first post-small-window hostile has a unit-coefficient relation of
    # h+m=17, but its physical components sit on defects +/-3 for that
    # presentation; their true zero-defect carrier is (2,5,7).
    boundary_relation = (1, 10, -7)
    check(sum(a * w for a, w in zip(boundary_relation, (5, 17, 25))) == 0,
          "boundary relation is not a relation")
    boundary_defects = sorted({sum(a * q for a, q in zip(boundary_relation, row[3]))
                               for row in hostile_rows[(5, 17, 25)][1]})
    check(boundary_defects == [-3, 3], ("boundary-defects", boundary_defects))

    # Presentation-overlap hostile: one physical comb, many algebraic labels.
    extra_patterns = MASTER_PATTERNS + defect_patterns
    extra_vectors = {p: coefficient_vectors(*p) for p in extra_patterns}
    multi_labels = tuple(p for p in extra_patterns
                         if any(sum(a*w for a,w in zip(v, (1, 5, 11))) == 0
                                for v in extra_vectors[p]))
    check(multi_labels == ((1, 2), (1, 6), (3, 4), (4, 9)),
          ("same-comb-multiple-presentations", multi_labels))
    check(physical_measure((1, 5, 11)) == (F(6, 77), 6),
          ("same-comb-hostile-measure", physical_measure((1, 5, 11))))
    check(physical_measure((1, 5, 13)) == (F(0), 0),
          ("master-overlap-hostile", physical_measure((1, 5, 13))))

    print("SMALL_ONE_COEFFICIENT_SPLIT")
    print("zero_defect", zero_patterns)
    print("nonzero_mod3_defect", defect_patterns)
    print("mod3_impossible", impossible_patterns)
    print("CANONICAL_CARRIER_L1_LE_14", shell)
    print("FIRST_MISSING_CARRIER_PATTERNS", missing_shell)
    print("H199_COUNTS", {
        "universe": len(records), "positive_comb": len(positive),
        "presentation_covered": len(covered),
        "presentation_uncovered": len(uncovered),
        "positive_covered": covered_positive,
        "positive_uncovered": uncovered_positive,
        "covered_zero": len(covered) - covered_positive,
    })
    print("H199_MULTIPLICITY", dict(sorted(multiplicity.items())))
    print("H199_SECTOR_COUNTS", sector_counts)
    print("HEIGHT_TABLE_all_covered_positive_positivecovered", height_rows)
    print("CROSS_SECTOR_OVERLAPS", cross)
    print("WITHIN_SECTOR_DUPLICATE_RAYS", within)
    print("DOUBLE_RELATION_RAYS", double_relation_rays)
    print("COVERED_ZERO_TRIPLES", covered_zero)
    print("EARLIEST_UNCOVERED", earliest_height, earliest)
    print("H199_UNCOVERED_PEAK", peak, peak_rows)
    print("MINIMAL_CARRIER_HOSTILES")
    for triple in hostile_expected:
        print(triple, hostile_rows[triple])
    print("DIV3_RELATION_DEFECTS", div3_relation, div3_defects)
    print("POST_WINDOW_RELATION_DEFECTS", boundary_relation, boundary_defects)
    print("SAME_COMB_MULTIPLE_PRESENTATIONS", (1, 5, 11), multi_labels, F(6, 77))
    print("PASS")


if __name__ == "__main__":
    main()
