#!/usr/bin/env python3
"""Independent exact audit for THM-4398 one-zero relation sectors through l1=14.

No repository or producer mathematics is imported.  Every proof gate is an
explicit RuntimeError check, so the audit remains live under ``python -O``.
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


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError("CHECK FAILED: " + label)


def gcd_all(values):
    ans = 0
    for value in values:
        ans = gcd(ans, abs(value))
    return ans


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def normalize(v):
    g = gcd_all(v)
    check(g > 0, "nonzero relation ray")
    v = tuple(x // g for x in v)
    if next(x for x in v if x) < 0:
        v = tuple(-x for x in v)
    return v


def signed_rays(pattern):
    out = set()
    for magnitudes in set(permutations(pattern)):
        for signs in product((-1, 1), repeat=3):
            c = tuple(a * s for a, s in zip(magnitudes, signs))
            if min(c) < 0 < max(c):
                out.add(normalize(c))
    return tuple(sorted(out))


def patterns_through(limit, zero_count):
    out = []
    for total in range(4, limit + 1, 2):
        for a in range(1, total):
            for b in range(a, total):
                c = total - a - b
                if c < b or c <= 0:
                    continue
                if gcd_all((a, b, c)) != 1:
                    continue
                if sum(x % 3 == 0 for x in (a, b, c)) != zero_count:
                    continue
                out.append((a, b, c))
    return tuple(out)


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
    ans = (
        x * s * (1 if a >= 0 else -1),
        y * s * (1 if b >= 0 else -1),
        t * (1 if c >= 0 else -1),
    )
    check(dot(ans, v) == g, "Bezout identity")
    return g, ans


def relation_defects(c):
    """All integer defects allowed by strictness and the mod-3 type."""
    norm = sum(abs(x) for x in c)
    bound = R * norm
    zeros = sum(x % 3 == 0 for x in c)
    residue = 0 if zeros == 0 else None
    out = []
    for delta in range(-ceil_q(bound), ceil_q(bound) + 1):
        if abs(delta) >= bound:
            continue
        if residue == 0 and delta % 3 == 0:
            out.append(delta)
        if residue is None and delta % 3 != 0:
            out.append(delta)
    return tuple(out)


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
    """Intrinsic C_delta+k*c enumeration for either residue type."""
    check(dot(w, c) == 0, "relation annihilates speed")
    g, v = bezout(c)
    check(g == 1, "primitive relation")
    expected_live_classes = 2 if all(x % 3 for x in c) else 1
    components = {}
    masses = {}
    for delta in relation_defects(c):
        C0 = cross(w, tuple(delta * x for x in v))
        check(dot(C0, w) == 0, "base carrier")
        residue_rows = tuple(
            tuple((C0[i] + k * c[i]) % 3 for i in range(3))
            for k in range(3)
        )
        live = sum(all(x != 0 for x in row) for row in residue_rows)
        check(live == expected_live_classes, "live k classes " + repr((w, c, delta, residue_rows)))

        lower = None
        upper = None
        for i in range(3):
            js = tuple(j for j in range(3) if j != i)
            bound = R * (w[js[0]] + w[js[1]])
            a = (-bound - C0[i]) / c[i]
            b = (bound - C0[i]) / c[i]
            if a > b:
                a, b = b, a
            lower = a if lower is None else max(lower, a)
            upper = b if upper is None else min(upper, b)

        mass = Q(0)
        for k in integer_open_interval(lower, upper):
            C = tuple(C0[i] + k * c[i] for i in range(3))
            if any(x % 3 == 0 for x in C):
                continue
            value = roof(w, C)
            if value <= 0:
                continue
            check(C not in components, "disjoint defect carrier fibres")
            components[C] = value
            mass += value
        masses[delta] = mass
    return components, masses


def presentation_universe(pattern, height):
    """Solve each labelled relation for the third sorted speed."""
    speeds = tuple(x for x in range(1, height + 1, 2) if x % 3)
    speed_set = set(speeds)
    out = defaultdict(set)
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
                out[w].add(c)
    return {w: tuple(sorted(rays)) for w, rays in sorted(out.items())}


def brute_universe(pattern, height):
    speeds = tuple(x for x in range(1, height + 1, 2) if x % 3)
    rays = signed_rays(pattern)
    out = {}
    for w in combinations(speeds, 3):
        if gcd_all(w) != 1:
            continue
        hits = tuple(c for c in rays if dot(c, w) == 0)
        if hits:
            out[w] = hits
    return out


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
    """Definition-level body-circle comb, with no quotient formula."""
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


PATTERNS = (
    (1, 2, 3),
    (1, 1, 6), (1, 3, 4),
    (2, 3, 5),
    (1, 2, 9), (1, 3, 8), (1, 5, 6), (2, 3, 7), (3, 4, 5),
    (1, 1, 12), (1, 3, 10), (1, 4, 9), (1, 6, 7), (3, 4, 7),
)


UNIT_PATTERNS = (
    (1, 1, 2), (1, 1, 4), (1, 2, 5),
    (1, 1, 8), (1, 2, 7), (1, 4, 5),
    (1, 1, 10), (1, 4, 7), (2, 5, 5),
    (1, 2, 11), (1, 5, 8), (2, 5, 7), (4, 5, 5),
)


EXPECTED = {
    (1, 2, 3): dict(height=30, triples=51, rays=55, positive=33,
                    maximum=Q(8, 245), winners=((1, 5, 7),),
                    masses={-1: Q(4, 245), 1: Q(4, 245)},
                    bulks={1: Q(1, 147)}, baseline=Q(2, 441), threshold=Q(945, 31)),
    (1, 1, 6): dict(height=14, triples=2, rays=2, positive=2,
                    maximum=Q(6, 77), winners=((1, 5, 11),),
                    masses={-1: Q(3, 77), 1: Q(3, 77)},
                    bulks={1: Q(17, 588)}, baseline=Q(17, 882), threshold=Q(8316, 569)),
    (1, 3, 4): dict(height=12, triples=3, rays=3, positive=2,
                    maximum=Q(6, 77), winners=((1, 5, 11),),
                    masses={-1: Q(3, 77), 1: Q(3, 77)},
                    bulks={1: Q(1, 56)}, baseline=Q(1, 84), threshold=Q(792, 61)),
    (2, 3, 5): dict(height=26, triples=20, rays=22, positive=15,
                    maximum=Q(6, 77), winners=((1, 5, 11),),
                    masses={-2: Q(0), -1: Q(3, 77), 1: Q(3, 77), 2: Q(0)},
                    bulks={1: Q(1, 49), 2: Q(1, 2940)}, baseline=Q(61, 4410), threshold=Q(83160, 3109)),
    (1, 2, 9): dict(height=34, triples=22, rays=22, positive=20,
                    maximum=Q(46, 665), winners=((5, 7, 19),),
                    masses={-2: Q(8, 665), -1: Q(3, 133), 1: Q(3, 133), 2: Q(8, 665)},
                    bulks={1: Q(1, 49), 2: Q(5, 588)}, baseline=Q(17, 882), threshold=Q(143640, 4181)),
    (1, 3, 8): dict(height=34, triples=27, rays=27, positive=24,
                    maximum=Q(58, 833), winners=((5, 7, 17),),
                    masses={-2: Q(8, 833), -1: Q(3, 119), 1: Q(3, 119), 2: Q(8, 833)},
                    bulks={1: Q(53, 2352), 2: Q(5, 784)}, baseline=Q(17, 882), threshold=Q(25704, 755)),
    (1, 5, 6): dict(height=32, triples=27, rays=28, positive=24,
                    maximum=Q(58, 833), winners=((5, 7, 17),),
                    masses={-2: Q(8, 833), -1: Q(3, 119), 1: Q(3, 119), 2: Q(8, 833)},
                    bulks={1: Q(19, 980), 2: Q(1, 196)}, baseline=Q(4, 245), threshold=Q(1190, 37)),
    (2, 3, 7): dict(height=28, triples=18, rays=19, positive=14,
                    maximum=Q(6, 77), winners=((1, 5, 11),),
                    masses={-2: Q(0), -1: Q(3, 77), 1: Q(3, 77), 2: Q(0)},
                    bulks={1: Q(23, 1029), 2: Q(4, 1029)}, baseline=Q(6, 343), threshold=Q(539, 19)),
    (3, 4, 5): dict(height=33, triples=38, rays=41, positive=32,
                    maximum=Q(6, 91), winners=((1, 7, 13),),
                    masses={-2: Q(0), -1: Q(3, 91), 1: Q(3, 91), 2: Q(0)},
                    bulks={1: Q(39, 1960), 2: Q(2, 735)}, baseline=Q(19, 1260), threshold=Q(28080, 833)),
    (1, 1, 12): dict(height=31, triples=8, rays=8, positive=7,
                    maximum=Q(12, 161), winners=((1, 11, 23),),
                    masses={-2: Q(3, 161), -1: Q(3, 161), 1: Q(3, 161), 2: Q(3, 161)},
                    bulks={1: Q(3, 196), 2: Q(3, 196)}, baseline=Q(1, 49), threshold=Q(1932, 61)),
    (1, 3, 10): dict(height=31, triples=20, rays=20, positive=17,
                    maximum=Q(12, 161), winners=((1, 11, 23),),
                    masses={-2: Q(3, 161), -1: Q(3, 161), 1: Q(3, 161), 2: Q(3, 161)},
                    bulks={1: Q(9, 490), 2: Q(11, 980)}, baseline=Q(29, 1470), threshold=Q(57960, 1853)),
    (1, 4, 9): dict(height=29, triples=17, rays=19, positive=16,
                    maximum=Q(6, 77), winners=((1, 5, 11),),
                    masses={-2: Q(0), -1: Q(3, 77), 1: Q(3, 77), 2: Q(0)},
                    bulks={1: Q(71, 3528), 2: Q(11, 1176)}, baseline=Q(26, 1323), threshold=Q(6237, 212)),
    (1, 6, 7): dict(height=33, triples=23, rays=23, positive=21,
                    maximum=Q(46, 665), winners=((5, 7, 19),),
                    masses={-2: Q(8, 665), -1: Q(3, 133), 1: Q(3, 133), 2: Q(8, 665)},
                    bulks={1: Q(25, 1372), 2: Q(11, 1372)}, baseline=Q(6, 343), threshold=Q(13965, 421)),
    (3, 4, 7): dict(height=34, triples=24, rays=24, positive=21,
                    maximum=Q(8, 119), winners=((5, 11, 17),),
                    masses={-2: Q(1, 119), -1: Q(3, 119), 1: Q(3, 119), 2: Q(1, 119)},
                    bulks={1: Q(167, 8232), 2: Q(1, 168)}, baseline=Q(6, 343), threshold=Q(4998, 145)),
}


WINNER_UNIT_RAY = {
    (1, 5, 7): (2, 1, -1),
    (1, 5, 11): (1, 2, -1),
    (5, 7, 19): (1, 2, -1),
    (5, 7, 17): (2, 1, -1),
    (1, 7, 13): (1, -2, 1),
    (1, 11, 23): (1, 2, -1),
    (5, 11, 17): (1, -2, 1),
}


def residue_dichotomy_audit():
    rows = defaultdict(int)
    for w in product((1, 2), repeat=3):
        live_carriers = {
            C for C in product(range(3), repeat=3)
            if dot(w, C) % 3 == 0 and all(x != 0 for x in C)
        }
        u = tuple(pow(x, -1, 3) for x in w)
        check(live_carriers == {u, tuple((-x) % 3 for x in u)}, "two global live carrier residues")
        for c in product(range(3), repeat=3):
            if c == (0, 0, 0) or dot(c, w) % 3:
                continue
            zero_count = sum(x == 0 for x in c)
            check(zero_count in (0, 1), "only unit or one-zero relation type")
            fibre_live = {}
            for delta in range(3):
                carriers = {
                    cross(w, n) for n in product(range(3), repeat=3)
                    if dot(c, n) % 3 == delta
                }
                carriers = {tuple(x % 3 for x in C) for C in carriers}
                fibre_live[delta] = carriers & live_carriers
            if zero_count == 0:
                check(tuple(len(fibre_live[d]) for d in range(3)) == (2, 0, 0), "unit 2/0/0 live split")
            else:
                check(tuple(len(fibre_live[d]) for d in range(3)) == (0, 1, 1), "one-zero 0/1/1 live split")
            for owners in permutations((0, 1, 2)):
                n = tuple((-w[i] * owners[i]) % 3 for i in range(3))
                delta = dot(c, n) % 3
                check((delta == 0) == (zero_count == 0), "physical defect dichotomy")
            rows[(zero_count, tuple(len(fibre_live[d]) for d in range(3)))] += 1
    return tuple(sorted(rows.items()))


def main():
    if "--tripwire" in sys.argv:
        check(False, "optimization-live tripwire")

    check(patterns_through(14, 1) == PATTERNS, "complete one-zero pattern list")
    check(patterns_through(14, 0) == UNIT_PATTERNS, "complete unit comparison list")
    residue_rows = residue_dichotomy_audit()

    result_rows = []
    total_components = 0
    unit_overlap_rows = []
    proof_union = set()
    for pattern in PATTERNS:
        expected = EXPECTED[pattern]
        height = expected["height"]
        presentations = presentation_universe(pattern, height)
        brute = brute_universe(pattern, height)
        check(presentations == brute, "two-route universe " + repr(pattern))
        check(len(presentations) == expected["triples"], "triple count " + repr(pattern))
        ray_count = sum(len(rays) for rays in presentations.values())
        check(ray_count == expected["rays"], "ray count " + repr(pattern))

        defects = relation_defects(signed_rays(pattern)[0])
        expected_defects = (-1, 1) if sum(pattern) <= 8 else (-2, -1, 1, 2)
        check(defects == expected_defects, "defect states " + repr(pattern))
        bulks = {}
        for delta in sorted(set(abs(d) for d in defects)):
            convolution = bulk_convolution(pattern, delta)
            polygon = bulk_polygon(pattern, delta)
            check(convolution == polygon, "two-route cube slice " + repr((pattern, delta)))
            check(convolution == expected["bulks"][delta], "cube slice value " + repr((pattern, delta)))
            bulks[delta] = convolution
        baseline = sum(bulk_convolution(pattern, d) for d in defects) / 3
        check(baseline == expected["baseline"], "bulk baseline " + repr(pattern))
        error_numerator = Q(3 * len(defects), 7)
        threshold = error_numerator / (expected["maximum"] - baseline)
        check(threshold == expected["threshold"], "tail threshold " + repr(pattern))
        check(floor_q(threshold) == height, "finite proof height " + repr(pattern))
        check(baseline + error_numerator / (height + 1) < expected["maximum"], "strict tail exclusion " + repr(pattern))

        measures = {}
        positive = 0
        components_here = 0
        for w, rays in presentations.items():
            proof_union.add(w)
            direct, pieces = physical_dict(w)
            check(pieces == len(direct), "physical component/carrier bijection " + repr(w))
            components_here += pieces
            value = sum(direct.values(), Q(0))
            measures[w] = value
            positive += value > 0
            for c in rays:
                quotient, masses = quotient_components(w, c)
                check(quotient == direct, "quotient/physical dictionary " + repr((pattern, w, c)))
                check(sum(masses.values(), Q(0)) == value, "defect mass sum " + repr((pattern, w, c)))

        maximum = max(measures.values())
        winners = tuple(w for w in sorted(measures) if measures[w] == maximum)
        check(positive == expected["positive"], "positive count " + repr(pattern))
        check(maximum == expected["maximum"], "sharp maximum " + repr(pattern))
        check(winners == expected["winners"], "unique unordered winner " + repr(pattern))
        winner_profiles = tuple(
            (c, quotient_components(winners[0], c)[1])
            for c in presentations[winners[0]]
        )
        check(winner_profiles[0][1] == expected["masses"], "winner defect masses " + repr(pattern))
        if pattern == (1, 4, 9):
            alternate = {-2: Q(3, 77), -1: Q(0), 1: Q(0), 2: Q(3, 77)}
            check(winner_profiles == (
                ((1, -9, 4), expected["masses"]),
                ((9, -4, 1), alternate),
            ), "chart-dependent 149 winner profiles")
        else:
            check(all(profile == expected["masses"] for _, profile in winner_profiles),
                  "all winning chart profiles " + repr(pattern))

        unit_hits = tuple(
            (unit_pattern, c)
            for unit_pattern in UNIT_PATTERNS
            for c in signed_rays(unit_pattern)
            if dot(c, winners[0]) == 0
        )
        check(unit_hits == (((1, 1, 2), WINNER_UNIT_RAY[winners[0]]),), "winner unit overlap " + repr(pattern))
        unit_dictionary, unit_masses = quotient_components(winners[0], unit_hits[0][1])
        direct, _ = physical_dict(winners[0])
        check(unit_dictionary == direct, "unit/nonunit equality dictionary " + repr(pattern))
        check(tuple(unit_masses) == (0,), "unit overlap is zero-defect " + repr(pattern))
        unit_overlap_rows.append((pattern, winners[0], unit_hits[0], maximum))

        total_components += components_here
        result_rows.append((pattern, sum(pattern), defects, height, len(presentations), ray_count,
                            positive, components_here, bulks, baseline, threshold,
                            maximum, winners, winner_profiles))

    # Common-height incidence is a sidecar only; no minimality filter is used.
    common = {pattern: presentation_universe(pattern, 34) for pattern in PATTERNS}
    common_union = defaultdict(set)
    common_rays = defaultdict(set)
    for pattern, records in common.items():
        for w, rays in records.items():
            common_union[w].add(pattern)
            common_rays[w].update(rays)
    common_unit_overlap = {}
    for w in common_union:
        hits = tuple(
            p for p in UNIT_PATTERNS
            if any(dot(c, w) == 0 for c in signed_rays(p))
        )
        if hits:
            common_unit_overlap[w] = hits
    common_summary = (
        sum(len(labels) for labels in common_union.values()),
        sum(len(rays) for rays in common_rays.values()),
        len(common_union),
        len(common_unit_overlap),
        sum(len(labels) > 1 for labels in common_union.values()),
    )

    overall = max(row[11] for row in result_rows)
    overall_rows = tuple((row[0], row[12]) for row in result_rows if row[11] == overall)
    check(overall == Q(6, 77), "overall presentation-sector maximum")
    check(set(w for _, winners in overall_rows for w in winners) == {(1, 5, 11)}, "overall unique physical maximizer")

    print("status=PASS")
    print("scope=all_one_zero_mod3_full_support_primitive_relation_patterns_l1_le_14; no_minimality_filter; LRC(14)=OPEN")
    print("residue_dichotomy=unit_relations_have_live_fibre_counts_2_0_0; one_zero_relations_have_0_1_1")
    print("residue_exhaustion=" + repr(residue_rows))
    print("patterns=" + repr(PATTERNS))
    print("result_rows=" + repr(tuple(result_rows)))
    proof_positive_union = sum(sum(physical_dict(w)[0].values(), Q(0)) > 0 for w in proof_union)
    print("proof_totals=" + repr((sum(row[4] for row in result_rows), sum(row[5] for row in result_rows), total_components,
                                  len(proof_union), proof_positive_union)))
    print("winner_unit_sector_equalities=" + repr(tuple(unit_overlap_rows)))
    print("common_height34_sidecar=" + repr(common_summary))
    print("overall_sharp=" + repr((overall, overall_rows)))
    print("physical_cache=" + repr(physical_components.cache_info()))
    print("optimization_safe_checks=yes")
    print("explicit_checks=" + str(CHECKS))


if __name__ == "__main__":
    main()
