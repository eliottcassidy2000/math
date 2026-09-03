#!/usr/bin/env python3
"""Exact verifier for the THM-4393 minimal ternary-unit l1=18 LRC shell.

No repository or producer implementation is imported.  All checks use exact
Fraction arithmetic and explicit RuntimeError guards that survive python -O.
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
    check(g > 0, f"nonzero direction {v}")
    q = tuple(x // g for x in v)
    for x in q:
        if x:
            return tuple(-z for z in q) if x < 0 else q
    raise RuntimeError("zero direction")


def extended_gcd(a, b):
    old_r, r = abs(a), abs(b)
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if a < 0:
        old_s = -old_s
    if b < 0:
        old_t = -old_t
    return old_r, old_s, old_t


def bezout_vector(c):
    g12, x, y = extended_gcd(c[0], c[1])
    g, s, t = extended_gcd(g12, c[2])
    check(g == 1, f"primitive relation {c}")
    v = (s * x, s * y, t)
    check(dot(c, v) == 1, f"relation Bezout {c}, {v}")
    return v


def shell_patterns(max_norm, exact=False):
    result = []
    for a in range(1, max_norm + 1):
        for b in range(a, max_norm + 1):
            for c in range(b, max_norm + 1):
                norm = a + b + c
                if norm > max_norm:
                    break
                if exact and norm != max_norm:
                    continue
                if norm % 2:
                    continue
                if any(x % 3 == 0 for x in (a, b, c)):
                    continue
                if gcd_many((a, b, c)) != 1:
                    continue
                result.append((a, b, c))
    return tuple(result)


@lru_cache(maxsize=None)
def relation_vectors(pattern):
    """Signed coordinate placements modulo overall sign."""
    result = set()
    for magnitude_order in set(permutations(pattern)):
        for signs in product((-1, 1), repeat=3):
            c = tuple(sign * magnitude for sign, magnitude in zip(signs, magnitude_order))
            if all(x > 0 for x in c) or all(x < 0 for x in c):
                continue
            result.add(primitive_direction(c))
    return tuple(sorted(result))


def admissible_speeds(height):
    return tuple(x for x in range(1, height + 1, 2) if x % 3)


def primitive_triples(height):
    for w in combinations(admissible_speeds(height), 3):
        if gcd_many(w) == 1:
            yield w


@lru_cache(maxsize=None)
def generated_relation_map(pattern, height):
    """Generate every sorted primitive speed triple on a coefficient ray."""
    values = admissible_speeds(height)
    value_set = set(values)
    result = {}
    for c in relation_vectors(pattern):
        c1, c2, c3 = c
        for w1 in values:
            for w2 in values:
                if w2 <= w1:
                    continue
                numerator = -(c1 * w1 + c2 * w2)
                if numerator % c3:
                    continue
                w3 = numerator // c3
                if w3 <= w2 or w3 > height or w3 not in value_set:
                    continue
                w = (w1, w2, w3)
                if gcd_many(w) != 1:
                    continue
                result.setdefault(w, set()).add(c)
    return result


def direct_relations(w, patterns):
    return {
        pattern: tuple(c for c in relation_vectors(pattern) if dot(c, w) == 0)
        for pattern in patterns
        if any(dot(c, w) == 0 for c in relation_vectors(pattern))
    }


def component_length(w, C):
    w1, w2, w3 = w
    c1, c2, c3 = C
    return max(Fraction(0), min(
        2 * R / w1,
        2 * R / w2,
        2 * R / w3,
        R / w1 + R / w2 - Fraction(abs(c3), w1 * w2),
        R / w1 + R / w3 - Fraction(abs(c2), w1 * w3),
        R / w2 + R / w3 - Fraction(abs(c1), w2 * w3),
    ))


def floor_fraction(x):
    return x.numerator // x.denominator


def ceil_fraction(x):
    return -floor_fraction(-x)


def real_line_support(w, c, C0):
    lower = []
    upper = []
    for i in range(3):
        j, k = tuple(index for index in range(3) if index != i)
        bound = R * (w[j] + w[k])
        left = Fraction(-bound - C0[i], c[i])
        right = Fraction(bound - C0[i], c[i])
        if left > right:
            left, right = right, left
        lower.append(left)
        upper.append(right)
    return max(lower), min(upper)


def integer_line_support(w, c, C0):
    left, right = real_line_support(w, c, C0)
    return floor_fraction(left) + 1, ceil_fraction(right) - 1


def raw_relation_components(w, c):
    """Exact raw-carrier components, retaining (delta,k) as a sidecar."""
    v = bezout_vector(c)
    components = {}
    metadata = {}
    for delta in (-3, 0, 3):
        n0 = tuple(delta * x for x in v)
        C0 = cross(w, n0)
        check(dot(C0, w) == 0, f"affine carrier orthogonality {w}, {c}, {delta}")
        lo, hi = integer_line_support(w, c, C0)
        for k in range(lo, hi + 1):
            C = tuple(C0[i] + k * c[i] for i in range(3))
            if any(x % 3 == 0 for x in C):
                continue
            length = component_length(w, C)
            if not length:
                continue
            check(C not in components, f"unique defect/carrier address {w}, {c}, {C}")
            components[C] = length
            metadata[C] = (delta, k)
    return components, metadata


def literal_physical_components(w):
    """Definition-level three-list intersection on the physical y circle."""
    interval_lists = []
    for speed in w:
        intervals = []
        inverse = pow(speed, -1, 3)
        for n in range(speed + 1):
            left = max(Fraction(0), (Fraction(n) - R) / speed)
            right = min(Fraction(1), (Fraction(n) + R) / speed)
            if left < right:
                intervals.append((left, right, n, (-inverse * n) % 3))
        interval_lists.append(intervals)
    components = {}
    representatives = {}
    indices = [0, 0, 0]
    while all(indices[i] < len(interval_lists[i]) for i in range(3)):
        current = tuple(interval_lists[i][indices[i]] for i in range(3))
        left = max(interval[0] for interval in current)
        right = min(interval[1] for interval in current)
        n = tuple(interval[2] for interval in current)
        owners = tuple(interval[3] for interval in current)
        if left < right and len(set(owners)) == 3:
            C = cross(w, n)
            for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
                expected = w[j] * w[k] * (owners[j] - owners[k])
                check((C[i] - expected) % 3 == 0,
                      f"literal owner/carrier congruence {w}, {C}, {i}")
            check(all(x % 3 for x in C), f"literal owner gate {w}, {C}")
            if C in representatives:
                old = representatives[C]
                difference = tuple(n[i] - old[i] for i in range(3))
                check(cross(w, difference) == (0, 0, 0),
                      f"same carrier lift class {w}, {C}")
            else:
                representatives[C] = n
            components[C] = components.get(C, Fraction(0)) + right - left
        earliest_right = min(interval[1] for interval in current)
        for i in range(3):
            if current[i][1] == earliest_right:
                indices[i] += 1
    return components, representatives


def slice_integral(pattern, delta):
    """Central/offset cube-slice integral by three-box convolution."""
    target = Fraction(abs(delta)) + R * sum(pattern)
    alternating_sum = Fraction(0)
    for mask in range(8):
        removed = 2 * R * sum(pattern[i] for i in range(3) if mask & (1 << i))
        positive = max(Fraction(0), target - removed)
        alternating_sum += (-1) ** mask.bit_count() * positive * positive
    return alternating_sum / (2 * pattern[0] * pattern[1] * pattern[2])


def affine_roof_integral(w, c, delta):
    """Integrate k -> L_w(C_delta+k c) by an unrelated lower-envelope route."""
    v = bezout_vector(c)
    C0 = cross(w, tuple(delta * x for x in v))
    left, right = real_line_support(w, c, C0)
    check(left < right, f"nonempty real defect slice {w}, {c}, {delta}")
    coarse = {left, right}
    for a, b in zip(C0, c):
        kink = Fraction(-a, b)
        if left < kink < right:
            coarse.add(kink)
    coarse = sorted(coarse)
    total = Fraction(0)
    pair_data = (
        (0, R / w[1] + R / w[2], w[1] * w[2]),
        (1, R / w[0] + R / w[2], w[0] * w[2]),
        (2, R / w[0] + R / w[1], w[0] * w[1]),
    )
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
        for a, b in zip(refined, refined[1:]):
            mid = (a + b) / 2
            alpha, beta = min(lines, key=lambda line: line[0] + line[1] * mid)
            check(alpha + beta * mid > 0,
                  f"positive affine roof {w}, {c}, {delta}, {a}, {b}")
            total += alpha * (b - a) + beta * (b * b - a * a) / 2
    return total


SHORT_EXPECTED = {
    (1, 1, 2), (1, 1, 4), (1, 1, 8), (1, 1, 10), (1, 1, 14),
    (1, 2, 5), (1, 2, 7), (1, 2, 11), (1, 2, 13),
    (1, 4, 5), (1, 4, 7), (1, 4, 11),
    (1, 5, 8), (1, 5, 10), (1, 7, 8),
    (2, 5, 5), (2, 5, 7), (2, 7, 7), (4, 5, 5), (4, 5, 7),
}

SHELL18 = (
    (1, 1, 16), (1, 4, 13), (1, 7, 10),
    (2, 5, 11), (4, 7, 7), (5, 5, 8),
)

CLAIMS = {
    (1, 1, 16): {
        "integrals": (Fraction(9, 784), Fraction(9, 784)),
        "bulk": Fraction(9, 392), "threshold": Fraction(400),
        "height": 400, "all": 501, "triples": 497, "rays": 497,
        "positive": 497, "maximum": Fraction(36, 1225),
        "maximizers": ((1, 11, 175),),
        "masses": (Fraction(12, 1225), Fraction(12, 1225), Fraction(12, 1225)),
    },
    (1, 4, 13): {
        "integrals": (Fraction(9, 637), Fraction(27, 5096)),
        "bulk": Fraction(3, 182), "threshold": Fraction(3692, 11),
        "height": 335, "all": 876, "triples": 866, "rays": 867,
        "positive": 864, "maximum": Fraction(12, 497),
        "maximizers": ((5, 7, 71),),
        "masses": (Fraction(3, 497), Fraction(6, 497), Fraction(3, 497)),
    },
    (1, 7, 10): {
        "integrals": (Fraction(9, 490), Fraction(27, 6860)),
        "bulk": Fraction(6, 343), "threshold": Fraction(2891, 13),
        "height": 222, "all": 495, "triples": 485, "rays": 486,
        "positive": 485, "maximum": Fraction(12, 413),
        "maximizers": ((1, 7, 59),),
        "masses": (Fraction(3, 413), Fraction(6, 413), Fraction(3, 413)),
    },
    (2, 5, 11): {
        "integrals": (Fraction(9, 539), Fraction(9, 2695)),
        "bulk": Fraction(6, 385), "threshold": Fraction(10285, 26),
        "height": 395, "all": 1439, "triples": 1429, "rays": 1430,
        "positive": 1428, "maximum": Fraction(318, 14399),
        "maximizers": ((11, 17, 121),),
        "masses": (Fraction(57, 14399), Fraction(12, 847), Fraction(57, 14399)),
    },
    (4, 7, 7): {
        "integrals": (Fraction(54, 2401), Fraction(9, 4802)),
        "bulk": Fraction(6, 343), "threshold": Fraction(357),
        "height": 357, "all": 789, "triples": 785, "rays": 785,
        "positive": 784, "maximum": Fraction(144, 5831),
        "maximizers": ((11, 17, 49),),
        "masses": (Fraction(22, 5831), Fraction(100, 5831), Fraction(22, 5831)),
    },
    (5, 5, 8): {
        "integrals": (Fraction(27, 1225), Fraction(9, 4900)),
        "bulk": Fraction(3, 175), "threshold": Fraction(18450, 49),
        "height": 376, "all": 855, "triples": 851, "rays": 851,
        "positive": 851, "maximum": Fraction(172, 7175),
        "maximizers": ((1, 25, 41),),
        "masses": (Fraction(22, 7175), Fraction(128, 7175), Fraction(22, 7175)),
    },
}

COMMON_COUNTS = {
    (1, 1, 16): (501, 497, 497),
    (1, 4, 13): (1237, 1227, 1228),
    (1, 7, 10): (1612, 1602, 1603),
    (2, 5, 11): (1463, 1453, 1454),
    (4, 7, 7): (981, 977, 977),
    (5, 5, 8): (972, 968, 968),
}

CROSS_OVERLAPS = {
    (1, 19, 35): ((1, 1, 16), (4, 7, 7)),
    (1, 25, 41): ((1, 1, 16), (5, 5, 8)),
    (1, 37, 53): ((1, 1, 16), (1, 7, 10)),
    (5, 13, 41): ((1, 4, 13), (1, 7, 10)),
    (7, 25, 29): ((1, 4, 13), (4, 7, 7)),
    (11, 29, 43): ((1, 4, 13), (1, 7, 10)),
    (13, 19, 47): ((1, 4, 13), (2, 5, 11)),
}

WITHIN_DUPLICATES = {
    (1, 17, 55): (1, 4, 13),
    (13, 23, 31): (1, 7, 10),
    (1, 17, 37): (2, 5, 11),
}

OVERLAP_MEASURES = {
    (1, 19, 35): Fraction(64, 4655),
    (1, 25, 41): Fraction(172, 7175),
    (1, 37, 53): Fraction(300, 13727),
    (5, 13, 41): Fraction(22, 3731),
    (7, 25, 29): Fraction(0),
    (11, 29, 43): Fraction(92, 8729),
    (13, 19, 47): Fraction(14, 893),
    (1, 17, 55): Fraction(2, 385),
    (13, 23, 31): Fraction(22, 4991),
    (1, 17, 37): Fraction(8, 4403),
}


def audit_shell_and_generator():
    short = set(shell_patterns(16))
    shell18 = shell_patterns(18, exact=True)
    check(short == SHORT_EXPECTED and len(short) == 20,
          f"complete shorter shell {short}")
    check(shell18 == SHELL18, f"complete l1=18 shell {shell18}")
    check(shell_patterns(17, exact=True) == (), "odd l1 shell impossible")
    # Independent small-box brute force checks the pair-generation route.
    brute_universe = tuple(primitive_triples(79))
    for pattern in tuple(sorted(short)) + shell18:
        generated = set(generated_relation_map(pattern, 79))
        brute = {w for w in brute_universe
                 if any(dot(c, w) == 0 for c in relation_vectors(pattern))}
        check(generated == brute, f"pair generator complete at H=79 for {pattern}")
    return tuple(sorted(short)), shell18


def audit_multiple_relations(short_patterns):
    labelled = [(pattern, c) for pattern in SHELL18 for c in relation_vectors(pattern)]
    candidates = set()
    for (_, c), (_, d) in combinations(labelled, 2):
        raw = cross(c, d)
        if raw == (0, 0, 0):
            continue
        if not (all(x > 0 for x in raw) or all(x < 0 for x in raw)):
            continue
        w = primitive_direction(raw)
        w = tuple(sorted(w))
        if any(x <= 0 or x % 2 == 0 or x % 3 == 0 for x in w):
            continue
        if len(set(w)) != 3 or gcd_many(w) != 1:
            continue
        if direct_relations(w, short_patterns):
            continue
        relations = direct_relations(w, SHELL18)
        ray_count = sum(len(cs) for cs in relations.values())
        if ray_count >= 2:
            candidates.add(w)
    expected = set(CROSS_OVERLAPS) | set(WITHIN_DUPLICATES)
    check(candidates == expected, f"complete multiple l18 relation rays {candidates}")
    for w in candidates:
        relations = direct_relations(w, SHELL18)
        labels = tuple(sorted(relations))
        if w in CROSS_OVERLAPS:
            check(labels == tuple(sorted(CROSS_OVERLAPS[w])),
                  f"cross-pattern labels {w}, {labels}")
            check(sum(len(cs) for cs in relations.values()) == 2,
                  f"two cross-pattern rays {w}")
        else:
            check(labels == (WITHIN_DUPLICATES[w],), f"within label {w}, {labels}")
            check(len(relations[labels[0]]) == 2, f"two within-pattern rays {w}")
    return tuple(sorted(candidates))


def audit_atlas(short_patterns):
    common_height = 400
    short_union = set()
    for pattern in short_patterns:
        short_union.update(generated_relation_map(pattern, common_height))

    all_rows = {}
    minimal_rows = {}
    for pattern in SHELL18:
        all_rows[pattern] = generated_relation_map(pattern, common_height)
        minimal_rows[pattern] = {
            w: cs for w, cs in all_rows[pattern].items() if w not in short_union
        }
        expected_all, expected_minimal, expected_rays = COMMON_COUNTS[pattern]
        check(len(all_rows[pattern]) == expected_all,
              f"H=400 all incidence {pattern}")
        check(len(minimal_rows[pattern]) == expected_minimal,
              f"H=400 minimal incidence {pattern}")
        check(sum(len(cs) for cs in minimal_rows[pattern].values()) == expected_rays,
              f"H=400 relation-ray count {pattern}")
        for w in all_rows[pattern]:
            brute_short = bool(direct_relations(w, short_patterns))
            check(brute_short == (w in short_union), f"two-route minimality filter {w}")

    physical_union = set().union(*(set(row) for row in minimal_rows.values()))
    sector_incidences = sum(len(row) for row in minimal_rows.values())
    relation_rays = sum(sum(len(cs) for cs in row.values()) for row in minimal_rows.values())
    multiplicities = {}
    cross_overlaps = {}
    for w in physical_union:
        labels = tuple(pattern for pattern in SHELL18 if w in minimal_rows[pattern])
        multiplicities[len(labels)] = multiplicities.get(len(labels), 0) + 1
        if len(labels) > 1:
            cross_overlaps[w] = labels
    check(len(physical_union) == 6717, "H=400 distinct physical triples")
    check(sector_incidences == 6724, "H=400 sector-pattern incidences")
    check(relation_rays == 6727, "H=400 signed relation rays")
    check(multiplicities == {1: 6710, 2: 7}, f"H=400 label multiplicities {multiplicities}")
    check(cross_overlaps == CROSS_OVERLAPS, f"H=400 cross overlaps {cross_overlaps}")

    empty_union = []
    for w in physical_union:
        first_pattern = next(pattern for pattern in SHELL18 if w in minimal_rows[pattern])
        first_relation = next(iter(minimal_rows[first_pattern][w]))
        raw_components, _ = raw_relation_components(w, first_relation)
        if not raw_components:
            empty_union.append(w)
    empty_union = tuple(sorted(empty_union))
    check(empty_union == ((7, 11, 43), (7, 11, 47), (7, 25, 29)),
          f"H=400 empty physical combs {empty_union}")

    physical_cache = {}
    row_summaries = []
    for pattern in SHELL18:
        claim = CLAIMS[pattern]
        height = claim["height"]
        all_proof = {w: cs for w, cs in all_rows[pattern].items() if max(w) <= height}
        proof = {w: cs for w, cs in minimal_rows[pattern].items() if max(w) <= height}
        check(len(all_proof) == claim["all"], f"finite all-presentations universe {pattern}")
        check(len(proof) == claim["triples"], f"finite minimal universe {pattern}")
        check(sum(len(cs) for cs in proof.values()) == claim["rays"],
              f"finite relation rays {pattern}")

        I0 = slice_integral(pattern, 0)
        I3 = slice_integral(pattern, 3)
        check((I0, I3) == claim["integrals"], f"slice integrals {pattern}, {I0}, {I3}")
        bulk = Fraction(2, 3) * (I0 + 2 * I3)
        check(bulk == claim["bulk"], f"bulk constant {pattern}, {bulk}")

        values = {}
        component_count = 0
        zero_triples = []
        for w, relations in proof.items():
            if w not in physical_cache:
                physical_cache[w] = literal_physical_components(w)
            direct_components, representatives = physical_cache[w]
            direct_measure = sum(direct_components.values(), Fraction(0))
            if not direct_measure:
                zero_triples.append(w)
            component_count += len(direct_components)
            for c in relations:
                raw_components, metadata = raw_relation_components(w, c)
                check(raw_components == direct_components,
                      f"raw/literal carrier dictionary {pattern}, {w}, {c}")
                for C, n in representatives.items():
                    delta = dot(c, n)
                    check(delta in (-3, 0, 3), f"literal defect range {pattern}, {w}, {c}, {delta}")
                    check(metadata[C][0] == delta,
                          f"raw/literal defect sidecar {pattern}, {w}, {c}, {C}")
            values[w] = direct_measure

        maximum = max(values.values())
        maximizers = tuple(w for w, value in values.items() if value == maximum)
        positive = sum(value > 0 for value in values.values())
        check(maximum == claim["maximum"], f"sharp maximum {pattern}, {maximum}")
        check(maximizers == claim["maximizers"], f"unique maximizer {pattern}, {maximizers}")
        check(positive == claim["positive"], f"positive count {pattern}, {positive}")

        maximizing_w = maximizers[0]
        maximizing_c = proof[maximizing_w]
        check(len(maximizing_c) == 1, f"unique maximizing relation ray {pattern}")
        maximizing_c = next(iter(maximizing_c))
        raw_components, metadata = raw_relation_components(maximizing_w, maximizing_c)
        masses = tuple(sum((raw_components[C] for C in raw_components
                            if metadata[C][0] == delta), Fraction(0))
                       for delta in (-3, 0, 3))
        check(masses == claim["masses"], f"maximizer defect masses {pattern}, {masses}")
        for delta, expected_integral in ((0, I0), (3, I3), (-3, I3)):
            actual_integral = affine_roof_integral(maximizing_w, maximizing_c, delta)
            check(actual_integral == expected_integral,
                  f"affine/convolution integral {pattern}, {delta}, {actual_integral}")

        threshold = Fraction(18, 7) / (maximum - bulk)
        check(threshold == claim["threshold"], f"tail threshold {pattern}, {threshold}")
        check(height == floor_fraction(threshold), f"finite cutoff floor {pattern}")
        check(bulk + Fraction(18, 7 * (height + 1)) < maximum,
              f"all larger integer heights excluded {pattern}")
        check(bulk + Fraction(18, 7 * height) >= maximum,
              f"declared broad cutoff minimal {pattern}")
        row_summaries.append({
            "pattern": pattern,
            "height": height,
            "all": len(all_proof),
            "triples": len(proof),
            "rays": sum(len(cs) for cs in proof.values()),
            "positive": positive,
            "components": component_count,
            "zeros": tuple(zero_triples),
            "I0": I0,
            "I3": I3,
            "bulk": bulk,
            "threshold": threshold,
            "maximum": maximum,
            "maximizers": maximizers,
            "masses": masses,
        })

    for w, expected_measure in OVERLAP_MEASURES.items():
        if w not in physical_cache:
            physical_cache[w] = literal_physical_components(w)
        components, representatives = physical_cache[w]
        measure = sum(components.values(), Fraction(0))
        check(measure == expected_measure, f"multiple-relation measure {w}, {measure}")
        relations = direct_relations(w, SHELL18)
        flat = [c for cs in relations.values() for c in cs]
        check(len(flat) == 2, f"exactly two minimal l18 rays {w}")
        for C, n in representatives.items():
            defects = tuple(dot(c, n) for c in flat)
            check(not (defects[0] == 0 and defects[1] == 0),
                  f"independent zero-defect obstruction survives {w}, {C}")

    global_maximum = max(row["maximum"] for row in row_summaries)
    global_maximizers = tuple(
        (row["pattern"], w)
        for row in row_summaries if row["maximum"] == global_maximum
        for w in row["maximizers"]
    )
    check(global_maximum == Fraction(36, 1225), "global minimal-l18 shell maximum")
    check(global_maximizers == (((1, 1, 16), (1, 11, 175)),),
          f"global minimal-l18 shell uniqueness {global_maximizers}")

    return row_summaries, {
        "height": common_height,
        "sector_incidences": sector_incidences,
        "relation_rays": relation_rays,
        "physical_triples": len(physical_union),
        "positive_combs": len(physical_union) - len(empty_union),
        "empty_combs": empty_union,
        "multiplicities": multiplicities,
    }, global_maximum, global_maximizers


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
    short_patterns, shell18 = audit_shell_and_generator()
    multiple_rays = audit_multiple_relations(short_patterns)
    rows, incidence, global_maximum, global_maximizers = audit_atlas(short_patterns)
    semantic = {
        "short_patterns": short_patterns,
        "shell18": shell18,
        "multiple_rays": multiple_rays,
        "rows": [
            {key: (ff(value) if isinstance(value, Fraction) else value)
             for key, value in row.items()}
            for row in rows
        ],
        "incidence": incidence,
        "global_maximum": global_maximum,
        "global_maximizers": global_maximizers,
        "overlap_measures": sorted((w, ff(value)) for w, value in OVERLAP_MEASURES.items()),
    }
    digest = sha256(json.dumps(json_ready(semantic), sort_keys=True).encode()).hexdigest()

    print("LRC14 MINIMAL TERNARY-UNIT L1=18 SHELL")
    print("implementation=standalone exact Fraction; raw carriers versus literal phase cells")
    print("minimality=no primitive full-support ternary-unit relation of l1<=16")
    print("patterns=" + ",".join("".join(map(str, pattern)) for pattern in shell18))
    for row in rows:
        print(
            "row=" + "".join(map(str, row["pattern"]))
            + f" I0={ff(row['I0'])} I3={ff(row['I3'])} bulk={ff(row['bulk'])}"
            + f" threshold={ff(row['threshold'])} height={row['height']}"
            + f" all={row['all']} minimal={row['triples']} rays={row['rays']}"
            + f" positive={row['positive']} components={row['components']}"
            + f" maximum={ff(row['maximum'])}"
            + " maximizers=" + ",".join("{" + ",".join(map(str, w)) + "}"
                                          for w in row["maximizers"])
            + " masses=" + ",".join(ff(value) for value in row["masses"])
            + " zeros=" + (",".join("{" + ",".join(map(str, w)) + "}"
                                       for w in row["zeros"]) or "none")
        )
    print(
        f"incidence_H={incidence['height']}"
        + f" sector_incidences={incidence['sector_incidences']}"
        + f" relation_rays={incidence['relation_rays']}"
        + f" physical_triples={incidence['physical_triples']}"
        + f" positive_combs={incidence['positive_combs']}"
        + " multiplicities=" + ",".join(f"{k}:{v}" for k, v in sorted(incidence["multiplicities"].items()))
    )
    print("cross_pattern_overlaps=" + str(len(CROSS_OVERLAPS)))
    print("within_pattern_duplicate_rays=" + str(len(WITHIN_DUPLICATES)))
    print("multiple_relation_nonempty=" + str(sum(value > 0 for value in OVERLAP_MEASURES.values())))
    print(f"global_maximum={ff(global_maximum)} at="
          + ",".join("".join(map(str, pattern)) + ":{" + ",".join(map(str, w)) + "}"
                     for pattern, w in global_maximizers))
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
