#!/usr/bin/env python3
"""Independent exact verifier for the norm-18 live/vanishing carrier gap.

This file imports no repository mathematics.  All geometry is computed with
fractions, and every theorem check remains active under Python optimization.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd
import json
import sys


R = Fraction(3, 14)
DEFECTS = (-3, 0, 3)
CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def scale_vector(k, a):
    return tuple(k * x for x in a)


def vector_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def primitive_patterns(norm):
    rows = []
    for a in range(1, norm):
        for b in range(a, norm):
            c = norm - a - b
            if c < b:
                continue
            if a % 3 and b % 3 and c % 3 and gcd(gcd(a, b), c) == 1:
                rows.append((a, b, c))
    return tuple(rows)


def canonical_sign(c):
    return c if c[0] > 0 else tuple(-x for x in c)


def coefficient_rays(pattern):
    rows = set()
    for magnitudes in set(permutations(pattern)):
        for signs in product((-1, 1), repeat=3):
            if len(set(signs)) == 1:
                continue
            c = tuple(signs[i] * magnitudes[i] for i in range(3))
            rows.add(canonical_sign(c))
    return tuple(sorted(rows))


def direct_unit_rays(norm):
    rows = set()
    for c0 in range(1, norm):
        for c1 in range(-norm, norm + 1):
            if c1 == 0:
                continue
            for c2 in range(-norm, norm + 1):
                if c2 == 0 or c0 + abs(c1) + abs(c2) != norm:
                    continue
                c = (c0, c1, c2)
                if any(x % 3 == 0 for x in c):
                    continue
                if gcd(gcd(abs(c0), abs(c1)), abs(c2)) != 1:
                    continue
                if not (any(x > 0 for x in c) and any(x < 0 for x in c)):
                    continue
                rows.add(c)
    return tuple(sorted(rows))


PATTERNS18 = primitive_patterns(18)
SHORT_PATTERNS = tuple(
    pattern
    for norm in range(2, 17, 2)
    for pattern in primitive_patterns(norm)
)
RAYS18_BY_PATTERN = {
    pattern: coefficient_rays(pattern) for pattern in PATTERNS18
}
RAYS18 = tuple(
    c for pattern in PATTERNS18 for c in RAYS18_BY_PATTERN[pattern]
)
SHORT_RAYS = tuple(
    c for pattern in SHORT_PATTERNS for c in coefficient_rays(pattern)
)


def pattern_of(c):
    return tuple(sorted(abs(x) for x in c))


def odd_three_units(height):
    return tuple(n for n in range(1, height + 1, 2) if n % 3)


def primitive_speed(w):
    return gcd(gcd(w[0], w[1]), w[2]) == 1


def brute_incidences(height):
    speeds = odd_three_units(height)
    all_rows = set()
    universe_count = 0
    for w in combinations(speeds, 3):
        if not primitive_speed(w):
            continue
        universe_count += 1
        for c in RAYS18:
            if dot(c, w) == 0:
                all_rows.add((w, c))
    return all_rows, universe_count


def solved_incidences(height):
    speeds = odd_three_units(height)
    speed_set = set(speeds)
    rows = set()
    for c in RAYS18:
        for w0 in speeds:
            for w1 in speeds:
                if w1 <= w0:
                    continue
                numerator = -(c[0] * w0 + c[1] * w1)
                if numerator % c[2]:
                    continue
                w2 = numerator // c[2]
                w = (w0, w1, w2)
                if w2 <= w1 or w2 not in speed_set or not primitive_speed(w):
                    continue
                rows.add((w, c))
    return rows


def is_minimal_norm18(w):
    return not any(dot(c, w) == 0 for c in SHORT_RAYS)


def egcd_nonnegative(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t


def bezout3(c):
    d1, x1, y1 = egcd_nonnegative(abs(c[0]), abs(c[1]))
    d, x2, z = egcd_nonnegative(d1, abs(c[2]))
    check(d == 1, f"primitive Bezout row {c}")
    v = (
        x2 * x1 * (1 if c[0] > 0 else -1),
        x2 * y1 * (1 if c[1] > 0 else -1),
        z * (1 if c[2] > 0 else -1),
    )
    check(dot(c, v) == 1, f"Bezout identity {c}, {v}")
    return v


def floor_fraction(x):
    return x.numerator // x.denominator


def ceil_fraction(x):
    return -((-x.numerator) // x.denominator)


def pair_margins(w, C):
    return tuple(
        3 * (w[(i + 1) % 3] + w[(i + 2) % 3]) - 14 * abs(C[i])
        for i in range(3)
    )


def integer_roof_gap(w, C):
    return -min(pair_margins(w, C))


def roof(w, C):
    candidates = [2 * R / wi for wi in w]
    candidates.extend((
        R / w[0] + R / w[1] - Fraction(abs(C[2]), w[0] * w[1]),
        R / w[0] + R / w[2] - Fraction(abs(C[1]), w[0] * w[2]),
        R / w[1] + R / w[2] - Fraction(abs(C[0]), w[1] * w[2]),
    ))
    return max(Fraction(0), min(candidates))


def fibre(w, c, delta):
    v = bezout3(c)
    base = cross(w, scale_vector(delta, v))
    check(cross(c, base) == scale_vector(delta, w),
          f"affine defect base {w}, {c}, {delta}")
    check(dot(base, w) == 0, f"base lies in carrier lattice {base}, {w}")

    forbidden = []
    allowed = []
    for k in range(3):
        residues = tuple((base[i] + k * c[i]) % 3 for i in range(3))
        if all(r == 0 for r in residues):
            forbidden.append(k)
        elif all(r != 0 for r in residues):
            allowed.append(k)
        else:
            raise RuntimeError(
                f"mixed carrier residue in unit defect fibre: {w}, {c}, {delta}, {k}, {residues}"
            )
    check(len(forbidden) == 1 and len(allowed) == 2,
          f"one deleted and two live residue classes {w}, {c}, {delta}")

    lower = None
    upper = None
    for i in range(3):
        q = R * (w[(i + 1) % 3] + w[(i + 2) % 3])
        endpoints = ((-q - base[i]) / c[i], (q - base[i]) / c[i])
        lo, hi = min(endpoints), max(endpoints)
        lower = lo if lower is None else max(lower, lo)
        upper = hi if upper is None else min(upper, hi)

    live = []
    geometric = []
    if lower < upper:
        first = floor_fraction(lower) + 1
        last = ceil_fraction(upper) - 1
        for k in range(first, last + 1):
            C = vector_add(base, scale_vector(k, c))
            value = roof(w, C)
            check(value > 0, f"strict support integer {w}, {c}, {delta}, {k}")
            geometric.append((k, C, value))
            if k % 3 in allowed:
                check(all(x % 3 for x in C), f"owner-live carrier {C}")
                live.append((k, C, value))
            else:
                check(all(x % 3 == 0 for x in C), f"deleted carrier {C}")
    return {
        "base": base,
        "allowed": tuple(allowed),
        "forbidden": tuple(forbidden),
        "interval": (lower, upper),
        "geometric": tuple(geometric),
        "live": tuple(live),
    }


def affine_components(w, c):
    result = {}
    for delta in DEFECTS:
        data = fibre(w, c, delta)
        for _, C, value in data["live"]:
            check(cross(c, C) == scale_vector(delta, w),
                  f"carrier determinant/defect {w}, {c}, {C}, {delta}")
            if C in result:
                check(result[C] == value, f"duplicate carrier length agrees {C}")
            result[C] = value
    return result


@lru_cache(maxsize=None)
def literal_components(w):
    interval_lists = []
    for speed in w:
        inverse = pow(speed, -1, 3)
        rows = []
        for nearest in range(speed + 1):
            left = max(Fraction(0), (Fraction(nearest) - R) / speed)
            right = min(Fraction(1), (Fraction(nearest) + R) / speed)
            if left < right:
                rows.append((left, right, nearest, (-inverse * nearest) % 3))
        interval_lists.append(rows)

    indices = [0, 0, 0]
    result = {}
    while all(indices[i] < len(interval_lists[i]) for i in range(3)):
        current = tuple(interval_lists[i][indices[i]] for i in range(3))
        left = max(row[0] for row in current)
        right = min(row[1] for row in current)
        if left < right and len({row[3] for row in current}) == 3:
            C = cross(w, tuple(row[2] for row in current))
            check(all(x % 3 for x in C), f"literal owner-live carrier {w}, {C}")
            result[C] = result.get(C, Fraction(0)) + right - left
        endpoint = min(row[1] for row in current)
        for i in range(3):
            if current[i][1] == endpoint:
                indices[i] += 1
    return result


def minimum_allowed_gap(w, c, delta):
    data = fibre(w, c, delta)
    seeds = []
    for k in (-2, -1, 0, 1, 2):
        if k % 3 in data["allowed"]:
            C = vector_add(data["base"], scale_vector(k, c))
            seeds.append((integer_roof_gap(w, C), k, C))
    initial = min(row[0] for row in seeds)
    i = max(range(3), key=lambda index: abs(c[index]))
    side = 3 * (w[(i + 1) % 3] + w[(i + 2) % 3])
    radius = (
        abs(data["base"][i]) + Fraction(initial + side, 14)
    ) / abs(c[i])
    bound = max(2, ceil_fraction(radius))
    rows = []
    for k in range(-bound, bound + 1):
        if k % 3 not in data["allowed"]:
            continue
        C = vector_add(data["base"], scale_vector(k, c))
        rows.append((integer_roof_gap(w, C), k, C))
    best = min(row[0] for row in rows)
    minimizers = tuple(row for row in rows if row[0] == best)
    check(all(abs(k) <= bound for _, k, _ in minimizers), "gap minimizers bounded")
    return best, minimizers


def shortest_basis_complements(w, c):
    base = fibre(w, c, 3)["base"]
    check(all(x % 3 == 0 for x in base), f"divisible defect-three base {base}")
    q0 = tuple(x // 3 for x in base)
    check(cross(c, q0) == w, f"basis complement {c}, {q0}, {w}")
    seed_norm = sum(abs(x) for x in q0)
    i = max(range(3), key=lambda index: abs(c[index]))
    bound = ceil_fraction(Fraction(abs(q0[i]) + seed_norm, abs(c[i])))
    rows = []
    for m in range(-bound, bound + 1):
        q = vector_add(q0, scale_vector(m, c))
        rows.append((sum(abs(x) for x in q), q))
    best = min(n for n, _ in rows)
    return tuple(sorted(q for n, q in rows if n == best))


def first_independent_unit_norm(w, c, maximum=40):
    for norm in range(2, maximum + 1, 2):
        found = []
        for pattern in primitive_patterns(norm):
            for d in coefficient_rays(pattern):
                if dot(d, w) == 0 and cross(c, d) != (0, 0, 0):
                    found.append(d)
        if found:
            return norm, tuple(sorted(set(found)))
    raise RuntimeError(f"no independent unit relation through {maximum}: {w}, {c}")


EXPECTED_CENTRAL_DEAD = {
    ((1, 17, 37), (2, -11, 5)): (40, -8, (1, 0, 1), Fraction(8, 4403)),
    ((1, 17, 55), (1, -13, 4)): (14, -34, (1, 0, 1), Fraction(2, 385)),
    ((1, 19, 35), (16, 1, -1)): (62, -4, (1, 0, 1), Fraction(64, 4655)),
    ((1, 25, 41), (16, 1, -1)): (26, -22, (2, 0, 2), Fraction(172, 7175)),
    ((5, 13, 41), (1, -13, 4)): (44, -22, (1, 0, 1), Fraction(22, 3731)),
    ((7, 11, 43), (5, -11, 2)): (4, 34, (0, 0, 0), Fraction(0)),
    ((7, 11, 47), (13, -4, -1)): (8, 20, (0, 0, 0), Fraction(0)),
    ((7, 25, 29), (13, 1, -4)): (20, 2, (0, 0, 0), Fraction(0)),
    ((7, 25, 29), (4, 7, -7)): (2, 20, (0, 0, 0), Fraction(0)),
    ((13, 23, 31), (1, -10, 7)): (8, -22, (1, 0, 1), Fraction(22, 4991)),
}


PATTERN_CUTOFFS = {
    (1, 1, 16): 73,
    (1, 4, 13): 59,
    (1, 7, 10): 43,
    (2, 5, 11): 49,
    (4, 7, 7): 31,
    (5, 5, 8): 37,
}


EXPECTED_PATTERN_CENSUS = {
    (1, 1, 16): (14, 14, 2),
    (1, 4, 13): (16, 17, 4),
    (1, 7, 10): (8, 9, 1),
    (2, 5, 11): (12, 13, 2),
    (4, 7, 7): (1, 1, 1),
    (5, 5, 8): (5, 5, 0),
}


def expected_family(t):
    m = 6 * t + 5
    W = 16 * m - 1
    K0 = (24 * m - 3) // 112
    A = (12 * m + 9) // 56
    U = 3 * (m - 1) // 14
    length = Fraction(3, 7 * W)
    rows = {}
    for k in range(-K0, K0 + 1):
        if k % 3:
            C = (k, -16 * k, k)
            check(C not in rows, f"family central carrier unique t={t}, k={k}")
            rows[C] = length
    for k in range(-A, U + 1):
        if k % 3:
            C = (3 * m + k, -3 - 16 * k, k)
            for signed in (C, tuple(-x for x in C)):
                check(signed not in rows, f"family off-central carrier unique t={t}, k={k}")
                rows[signed] = length
    return (1, m, W), (K0, A, U), rows


def count_nonmultiples_of_three(lo, hi):
    return sum(1 for k in range(lo, hi + 1) if k % 3)


def family_count(t):
    w, (K0, A, U), _ = expected_family(t)
    N0 = count_nonmultiples_of_three(-K0, K0)
    N3 = count_nonmultiples_of_three(-A, U)
    return w, K0, A, U, N0, N3, N0 + 2 * N3


def audit_family():
    offsets = (4, 12, 12, 20, 24, 30, 36)
    c = (1, -16, 1)
    selected = tuple(range(14)) + (20, 50, 99, 100, 101, 106)
    literal_t = tuple(range(14)) + (20, 50)
    for t in selected:
        w, (K0, A, U), expected = expected_family(t)
        m = w[1]
        check(dot(c, w) == 0, f"family relation t={t}")
        check(K0 == (3 * m) // 14, f"central endpoint floor t={t}")
        check(A == (24 * m + 21) // 112, f"negative defect endpoint floor t={t}")
        check(U == (24 * m - 21) // 112, f"positive defect endpoint floor t={t}")
        actual = affine_components(w, c)
        check(actual == expected, f"family affine formula complete t={t}")
        expected_length = Fraction(3, 7 * w[2])
        check(set(actual.values()) == {expected_length}, f"constant family roof t={t}")
        if t in literal_t:
            check(literal_components(w) == expected, f"family literal control t={t}")
        if t <= 50:
            check(is_minimal_norm18(w), f"family finite minimality control t={t}")

    for t in range(501):
        w, K0, A, U, N0, N3, N = family_count(t)
        check(N == 36 * (t // 7) + offsets[t % 7], f"family quasipolynomial t={t}")
        if t + 7 <= 500:
            later = family_count(t + 7)
            check(later[1] == K0 + 9 and later[2] == A + 9 and later[3] == U + 9,
                  f"family endpoint recurrence t={t}")
            check(later[-1] == N + 36, f"family carrier recurrence t={t}")
    return offsets


def main():
    sys.stdout.reconfigure(newline="\n")

    expected_patterns = (
        (1, 1, 16), (1, 4, 13), (1, 7, 10),
        (2, 5, 11), (4, 7, 7), (5, 5, 8),
    )
    check(PATTERNS18 == expected_patterns, f"six norm-18 patterns {PATTERNS18}")
    check(len(SHORT_PATTERNS) == 20, f"twenty shorter patterns {SHORT_PATTERNS}")
    check(tuple(sorted(RAYS18)) == direct_unit_rays(18), "direct norm-18 coefficient universe")
    check(tuple(sorted(SHORT_RAYS)) == tuple(sorted(
        c for norm in range(2, 17, 2) for c in direct_unit_rays(norm)
    )), "direct shorter coefficient universe")
    check(len(RAYS18) == 81 and len(SHORT_RAYS) == 288, "coefficient ray counts")

    brute, universe_count = brute_incidences(73)
    solved = solved_incidences(73)
    check(brute == solved, "coefficient-first/triple-first incidence equality")
    check(universe_count == 2289, f"primitive speed universe {universe_count}")
    check(len(brute) == 232 and len({w for w, _ in brute}) == 209,
          "unfiltered norm-18 relation universe")
    minimal = {(w, c) for w, c in brute if is_minimal_norm18(w)}
    check(len(minimal) == 190 and len({w for w, _ in minimal}) == 180,
          "minimal norm-18 finite core")

    physical_cache = {}
    for w, c in sorted(minimal):
        physical = physical_cache.setdefault(w, literal_components(w))
        chart = affine_components(w, c)
        check(chart == physical, f"affine/literal dictionary {w}, {c}")
        for C, value in chart.items():
            margins = pair_margins(w, C)
            check(all(margin > 0 and margin % 2 == 0 for margin in margins),
                  f"positive even nonzero roof margins {w}, {C}, {margins}")
            check(value >= Fraction(1, 7 * w[1] * w[2]),
                  f"quantized component lower bound {w}, {C}, {value}")

    central_dead = {}
    for w, c in minimal:
        if roof(w, c) > 0:
            continue
        fibres = {delta: fibre(w, c, delta) for delta in DEFECTS}
        profile = tuple(len(fibres[delta]["live"]) for delta in DEFECTS)
        central_gap = integer_roof_gap(w, c)
        off_gap, _ = minimum_allowed_gap(w, c, 3)
        mass = sum(literal_components(w).values(), Fraction(0))
        central_dead[(w, c)] = (central_gap, off_gap, profile, mass)
    check(central_dead == EXPECTED_CENTRAL_DEAD, f"complete central-dead atlas {central_dead}")

    for pattern in PATTERNS18:
        cutoff = PATTERN_CUTOFFS[pattern]
        rows = [(w, c) for w, c in minimal
                if pattern_of(c) == pattern and max(w) <= cutoff]
        triple_count = len({w for w, _ in rows})
        dead_count = sum((w, c) in central_dead for w, c in rows)
        check((triple_count, len(rows), dead_count) == EXPECTED_PATTERN_CENSUS[pattern],
              f"pattern cutoff census {pattern}")
        check(cutoff <= Fraction(14 * max(pattern), 3), f"valid pattern cutoff {pattern}")
        next_unit = cutoff + 1
        while next_unit % 2 == 0 or next_unit % 3 == 0:
            next_unit += 1
        check(next_unit > Fraction(14 * max(pattern), 3),
              f"first omitted speed beyond analytic cutoff {pattern}")

    physical_empty = tuple(sorted(w for w in {w for w, _ in minimal}
                                  if not literal_components(w)))
    expected_empty = ((7, 11, 43), (7, 11, 47), (7, 25, 29))
    check(physical_empty == expected_empty, f"complete empty physical atlas {physical_empty}")
    check(min(max(w) for w in physical_empty) == 29, "minimal empty height is 29")
    rescued = tuple(sorted({w for (w, c), row in central_dead.items() if row[3] > 0}))
    check(rescued == ((1, 17, 37), (1, 17, 55), (1, 19, 35),
                      (1, 25, 41), (5, 13, 41), (13, 23, 31)),
          f"central-dead defect-three rescues {rescued}")
    check(min(rescued, key=lambda w: (max(w), w)) == (13, 23, 31),
          "smallest rescue by height")

    empty_q = {}
    deleted_roofs = {}
    independent_norms = {}
    expected_deleted_roofs = {
        (7, 11, 43): Fraction(18, 3311),
        (7, 11, 47): Fraction(18, 2303),
        (7, 25, 29): Fraction(18, 5075),
    }
    for w in expected_empty:
        for c in sorted(c for ww, c in minimal if ww == w):
            complements = shortest_basis_complements(w, c)
            check(all(pattern_of(q) == (1, 2, 3) for q in complements),
                  f"empty ray has 123 complement {w}, {c}, {complements}")
            roof_values = tuple(roof(w, scale_vector(3, q)) for q in complements)
            check(set(roof_values) == {expected_deleted_roofs[w]},
                  f"deleted three-complement roof {w}, {c}, {roof_values}")
            for q in complements:
                deleted = scale_vector(3, q)
                check(cross(c, deleted) == scale_vector(3, w),
                      f"deleted carrier occupies defect-three fibre {w}, {c}, {q}")
                check(all(x % 3 == 0 for x in deleted),
                      f"three-complement fails owner gate {w}, {c}, {q}")
            empty_q[(w, c)] = complements
            deleted_roofs[(w, c)] = roof_values
            independent_norms[(w, c)] = first_independent_unit_norm(w, c)
    check(independent_norms[((7, 25, 29), (13, 1, -4))][0] == 18,
          "multiply-related empty next unit norm 18")
    check(independent_norms[((7, 25, 29), (4, 7, -7))][0] == 18,
          "second multiply-related chart next unit norm 18")
    check(independent_norms[((7, 11, 43), (5, -11, 2))][0] == 20,
          "single norm-18 empty next independent norm 20")
    check(independent_norms[((7, 11, 47), (13, -4, -1))][0] == 22,
          "single norm-18 empty next independent norm 22")

    dead_w, dead_c, dead_C = (7, 25, 29), (13, 1, -4), (4, 7, -7)
    check(cross(dead_c, dead_C) == scale_vector(3, dead_w), "sharp dead defect-three carrier")
    check(all(x % 3 for x in dead_C), "sharp dead carrier owner-admissible")
    check(integer_roof_gap(dead_w, dead_C) == 2 and roof(dead_w, dead_C) == 0,
          "sharp positive-side gap two")

    live_w, live_c, live_C = (29, 73, 77), (7, -7, 4), (32, 1, -13)
    check(dot(live_c, live_w) == 0 and is_minimal_norm18(live_w),
          "sharp live witness lies in minimal norm-18 shell")
    check(cross(live_c, live_C) == scale_vector(3, live_w), "sharp live defect-three carrier")
    check(integer_roof_gap(live_w, live_C) == -2, "sharp negative-side gap two")
    check(roof(live_w, live_C) == Fraction(1, 39347), "sharp live component length")
    check(roof(live_w, live_C) == Fraction(1, 7 * live_w[1] * live_w[2]),
          "sharp general component lower bound")
    check(literal_components(live_w).get(live_C) == roof(live_w, live_C),
          "sharp live literal control")

    family_offsets = audit_family()

    w18 = (7, 25, 29)
    c18a, c18b, q18 = (13, 1, -4), (4, 7, -7), (-3, 2, -1)
    check(cross(c18a, q18) == w18 and cross(c18b, q18) == w18,
          "norm-18 empty common complementary basis")
    check(vector_add(c18a, scale_vector(-1, c18b)) == scale_vector(-3, q18),
          "norm-18 unit rows differ by three times one-zero row")
    check(cross(c18a, c18b) == scale_vector(3, w18), "norm-18 index-three unit pair")

    w20 = (1, 19, 41)
    c20a, c20b, q20 = (13, -5, 2), (4, -11, 5), (3, 2, -1)
    check(cross(c20a, q20) == w20 and cross(c20b, q20) == w20,
          "norm-20 empty common complementary basis")
    check(vector_add(c20a, scale_vector(-1, c20b)) == scale_vector(3, q20),
          "norm-20 unit rows differ by three times one-zero row")
    check(cross(c20a, c20b) == scale_vector(-3, w20), "norm-20 index-three unit pair")
    check(not any(dot(c, w20) == 0 for c in SHORT_RAYS + RAYS18),
          "norm-20 comparison has no unit relation through norm 18")
    check(affine_components(w20, c20a) == {} and affine_components(w20, c20b) == {},
          "norm-20 comparison empty in both charts")
    check(literal_components(w20) == {}, "norm-20 comparison literal empty")
    check((integer_roof_gap(w20, c20a), integer_roof_gap(w20, c20b)) == (2, 28),
          "norm-20 central dead gaps")
    check(roof(w20, scale_vector(3, q20)) == Fraction(27, 5453),
          "norm-20 deleted geometric carrier is positive")

    hybrid_w, hybrid_C = (11, 13, 17), (-5, -1, 4)
    check(integer_roof_gap(hybrid_w, hybrid_C) == -16, "THM-4396 carrier gap control")
    check(roof(hybrid_w, hybrid_C) == Fraction(10, 1547),
          "THM-4396 component control")
    check(sum(literal_components(hybrid_w).values(), Fraction(0)) == Fraction(20, 1547),
          "THM-4396 physical mass control")

    semantic = {
        "patterns18": PATTERNS18,
        "short_pattern_count": len(SHORT_PATTERNS),
        "coefficient_counts": (len(RAYS18), len(SHORT_RAYS)),
        "finite_core": (universe_count, len(brute), len(minimal), len({w for w, _ in minimal})),
        "pattern_cutoffs": tuple(sorted(PATTERN_CUTOFFS.items())),
        "pattern_census": tuple(sorted(EXPECTED_PATTERN_CENSUS.items())),
        "central_dead": tuple(sorted(
            (w, c, row) for (w, c), row in EXPECTED_CENTRAL_DEAD.items()
        )),
        "empty": expected_empty,
        "rescued": rescued,
        "empty_q": tuple(sorted((w, c, rows) for (w, c), rows in empty_q.items())),
        "deleted_roofs": tuple(sorted(
            (w, c, tuple(str(value) for value in rows))
            for (w, c), rows in deleted_roofs.items()
        )),
        "independent_norms": tuple(sorted(
            (w, c, row) for (w, c), row in independent_norms.items()
        )),
        "sharp_gaps": {"dead": 2, "live": -2, "length": "1/39347"},
        "family_offsets": family_offsets,
        "norm20": {"w": w20, "central_gaps": (2, 28), "mass": "0"},
        "hybrid_control": (hybrid_w, hybrid_C, "10/1547", "20/1547"),
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, default=str).encode()).hexdigest()

    print("LRC14 NORM-18 VANISHING/LIVE RAW-CARRIER GAP -- INDEPENDENT AUDIT")
    print("status=PASS exact_local_shell_theorem; LRC14=OPEN")
    print("patterns=1-1-16,1-4-13,1-7-10,2-5-11,4-7-7,5-5-8")
    print("carrier_self_duality=C_in_w_perp; c_cross_C=delta*w; delta=-3,0,3")
    print("roof_gate=all_i(3*(w_j+w_k)-14*abs(C_i)>0)")
    print("owner_gap=no_zero_margin; live_D<=-2; dead_D>=2; sharp=yes")
    print("live_component_floor=1/(7*w_2*w_3)_for_sorted_w; sharp_at=(29,73,77)")
    print("central_failure_cutoffs=116:73,1413:59,1710:43,2511:49,477:31,558:37")
    print("finite_core=2289_speed_triples,232_raw_incidences,190_minimal_incidences,180_minimal_triples")
    print("central_dead_presentations=10; rescued_by_defect_three=6; dead_presentations=4")
    print("physical_empty={(7,11,43),(7,11,47),(7,25,29)}; minimal_empty_height=29")
    print("empty_next_independent_unit_norms=(7,25,29):18;(7,11,43):20;(7,11,47):22")
    print("family=w=(1,m,16m-1),m=6t+5,t>=0; all_components=3/[7*(16m-1)]")
    print("family_ranges=K0=floor((24m-3)/112),A=floor((12m+9)/56),U=floor(3*(m-1)/14)")
    print("family_N_offsets_mod7=4,12,12,20,24,30,36; N(t+7)=N(t)+36; limit=9/392")
    print("norm20_empty_comparison=(1,19,41); same_index3_plus_123_complement; gaps=2,28")
    print("composition=THM4396_control_only; THM4398_complementary_123_bridge_only; no_LRC14_step")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
