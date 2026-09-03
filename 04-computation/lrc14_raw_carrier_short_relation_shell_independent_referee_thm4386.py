#!/usr/bin/env python3
"""Clean-room exact audit for the THM-4386 carrier theorem.

This file deliberately imports no repository mathematics and shares no code
with the producer packet.  All rational arithmetic is exact.  Every check is
an explicit exception and therefore remains active under ``python -O``.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd
from functools import lru_cache
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


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def primitive_direction(v):
    g = gcd_many(v)
    check(g > 0, f"nonzero primitive direction {v}")
    q = tuple(x // g for x in v)
    for x in q:
        if x:
            return tuple(-z for z in q) if x < 0 else q
    raise RuntimeError("unreachable zero direction")


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


def bezout_vector(w):
    g12, a, b = extended_gcd(w[0], w[1])
    g, c, d = extended_gcd(g12, w[2])
    check(g == 1, f"primitive Bezout input {w}")
    u = (c * a, c * b, d)
    check(dot(u, w) == 1, f"Bezout identity {w}, {u}")
    return u


def strict_integer_bound(x):
    """Largest nonnegative integer k satisfying k < positive Fraction x."""
    check(x > 0, f"positive strict bound {x}")
    return (x.numerator - 1) // x.denominator


def component_length(w, C):
    w1, w2, w3 = w
    c1, c2, c3 = C
    candidates = (
        2 * R / w1,
        2 * R / w2,
        2 * R / w3,
        R / w1 + R / w2 - Fraction(abs(c3), w1 * w2),
        R / w1 + R / w3 - Fraction(abs(c2), w1 * w3),
        R / w2 + R / w3 - Fraction(abs(c1), w2 * w3),
    )
    return max(Fraction(0), min(candidates))


def lattice_components(w):
    """Enumerate Lambda_w directly from C.w=0 inside the strict support box."""
    b1 = strict_integer_bound(R * (w[1] + w[2]))
    b2 = strict_integer_bound(R * (w[0] + w[2]))
    b3 = strict_integer_bound(R * (w[0] + w[1]))
    result = {}
    for c1 in range(-b1, b1 + 1):
        for c2 in range(-b2, b2 + 1):
            numerator = -(w[0] * c1 + w[1] * c2)
            if numerator % w[2]:
                continue
            c3 = numerator // w[2]
            if abs(c3) > b3:
                continue
            C = (c1, c2, c3)
            if any(x % 3 == 0 for x in C):
                continue
            length = component_length(w, C)
            if length:
                result[C] = length
    return result


def nearest_integer(x):
    return (2 * x.numerator + x.denominator) // (2 * x.denominator)


def physical_components(w):
    """Literal endpoint-cell implementation on y in [0,1)."""
    endpoints = {Fraction(0), Fraction(1)}
    for speed in w:
        for n in range(speed + 1):
            for sign in (-1, 1):
                endpoint = (Fraction(n) + sign * R) / speed
                if 0 < endpoint < 1:
                    endpoints.add(endpoint)
    endpoints = sorted(endpoints)
    result = {}
    for left, right in zip(endpoints, endpoints[1:]):
        y = (left + right) / 2
        n = tuple(nearest_integer(speed * y) for speed in w)
        errors = tuple(speed * y - integer for speed, integer in zip(w, n))
        if not all(abs(error) < R for error in errors):
            continue
        owners = tuple((-pow(speed, -1, 3) * integer) % 3
                       for speed, integer in zip(w, n))
        if len(set(owners)) != 3:
            continue
        C = cross(w, n)
        result[C] = result.get(C, Fraction(0)) + right - left
    return result


@lru_cache(maxsize=None)
def relation_vectors(pattern):
    """All signed coordinate placements, modulo simultaneous negation."""
    vectors = set()
    for magnitude_order in set(permutations(pattern)):
        for signs in product((-1, 1), repeat=3):
            v = tuple(s * a for s, a in zip(signs, magnitude_order))
            if all(x > 0 for x in v) or all(x < 0 for x in v):
                continue
            vectors.add(primitive_direction(v))
    return tuple(sorted(vectors))


def relations_for(w, pattern):
    return tuple(c for c in relation_vectors(pattern) if dot(c, w) == 0)


def admissible_speeds(height):
    return tuple(x for x in range(1, height + 1, 2) if x % 3)


def primitive_triples(height):
    for w in combinations(admissible_speeds(height), 3):
        if gcd_many(w) == 1:
            yield w


def ray_measure(w, c):
    """Sum the nonzero ternary-unit multiples of one primitive carrier ray."""
    coordinate_bounds = (
        R * (w[1] + w[2]) / abs(c[0]),
        R * (w[0] + w[2]) / abs(c[1]),
        R * (w[0] + w[1]) / abs(c[2]),
    )
    kmax = min(strict_integer_bound(bound) for bound in coordinate_bounds)
    total = Fraction(0)
    for k in range(1, kmax + 1):
        if k % 3:
            total += 2 * component_length(w, tuple(k * x for x in c))
    return total


def box_spline_integral(pattern):
    """Three-box convolution evaluation of the central slice volume."""
    a, b, c = pattern
    values = (a, b, c)
    alternating_sum = 0
    for mask in range(8):
        removed = 2 * sum(values[i] for i in range(3) if mask & (1 << i))
        positive = max(0, sum(values) - removed)
        alternating_sum += (-1) ** mask.bit_count() * positive * positive
    return R * R * alternating_sum / (2 * a * b * c)


def lower_envelope_integral(w, c):
    """Independently integrate min_j(alpha_j-beta_j*t)_+ on t>=0."""
    w1, w2, w3 = w
    c1, c2, c3 = map(abs, c)
    lines = (
        (2 * R / w1, Fraction(0)),
        (2 * R / w2, Fraction(0)),
        (2 * R / w3, Fraction(0)),
        (R / w1 + R / w2, Fraction(c3, w1 * w2)),
        (R / w1 + R / w3, Fraction(c2, w1 * w3)),
        (R / w2 + R / w3, Fraction(c1, w2 * w3)),
    )
    stop = min(alpha / beta for alpha, beta in lines if beta)
    breaks = {Fraction(0), stop}
    for i, (a1, b1) in enumerate(lines):
        for a2, b2 in lines[i + 1:]:
            if b1 == b2:
                continue
            t = (a1 - a2) / (b1 - b2)
            if 0 < t < stop:
                breaks.add(t)
    breaks = sorted(breaks)
    half = Fraction(0)
    for left, right in zip(breaks, breaks[1:]):
        midpoint = (left + right) / 2
        alpha, beta = min(lines, key=lambda line: line[0] - line[1] * midpoint)
        check(alpha - beta * midpoint >= 0,
              f"positive lower envelope {w}, {c}, {left}, {right}")
        half += alpha * (right - left) - beta * (right * right - left * left) / 2
    return 2 * half


def shell_patterns(max_norm):
    found = []
    for a in range(1, max_norm + 1):
        for b in range(a, max_norm + 1):
            for c in range(b, max_norm + 1):
                if a + b + c > max_norm:
                    break
                if any(x % 3 == 0 for x in (a, b, c)):
                    continue
                if gcd_many((a, b, c)) != 1:
                    continue
                if (a + b + c) % 2:
                    continue
                found.append((a, b, c))
    return tuple(found)


def audit_constructive_isomorphism_and_lattice_formula():
    checked_triples = 0
    checked_components = 0
    # These six are replayed after the independent H=79 universe so the
    # aggregate can be compared exactly with the producer transcript.
    high_controls = ((1, 11, 121), (1, 11, 175), (5, 53, 55),
                     (35, 55, 77), (67, 131, 199), (191, 193, 199))
    universe = list(primitive_triples(79)) + list(high_controls)
    check(len(list(primitive_triples(79))) == 2910, "H=79 universe size")
    for w in universe:
        direct = physical_components(w)
        lattice = lattice_components(w)
        check(direct == lattice, f"literal/lattice carrier dictionary {w}")
        u = bezout_vector(w)
        for C in lattice:
            check(dot(C, w) == 0, f"carrier orthogonality {w}, {C}")
            n = cross(C, u)
            check(cross(w, n) == C, f"constructive carrier inverse {w}, {C}")
            residues = tuple(C[i] % 3 for i in range(3))
            check(all(residues), f"carrier owner gate residue {w}, {C}")
        checked_triples += 1
        checked_components += len(direct)
    return checked_triples, checked_components


def audit_shell_enumeration():
    expected = (
        (1, 1, 2),
        (1, 1, 4),
        (1, 2, 5),
        (1, 1, 8), (1, 2, 7), (1, 4, 5),
        (1, 1, 10), (1, 4, 7), (2, 5, 5),
        (1, 2, 11), (1, 5, 8), (2, 5, 7), (4, 5, 5),
    )
    actual = shell_patterns(14)
    check(set(actual) == set(expected) and len(actual) == len(expected),
          f"complete l1<=14 shell: {actual}")
    return actual


def sector_triples(pattern, height):
    return tuple(w for w in primitive_triples(height) if relations_for(w, pattern))


def audit_missing_sharp_sectors():
    claimed = {
        (2, 5, 5): (Fraction(564, 20405), ((5, 53, 55),),
                    Fraction(204050, 1333), 208, Fraction(81, 2450)),
        (2, 5, 7): (Fraction(12, 539), ((7, 17, 77), (7, 53, 77)),
                    Fraction(539, 3), 442, Fraction(9, 343)),
        (4, 5, 5): (Fraction(444, 18179), ((5, 49, 53),),
                    Fraction(64925, 366), 252, Fraction(36, 1225)),
    }
    rows = []
    direct_cache = {}
    for pattern, (expected_max, expected_argmax, expected_threshold,
                  expected_count, expected_integral) in claimed.items():
        discovery = sector_triples(pattern, 300)
        ray_values = {}
        for w in discovery:
            relations = relations_for(w, pattern)
            value = ray_measure(w, relations[0])
            for relation in relations[1:]:
                check(ray_measure(w, relation) == value,
                      f"presentation-independent ray measure {pattern}, {w}")
            ray_values[w] = value
        discovered_max = max(ray_values.values())
        discovered_argmax = tuple(w for w, value in ray_values.items()
                                  if value == discovered_max)
        check(discovered_max == expected_max,
              f"discovery maximum {pattern}: {discovered_max}")
        check(discovered_argmax == expected_argmax,
              f"discovery maximizers {pattern}: {discovered_argmax}")

        integral = box_spline_integral(pattern)
        check(integral == expected_integral,
              f"central slice integral {pattern}: {integral}")
        baseline = 2 * integral / 3
        threshold = Fraction(6, 7) / (discovered_max - baseline)
        check(threshold == expected_threshold,
              f"analytic cutoff {pattern}: {threshold}")
        proof_height = threshold.numerator // threshold.denominator
        proof_universe = sector_triples(pattern, proof_height)
        check(len(proof_universe) == expected_count,
              f"finite proof universe {pattern}: {len(proof_universe)}")

        values = {}
        component_count = 0
        for w in proof_universe:
            relations = relations_for(w, pattern)
            if w not in direct_cache:
                direct_cache[w] = physical_components(w)
            direct_components = direct_cache[w]
            direct_value = sum(direct_components.values(), Fraction(0))
            for relation in relations:
                check(lower_envelope_integral(w, relation) == integral,
                      f"two-route integral {pattern}, {w}, {relation}")
                check(ray_measure(w, relation) == direct_value,
                      f"literal/ray measure {pattern}, {w}, {relation}")
                for C in direct_components:
                    check(cross(C, relation) == (0, 0, 0),
                          f"one-ray carrier support {pattern}, {w}, {C}")
                    multiplier = next(Fraction(C[i], relation[i])
                                      for i in range(3) if relation[i])
                    check(multiplier.denominator == 1 and multiplier.numerator % 3,
                          f"integer ternary-unit ray address {pattern}, {w}, {C}")
            values[w] = direct_value
            component_count += len(direct_components)
        exact_max = max(values.values())
        exact_argmax = tuple(w for w, value in values.items() if value == exact_max)
        check(exact_max == expected_max, f"certified maximum {pattern}")
        check(exact_argmax == expected_argmax,
              f"certified maximizers {pattern}: {exact_argmax}")
        check(baseline + Fraction(6, 7 * (proof_height + 1)) < exact_max,
              f"cutoff excludes all integer heights above proof window {pattern}")
        rows.append({
            "pattern": pattern,
            "integral": integral,
            "baseline": baseline,
            "threshold": threshold,
            "proof_height": proof_height,
            "triples": len(proof_universe),
            "components": component_count,
            "maximum": exact_max,
            "maximizers": exact_argmax,
        })
    return rows


def audit_double_relation_obstruction():
    ten_patterns = (
        (1, 1, 2), (1, 1, 4), (1, 1, 8), (1, 1, 10),
        (1, 2, 5), (1, 2, 7), (1, 2, 11),
        (1, 4, 5), (1, 4, 7), (1, 5, 8),
    )
    labelled_vectors = []
    for pattern in ten_patterns:
        for c in relation_vectors(pattern):
            labelled_vectors.append((pattern, c))
    rays = {}
    for (pattern1, c), (pattern2, d) in combinations(labelled_vectors, 2):
        if cross(c, d) == (0, 0, 0):
            continue
        raw_w = cross(c, d)
        if not (all(x > 0 for x in raw_w) or all(x < 0 for x in raw_w)):
            continue
        w = primitive_direction(raw_w)
        if any(x <= 0 or x % 2 == 0 or x % 3 == 0 for x in w):
            continue
        if len(set(w)) != 3:
            continue
        speed_set = tuple(sorted(w))
        rays.setdefault(speed_set, set()).update((pattern1, pattern2))
    expected = {
        (1, 5, 13), (1, 5, 17), (1, 7, 11), (1, 7, 17), (1, 7, 25),
        (1, 11, 19), (1, 13, 23), (5, 7, 11), (5, 7, 31),
        (5, 7, 41), (7, 11, 19), (7, 19, 29), (11, 13, 23),
    }
    check(set(rays) == expected, f"double-relation rays {set(rays)}")
    for w in sorted(rays):
        direct = physical_components(w)
        check(not direct, f"double-relation physical emptiness {w}")
    return tuple(sorted(rays))


def audit_h199_incidence():
    labels = (
        (1, 2), (1, 4), (1, 8), (1, 10),
        (2, 5), (2, 7), (2, 11), (4, 5), (4, 7), (5, 8),
    )
    patterns = {label: tuple(sorted((1,) + label)) for label in labels}
    per_label = {label: 0 for label in labels}
    multiplicities = {}
    universe = 0
    covered = 0
    first_uncovered_height = None
    first_uncovered = []
    for w in primitive_triples(199):
        universe += 1
        incident = []
        for label in labels:
            if relations_for(w, patterns[label]):
                incident.append(label)
                per_label[label] += 1
        multiplicities[len(incident)] = multiplicities.get(len(incident), 0) + 1
        if incident:
            covered += 1
        else:
            height = max(w)
            if first_uncovered_height is None or height < first_uncovered_height:
                first_uncovered_height = height
                first_uncovered = [w]
            elif height == first_uncovered_height:
                first_uncovered.append(w)
    check(universe == 47499, f"H=199 incidence universe {universe}")
    check(covered == 5674, f"H=199 ten-sector union {covered}")
    check(multiplicities == {0: 41825, 1: 5663, 2: 10, 3: 1},
          f"H=199 incidence multiplicities {multiplicities}")
    expected_per_label = {
        (1, 2): 1023, (1, 4): 515, (1, 8): 257, (1, 10): 203,
        (2, 5): 825, (2, 7): 588, (2, 11): 370,
        (4, 5): 800, (4, 7): 589, (5, 8): 516,
    }
    check(per_label == expected_per_label, f"H=199 label counts {per_label}")
    expected_first = ((5, 19, 23), (7, 19, 23), (11, 19, 23), (17, 19, 23))
    check(first_uncovered_height == 23, "first uncovered height")
    check(tuple(first_uncovered) == expected_first,
          f"first uncovered triples {first_uncovered}")
    return universe, covered, multiplicities, per_label, tuple(first_uncovered)


def format_fraction(value):
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main():
    shell = audit_shell_enumeration()
    carrier_triples, carrier_components = audit_constructive_isomorphism_and_lattice_formula()
    sharp_rows = audit_missing_sharp_sectors()
    overlap_rays = audit_double_relation_obstruction()
    universe, covered, multiplicities, per_label, first_uncovered = audit_h199_incidence()

    semantic = {
        "shell": shell,
        "carrier_test": [carrier_triples, carrier_components],
        "sharp": [
            {
                "pattern": row["pattern"],
                "integral": format_fraction(row["integral"]),
                "baseline": format_fraction(row["baseline"]),
                "threshold": format_fraction(row["threshold"]),
                "proof_height": row["proof_height"],
                "triples": row["triples"],
                "components": row["components"],
                "maximum": format_fraction(row["maximum"]),
                "maximizers": row["maximizers"],
            }
            for row in sharp_rows
        ],
        "overlap_rays": overlap_rays,
        "incidence": {
            "universe": universe,
            "covered": covered,
            "multiplicities": sorted(multiplicities.items()),
            "per_label": sorted((str(k), v) for k, v in per_label.items()),
            "first_uncovered": first_uncovered,
        },
    }
    digest = sha256(json.dumps(semantic, sort_keys=True).encode()).hexdigest()

    print("THM4386 CLEAN-ROOM CARRIER AUDIT")
    print("implementation=standalone exact Fraction; no repository imports; explicit live checks")
    print("constructive_inverse=n=C_cross_u with u_dot_w=1")
    print(f"shell_count={len(shell)}")
    print("shell=" + ",".join("".join(map(str, pattern)) for pattern in shell))
    print(f"literal_lattice_triples={carrier_triples}")
    print(f"literal_lattice_components={carrier_components}")
    for row in sharp_rows:
        print(
            "sharp=" + "".join(map(str, row["pattern"]))
            + f" integral={format_fraction(row['integral'])}"
            + f" baseline={format_fraction(row['baseline'])}"
            + f" cutoff={format_fraction(row['threshold'])}"
            + f" proof_height={row['proof_height']}"
            + f" triples={row['triples']}"
            + f" components={row['components']}"
            + f" maximum={format_fraction(row['maximum'])}"
            + " maximizers=" + ",".join("{" + ",".join(map(str, w)) + "}"
                                          for w in row["maximizers"])
        )
    print(f"double_relation_rays={len(overlap_rays)}")
    print(f"h199_universe={universe}")
    print(f"h199_covered={covered}")
    print("h199_multiplicities=" + ",".join(f"{k}:{v}" for k, v in sorted(multiplicities.items())))
    print("h199_first_uncovered=" + ",".join("{" + ",".join(map(str, w)) + "}"
                                               for w in first_uncovered))
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
