#!/usr/bin/env python3
"""Clean-room exact audit of the THM-4387 defect-three LRC14 slice.

This program intentionally does not import the producer's verifier.  Relations
are found by scanning unordered speed triples and assigning their roles; the
definition-level check intersects owner-labelled rational intervals on the
circle.  All arithmetic is exact and all checks remain active under ``-O``.
"""

from fractions import Fraction
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import sys


Q = Fraction
RADIUS = Q(3, 14)
DEFECTS = (-3, 0, 3)
NEW_PATTERNS = ((1, 14), (2, 13), (4, 11), (5, 10), (7, 8))
OLD_PATTERNS = (
    (1, 2), (1, 4), (1, 8), (1, 10), (2, 5),
    (2, 7), (2, 11), (4, 5), (4, 7), (5, 8),
)

CLAIMS = {
    (1, 14): (Q(4, 133), (1, 5, 19), 205, 298, 149,
              ((-3, Q(2, 133)), (0, Q(0)), (3, Q(2, 133)))),
    (2, 13): (Q(228, 11165), (1, 11, 145), 409, 1295, 1295,
              ((-3, Q(48, 11165)), (0, Q(12, 1015)), (3, Q(48, 11165)))),
    (4, 11): (Q(218, 10465), (7, 13, 115), 372, 1270, 1270,
              ((-3, Q(31, 10465)), (0, Q(12, 805)), (3, Q(31, 10465)))),
    (5, 10): (Q(16, 715), (1, 11, 65), 335, 1105, 1104,
              ((-3, Q(23, 5005)), (0, Q(6, 455)), (3, Q(23, 5005)))),
    (7, 8): (Q(36, 1309), (1, 11, 85), 257, 792, 792,
              ((-3, Q(5, 1309)), (0, Q(26, 1309)), (3, Q(5, 1309)))),
}

EXPECTED_OLD_NEW = {
    ((1, 14), (1, 4)): ((1, 5, 19),),
    ((1, 14), (5, 8)): ((1, 23, 37),),
    ((2, 13), (1, 10)): ((5, 7, 43),),
    ((2, 13), (2, 5)): ((1, 5, 23),),
    ((2, 13), (4, 7)): ((5, 7, 29), (5, 17, 31)),
    ((4, 11), (1, 8)): ((5, 11, 29),),
    ((4, 11), (1, 10)): ((1, 7, 17), (5, 13, 37)),
    ((4, 11), (2, 5)): ((1, 7, 17), (1, 7, 19)),
    ((4, 11), (4, 7)): ((1, 5, 31), (5, 13, 17)),
    ((5, 10), (1, 4)): ((5, 7, 13),),
    ((5, 10), (2, 7)): ((1, 7, 25), (5, 11, 23)),
    ((5, 10), (2, 11)): ((1, 7, 25), (5, 19, 37)),
    ((5, 10), (5, 8)): ((1, 5, 35), (5, 13, 25), (11, 17, 25)),
}

EXPECTED_NEW_NEW = {
    ((1, 14), (5, 10)): ((5, 23, 47),),
    ((2, 13), (4, 11)): ((1, 17, 47),),
}


class AuditFailure(RuntimeError):
    pass


CHECK_COUNT = 0


def require(condition, label):
    """Optimization-safe replacement for assert."""
    global CHECK_COUNT
    CHECK_COUNT += 1
    if not condition:
        raise AuditFailure(label)


def ceil_q(value):
    return -((-value.numerator) // value.denominator)


def nearest_integer(value):
    shifted = value + Q(1, 2)
    return shifted.numerator // shifted.denominator


def bezout(a, b):
    """Return positive gcd and x,y with ax+by=gcd(a,b)."""
    old_r, r = a, b
    old_x, x = 1, 0
    old_y, y = 0, 1
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_x, x = x, old_x - quotient * x
        old_y, y = y, old_y - quotient * y
    if old_r < 0:
        return -old_r, -old_x, -old_y
    return old_r, old_x, old_y


def owner(speed, lift):
    return (-pow(speed, -1, 3) * lift) % 3


def section_lift(p, b, q, h, m, s, t, delta, k):
    """Choose one integral lift in the affine (delta,k) class."""
    divisor, u, v = bezout(b, p)
    require(divisor == 1, ("nonprimitive-pb", p, b, q))
    n_p = u * k
    n_b = -v * k
    n_q = t * (m * n_b - s * h * n_p - delta)

    require(b * n_p - p * n_b == k, ("section-k", p, b, delta, k))
    require(m * n_b - s * h * n_p - t * n_q == delta,
            ("section-delta", p, b, q, delta, k))
    require(q * n_p - p * n_q == t * (m * k + p * delta),
            ("endpoint-determinant", p, b, q, delta, k))
    require(q * n_b - b * n_q == t * (s * h * k + b * delta),
            ("bq-determinant", p, b, q, delta, k))

    owners = (owner(p, n_p), owner(b, n_b), owner(q, n_q))
    require((len(set(owners)) == 3) == (k % 3 != 0),
            ("owner-gate", p, b, q, delta, k, owners))
    require(sum(owners) % 3 == 0,
            ("owner-sum", p, b, q, delta, k, owners))

    # Two arbitrary translations must leave both affine coordinates fixed.
    for shift in (-2, 3):
        translated = (n_p + shift * p, n_b + shift * b, n_q + shift * q)
        require(b * translated[0] - p * translated[1] == k,
                ("translation-k", shift, p, b, q))
        require(m * translated[1] - s * h * translated[0] - t * translated[2] == delta,
                ("translation-delta", shift, p, b, q))
    return n_p, n_b, n_q, owners


def interval_roof(p, b, q, m, t, delta, k):
    """Intersection length from the three independently written intervals."""
    centers = (
        Q(0),
        -Q(k, p * b),
        -Q(t * (m * k + p * delta), p * q),
    )
    radii = (RADIUS / p, RADIUS / b, RADIUS / q)
    left = max(center - rad for center, rad in zip(centers, radii))
    right = min(center + rad for center, rad in zip(centers, radii))
    return max(Q(0), right - left)


def six_term_roof(p, b, q, h, m, s, delta, k):
    """The proposed six-term closed form, evaluated independently."""
    candidates = (
        2 * RADIUS / p,
        2 * RADIUS / b,
        2 * RADIUS / q,
        RADIUS / p + RADIUS / b - Q(abs(k), p * b),
        RADIUS / p + RADIUS / q - Q(abs(m * k + p * delta), p * q),
        RADIUS / b + RADIUS / q - Q(abs(s * h * k + b * delta), b * q),
    )
    return max(Q(0), min(candidates))


def presentation_measure(p, b, q, h, m, s, t):
    """Sum all allowed affine classes using a p-b support bound."""
    support_limit = ceil_q(RADIUS * (p + b)) + 1
    pieces = {delta: [] for delta in DEFECTS}
    for delta in DEFECTS:
        for k in range(-support_limit, support_limit + 1):
            section_lift(p, b, q, h, m, s, t, delta, k)
            geometric = interval_roof(p, b, q, m, t, delta, k)
            closed = six_term_roof(p, b, q, h, m, s, delta, k)
            require(geometric == closed,
                    ("six-term-roof", p, b, q, h, m, s, t, delta, k,
                     geometric, closed))
            reflected = interval_roof(p, b, q, m, t, -delta, -k)
            require(geometric == reflected,
                    ("roof-reflection", p, b, q, delta, k))
            if k % 3 and geometric:
                pieces[delta].append((k, geometric))

        # The p-b overlap makes every omitted k identically zero.
        for k in (-support_limit - 1, support_limit + 1):
            require(interval_roof(p, b, q, m, t, delta, k) == 0,
                    ("support-bound", p, b, q, delta, k))

    masses = tuple((delta, sum((value for _, value in pieces[delta]), Q(0)))
                   for delta in DEFECTS)
    require(dict(masses)[-3] == dict(masses)[3],
            ("defect-reflection-total", p, b, q, masses))
    return sum((mass for _, mass in masses), Q(0)), masses, pieces


def line_intersection(first, second):
    slope_a, intercept_a = first
    slope_b, intercept_b = second
    if slope_a == slope_b:
        return None
    return (intercept_b - intercept_a) / (slope_a - slope_b)


def line_value(line, x):
    return line[0] * x + line[1]


def exact_roof_integral(p, b, q, m, t, delta):
    """Integrate the piecewise-affine interval roof by all line crossings."""
    centers = (
        (Q(0), Q(0)),
        (-Q(1, p * b), Q(0)),
        (-Q(t * m, p * q), -Q(t * delta, q)),
    )
    radii = (RADIUS / p, RADIUS / b, RADIUS / q)
    lowers = tuple((slope, intercept - rad)
                   for (slope, intercept), rad in zip(centers, radii))
    uppers = tuple((slope, intercept + rad)
                   for (slope, intercept), rad in zip(centers, radii))
    lines = lowers + uppers
    crossings = sorted({point for a, b_line in combinations(lines, 2)
                        if (point := line_intersection(a, b_line)) is not None})
    require(len(crossings) >= 2, ("integral-crossings", p, b, q, delta))

    area = Q(0)
    for left, right in zip(crossings, crossings[1:]):
        mid = (left + right) / 2
        lower_line = max(lowers, key=lambda line: line_value(line, mid))
        upper_line = min(uppers, key=lambda line: line_value(line, mid))
        mid_gap = line_value(upper_line, mid) - line_value(lower_line, mid)
        if mid_gap <= 0:
            continue
        gap_left = line_value(upper_line, left) - line_value(lower_line, left)
        gap_right = line_value(upper_line, right) - line_value(lower_line, right)
        require(gap_left >= 0 and gap_right >= 0,
                ("integral-sign", p, b, q, delta, left, right))
        area += (gap_left + gap_right) * (right - left) / 2
    return area


OWNER_INTERVAL_CACHE = {}


def owner_intervals(speed):
    """Definition-level eligible intervals in [0,1], split by actual lift owner."""
    if speed in OWNER_INTERVAL_CACHE:
        return OWNER_INTERVAL_CACHE[speed]
    groups = [[], [], []]
    for lift in range(speed + 1):
        left = max(Q(0), (Q(lift) - RADIUS) / speed)
        right = min(Q(1), (Q(lift) + RADIUS) / speed)
        if left < right:
            groups[owner(speed, lift)].append((left, right))
    packed = tuple(tuple(sorted(group)) for group in groups)
    OWNER_INTERVAL_CACHE[speed] = packed
    return packed


def intersect_interval_lists(left_list, right_list):
    result = []
    i = j = 0
    while i < len(left_list) and j < len(right_list):
        left = max(left_list[i][0], right_list[j][0])
        right = min(left_list[i][1], right_list[j][1])
        if left < right:
            result.append((left, right))
        end_left = left_list[i][1]
        end_right = right_list[j][1]
        if end_left <= end_right:
            i += 1
        if end_right <= end_left:
            j += 1
    return tuple(result)


PHYSICAL_CACHE = {}


def definition_measure(speed_set):
    """Measure the owner-distinct triple comb without any relation coordinates."""
    key = tuple(sorted(speed_set))
    if key in PHYSICAL_CACHE:
        return PHYSICAL_CACHE[key]
    labelled = tuple(owner_intervals(speed) for speed in key)
    total = Q(0)
    components = 0
    for owner_assignment in permutations((0, 1, 2)):
        overlap = labelled[0][owner_assignment[0]]
        overlap = intersect_interval_lists(overlap, labelled[1][owner_assignment[1]])
        overlap = intersect_interval_lists(overlap, labelled[2][owner_assignment[2]])
        total += sum((right - left for left, right in overlap), Q(0))
        components += len(overlap)
    PHYSICAL_CACHE[key] = total, components
    return total, components


def full_boundary_shapes():
    shapes = []
    for a in range(1, 17):
        for b in range(a, 17):
            c = 16 - a - b
            if c < b:
                continue
            if c <= 0:
                continue
            if gcd(gcd(a, b), c) != 1:
                continue
            if any(value % 3 == 0 for value in (a, b, c)):
                continue
            shapes.append((a, b, c))
    return tuple(shapes)


def relation_presentations(h, m, width):
    """Independent universe: unordered triples first, then all role assignments."""
    speeds = tuple(value for value in range(1, width, 2) if value % 3)
    for speed_set in combinations(speeds, 3):
        if gcd(gcd(speed_set[0], speed_set[1]), speed_set[2]) != 1:
            continue
        for p, b, q in permutations(speed_set):
            for s, t in product((-1, 1), repeat=2):
                if m * b == s * h * p + t * q:
                    yield p, b, q, s, t


def continuous_bulk(h, m):
    return Q(6 * (h + 1), 49 * m * h)


def analytic_cutoff(h, m, target):
    width = 1
    while continuous_bulk(h, m) + Q(18, 7 * width) >= target:
        width += 1
    return width


def signed_coefficient_rows(pattern):
    h, m = pattern
    rows = set()
    for magnitudes in sorted(set(permutations((1, h, m)))):
        for signs in product((-1, 1), repeat=3):
            row = tuple(magnitude * sign for magnitude, sign in zip(magnitudes, signs))
            if row[0] < 0:
                row = tuple(-entry for entry in row)
            rows.add(row)
    return tuple(sorted(rows))


def cross_product(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def common_positive_rays(left_pattern, right_pattern):
    rays = set()
    for left in signed_coefficient_rows(left_pattern):
        for right in signed_coefficient_rows(right_pattern):
            if left == right:
                continue
            ray = cross_product(left, right)
            if 0 in ray:
                continue
            if all(entry < 0 for entry in ray):
                ray = tuple(-entry for entry in ray)
            if not all(entry > 0 for entry in ray):
                continue
            divisor = gcd(gcd(ray[0], ray[1]), ray[2])
            primitive = tuple(entry // divisor for entry in ray)
            if len(set(primitive)) != 3:
                continue
            if any(entry % 2 == 0 or entry % 3 == 0 for entry in primitive):
                continue
            rays.add(tuple(sorted(primitive)))
    return tuple(sorted(rays))


def chart_incidence_audit():
    old_new = {}
    for fresh in NEW_PATTERNS:
        for old in OLD_PATTERNS:
            rays = common_positive_rays(fresh, old)
            if rays:
                old_new[(fresh, old)] = rays
    require(old_new == EXPECTED_OLD_NEW, ("old-new-incidence", old_new))

    new_new = {}
    for index, left in enumerate(NEW_PATTERNS):
        for right in NEW_PATTERNS[index + 1:]:
            rays = common_positive_rays(left, right)
            if rays:
                new_new[(left, right)] = rays
    require(new_new == EXPECTED_NEW_NEW, ("new-new-incidence", new_new))

    # Derive, rather than assume, the address transition at the strongest ray.
    require(14 * 1 == 19 - 5, "new winner relation")
    require(4 * 5 == 1 + 19, "old winner relation")
    total, masses, pieces = presentation_measure(19, 1, 5, 1, 14, 1, -1)
    require(total == Q(4, 133) and dict(masses)[0] == 0,
            ("winner-new-chart", total, masses))
    transport = []
    for delta in (-3, 3):
        for k, length in pieces[delta]:
            old_k = -k - delta
            old_length = interval_roof(1, 5, 19, 4, 1, 0, old_k)
            require(length == old_length,
                    ("chart-length", delta, k, old_k, length, old_length))
            transport.append((delta, k, old_k, length))
    transport = tuple(sorted(transport))
    require(transport == ((-3, 4, -1, Q(2, 133)),
                          (3, -4, 1, Q(2, 133))),
            ("chart-address", transport))
    return tuple(old_new.items()), tuple(new_new.items()), transport


def endpoint_packet(p, b, q, h, m, s, t, delta, k):
    n_p, n_b, n_q, owners = section_lift(p, b, q, h, m, s, t, delta, k)
    g = gcd(p, q)
    P, Qred, M = p // g, q // g, m // g
    divisor, alpha, beta = bezout(P, Qred)
    require(divisor == 1, ("endpoint-bezout", p, q))
    determinant = q * n_p - p * n_q
    require(determinant % g == 0, ("endpoint-divisibility", p, q, determinant))
    j = determinant // g
    torsion = (alpha * n_p + beta * n_q) % g
    require(j == t * (M * k + P * delta),
            ("normalized-j", p, q, delta, k, j))
    require((j % 3 != 0) == (k % 3 != 0),
            ("normalized-owner-gate", p, q, delta, k, j))
    return j, torsion, (n_p, n_b, n_q), owners


def seam_audit():
    examples = (
        # h,m,p,b,q,s,t,states
        (2, 13, 13, 7, 65, 1, 1, (-3, 0, 3)),
        (4, 11, 11, 1, 55, -1, 1, (-3, 0, 3)),
        (1, 14, 7, 5, 77, -1, 1, (-3, 3)),
        (5, 10, 25, 13, 5, 1, 1, (-3, 3)),
    )
    output = []
    for h, m, p, b, q, s, t, states in examples:
        require(m * b == s * h * p + t * q,
                ("seam-relation", h, m, p, b, q))
        require(gcd(gcd(p, b), q) == 1, ("seam-primitivity", p, b, q))
        g = gcd(p, q)
        P, M = p // g, m // g
        packets = []
        for delta in states:
            numerator = t - P * delta
            require(numerator % M == 0,
                    ("seam-integral-k", h, m, delta, numerator, M))
            k = numerator // M
            require(k % 3 != 0, ("seam-owner-k", h, m, delta, k))
            packet = endpoint_packet(p, b, q, h, m, s, t, delta, k)
            require(packet[0] == 1, ("seam-j-one", h, m, delta, packet))
            packets.append((delta, k, packet[0], packet[1]))
        require(len({row[2] for row in packets}) == 1,
                ("seam-j-collision", h, m, packets))
        require(len({row[3] for row in packets}) == len(states),
                ("seam-torsion-separation", h, m, packets))
        output.append(((h, m), g, M, tuple(packets)))
    return tuple(output)


def hostile_witness_audit():
    witnesses = (
        (1, 14, 19, 1, 5, 1, -1, Q(211, 266), 3),
        (2, 13, 5, 1, 23, -1, 1, Q(34, 161), -3),
        (4, 11, 19, 7, 1, 1, 1, Q(213, 1862), 3),
        (5, 10, 13, 7, 5, 1, 1, Q(509, 1274), 3),
        (7, 8, 23, 17, 25, 1, -1, Q(2189, 5474), 3),
    )
    output = []
    for h, m, p, b, q, s, t, y, claimed_delta in witnesses:
        require(0 < y < 1, ("witness-circle", h, m, y))
        require(m * b == s * h * p + t * q,
                ("witness-relation", h, m, p, b, q))
        lifts = tuple(nearest_integer(speed * y) for speed in (p, b, q))
        errors = tuple(speed * y - lift for speed, lift in zip((p, b, q), lifts))
        owners = tuple(owner(speed, lift) for speed, lift in zip((p, b, q), lifts))
        delta = m * lifts[1] - s * h * lifts[0] - t * lifts[2]
        k = b * lifts[0] - p * lifts[1]
        require(all(abs(error) < RADIUS for error in errors),
                ("witness-eligibility", h, m, errors))
        require(set(owners) == {0, 1, 2},
                ("witness-owners", h, m, owners))
        require(delta == claimed_delta, ("witness-defect", h, m, delta))
        require(k % 3 != 0, ("witness-owner-address", h, m, k))
        require(interval_roof(p, b, q, m, t, delta, k) > 0,
                ("witness-component", h, m, delta, k))

        # Reflection y -> 1-y realizes the opposite defect with distinct owners.
        reflected_y = 1 - y
        reflected_lifts = tuple(nearest_integer(speed * reflected_y)
                                for speed in (p, b, q))
        reflected_delta = (m * reflected_lifts[1] - s * h * reflected_lifts[0]
                           - t * reflected_lifts[2])
        reflected_owners = tuple(owner(speed, lift)
                                 for speed, lift in zip((p, b, q), reflected_lifts))
        require(reflected_delta == -delta,
                ("witness-reflected-defect", h, m, reflected_delta))
        require(set(reflected_owners) == {0, 1, 2},
                ("witness-reflected-owners", h, m, reflected_owners))
        output.append(((h, m), (p, b, q), y, lifts, errors, owners, delta, k))
    return tuple(output)


def main():
    global CHECK_COUNT
    if "--check-tripwire" in sys.argv:
        try:
            require(False, "deliberate optimization tripwire")
        except AuditFailure:
            print("optimization_tripwire=LIVE")
            print(f"explicit_checks={CHECK_COUNT}")
            return
        raise AuditFailure("optimization removed explicit check")

    require(Path(__file__).name ==
            "lrc14_defect_three_coefficient_one_boundary_independent_referee_thm4387.py",
            "wrong canonical referee filename")

    shapes = full_boundary_shapes()
    expected_shapes = (
        (1, 1, 14), (1, 2, 13), (1, 4, 11), (1, 5, 10),
        (1, 7, 8), (2, 7, 7), (4, 5, 7),
    )
    require(shapes == expected_shapes, ("coefficient-shapes", shapes))
    included_shapes = tuple(sorted(tuple(sorted((1, h, m))) for h, m in NEW_PATTERNS))
    omitted_shapes = tuple(shape for shape in shapes if shape not in included_shapes)
    require(omitted_shapes == ((2, 7, 7), (4, 5, 7)),
            ("omitted-coefficient-shapes", omitted_shapes))
    require(all(1 not in shape for shape in omitted_shapes),
            ("omitted-not-coefficient-one", omitted_shapes))

    expected_bulk = {
        (1, 14): Q(6, 343),
        (2, 13): Q(9, 637),
        (4, 11): Q(15, 1078),
        (5, 10): Q(18, 1225),
        (7, 8): Q(6, 343),
    }

    summaries = []
    sector_keys = {}
    gcd_census = {}
    duplicate_rows = []
    duplicate_state_rows = []
    observed_defects = {}
    for h, m in NEW_PATTERNS:
        target, winner, expected_width, expected_presentations, expected_triples, expected_parts = CLAIMS[(h, m)]
        bulk = continuous_bulk(h, m)
        require(bulk == expected_bulk[(h, m)], ("bulk", h, m, bulk))
        width = analytic_cutoff(h, m, target)
        require(width == expected_width, ("cutoff", h, m, width))
        require(bulk + Q(18, 7 * width) < target,
                ("tail-strict", h, m, width))
        require(bulk + Q(18, 7 * (width - 1)) >= target,
                ("tail-minimal", h, m, width))

        # Exact continuous integrals are checked at the claimed maximizing chart.
        candidate_rows = list(relation_presentations(h, m, width))
        require(len(candidate_rows) == expected_presentations,
                ("presentation-universe", h, m, len(candidate_rows)))

        fibres = {}
        gcd_values = set()
        state_values = set()
        for p, b, q, s, t in candidate_rows:
            require(max(p, b, q) < width, ("finite-width", h, m, p, b, q))
            require(m * b == s * h * p + t * q,
                    ("relation-enumerator", h, m, p, b, q, s, t))
            require(gcd(gcd(p, b), q) == 1,
                    ("primitive-enumerator", h, m, p, b, q))
            require(gcd(p, b) == 1,
                    ("primitive-pb-consequence", h, m, p, b, q))
            g = gcd(p, q)
            gcd_values.add(g)
            require(g == gcd(p, m), ("gcd-identity", h, m, p, b, q, g))
            require(m % g == 0 and gcd(g, b) == 1,
                    ("gcd-seam", h, m, p, b, q, g))

            measure, masses, pieces = presentation_measure(p, b, q, h, m, s, t)
            W = max(p, b, q)
            require(measure <= bulk + Q(18, 7 * W),
                    ("bulk-envelope", h, m, p, b, q, measure))
            state_values.update(delta for delta in DEFECTS if dict(masses)[delta] > 0)

            speed_set = tuple(sorted((p, b, q)))
            physical, component_count = definition_measure(speed_set)
            require(measure == physical,
                    ("definition-comparison", h, m, p, b, q, measure, physical,
                     component_count))
            if speed_set in fibres:
                old_measure, old_masses, old_rows = fibres[speed_set]
                require(old_measure == measure,
                        ("multi-chart-measure", h, m, speed_set))
                old_rows.append((p, b, q, s, t))
            else:
                fibres[speed_set] = [measure, masses, [(p, b, q, s, t)]]

            if (h, m) == (5, 10) and speed_set == (5, 17, 35):
                duplicate_state_rows.append(((p, b, q, s, t), masses))

        require(len(fibres) == expected_triples,
                ("triple-universe", h, m, len(fibres)))
        maximizers = tuple(sorted(speed_set for speed_set, row in fibres.items()
                                  if row[0] == max(item[0] for item in fibres.values())))
        require(maximizers == (winner,), ("sharp-maximizer", h, m, maximizers))
        require(fibres[winner][0] == target,
                ("sharp-value", h, m, fibres[winner][0], target))
        require(fibres[winner][1] == expected_parts,
                ("sharp-parts", h, m, fibres[winner][1], expected_parts))

        # Pick its first role assignment only after the independent maximum scan.
        winner_presentation = fibres[winner][2][0]
        wp, wb, wq, ws, wt = winner_presentation
        expected_integrals = {
            -3: Q(2) * RADIUS * RADIUS / (m * h),
            0: Q(4) * RADIUS * RADIUS / m,
            3: Q(2) * RADIUS * RADIUS / (m * h),
        }
        actual_integrals = tuple(
            (delta, exact_roof_integral(wp, wb, wq, m, wt, delta))
            for delta in DEFECTS
        )
        require(dict(actual_integrals) == expected_integrals,
                ("continuous-roof-integrals", h, m, actual_integrals, expected_integrals))

        if (h, m) == (5, 10):
            duplicate_rows = sorted(fibres[(5, 17, 35)][2])

        sector_keys[(h, m)] = set(fibres)
        gcd_census[(h, m)] = tuple(sorted(gcd_values))
        observed_defects[(h, m)] = tuple(sorted(state_values))
        summaries.append(((h, m), target, winner, bulk, width,
                          len(candidate_rows), len(fibres), expected_parts,
                          actual_integrals))

    expected_gcds = {
        (1, 14): (1, 7),
        (2, 13): (1, 13),
        (4, 11): (1, 11),
        (5, 10): (1, 5),
        (7, 8): (1,),
    }
    require(gcd_census == expected_gcds, ("gcd-census", gcd_census))
    require(all(states == DEFECTS for states in observed_defects.values()),
            ("all-defect-states-observed", observed_defects))
    require(duplicate_rows == [(17, 5, 35, 1, -1), (35, 17, 5, 1, -1)],
            ("duplicate-presentation", duplicate_rows))
    expected_duplicate_states = [
        ((17, 5, 35, 1, -1),
         ((-3, Q(1, 4165)), (0, Q(16, 4165)), (3, Q(1, 4165)))),
        ((35, 17, 5, 1, -1),
         ((-3, Q(9, 4165)), (0, Q(0)), (3, Q(9, 4165)))),
    ]
    require(duplicate_state_rows == expected_duplicate_states,
            ("duplicate-state-chart", duplicate_state_rows))
    duplicate_measure = definition_measure((5, 17, 35))[0]
    require(duplicate_measure == Q(18, 4165),
            ("duplicate-measure", duplicate_measure))

    finite_intersections = {}
    for index, left in enumerate(NEW_PATTERNS):
        for right in NEW_PATTERNS[index + 1:]:
            overlap = tuple(sorted(sector_keys[left] & sector_keys[right]))
            if overlap:
                finite_intersections[(left, right)] = overlap
    require(finite_intersections == EXPECTED_NEW_NEW,
            ("finite-sector-intersections", finite_intersections))
    all_physical = set().union(*sector_keys.values())
    require(sum(len(keys) for keys in sector_keys.values()) == 4610,
            "sector-labelled triple total")
    require(len(all_physical) == 4608, ("physical-union", len(all_physical)))

    seam_rows = seam_audit()
    hostile_rows = hostile_witness_audit()
    old_new, new_new, transport = chart_incidence_audit()

    print("status=PASS")
    print("audit=clean_room_no_producer_import")
    print(f"all_l1_16_shapes={shapes}")
    print(f"included_coefficient_one_shapes={included_shapes}")
    print(f"explicitly_omitted_shapes={omitted_shapes}")
    for row in summaries:
        print("pattern=%s max=%s winner=%s bulk=%s cutoff=%s presentations=%s triples=%s parts=%s integrals=%s" % row)
    print(f"observed_defect_states={tuple(observed_defects.items())}")
    print(f"gcd_census={tuple(gcd_census.items())}")
    print(f"endpoint_torsion_seams={seam_rows}")
    print(f"hostile_boundary_witnesses={hostile_rows}")
    print(f"duplicate_fibre_5_10={tuple(duplicate_rows)} measure={duplicate_measure} state_charts={tuple(duplicate_state_rows)}")
    print(f"finite_new_sector_intersections={tuple(finite_intersections.items())}")
    print(f"old_new_relation_rays={old_new}")
    print(f"new_new_relation_rays={new_new}")
    print(f"winner_chart_transport={transport}")
    print(f"sector_labelled_triples={sum(len(keys) for keys in sector_keys.values())}")
    print(f"physical_union_triples={len(all_physical)}")
    print(f"definition_level_measures={len(PHYSICAL_CACHE)}")
    print("optimization_safe_checks=yes")
    print(f"explicit_checks={CHECK_COUNT}")
    print("scope=coefficient_one_l1_16_only; LRC(14)=OPEN")


if __name__ == "__main__":
    main()
