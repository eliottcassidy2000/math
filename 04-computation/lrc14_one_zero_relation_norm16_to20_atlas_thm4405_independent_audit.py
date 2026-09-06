#!/usr/bin/env python3
"""Small dependency-free clean-room audit of one-zero relation norms 16--20.

Nothing is imported from repository mathematics or from the producer packet.
The verifier reconstructs coefficient shapes, signed labelled relation rays,
two finite presentation universes, exact cube slices, literal physical-circle
components, and affine carrier roofs from their definitions.  All validation
uses explicit checks which remain live under ``python -O``.
"""

from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
from json import dumps
from math import gcd
import sys


RADIUS = Q(3, 14)
TARGET = Q(6, 77)
NORMS = (16, 18, 20)
SAFE_HEIGHT = 65
CHECKS = 0


class VerificationError(RuntimeError):
    pass


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise VerificationError(message)


def guard(condition, message):
    if not condition:
        raise VerificationError(message)


def ff(value):
    return f"{value.numerator}/{value.denominator}"


def dot(first, second):
    return sum(a * b for a, b in zip(first, second))


def cross(first, second):
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def primitive_triple(values):
    return gcd(gcd(abs(values[0]), abs(values[1])), abs(values[2])) == 1


def coefficient_shapes():
    rows = []
    for norm in NORMS:
        for first in range(1, norm + 1):
            for second in range(first, norm + 1):
                third = norm - first - second
                if third < second:
                    continue
                shape = (first, second, third)
                if primitive_triple(shape) and sum(value % 3 == 0 for value in shape) == 1:
                    rows.append(shape)
    return tuple(rows)


def normalized_ray(vector):
    guard(all(vector), ("full-support-ray", vector))
    if vector[0] < 0:
        vector = tuple(-value for value in vector)
    return tuple(vector)


@lru_cache(maxsize=None)
def signed_labelled_rays(shape):
    rays = set()
    for placement in set(permutations(shape)):
        for signs in product((-1, 1), repeat=3):
            ray = normalized_ray(tuple(sign * value for sign, value in zip(signs, placement)))
            if min(ray) < 0 < max(ray):
                rays.add(ray)
    return tuple(sorted(rays))


def eligible_speed(value):
    return value > 0 and value % 2 == 1 and value % 3 != 0


SPEEDS = tuple(value for value in range(1, SAFE_HEIGHT + 1) if eligible_speed(value))
TRIPLES = tuple(
    triple for triple in combinations(SPEEDS, 3)
    if primitive_triple(triple)
)


def brute_presentations(shape):
    rays = signed_labelled_rays(shape)
    return {
        w: tuple(ray for ray in rays if dot(ray, w) == 0)
        for w in TRIPLES
        if any(dot(ray, w) == 0 for ray in rays)
    }


def solved_presentations(shape):
    rows = {}
    for ray in signed_labelled_rays(shape):
        for first in SPEEDS:
            for second in SPEEDS:
                if second <= first:
                    continue
                numerator = -(ray[0] * first + ray[1] * second)
                if numerator % ray[2]:
                    continue
                third = numerator // ray[2]
                w = (first, second, third)
                if (
                    third <= second
                    or third > SAFE_HEIGHT
                    or not eligible_speed(third)
                    or not primitive_triple(w)
                ):
                    continue
                rows.setdefault(w, set()).add(ray)
    return {w: tuple(sorted(rays)) for w, rays in sorted(rows.items())}


def clip_polygon(polygon, coefficient_x, coefficient_y, bound):
    result = []
    for index, current in enumerate(polygon):
        previous = polygon[index - 1]
        previous_value = coefficient_x * previous[0] + coefficient_y * previous[1] - bound
        current_value = coefficient_x * current[0] + coefficient_y * current[1] - bound
        previous_inside = previous_value <= 0
        current_inside = current_value <= 0
        if previous_inside != current_inside:
            direction = (current[0] - previous[0], current[1] - previous[1])
            denominator = coefficient_x * direction[0] + coefficient_y * direction[1]
            guard(denominator != 0, ("clip-crossing-denominator", polygon, bound))
            parameter = -previous_value / denominator
            result.append((
                previous[0] + parameter * direction[0],
                previous[1] + parameter * direction[1],
            ))
        if current_inside:
            result.append(current)
    return tuple(result)


def polygon_area(polygon):
    if len(polygon) < 3:
        return Q(0)
    twice = sum(
        polygon[index - 1][0] * polygon[index][1]
        - polygon[index][0] * polygon[index - 1][1]
        for index in range(len(polygon))
    )
    return abs(twice) / 2


def slice_by_clipping(shape, defect):
    first, second, third = shape
    square = (
        (-RADIUS, -RADIUS),
        (RADIUS, -RADIUS),
        (RADIUS, RADIUS),
        (-RADIUS, RADIUS),
    )
    # |defect + first*x + second*y| <= third*RADIUS.
    polygon = clip_polygon(square, first, second, third * RADIUS - defect)
    polygon = clip_polygon(polygon, -first, -second, third * RADIUS + defect)
    return polygon_area(polygon) / third


def positive_part(value):
    return max(Q(0), value)


def slice_by_convolution(shape, defect):
    norm = sum(shape)
    total = Q(0)
    for mask in range(8):
        subset_sum = sum(shape[index] for index in range(3) if mask & (1 << index))
        argument = abs(defect) + RADIUS * norm - 2 * RADIUS * subset_sum
        term = positive_part(argument) ** 2
        total += -term if (mask.bit_count() % 2) else term
    return total / (2 * shape[0] * shape[1] * shape[2])


def allowed_defects(shape):
    bound = RADIUS * sum(shape)
    return tuple(
        defect for defect in range(-int(bound) - 1, int(bound) + 2)
        if defect and defect % 3 and abs(defect) < bound
    )


@lru_cache(maxsize=None)
def nearest_intervals(speed):
    inverse = pow(speed, -1, 3)
    rows = []
    for nearest in range(speed + 1):
        left = max(Q(0), Q(nearest, speed) - RADIUS / speed)
        right = min(Q(1), Q(nearest, speed) + RADIUS / speed)
        if left < right:
            owner = (-inverse * nearest) % 3
            rows.append((left, right, nearest, owner))
    rows.sort()
    return tuple(rows)


def intersect_tagged_lists(first, second):
    """Two-pointer tagged intersection without relying on mutable cursors."""
    first_cursor = second_cursor = 0
    rows = []
    while first_cursor < len(first) and second_cursor < len(second):
        first_row = first[first_cursor]
        second_row = second[second_cursor]
        left = max(first_row[0], second_row[0])
        right = min(first_row[1], second_row[1])
        if left < right:
            first_nearest = first_row[2] if isinstance(first_row[2], tuple) else (first_row[2],)
            first_owner = first_row[3] if isinstance(first_row[3], tuple) else (first_row[3],)
            rows.append((left, right, first_nearest + (second_row[2],), first_owner + (second_row[3],)))
        if first_row[1] <= second_row[1]:
            first_cursor += 1
        if second_row[1] <= first_row[1]:
            second_cursor += 1
    return tuple(rows)


def carrier_roof(w, carrier):
    return max(Q(0), min(
        2 * RADIUS / w[0],
        2 * RADIUS / w[1],
        2 * RADIUS / w[2],
        RADIUS / w[0] + RADIUS / w[1] - Q(abs(carrier[2]), w[0] * w[1]),
        RADIUS / w[0] + RADIUS / w[2] - Q(abs(carrier[1]), w[0] * w[2]),
        RADIUS / w[1] + RADIUS / w[2] - Q(abs(carrier[0]), w[1] * w[2]),
    ))


@lru_cache(maxsize=None)
def physical_carriers(w):
    pair = intersect_tagged_lists(nearest_intervals(w[0]), nearest_intervals(w[1]))
    triple = intersect_tagged_lists(pair, nearest_intervals(w[2]))
    dictionary = {}
    for left, right, nearest, owners in triple:
        if len(set(owners)) != 3:
            continue
        carrier = cross(w, nearest)
        dictionary[carrier] = dictionary.get(carrier, Q(0)) + right - left
    for carrier, length in dictionary.items():
        require(dot(w, carrier) == 0, ("physical-carrier-orthogonal", w, carrier))
        require(all(value % 3 for value in carrier), ("physical-owner-gate", w, carrier))
        require(length == carrier_roof(w, carrier),
                ("physical-roof", w, carrier, length, carrier_roof(w, carrier)))
    return tuple(sorted(dictionary.items()))


def extended_gcd(first, second):
    old_r, r = abs(first), abs(second)
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    coefficient_first = old_s if first >= 0 else -old_s
    coefficient_second = old_t if second >= 0 else -old_t
    return old_r, coefficient_first, coefficient_second


def bezout_vector(ray):
    common, first_coefficient, second_coefficient = extended_gcd(ray[0], ray[1])
    final_gcd, common_coefficient, third_coefficient = extended_gcd(common, ray[2])
    guard(final_gcd == 1, ("primitive-ray", ray, final_gcd))
    vector = (
        common_coefficient * first_coefficient,
        common_coefficient * second_coefficient,
        third_coefficient,
    )
    guard(dot(ray, vector) == 1, ("bezout", ray, vector, dot(ray, vector)))
    return vector


def affine_chart_carriers(w, ray, shape):
    vector = bezout_vector(ray)
    dictionary = {}
    profile = {}
    for defect in allowed_defects(shape):
        base = cross(w, tuple(defect * value for value in vector))
        bounds = tuple(
            RADIUS * (w[(index + 1) % 3] + w[(index + 2) % 3])
            for index in range(3)
        )
        limit = max(
            (bounds[index] + abs(base[index])) / abs(ray[index])
            for index in range(3)
        )
        integer_limit = (limit.numerator + limit.denominator - 1) // limit.denominator + 1
        for step in range(-integer_limit, integer_limit + 1):
            carrier = tuple(base[index] + step * ray[index] for index in range(3))
            if not all(value % 3 for value in carrier):
                continue
            length = carrier_roof(w, carrier)
            if not length:
                continue
            require(dot(w, carrier) == 0, ("chart-carrier-orthogonal", w, ray, defect, carrier))
            require(carrier not in dictionary, ("chart-duplicate-carrier", w, ray, carrier))
            dictionary[carrier] = length
            profile[defect] = profile.get(defect, Q(0)) + length
    return tuple(sorted(dictionary.items())), tuple(sorted(profile.items()))


def next_admissible(value):
    candidate = value + 1
    while not eligible_speed(candidate):
        candidate += 1
    return candidate


EXPECTED = {
    (2, 3, 11): (Q(191, 9702), Q(282744, 8237), 31, (19, 19, 15, 38), Q(58, 833), (5, 7, 17)),
    (2, 5, 9): (Q(124, 6615), Q(65205, 2122), 29, (17, 19, 15, 36), Q(12, 161), (1, 11, 23)),
    (3, 5, 8): (Q(8, 441), Q(4158, 145), 25, (11, 11, 8, 18), Q(6, 77), (1, 5, 11)),
    (5, 5, 6): (Q(127, 7350), Q(214200, 6241), 31, (13, 13, 12, 32), Q(8, 119), (5, 11, 17)),
    (1, 2, 15): (Q(4, 245), Q(15015, 382), 37, (15, 15, 14, 38), Q(60, 1001), (1, 11, 13)),
    (1, 3, 14): (Q(6, 343), Q(2450, 73), 31, (13, 13, 9, 28), Q(12, 175), (1, 13, 25)),
    (1, 5, 12): (Q(53, 2940), Q(95760, 2833), 31, (15, 16, 12, 32), Q(64, 931), (7, 13, 19)),
    (1, 6, 11): (Q(8, 441), Q(4914, 137), 35, (19, 20, 17, 38), Q(6, 91), (1, 7, 13)),
    (1, 8, 9): (Q(5, 294), Q(138600, 3581), 37, (25, 25, 19, 56), Q(118, 1925), (7, 11, 25)),
    (2, 3, 13): (Q(106, 5733), Q(113022, 3695), 29, (8, 8, 7, 22), Q(12, 161), (1, 11, 23)),
    (2, 7, 9): (Q(6, 343), Q(3381, 88), 37, (25, 25, 23, 56), Q(10, 161), (5, 13, 23)),
    (3, 4, 11): (Q(367, 19404), Q(3160080, 92647), 31, (17, 18, 14, 36), Q(46, 665), (5, 7, 19)),
    (3, 5, 10): (Q(407, 22050), Q(642600, 19181), 31, (16, 16, 12, 32), Q(58, 833), (5, 7, 17)),
    (3, 7, 8): (Q(6, 343), Q(539, 19), 25, (14, 15, 13, 36), Q(6, 77), (1, 5, 11)),
    (4, 5, 9): (Q(479, 26460), Q(1043280, 26783), 37, (25, 25, 22, 56), Q(10, 161), (5, 13, 23)),
    (5, 6, 7): (Q(6, 343), Q(5586, 167), 31, (23, 25, 19, 52), Q(64, 931), (7, 13, 19)),
    (1, 1, 18): (Q(20, 1323), Q(62937, 953), 65, (23, 23, 20, 68), Q(2, 37), (1, 19, 37)),
    (1, 3, 16): (Q(1, 63), Q(4698, 79), 59, (41, 42, 37, 116), Q(12, 203), (5, 17, 29)),
    (1, 4, 15): (Q(37, 2205), Q(13770, 283), 47, (29, 29, 25, 72), Q(58, 833), (5, 7, 17)),
    (1, 6, 13): (Q(199, 11466), Q(324324, 7639), 41, (24, 25, 21, 58), Q(6, 77), (1, 5, 11)),
    (1, 7, 12): (Q(6, 343), Q(7497, 145), 49, (42, 42, 40, 136), Q(8, 119), (5, 11, 17)),
    (1, 9, 10): (Q(22, 1323), Q(49329, 815), 59, (59, 59, 57, 200), Q(12, 203), (7, 11, 29)),
    (3, 4, 13): (Q(137, 7644), Q(216216, 5045), 41, (27, 27, 21, 54), Q(6, 77), (1, 5, 11)),
    (3, 7, 10): (Q(6, 343), Q(637, 12), 53, (47, 47, 45, 138), Q(6, 91), (1, 7, 13)),
    (4, 7, 9): (Q(6, 343), Q(1617, 38), 41, (37, 37, 31, 90), Q(6, 77), (1, 5, 11)),
    (6, 7, 7): (Q(6, 343), Q(8379, 167), 49, (32, 32, 31, 100), Q(64, 931), (7, 13, 19)),
}


def main():
    sys.stdout.reconfigure(newline="\n")
    tripwire = "--tripwire" in sys.argv[1:]
    require(all(argument == "--tripwire" for argument in sys.argv[1:]),
            ("unknown-arguments", sys.argv[1:]))
    shapes = coefficient_shapes()
    expected_shapes = tuple(EXPECTED)
    require(shapes == expected_shapes, ("shape-list", shapes, expected_shapes))
    require(len(shapes) == 26, ("shape-count", len(shapes)))
    require(len(SPEEDS) == 22, ("eligible-speed-count", len(SPEEDS)))

    presentations = {}
    for shape in shapes:
        brute = brute_presentations(shape)
        solved = solved_presentations(shape)
        require(brute == solved, ("two-route-presentations", shape, brute, solved))
        presentations[shape] = brute

    # Two independent exact cube-slice constructions.
    bulks = {}
    for shape in shapes:
        defects = allowed_defects(shape)
        require(defects == (
            (-2, -1, 1, 2) if sum(shape) in (16, 18) else (-4, -2, -1, 1, 2, 4)
        ), ("defect-set", shape, defects))
        slices = {}
        for defect in defects:
            convolution = slice_by_convolution(shape, defect)
            clipping = slice_by_clipping(shape, defect)
            require(convolution == clipping, ("slice-two-route", shape, defect, convolution, clipping))
            slices[defect] = convolution
        bulk = sum(slices.values(), Q(0)) / 3
        require(bulk == EXPECTED[shape][0], ("bulk", shape, bulk, EXPECTED[shape][0]))
        bulks[shape] = bulk

    row_summaries = []
    all_row_triples = set()
    all_positive_triples = set()
    equality_shapes = []
    equality_profiles = {}
    unit_rays = signed_labelled_rays((1, 1, 2))
    total_chart_checks = 0
    repeated_shape_triples = 0
    equality_chart_checks = 0
    shape_memberships = {}

    for shape in shapes:
        bulk, expected_threshold, expected_height, expected_counts, expected_mass, expected_winner = EXPECTED[shape]
        defects = allowed_defects(shape)
        error = Q(3 * len(defects), 7)
        all_rows = presentations[shape]

        masses = {
            w: sum((length for _, length in physical_carriers(w)), Q(0))
            for w in all_rows
        }
        candidate_mass = max(masses.values())
        winners = tuple(sorted(w for w, mass in masses.items() if mass == candidate_mass))
        require(candidate_mass == expected_mass, ("row-mass", shape, candidate_mass, expected_mass))
        require(winners == (expected_winner,), ("unique-winner", shape, winners, expected_winner))
        require(any(dot(ray, expected_winner) == 0 for ray in unit_rays),
                ("winner-unit-112", shape, expected_winner))

        threshold = error / (candidate_mass - bulk)
        require(threshold == expected_threshold, ("threshold", shape, threshold, expected_threshold))
        admissible_at_or_below = tuple(value for value in SPEEDS if Q(value) <= threshold)
        height = admissible_at_or_below[-1]
        require(height == expected_height, ("height", shape, height, expected_height))
        following = next_admissible(height)
        require(Q(height) <= threshold < Q(following),
                ("admissible-cutoff", shape, height, threshold, following))
        require(bulk + error / height >= candidate_mass, ("cutoff-not-strict", shape))
        require(bulk + error / following < candidate_mass, ("next-cutoff-strict", shape))

        row = {w: rays for w, rays in all_rows.items() if max(w) <= height}
        require(any(max(w) == height for w in row), ("cutoff-occupied", shape, height))
        require(all(masses[w] < candidate_mass for w in all_rows if max(w) > height),
                ("scanned-tail-control", shape, height))
        positive = tuple(w for w in row if masses[w] > 0)
        relation_count = sum(len(rays) for rays in row.values())
        component_count = sum(len(physical_carriers(w)) for w in row)
        counts = (len(row), relation_count, len(positive), component_count)
        require(counts == expected_counts, ("row-counts", shape, counts, expected_counts))

        for w, rays in row.items():
            physical = physical_carriers(w)
            if len(rays) > 1:
                repeated_shape_triples += 1
            for ray in rays:
                chart, profile = affine_chart_carriers(w, ray, shape)
                require(chart == physical, ("chart-vs-physical", shape, w, ray, chart, physical))
                total_chart_checks += 1
                if w == (1, 5, 11) and candidate_mass == TARGET:
                    equality_profiles.setdefault(shape, set()).add(profile)
                    equality_chart_checks += 1

        all_row_triples.update(row)
        all_positive_triples.update(positive)
        for w in row:
            shape_memberships.setdefault(w, []).append(shape)
        if candidate_mass == TARGET:
            equality_shapes.append(shape)
        row_summaries.append((
            shape,
            ff(bulk),
            ff(threshold),
            height,
            counts,
            ff(candidate_mass),
            expected_winner,
        ))

    require(sum(row[4][0] for row in row_summaries) == 636, "total-row-triples")
    require(sum(row[4][1] for row in row_summaries) == 646, "total-relation-presentations")
    require(sum(row[4][2] for row in row_summaries) == 559, "total-positive-row-triples")
    require(sum(row[4][3] for row in row_summaries) == 1638, "total-components")
    require(len(all_row_triples) == 370, ("union-triples", len(all_row_triples)))
    require(len(all_positive_triples) == 340, ("union-positive", len(all_positive_triples)))
    require(total_chart_checks == 646, ("chart-check-count", total_chart_checks))
    require(repeated_shape_triples == 10,
            ("same-shape-multiple-ray-triples", repeated_shape_triples))
    require(equality_chart_checks == 5, ("equality-chart-count", equality_chart_checks))
    expected_equality_shapes = ((3, 5, 8), (3, 7, 8), (1, 6, 13), (3, 4, 13), (4, 7, 9))
    require(tuple(equality_shapes) == expected_equality_shapes,
            ("equality-shapes", equality_shapes, expected_equality_shapes))

    profile_a = ((-1, Q(3, 77)), (1, Q(3, 77)))
    profile_b = ((-2, Q(3, 77)), (2, Q(3, 77)))
    for shape in expected_equality_shapes:
        expected_profile = profile_a if shape in ((3, 5, 8), (1, 6, 13), (4, 7, 9)) else profile_b
        require(equality_profiles[shape] == {expected_profile},
                ("equality-defect-profile", shape, equality_profiles[shape], expected_profile))

    physical_equality = physical_carriers((1, 5, 11))
    require(sum((length for _, length in physical_equality), Q(0)) == TARGET,
            ("physical-equality-mass", physical_equality))
    require(len(physical_equality) == 2, ("physical-equality-components", physical_equality))

    membership_histogram = {}
    for memberships in shape_memberships.values():
        membership_histogram[len(memberships)] = membership_histogram.get(len(memberships), 0) + 1
    membership_histogram = tuple(sorted(membership_histogram.items()))
    membership_excess = sum(len(memberships) - 1 for memberships in shape_memberships.values())
    expected_membership_histogram = (
        (1, 219), (2, 88), (3, 31), (4, 20), (5, 7), (6, 3), (7, 1), (8, 1)
    )
    require(membership_histogram == expected_membership_histogram,
            ("cross-shape-membership-histogram", membership_histogram))
    require(membership_excess == 266, ("row-membership-excess", membership_excess))
    maximally_overlapping = tuple(
        sorted(w for w, memberships in shape_memberships.items() if len(memberships) == 8)
    )
    require(maximally_overlapping == ((1, 17, 19),),
            ("maximally-overlapping-triple", maximally_overlapping))

    semantic = {
        "shapes": shapes,
        "rows": row_summaries,
        "totals": (636, 646, 559, 1638, 370, 340),
        "equality_shapes": expected_equality_shapes,
        "equality_profiles": {
            str(shape): [
                [(defect, ff(mass)) for defect, mass in profile]
                for profile in sorted(equality_profiles[shape])
            ]
            for shape in expected_equality_shapes
        },
        "repeated_shape_triples": repeated_shape_triples,
        "membership_histogram": membership_histogram,
        "scope": "one-zero presentation sectors norms 16,18,20; no minimality; LRC14 OPEN",
    }
    digest = sha256(dumps(semantic, sort_keys=True).encode()).hexdigest()

    if tripwire:
        require(False, "intentional optimization-safe tripwire")

    print("LRC ONE-ZERO NORM-16/18/20 CLEAN-ROOM AUDIT")
    print("status=PASS FINITE_EXACT_RELATIVE_TO_TAIL_LEMMA; LRC14=OPEN")
    print("shape_counts=norm16:4,norm18:12,norm20:10,total:26")
    print("shapes=" + str(shapes))
    print("safe_scan_height=65; eligible_speeds=22; primitive_triples=" + str(len(TRIPLES)))
    print("presentation_universes=relation_solve_equals_brute_triples:yes")
    print("cube_slices=signed_convolution_equals_rational_polygon_clipping:yes")
    print("row_summaries=(shape,B,threshold,H,U/R/P/C,mass,winner)")
    for row in row_summaries:
        print(str(row))
    print("totals=U:636,R:646,P:559,C:1638,union_U:370,union_positive:340")
    print("chart_vs_literal_physical_dictionary_checks=" + str(total_chart_checks))
    print("same_shape_multi_ray_triples=" + str(repeated_shape_triples))
    print("cross_shape_membership_histogram=" + str(membership_histogram)
          + "; excess_row_memberships=" + str(membership_excess))
    print("maximum_new_shape_multiplicity=8 at triple=(1,17,19)")
    print("global_max=6/77; equality_shapes=" + str(expected_equality_shapes))
    print("common_equality_comb=(1,5,11); physical_components=2")
    print("equality_profiles=A[-1,1]x3/77 for 358,1-6-13,479; B[-2,2]x3/77 for 378,3-4-13")
    print("all_row_winners_have_signed_(1,1,2)_relation=yes")
    print("physical_mass_used_to_find_finite_row_maxima_and_tail_cutoffs=yes")
    print("producer_tables_used_as_expected-value assertions_only=yes")
    print("optimization_safe_checks=yes")
    print("checks=" + str(CHECKS))
    print("semantic_sha256=" + digest)
    print("verdict=PASS")


if __name__ == "__main__":
    main()
