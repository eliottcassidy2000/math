#!/usr/bin/env python3
"""Exact audit of the two semantically different 72/539 occurrences.

This is self-contained Fraction arithmetic.  It reconstructs the THM-4437
coefficient slice, the THM-4449 dyadic owner cut, and exact hostiles to an
affine/projective identification.  Runtime gates survive python -O.
"""

from fractions import Fraction as Q
from itertools import combinations_with_replacement, permutations
from math import gcd


R = Q(3, 14)
TARGET = Q(72, 539)
CHECKS = 0


def need(test: bool, detail: object) -> None:
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def merge(intervals: list[tuple[Q, Q]]) -> list[tuple[Q, Q]]:
    answer: list[list[Q]] = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if not answer or left > answer[-1][1]:
            answer.append([left, right])
        else:
            answer[-1][1] = max(answer[-1][1], right)
    return [(left, right) for left, right in answer]


def intersection(
    first: list[tuple[Q, Q]], second: list[tuple[Q, Q]]
) -> list[tuple[Q, Q]]:
    """Two-cursor intersection with unambiguous equal-end handling."""
    answer: list[tuple[Q, Q]] = []
    i = j = 0
    while i < len(first) and j < len(second):
        left, right = max(first[i][0], second[j][0]), min(first[i][1], second[j][1])
        if left < right:
            answer.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        elif second[j][1] < first[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return merge(answer)


def danger(speed: int) -> list[tuple[Q, Q]]:
    arcs: list[tuple[Q, Q]] = []
    radius = Q(1, 14 * speed)
    for k in range(speed):
        left, right = Q(k, speed) - radius, Q(k, speed) + radius
        if left < 0:
            arcs.extend(((Q(0), right), (1 + left, Q(1))))
        elif right > 1:
            arcs.extend(((left, Q(1)), (Q(0), right - 1)))
        else:
            arcs.append((left, right))
    return merge(arcs)


def shift_half(intervals: list[tuple[Q, Q]]) -> list[tuple[Q, Q]]:
    shifted: list[tuple[Q, Q]] = []
    for left, right in intervals:
        left, right = left + Q(1, 2), right + Q(1, 2)
        if left >= 1:
            shifted.append((left - 1, right - 1))
        elif right > 1:
            shifted.extend(((left, Q(1)), (Q(0), right - 1)))
        else:
            shifted.append((left, right))
    return merge(shifted)


def measure(intervals: list[tuple[Q, Q]]) -> Q:
    return sum((right - left for left, right in intervals), Q())


def physical_failure(tails: tuple[int, int, int]) -> tuple[Q, list[tuple[Q, Q]]]:
    union = merge([arc for speed in tails for arc in danger(speed)])
    cells = intersection(union, shift_half(union))
    return measure(cells), cells


def circle_distance(value: Q) -> Q:
    value %= 1
    return min(value, 1 - value)


def bad(speed: int, phase: Q) -> bool:
    return circle_distance(speed * phase) < Q(1, 14)


def both_sheets_bad(p: int, q: int, y: Q) -> bool:
    return all(bad(p, (y + j) / 2) or bad(q, (y + j) / 2) for j in (0, 1))


def cross_components(p: int, q: int) -> tuple[tuple[Q, Q], ...]:
    walls = {Q(0), Q(1)}
    for shift in (Q(0), Q(1, 2)):
        for speed in (p, q):
            for integer in range(speed):
                for sign in (-1, 1):
                    walls.add((2 * ((Q(integer) + Q(sign, 14)) / speed - shift)) % 1)
    ordered = sorted(walls)
    raw = [
        (left, right) for left, right in zip(ordered, ordered[1:])
        if both_sheets_bad(p, q, (left + right) / 2)
    ]
    answer: list[tuple[Q, Q]] = []
    for left, right in raw:
        if answer and answer[-1][1] == left and both_sheets_bad(p, q, left):
            answer[-1] = (answer[-1][0], right)
        else:
            answer.append((left, right))
    return tuple(answer)


def keep_halfplane(
    polygon: list[tuple[Q, Q]], A: int, B: int, C: Q
) -> list[tuple[Q, Q]]:
    if not polygon:
        return []
    result: list[tuple[Q, Q]] = []
    for u, v in zip(polygon, polygon[1:] + polygon[:1]):
        fu, fv = A * u[0] + B * u[1] - C, A * v[0] + B * v[1] - C
        if fu <= 0:
            result.append(u)
        if fu * fv < 0:
            t = fu / (fu - fv)
            result.append((u[0] + t * (v[0] - u[0]), u[1] + t * (v[1] - u[1])))
    return list(dict.fromkeys(result))


def plane_polygon(
    v: tuple[int, int, int], delta: Q, lo: Q, hi: Q
) -> tuple[tuple[Q, Q, Q], ...]:
    need(v[0] != 0, ("pivot", v))
    polygon = [(lo, lo), (hi, lo), (hi, hi), (lo, hi)]
    lower, upper = sorted((v[0] * lo, v[0] * hi))
    polygon = keep_halfplane(polygon, -v[1], -v[2], upper - delta)
    polygon = keep_halfplane(polygon, v[1], v[2], delta - lower)
    vertices = tuple(((delta - v[1] * y - v[2] * z) / v[0], y, z) for y, z in polygon)
    need(bool(vertices), ("nonempty slice", v, delta))
    need(all(sum(v[i] * e[i] for i in range(3)) == delta for e in vertices),
         ("slice equations", v, delta))
    return vertices


def sectors(pattern: tuple[int, int, int]) -> tuple[tuple[int, int, int], ...]:
    return tuple(sorted({
        tuple(-a[i] if i == k else a[i] for i in range(3))
        for a in permutations(pattern)
        for k in range(3)
        if a[k] != 0
        if any((-a[i] if i == k else a[i]) > 0 for i in range(3))
        and any((-a[i] if i == k else a[i]) < 0 for i in range(3))
    }))


def compile_pattern(
    pattern: tuple[int, int, int]
) -> tuple[tuple[int, ...], Q, Q, tuple[tuple[int, int, int], tuple[Q, Q, Q]]]:
    unit = all(x % 3 for x in pattern)
    bound = (3 * sum(pattern) - 1) // 14
    defects = tuple(d for d in range(-bound, bound + 1) if (d % 3 == 0) == unit)
    rho = Q(2, 3) if unit else Q(1, 3)
    intercept = (Q(4, 3) if unit else Q(1)) * len(defects)
    best = Q(-1)
    witness = None
    for original in sectors(pattern):
        pivot = next(i for i in range(3) if original[i])
        v = original[pivot:] + original[:pivot]
        i = next(i for i in range(1, 3) if v[i])
        j, k = (i + 1) % 3, (i + 2) % 3
        error_polygons = [plane_polygon(v, Q(d), -R, R) for d in defects]
        speed_polygon = plane_polygon(v, Q(0), Q(0), Q(1))
        for w in speed_polygon:
            value = Q()
            for polygon in error_polygons:
                scalars = [(w[j] * e[k] - w[k] * e[j]) / v[i] for e in polygon]
                value += rho * (max(scalars) - min(scalars))
            if value > best:
                best, witness = value, (v, w)
    need(witness is not None, ("compiled witness", pattern))
    return defects, best, intercept, witness


def additive_witness_parts() -> tuple[tuple[int, Q], ...]:
    v = (-11, 4, 7)
    w = (Q(1), Q(1), Q(1))
    i, j, k = 1, 2, 0
    parts = []
    for d in (-3, 0, 3):
        polygon = plane_polygon(v, Q(d), -R, R)
        scalars = [(w[j] * e[k] - w[k] * e[j]) / v[i] for e in polygon]
        parts.append((d, Q(2, 3) * (max(scalars) - min(scalars))))
    return tuple(parts)


def affine_label_maps(
    source: tuple[int, int, int], target: tuple[int, int, int]
) -> tuple[tuple[Q, Q, tuple[int, int, int]], ...]:
    maps = []
    for image in permutations(target):
        a = Q(image[1] - image[0], source[1] - source[0])
        b = Q(image[0]) - a * source[0]
        if a * source[2] + b == image[2]:
            maps.append((a, b, image))
    return tuple(maps)


def projectively_equivalent(
    source: tuple[int, int, int], target: tuple[int, int, int]
) -> bool:
    return any(
        all(Q(image[i], source[i]) == Q(image[0], source[0]) for i in range(3))
        for image in permutations(target)
    )


def eligible(pattern: tuple[int, int, int]) -> bool:
    return (
        sum(x != 0 for x in pattern) >= 2
        and gcd(gcd(pattern[0], pattern[1]), pattern[2]) == 1
        and sum(x % 3 == 0 for x in pattern) <= 1
        and pattern != (0, 1, 1)
    )


def main() -> None:
    coefficient = compile_pattern((4, 7, 11))
    need(coefficient[0] == (-3, 0, 3) and coefficient[1] == TARGET,
         ("coefficient occurrence", coefficient))
    need(coefficient[3] == ((-11, 4, 7), (Q(1), Q(1), Q(1))),
         ("equal-speed boundary witness", coefficient[3]))
    parts = additive_witness_parts()
    expected_parts = ((-3, Q(2, 49)), (0, Q(4, 77)), (3, Q(2, 49)))
    need(parts == expected_parts and sum((value for _, value in parts), Q()) == TARGET,
         ("coefficient defect decomposition", parts))

    pair_cells = {
        (1, 7): cross_components(1, 7),
        (1, 11): cross_components(1, 11),
        (7, 11): cross_components(7, 11),
    }
    expected_cells = {
        (1, 7): ((Q(6, 49), Q(1, 7)), (Q(6, 7), Q(43, 49))),
        (1, 11): ((Q(6, 77), Q(8, 77)), (Q(69, 77), Q(71, 77))),
        (7, 11): ((Q(13, 49), Q(2, 7)), (Q(5, 7), Q(36, 49))),
    }
    need(pair_cells == expected_cells, ("pair component atlas", pair_cells))
    pair_masses = tuple(
        (pair, sum((right - left for left, right in pair_cells[pair]), Q()))
        for pair in ((1, 7), (1, 11), (7, 11))
    )
    need(tuple(value for _, value in pair_masses) == (Q(2, 49), Q(4, 77), Q(2, 49)),
         ("pair energy decomposition", pair_masses))
    all_pair_cells = [cell for cells in pair_cells.values() for cell in cells]
    need(sum(len(cells) for cells in pair_cells.values()) == len(merge(all_pair_cells)),
         "pair cross-combs disjoint")

    dyadic_mass, dyadic_cells = physical_failure((1, 7, 11))
    need(dyadic_mass == TARGET and len(dyadic_cells) == 12,
         ("dyadic physical equality", dyadic_mass, dyadic_cells))
    energy = sum((value for _, value in pair_masses), Q())
    need(energy == dyadic_mass, ("zero mixed-owner correction", energy, dyadic_mass))

    # Category-error hostile: the coefficient tuple used as a speed tuple.
    hostile_mass, hostile_cells = physical_failure((4, 7, 11))
    need(danger(4) == shift_half(danger(4)), "even speed owns both half-lifts")
    need(hostile_mass == Q(9, 49) and hostile_mass != TARGET and len(hostile_cells) == 9,
         ("coefficient-as-speed hostile", hostile_mass, hostile_cells))

    need(not projectively_equivalent((1, 7, 11), (4, 7, 11)),
         "no common-dilation/permutation equivalence")
    affine = affine_label_maps((1, 7, 11), (4, 7, 11))
    need(affine == (), ("no affine label bijection", affine))

    # The tempting transformation (1,p,q) -> (q-p,p,q) is not a theorem.
    nearby_coefficient = compile_pattern((6, 7, 13))[1]
    nearby_physical = physical_failure((1, 7, 13))[0]
    need(nearby_coefficient == Q(233, 1911), nearby_coefficient)
    need(nearby_physical == Q(80, 637) == Q(240, 1911), nearby_physical)
    need(nearby_coefficient != nearby_physical, "nearby transformation hostile")

    # Exact uniqueness inside THM-4437's complete max-coefficient-18 box.
    hits = []
    count = 0
    for pattern in combinations_with_replacement(range(19), 3):
        if not eligible(pattern):
            continue
        count += 1
        if compile_pattern(pattern)[1] == TARGET:
            hits.append(pattern)
    need(count == 750 and hits == [(4, 7, 11)], ("coefficient-box uniqueness", count, hits))

    print("coefficient_pattern=(4,7,11);defects=-3,0,3;boundary_w=(1,1,1);slope=72/539")
    print("coefficient_weighted_parts=-3:2/49,0:4/77,3:2/49")
    print("dyadic_speed_shape=(1,7,11);pair_masses=1:7->2/49,1:11->4/77,7:11->2/49")
    print("dyadic_owner_correction=0;physical_mass=pair_energy=72/539;quotient_pair_components=6")
    print("structural_bridge=exact_termwise_scalar_decomposition;no_canonical_outer-sign_owner_label")
    print("projective_common_scale=False;affine_label_maps=0")
    print("coefficient_as_speed_hostile=(4,7,11):9/49;even_4_half_shift_invariant=True")
    print("nearby_hostile=(6,7,13)_coefficient:233/1911;(1,7,13)_physical:240/1911")
    print("coefficient_box_size=750;72/539_patterns=(4,7,11)")
    print(f"status=PASS;checks={CHECKS};typed=SCALAR_BRIDGE_NOT_AFFINE_PROJECTIVE_OWNER_CONJUGACY")


if __name__ == "__main__":
    main()
