#!/usr/bin/env python3
"""Independent exact audit of THM-3391 weighted common-source capacity.

This companion deliberately does not import the frozen probe.  It obtains cyclic
band traces by sweeping the exact phase-boundary arrangement, and obtains q=2
cover measures by counting cells on a common rational grid.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import gcd, lcm


RHO = Fraction(1, 14)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def mod_one(value):
    return value - value.numerator // value.denominator


def cyclic_distance(value):
    residue = mod_one(value)
    return min(residue, 1 - residue)


def phase_probes(n, residues):
    boundaries = {Fraction(0)}
    for residue in residues:
        centre = -Fraction(residue, n)
        boundaries.add(mod_one(centre - RHO))
        boundaries.add(mod_one(centre + RHO))
    ordered = sorted(boundaries)
    probes = set(ordered)
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += 1
        probes.add(mod_one((left + right) / 2))
    return tuple(sorted(probes))


def exact_edges(n, images):
    """All exact source masks attained at boundaries or in open phase cells."""
    probes = phase_probes(n, images)
    edges = set()
    for theta in probes:
        mask = 0
        for index, residue in enumerate(images):
            if cyclic_distance(theta + Fraction(residue, n)) < RHO:
                mask |= 1 << index
        edges.add(mask)
    return tuple(sorted(edges))


def exact_lambda(n, weights, multiplier):
    images = tuple((multiplier * residue) % n for residue in range(n))
    edges = exact_edges(n, images)
    return max(sum(weight for i, weight in enumerate(weights) if edge >> i & 1) for edge in edges)


def block_lambda(n, weights, multiplier):
    pushed = [0] * n
    for residue, weight in enumerate(weights):
        pushed[(multiplier * residue) % n] += weight
    width = (n + 6) // 7
    return max(
        sum(pushed[(start + offset) % n] for offset in range(width))
        for start in range(n)
    )


def edge_weight(edge, weights):
    return sum(weight for index, weight in enumerate(weights) if edge >> index & 1)


def q2_danger(speed, sheet, cell_index, grid):
    # y=(2j+1)/(2 grid), x=(y+sheet)/2.
    denominator = 4 * grid
    numerator = speed * (2 * cell_index + 1 + 2 * grid * sheet)
    residue = numerator % denominator
    return 14 * min(residue, denominator - residue) < denominator


def q2_grid_measure(first, second):
    grid = 7 * lcm(first, second)
    covered = 0
    for cell in range(grid):
        if all(
            q2_danger(first, sheet, cell, grid)
            or q2_danger(second, sheet, cell, grid)
            for sheet in (0, 1)
        ):
            covered += 1
    return Fraction(covered, grid)


def danger_at_point(speed, numerator, denominator):
    residue = (speed * numerator) % denominator
    return 14 * min(residue, denominator - residue) < denominator


def endpoint_witness():
    core = (2, 4, 5, 9)
    transverse = (5, 9)
    degree_grid = 2520
    bad_open_cells = []
    for cell in range(degree_grid):
        # y=(2j+1)/(2D); sheet preimages have denominator 4D.
        y_num = 2 * cell + 1
        core_danger = any(danger_at_point(c, y_num, 2 * degree_grid) for c in core)
        transverse_cover = all(
            any(
                danger_at_point(u, y_num + 2 * degree_grid * sheet, 4 * degree_grid)
                for u in transverse
            )
            for sheet in (0, 1)
        )
        if transverse_cover and not core_danger:
            bad_open_cells.append(cell)
    residual_points = []
    for point in range(degree_grid):
        core_danger = any(danger_at_point(c, point, degree_grid) for c in core)
        transverse_cover = all(
            any(
                danger_at_point(u, point + degree_grid * sheet, 2 * degree_grid)
                for u in transverse
            )
            for sheet in (0, 1)
        )
        if transverse_cover and not core_danger:
            residual_points.append(Fraction(point, degree_grid))
    return tuple(bad_open_cells), tuple(residual_points)


def main():
    formula_cases = 0
    # Exhaust all 0/1 weights through n=10, including the strict n=7 case.
    for n in range(2, 11):
        for bits in range(1 << n):
            weights = tuple((bits >> residue) & 1 for residue in range(n))
            for multiplier in range(n):
                observed = exact_lambda(n, weights, multiplier)
                predicted = block_lambda(n, weights, multiplier)
                require(observed == predicted, ("window", n, bits, multiplier, observed, predicted))
                formula_cases += 1
    # Weighted, collision-heavy controls at two further multiples of seven.
    for n in (14, 21):
        for seed in range(7):
            weights = tuple((seed + 3 * j + j * j) % 5 for j in range(n))
            for multiplier in range(n):
                observed = exact_lambda(n, weights, multiplier)
                predicted = block_lambda(n, weights, multiplier)
                require(observed == predicted, ("weighted window", n, seed, multiplier))
                formula_cases += 1

    equality_cases = 0
    # Independently check the equality-partition iff with arbitrary positive
    # weights (zeros are discarded) for every pair of multipliers through n=6.
    for n in range(2, 7):
        edge_cache = {
            a: exact_edges(n, tuple((a * residue) % n for residue in range(n)))
            for a in range(n)
        }
        for weights in product(range(3), repeat=n):
            positive = sum((1 << i) for i, weight in enumerate(weights) if weight)
            total = sum(weights)
            if not total:
                continue
            for first in range(n):
                first_edges = edge_cache[first]
                first_lambda = max(edge_weight(edge, weights) for edge in first_edges)
                for second in range(n):
                    second_edges = edge_cache[second]
                    second_lambda = max(edge_weight(edge, weights) for edge in second_edges)
                    if first_lambda + second_lambda != total:
                        continue
                    cover = any(
                        ((left | right) & positive) == positive
                        for left in first_edges
                        for right in second_edges
                    )
                    partition = any(
                        edge_weight(left, weights) == first_lambda
                        and edge_weight(right, weights) == second_lambda
                        and not (left & right & positive)
                        and ((left | right) & positive) == positive
                        for left in first_edges
                        for right in second_edges
                    )
                    require(cover == partition, ("partition", n, weights, first, second))
                    equality_cases += 1

    mixed_map_cases = 0
    # Ground-set identity is retained while the two maps land in different
    # quotients.  Exhaust all maps C_3 -> X_2 and C_3 -> X_3.
    weights = (1, 2, 3)
    total = sum(weights)
    for map_two in product(range(2), repeat=3):
        edges_two = exact_edges(2, map_two)
        lambda_two = max(edge_weight(edge, weights) for edge in edges_two)
        for map_three in product(range(3), repeat=3):
            edges_three = exact_edges(3, map_three)
            lambda_three = max(edge_weight(edge, weights) for edge in edges_three)
            cover = any((left | right) == 0b111 for left in edges_two for right in edges_three)
            if lambda_two + lambda_three < total:
                require(not cover, ("mixed quotient union bound", map_two, map_three))
            if lambda_two + lambda_three == total and cover:
                require(
                    any(
                        edge_weight(left, weights) == lambda_two
                        and edge_weight(right, weights) == lambda_three
                        and not (left & right)
                        and (left | right) == 0b111
                        for left in edges_two
                        for right in edges_three
                    ),
                    ("mixed quotient equality", map_two, map_three),
                )
            mixed_map_cases += 1

    unit_support_cases = 0
    for n in range(2, 15):
        units = tuple(a for a in range(n) if gcd(a, n) == 1)
        edge_cache = {
            a: exact_edges(n, tuple((a * residue) % n for residue in range(n)))
            for a in units
        }
        for support in range(1, 1 << n):
            weights = tuple((support >> residue) & 1 for residue in range(n))
            universal_survival = all(block_lambda(n, weights, a) < support.bit_count() for a in units)
            no_full_edge = not any(
                (edge & support) == support for a in units for edge in edge_cache[a]
            )
            require(universal_survival == no_full_edge, ("unit maximum", n, support))
            unit_support_cases += 1

    coset_cases = 0
    singleton_cases = 0
    for d in range(2, 31):
        for order in range(2, d + 1):
            if d % order:
                continue
            step = d // order
            for offset in range(step):
                source = tuple(offset + j * step for j in range(order))
                for multiplier in range(d):
                    images = tuple((multiplier * residue) % d for residue in source)
                    edges = exact_edges(d, images)
                    observed = max(edge.bit_count() for edge in edges)
                    common = gcd(multiplier, order)
                    predicted = common * ((order // common + 6) // 7)
                    require(observed == predicted, ("coset", d, order, offset, multiplier))
                    coset_cases += 1
                    if order <= 7 and gcd(multiplier, d) == 1:
                        require(observed == 1, ("unit coset load", d, order, multiplier))
                        attained = {edge for edge in edges if edge.bit_count() == 1}
                        require(len(attained) == order, ("singleton availability", d, order, multiplier))
                        singleton_cases += 1

    q2_cases = 0
    max_measure = Fraction(0)
    for first in range(1, 40, 2):
        for second in range(1, 40, 2):
            measure = q2_grid_measure(first, second)
            edge = first + second > 7 * gcd(first, second)
            require((measure > 0) == edge, ("q2 edge", first, second, measure))
            require(measure <= Fraction(4, 7), ("q2 measure", first, second, measure))
            max_measure = max(max_measure, measure)
            q2_cases += 1
    require(q2_grid_measure(1, 3) == 0, "q2 nonedge")
    require(q2_grid_measure(1, 9) == Fraction(4, 63), "q2 edge measure")

    bad_cells, residual_points = endpoint_witness()
    require(not bad_cells, bad_cells)
    require(residual_points == (Fraction(3, 14), Fraction(11, 14)), residual_points)

    packet = (
        formula_cases,
        equality_cases,
        mixed_map_cases,
        unit_support_cases,
        coset_cases,
        singleton_cases,
        q2_cases,
        max_measure,
        bad_cells,
        residual_points,
    )
    print("INDEPENDENT WEIGHTED SUPPORT AUDIT")
    print("method=exact_phase_arrangement_plus_common_grid_cell_count")
    print(
        f"cases=window:{formula_cases};equality:{equality_cases};mixed_maps:{mixed_map_cases};"
        f"unit_support:{unit_support_cases};coset:{coset_cases};q2:{q2_cases}"
    )
    print(f"q2_1_9={q2_grid_measure(1, 9)};q2_max={max_measure}")
    print(f"endpoint_residual={residual_points};open_cell_residual={len(bad_cells)}")
    print(f"audit_sha256={sha256(repr(packet).encode('ascii')).hexdigest()}")
    print("all_independent_controls=PASS")


if __name__ == "__main__":
    main()
