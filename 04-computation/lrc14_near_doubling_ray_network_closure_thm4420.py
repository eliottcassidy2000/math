#!/usr/bin/env python3
"""Standalone exact verifier for the THM-4420 near-doubling ray theorem.

Analytic family:
    w=(1,m,2m+sigma), sigma in {+1,-1}, all speeds distinct odd 3-units.

The main loop uses only the proved carrier formula.  A separate finite control
path reconstructs carriers from the relation lattice and reconstructs the six
literal sheet contact networks, including exact rational max flow.
"""

from collections import deque
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from json import dumps
import argparse
import sys


TARGET = Q(6, 77)
R = Q(3, 14)
CHECKS = 0


class VerificationError(RuntimeError):
    pass


def require(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise VerificationError(payload)


def strict_floor(numerator, denominator):
    return (numerator - 1) // denominator


def ff(value):
    return f"{value.numerator}/{value.denominator}"


def family_parameters(m, sigma):
    c = 2 * m + sigma
    require(sigma in (-1, 1), ("sigma", sigma))
    require(m >= (5 if sigma == 1 else 7), ("base-m", m, sigma))
    require(m % 2 == 1 and c % 2 == 1, ("odd", m, c))
    require(m % 3 and c % 3, ("ternary-unit", m, c))
    require(1 < m < c, ("sorted-distinct", m, c))
    K = strict_floor(3 * (m + 1 if sigma == 1 else m), 14)
    carriers = tuple(
        (-sigma * k, -2 * k, k)
        for k in range(-K, K + 1)
        if k % 3
    )
    return (1, m, c), K, carriers


def pair_roofs(w, carrier):
    a, b, c = w
    first, second, third = carrier
    return (
        Q(3 * (a + b) - 14 * abs(third), 14 * a * b),
        Q(3 * (a + c) - 14 * abs(second), 14 * a * c),
        Q(3 * (b + c) - 14 * abs(first), 14 * b * c),
    )


def raw_values(w, carriers):
    q = Q(3, 7 * w[2])
    capacities = [Q(0), Q(0), Q(0)]
    physical = Q(0)
    for carrier in carriers:
        require(sum(w[i] * carrier[i] for i in range(3)) == 0,
                ("carrier-relation", w, carrier))
        require(all(value % 3 for value in carrier), ("carrier-live", w, carrier))
        roofs = pair_roofs(w, carrier)
        require(all(value > 0 for value in roofs), ("carrier-roof", w, carrier, roofs))
        for index in range(3):
            capacities[index] += min(q, roofs[index])
        physical += min((q,) + roofs)
    return tuple(capacities), physical


def predicted_outside_roof(w, K, sigma):
    carrier = (-sigma * (K + 1), -2 * (K + 1), K + 1)
    return carrier, pair_roofs(w, carrier)


def lattice_box_carriers(w):
    """Independent ker(w) enumeration, not using the one-ray formula."""
    bounds = (
        strict_floor(3 * (w[1] + w[2]), 14),
        strict_floor(3 * (w[0] + w[2]), 14),
        strict_floor(3 * (w[0] + w[1]), 14),
    )
    result = set()
    for first in range(-bounds[0], bounds[0] + 1):
        for second in range(-bounds[1], bounds[1] + 1):
            numerator = -(w[0] * first + w[1] * second)
            if numerator % w[2]:
                continue
            third = numerator // w[2]
            carrier = (first, second, third)
            if abs(third) > bounds[2] or not all(value % 3 for value in carrier):
                continue
            if all(value > 0 for value in pair_roofs(w, carrier)):
                result.add(carrier)
    return result


@lru_cache(maxsize=None)
def sheet_intervals(speed, sheet):
    radius = Q(1, 14 * speed)
    pieces = []
    for owner in range(speed):
        center = Q(owner, speed) - Q(sheet, 3)
        center -= center.numerator // center.denominator
        left, right = center - radius, center + radius
        if left < 0:
            pieces.extend(((Q(0), right), (left + 1, Q(1))))
        elif right > 1:
            pieces.extend(((left, Q(1)), (Q(0), right - 1)))
        else:
            pieces.append((left, right))
    pieces.sort()
    require(sum((right - left for left, right in pieces), Q(0)) == Q(1, 7),
            ("sheet-measure", speed, sheet))
    for index in range(1, len(pieces)):
        require(pieces[index - 1][1] <= pieces[index][0],
                ("sheet-disjoint", speed, sheet, index))
    return tuple(pieces)


def intersect(first, second):
    i = j = 0
    result = []
    while i < len(first) and j < len(second):
        first_right = first[i][1]
        second_right = second[j][1]
        left = max(first[i][0], second[j][0])
        right = min(first_right, second_right)
        if left < right:
            result.append((left, right))
        if first_right <= second_right:
            i += 1
        if second_right <= first_right:
            j += 1
    return tuple(result)


def exact_maxflow(left, right, edges):
    left_count = len(left)
    right_count = len(right)
    source = left_count + right_count
    sink = source + 1
    adjacency = [[] for _ in range(sink + 1)]
    residual = {}

    def add_edge(u, v, capacity):
        adjacency[u].append(v)
        adjacency[v].append(u)
        residual[u, v] = residual.get((u, v), Q(0)) + capacity
        residual.setdefault((v, u), Q(0))

    left_lengths = tuple(b - a for a, b in left)
    right_lengths = tuple(b - a for a, b in right)
    infinity = sum(left_lengths, Q(0)) + sum(right_lengths, Q(0))
    for index, capacity in enumerate(left_lengths):
        add_edge(source, index, capacity)
    for index, capacity in enumerate(right_lengths):
        add_edge(left_count + index, sink, capacity)
    for i, j in edges:
        add_edge(i, left_count + j, infinity)
    flow = Q(0)
    while True:
        parent = {source: source}
        queue = deque((source,))
        while queue and sink not in parent:
            vertex = queue.popleft()
            for neighbor in adjacency[vertex]:
                if neighbor not in parent and residual[vertex, neighbor] > 0:
                    parent[neighbor] = vertex
                    queue.append(neighbor)
        if sink not in parent:
            return flow
        increment = None
        vertex = sink
        while vertex != source:
            previous = parent[vertex]
            increment = (residual[previous, vertex] if increment is None
                         else min(increment, residual[previous, vertex]))
            vertex = previous
        vertex = sink
        while vertex != source:
            previous = parent[vertex]
            residual[previous, vertex] -= increment
            residual[vertex, previous] += increment
            vertex = previous
        flow += increment


def literal_network_values(w):
    capacities = [Q(0), Q(0), Q(0)]
    physical_by_pair = [Q(0), Q(0), Q(0)]
    pair_choices = ((0, 1), (0, 2), (1, 2))
    total_edges = 0
    for pair_index, pair in enumerate(pair_choices):
        first, second = pair
        third = 3 - first - second
        for assignment in permutations((0, 1, 2)):
            pair_pieces = intersect(
                sheet_intervals(w[first], assignment[first]),
                sheet_intervals(w[second], assignment[second]),
            )
            third_pieces = sheet_intervals(w[third], assignment[third])
            edges = []
            edge_sum = exact = Q(0)
            left_load = [Q(0)] * len(pair_pieces)
            right_load = [Q(0)] * len(third_pieces)
            i = j = 0
            while i < len(pair_pieces) and j < len(third_pieces):
                pair_right = pair_pieces[i][1]
                third_right = third_pieces[j][1]
                overlap_left = max(pair_pieces[i][0], third_pieces[j][0])
                overlap_right = min(pair_right, third_right)
                if overlap_left < overlap_right:
                    edges.append((i, j))
                    exact += overlap_right - overlap_left
                    load = min(pair_right - pair_pieces[i][0],
                               third_right - third_pieces[j][0])
                    edge_sum += load
                    left_load[i] += load
                    right_load[j] += load
                if pair_right <= third_right:
                    i += 1
                if third_right <= pair_right:
                    j += 1
            for index, load in enumerate(left_load):
                require(load <= pair_pieces[index][1] - pair_pieces[index][0],
                        ("left-edgewise-feasible", w, pair, assignment, index, load))
            for index, load in enumerate(right_load):
                require(load <= third_pieces[index][1] - third_pieces[index][0],
                        ("right-edgewise-feasible", w, pair, assignment, index, load))
            flow = exact_maxflow(pair_pieces, third_pieces, edges)
            require(flow == edge_sum, ("flow-vs-edge", w, pair, assignment, flow, edge_sum))
            capacities[pair_index] += flow
            physical_by_pair[pair_index] += exact
            total_edges += len(edges)
    require(physical_by_pair[0] == physical_by_pair[1] == physical_by_pair[2],
            ("literal-mass-pair-independence", w, physical_by_pair))
    return tuple(capacities), physical_by_pair[0], total_edges


def admissible_rows(height):
    for sigma, start in ((1, 5), (-1, 7)):
        for m in range(start, (height - sigma) // 2 + 1, 6):
            c = 2 * m + sigma
            if c <= height:
                yield m, sigma


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=4999)
    args = parser.parse_args()
    require(args.height >= 11, ("height", args.height))
    sys.stdout.reconfigure(newline="\n")

    controls = {
        (5, 1), (11, 1), (17, 1), (23, 1), (41, 1), (83, 1), (2495, 1),
        (7, -1), (13, -1), (19, -1), (25, -1), (61, -1), (103, -1),
        (145, -1), (2497, -1),
    }
    rows = plus = minus = base_checks = large_checks = 0
    literal_controls = []
    maximum = Q(0)
    maximizers = []
    semantic = []
    for m, sigma in admissible_rows(args.height):
        w, K, carriers = family_parameters(m, sigma)
        plus += sigma == 1
        minus += sigma == -1
        rows += 1
        N = 2 * (K - K // 3)
        require(len(carriers) == N, ("carrier-count", w, len(carriers), N))
        outside, outside_roofs = predicted_outside_roof(w, K, sigma)
        require(min(outside_roofs) <= 0,
                ("strict-endpoint-cutoff", w, K, outside, outside_roofs))
        capacities, physical = raw_values(w, carriers)
        q_count = Q(3 * N, 7 * w[2])
        require(all(value <= q_count for value in capacities),
                ("q-count-envelope", w, capacities, q_count))

        if (sigma == 1 and m >= 23) or (sigma == -1 and m >= 25):
            # K-floor(K/3) <= (2K+2)/3, followed by the displayed threshold
            # comparison in the report.
            require(3 * (K - K // 3) <= 2 * K + 2,
                    ("nonmultiple-count-bound", w, K))
            comparison = (Q(2 * (m + 1), 7) + Q(4, 3) if sigma == 1
                          else Q(2 * m, 7) + Q(4, 3))
            threshold = Q(2 * w[2], 11)
            require(comparison <= threshold, ("large-threshold", w, comparison, threshold))
            require(Q(N) < comparison, ("strict-large-count", w, N, comparison))
            large_checks += 1
        else:
            require((m, sigma) in {(5, 1), (11, 1), (17, 1),
                                            (7, -1), (13, -1), (19, -1)},
                    ("finite-base-list", m, sigma))
            require(q_count <= TARGET, ("base-target", w, q_count))
            base_checks += 1
        require(q_count <= TARGET, ("family-target", w, q_count))

        best = min(capacities)
        row = (w, tuple(ff(value) for value in capacities), ff(physical), N, K, ff(q_count))
        semantic.append(row)
        if best > maximum:
            maximum = best
            maximizers = [row]
        elif best == maximum:
            maximizers.append(row)

        if (m, sigma) in controls:
            independent = lattice_box_carriers(w)
            require(independent == set(carriers),
                    ("lattice-box-vs-formula", w, independent ^ set(carriers)))
            literal_caps, literal_mass, edge_count = literal_network_values(w)
            require(literal_caps == capacities,
                    ("literal-network-vs-raw", w, literal_caps, capacities))
            require(literal_mass == physical,
                    ("literal-mass-vs-raw", w, literal_mass, physical))
            literal_controls.append((w, tuple(ff(v) for v in literal_caps),
                                     ff(literal_mass), edge_count, N, K))

    expected_max = ((1, 5, 11), ("6/77", "6/77", "6/77"), "6/77", 2, 1, "6/77")
    require(maximizers == [expected_max], ("unique-equality", maximizers))
    require(rows == plus + minus and base_checks == 6, ("row-accounting", rows, plus, minus))
    digest = sha256(dumps(semantic, sort_keys=True).encode()).hexdigest()
    print("LRC14 NEAR-DOUBLING RAY NETWORK CLOSURE THM-4420 VERIFIER")
    print("status=PASS; theorem=PROVED_RELATIVE_TO_THM_4409_THM_4414; LRC14=OPEN")
    print("family=w=(1,m,2m+sigma), sigma=+/-1, sorted distinct odd ternary units")
    print("height=%d; rows=%d; plus=%d; minus=%d" % (args.height, rows, plus, minus))
    print("carrier_formula=C=k*(-sigma,-2,1), 0<abs(k)<=K, 3_not_divide_k")
    print("K_plus=strict_floor(3(m+1)/14); K_minus=strict_floor(3m/14)")
    print("N=2(K-floor(K/3)); capacity_each_pair<=3N/(7(2m+sigma))<=6/77")
    print("large_threshold_rows=%d; finite_base_rows=%d" % (large_checks, base_checks))
    print("unique_equality=" + str(expected_max))
    print("literal_network_control_count=%d" % len(literal_controls))
    for control in literal_controls:
        print("literal_control=" + str(control))
    print("semantic_sha256=" + digest)
    print("checks=%d" % CHECKS)
    print("verdict=PASS")


if __name__ == "__main__":
    main()
