#!/usr/bin/env python3
"""Exact rational census for degree-zero third-sheet component sidecars.

The certificate is selected without consulting the exact physical mass.  That
mass is computed afterward by literal three-list intersection and is used only
as a validity control.  Besides the first contact/disjoint-bit certificate,
the script computes the strongest separated transportation upper bound that
retains only the bipartite contact graph and component lengths.
"""

from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
from json import dumps
from math import gcd
from collections import deque
import sys


LAMBDA = Q(1, 14)
TARGET = Q(6, 77)
HEIGHT = 79
CHECKS = 0


class VerificationError(RuntimeError):
    pass


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise VerificationError(message)


def ff(value):
    return f"{value.numerator}/{value.denominator}"


@lru_cache(maxsize=None)
def sheet_intervals(speed, sheet):
    radius = LAMBDA / speed
    result = []
    for owner in range(speed):
        center = Q(owner, speed) - Q(sheet, 3)
        center -= center.numerator // center.denominator
        left, right = center - radius, center + radius
        if left < 0:
            result.extend(((Q(0), right), (left + 1, Q(1))))
        elif right > 1:
            result.extend(((left, Q(1)), (Q(0), right - 1)))
        else:
            result.append((left, right))
    result.sort()
    for index, (left, right) in enumerate(result):
        require(Q(0) <= left < right <= Q(1), ("sheet-piece", speed, sheet, left, right))
        if index:
            require(result[index - 1][1] <= left, ("sheet-disjoint", speed, sheet, index))
    require(sum((right - left for left, right in result), Q(0)) == Q(1, 7),
            ("sheet-measure", speed, sheet))
    return tuple(result)


def intersect(first, second):
    i = j = 0
    result = []
    while i < len(first) and j < len(second):
        first_right, second_right = first[i][1], second[j][1]
        left, right = max(first[i][0], second[j][0]), min(first_right, second_right)
        if left < right:
            result.append((left, right))
        if first_right <= second_right:
            i += 1
        if second_right <= first_right:
            j += 1
    return tuple(result)


def bipartite_capacity(left_pieces, right_pieces, edges):
    """Exact maximum transportation flow using only contact and lengths."""
    left_count = len(left_pieces)
    right_count = len(right_pieces)
    source = left_count + right_count
    sink = source + 1
    adjacency = [set() for _ in range(sink + 1)]
    residual = {}

    def add_edge(first, second, capacity):
        adjacency[first].add(second)
        adjacency[second].add(first)
        residual[first, second] = residual.get((first, second), Q(0)) + capacity
        residual.setdefault((second, first), Q(0))

    left_lengths = tuple(right - left for left, right in left_pieces)
    right_lengths = tuple(right - left for left, right in right_pieces)
    infinite = sum(left_lengths, Q(0)) + sum(right_lengths, Q(0))
    for index, capacity in enumerate(left_lengths):
        add_edge(source, index, capacity)
    for index, capacity in enumerate(right_lengths):
        add_edge(left_count + index, sink, capacity)
    for left_index, right_index in edges:
        add_edge(left_index, left_count + right_index, infinite)

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
            capacity = residual[previous, vertex]
            increment = capacity if increment is None else min(increment, capacity)
            vertex = previous
        vertex = sink
        while vertex != source:
            previous = parent[vertex]
            residual[previous, vertex] -= increment
            residual[vertex, previous] += increment
            vertex = previous
        flow += increment


def contact_bit_certificate(w, pair, compute_control=False):
    """Degree-zero bound and optional exact physical mass control."""
    first, second = pair
    exact_index = 3 - first - second
    touch_length = Q(0)
    disjoint_third_length = Q(0)
    edge_capacity = Q(0)
    network_capacity = Q(0)
    exact_mass = Q(0)
    pair_piece_count = touch_count = third_piece_count = disjoint_count = 0
    edge_count = nested_edge_count = crossing_edge_count = 0
    for assignment in permutations((0, 1, 2)):
        pair_pieces = intersect(
            sheet_intervals(w[first], assignment[first]),
            sheet_intervals(w[second], assignment[second]),
        )
        third_pieces = sheet_intervals(w[exact_index], assignment[exact_index])
        pair_piece_count += len(pair_pieces)
        third_piece_count += len(third_pieces)
        # One independent sweep marks only component incidence.  In
        # particular, the upper certificate never records overlap lengths.
        pair_touched = [False] * len(pair_pieces)
        third_touched = [False] * len(third_pieces)
        contact_edges = []
        pair_cursor = third_cursor = 0
        while pair_cursor < len(pair_pieces) and third_cursor < len(third_pieces):
            pair_piece = pair_pieces[pair_cursor]
            third_piece = third_pieces[third_cursor]
            left = max(pair_piece[0], third_piece[0])
            right = min(pair_piece[1], third_piece[1])
            if left < right:
                pair_touched[pair_cursor] = True
                third_touched[third_cursor] = True
                contact_edges.append((pair_cursor, third_cursor))
                edge_count += 1
                pair_length = pair_piece[1] - pair_piece[0]
                third_length = third_piece[1] - third_piece[0]
                edge_capacity += min(pair_length, third_length)
                if (
                    (pair_piece[0] >= third_piece[0] and pair_piece[1] <= third_piece[1])
                    or (third_piece[0] >= pair_piece[0] and third_piece[1] <= pair_piece[1])
                ):
                    nested_edge_count += 1
                else:
                    crossing_edge_count += 1
                if compute_control:
                    exact_mass += right - left
            if pair_piece[1] <= third_piece[1]:
                pair_cursor += 1
            if third_piece[1] <= pair_piece[1]:
                third_cursor += 1

        touch_count += sum(pair_touched)
        touch_length += sum(
            (right - left for touched, (left, right) in zip(pair_touched, pair_pieces) if touched),
            Q(0),
        )
        disjoint_count += len(third_pieces) - sum(third_touched)
        disjoint_third_length += sum(
            (right - left for touched, (left, right) in zip(third_touched, third_pieces) if not touched),
            Q(0),
        )
        network_capacity += bipartite_capacity(pair_pieces, third_pieces, contact_edges)

    # g_i=g_j=1/7, finite dual term is the zero frequency 6/343.
    upper = Q(6, 343) + Q(48, 49) * touch_length - Q(1, 49) * disjoint_third_length
    summary = (
        pair,
        upper,
        network_capacity,
        edge_capacity,
        touch_length,
        disjoint_third_length,
        pair_piece_count,
        touch_count,
        third_piece_count,
        disjoint_count,
        edge_count,
        nested_edge_count,
        crossing_edge_count,
    )
    if compute_control:
        require(exact_mass <= upper, ("certificate-validity", w, summary, exact_mass))
        return summary, exact_mass
    return summary


def strict_integer_below(value):
    """Largest integer strictly below a rational value."""
    return (value.numerator - 1) // value.denominator


def raw_carrier_measure(w):
    """Independent THM-4392 raw-carrier reconstruction of the comb mass."""
    radius = Q(3, 14)
    bounds = (
        strict_integer_below(radius * (w[1] + w[2])),
        strict_integer_below(radius * (w[0] + w[2])),
        strict_integer_below(radius * (w[0] + w[1])),
    )
    mass = Q(0)
    carriers = 0
    for first in range(-bounds[0], bounds[0] + 1):
        for second in range(-bounds[1], bounds[1] + 1):
            numerator = -(w[0] * first + w[1] * second)
            if numerator % w[2]:
                continue
            third = numerator // w[2]
            coefficient = (first, second, third)
            if abs(third) > bounds[2] or any(value % 3 == 0 for value in coefficient):
                continue
            length = max(Q(0), min(
                2 * radius / w[0],
                2 * radius / w[1],
                2 * radius / w[2],
                radius / w[0] + radius / w[1] - Q(abs(third), w[0] * w[1]),
                radius / w[0] + radius / w[2] - Q(abs(second), w[0] * w[2]),
                radius / w[1] + radius / w[2] - Q(abs(first), w[1] * w[2]),
            ))
            if length:
                carriers += 1
                mass += length
    return mass, carriers


def main():
    sys.stdout.reconfigure(newline="\n")
    values = tuple(value for value in range(1, HEIGHT + 1, 2) if value % 3)
    triples = tuple(
        triple for triple in combinations(values, 3)
        if gcd(gcd(triple[0], triple[1]), triple[2]) == 1
    )
    require(values == tuple(value for value in range(1, 80, 2) if value % 3), "value universe")

    closed = 0
    capacity_closed = 0
    capacity_exact = 0
    strict_physical = equality_physical = above_physical = 0
    pair_successes = {(0, 1): 0, (0, 2): 0, (1, 2): 0}
    capacity_pair_successes = {(0, 1): 0, (0, 2): 0, (1, 2): 0}
    best_pair_counts = {(0, 1): 0, (0, 2): 0, (1, 2): 0}
    capacity_best_pair_counts = {(0, 1): 0, (0, 2): 0, (1, 2): 0}
    failures = []
    capacity_failures = []
    equality_rows = []
    carrier_counts = []
    strongest = None
    weakest_closed = None
    capacity_weakest_strict = None
    capacity_target_rows = []
    network_strict_improvements = 0
    semantic_rows = []

    for w in triples:
        # Certificate selection sees only the contact/disjoint sidecar.
        candidates = tuple(
            contact_bit_certificate(w, pair, compute_control=False)
            for pair in ((0, 1), (0, 2), (1, 2))
        )
        best = min(candidates, key=lambda row: (row[1], row[0]))
        capacity_best = min(candidates, key=lambda row: (row[2], row[0]))
        best_pair_counts[best[0]] += 1
        capacity_best_pair_counts[capacity_best[0]] += 1
        for row in candidates:
            if row[1] <= TARGET:
                pair_successes[row[0]] += 1
            if row[2] <= TARGET:
                capacity_pair_successes[row[0]] += 1

        # Only after selection, reconstruct the literal physical mass from a
        # separate three-list intersection control.
        checked_best, exact_mass = contact_bit_certificate(w, best[0], compute_control=True)
        require(checked_best == best, ("stable-certificate-selection", w, best, checked_best))
        for row in candidates:
            require(exact_mass <= row[1], ("all-pair-upper", w, row, exact_mass))
            require(exact_mass <= row[2], ("all-pair-capacity-upper", w, row, exact_mass))
            require(row[2] <= row[3], ("network-no-worse-than-edge-sum", w, row))
            if row[2] < row[3]:
                network_strict_improvements += 1
            if row[12] == 0:
                require(row[2] == exact_mass,
                        ("nested-contact-network-exact", w, row, exact_mass))
        raw_mass, carrier_count = raw_carrier_measure(w)
        require(raw_mass == exact_mass, ("raw-vs-physical", w, raw_mass, exact_mass))
        carrier_counts.append(carrier_count)

        if exact_mass < TARGET:
            strict_physical += 1
        elif exact_mass == TARGET:
            equality_physical += 1
            equality_rows.append((w, best))
        else:
            above_physical += 1

        if best[1] <= TARGET:
            closed += 1
            margin = TARGET - best[1]
            row = (margin, w, best)
            if strongest is None or row > strongest:
                strongest = row
            if weakest_closed is None or row < weakest_closed:
                weakest_closed = row
        else:
            failures.append((w, best, exact_mass))

        if capacity_best[2] <= TARGET:
            capacity_closed += 1
        else:
            capacity_failures.append((w, capacity_best, exact_mass))
        if capacity_best[2] == exact_mass:
            capacity_exact += 1
        if capacity_best[2] == TARGET:
            capacity_target_rows.append((w, capacity_best, exact_mass))
        if exact_mass < TARGET:
            capacity_margin_row = (TARGET - capacity_best[2], w, capacity_best)
            if capacity_weakest_strict is None or capacity_margin_row < capacity_weakest_strict:
                capacity_weakest_strict = capacity_margin_row

        semantic_rows.append(
            (w, ff(best[1]), best[0], ff(capacity_best[2]), capacity_best[0], ff(exact_mass))
        )

    require(len(triples) == 2910, ("triple-count", len(triples)))
    require(strict_physical + equality_physical + above_physical == len(triples), "physical trichotomy")
    require(above_physical == 0, ("finite physical above target", above_physical))
    require(equality_rows == [((1, 5, 11), min(
        (contact_bit_certificate((1, 5, 11), pair) for pair in ((0, 1), (0, 2), (1, 2))),
        key=lambda row: (row[1], row[0]),
    ))], ("unique-equality-control", equality_rows))
    require(len(failures) == 1 and failures[0][0] == (1, 5, 11), ("unique-certificate-failure", failures))
    require(closed == len(triples) - 1, ("all-strict-closed", closed))
    require(not capacity_failures, ("capacity-failures", capacity_failures))
    require(capacity_closed == len(triples), ("capacity-closes-entire-universe", capacity_closed))
    equality_capacities = tuple(
        contact_bit_certificate((1, 5, 11), pair)
        for pair in ((0, 1), (0, 2), (1, 2))
    )
    require(all(row[2] == TARGET for row in equality_capacities),
            ("equality-capacity", equality_capacities))
    require(all(row[10:] == (6, 6, 0) for row in equality_capacities),
            ("equality-all-contacts-nested", equality_capacities))
    require(all(row[7] == row[10] and row[8] - row[9] == row[10]
                for row in equality_capacities),
            ("equality-contact-graphs-are-matchings", equality_capacities))
    require(len(capacity_target_rows) == 1 and capacity_target_rows[0][0] == (1, 5, 11),
            ("unique-capacity-target-row", capacity_target_rows))
    require((strict_physical, equality_physical, above_physical) == (2909, 1, 0),
            ("physical-counts", strict_physical, equality_physical, above_physical))
    require(capacity_exact == 1747, ("capacity-exact-count", capacity_exact))
    require(network_strict_improvements == 0,
            ("network-vs-edge-count", network_strict_improvements))
    require(pair_successes == {(0, 1): 2535, (0, 2): 2854, (1, 2): 2833},
            ("contact-pair-success-counts", pair_successes))
    require(capacity_pair_successes == {(0, 1): 2818, (0, 2): 2855, (1, 2): 2859},
            ("capacity-pair-success-counts", capacity_pair_successes))
    require(best_pair_counts == {(0, 1): 198, (0, 2): 952, (1, 2): 1760},
            ("contact-best-pair-counts", best_pair_counts))
    require(capacity_best_pair_counts == {(0, 1): 400, (0, 2): 533, (1, 2): 1977},
            ("capacity-best-pair-counts", capacity_best_pair_counts))
    require((min(carrier_counts), max(carrier_counts)) == (0, 22),
            ("carrier-count-range", min(carrier_counts), max(carrier_counts)))
    require(
        capacity_weakest_strict[0] == Q(6, 1771)
        and capacity_weakest_strict[1] == (1, 11, 23)
        and capacity_weakest_strict[2][0] == (0, 1)
        and capacity_weakest_strict[2][2] == Q(12, 161),
        ("capacity-weakest-strict", capacity_weakest_strict),
    )

    semantic = {
        "height": HEIGHT,
        "values": values,
        "triple_count": len(triples),
        "closed": closed,
        "capacity_closed": capacity_closed,
        "capacity_exact": capacity_exact,
        "capacity_target_count": len(capacity_target_rows),
        "network_strict_improvements": network_strict_improvements,
        "failures": [
            (w, tuple(ff(value) if isinstance(value, Q) else value for value in best), ff(mass))
            for w, best, mass in failures
        ],
        "physical_counts": (strict_physical, equality_physical, above_physical),
        "carrier_count_range": (min(carrier_counts), max(carrier_counts)),
        "pair_successes": sorted(pair_successes.items()),
        "capacity_pair_successes": sorted(capacity_pair_successes.items()),
        "best_pair_counts": sorted(best_pair_counts.items()),
        "capacity_best_pair_counts": sorted(capacity_best_pair_counts.items()),
        "strongest": (ff(strongest[0]), strongest[1], strongest[2][0], ff(strongest[2][1])),
        "weakest_closed": (
            ff(weakest_closed[0]), weakest_closed[1], weakest_closed[2][0], ff(weakest_closed[2][1])
        ),
        "capacity_weakest_strict": (
            ff(capacity_weakest_strict[0]),
            capacity_weakest_strict[1],
            capacity_weakest_strict[2][0],
            ff(capacity_weakest_strict[2][2]),
        ),
        "rows_sha256": sha256(dumps(semantic_rows, sort_keys=True).encode()).hexdigest(),
        "scope": "primitive sorted distinct positive odd ternary units height<=79; LRC14 OPEN",
    }
    digest = sha256(dumps(semantic, sort_keys=True).encode()).hexdigest()

    failure_w, failure_best, failure_mass = failures[0]
    print("LRC THIRD-SHEET DEGREE-ZERO COMPONENT-BIT CENSUS")
    print("status=PASS FINITE_EXACT; LRC14=OPEN")
    print("universe=sorted primitive distinct positive odd speeds prime_to_3; max_speed<=79")
    print("eligible_speed_count=" + str(len(values)) + "; triples=" + str(len(triples)))
    print("certificate=6/343+(48/49)*touch_pair_length-(1/49)*disjoint_third_length")
    print("closed=" + str(closed) + "; failed=" + str(len(failures)))
    print("network_capacity=max_transport_flow(contact_graph; component_lengths)")
    print("edge_capacity=sum_over_contact_edges(min(pair_component_length,third_component_length))")
    print("network_strict_improvements_over_edge_capacity=" + str(network_strict_improvements))
    print("capacity_closed=" + str(capacity_closed) + "; capacity_failed="
          + str(len(capacity_failures)) + "; capacity_exact=" + str(capacity_exact)
          + "; capacity_target_equalities=" + str(len(capacity_target_rows)))
    print("physical_controls=strict:" + str(strict_physical) + ",equality:" + str(equality_physical)
          + ",above:" + str(above_physical))
    print("raw_carrier_controls=all_triples; carrier_count_range="
          + str((min(carrier_counts), max(carrier_counts))))
    print("unique_failure=" + str(failure_w) + "; best_pair=" + str(failure_best[0])
          + "; upper=" + ff(failure_best[1]) + "; exact_control=" + ff(failure_mass))
    print("failure_sidecar=network_capacity:" + ff(failure_best[2])
          + ",edge_capacity:" + ff(failure_best[3])
          + ",touch_length:" + ff(failure_best[4])
          + ",disjoint_third_length:" + ff(failure_best[5])
          + ",pair_pieces:" + str(failure_best[6])
          + ",touches:" + str(failure_best[7])
          + ",third_pieces:" + str(failure_best[8])
          + ",disjoint:" + str(failure_best[9])
          + ",edges/nested/crossing:" + str(failure_best[10:]))
    print("equality_network_rows=" + str(tuple(
        (row[0], ff(row[2]), row[10:]) for row in equality_capacities
    )))
    print("contact_bit_equality_slack=" + ff(failure_best[1] - failure_mass))
    print("pair_success_counts=" + str(tuple(sorted(pair_successes.items()))))
    print("capacity_pair_success_counts=" + str(tuple(sorted(capacity_pair_successes.items()))))
    print("best_pair_counts=" + str(tuple(sorted(best_pair_counts.items()))))
    print("capacity_best_pair_counts=" + str(tuple(sorted(capacity_best_pair_counts.items()))))
    print("capacity_weakest_strict=triple:" + str(capacity_weakest_strict[1])
          + ",pair:" + str(capacity_weakest_strict[2][0])
          + ",upper:" + ff(capacity_weakest_strict[2][2])
          + ",margin:" + ff(capacity_weakest_strict[0]))
    print("weakest_closed=triple:" + str(weakest_closed[1]) + ",pair:" + str(weakest_closed[2][0])
          + ",upper:" + ff(weakest_closed[2][1]) + ",margin:" + ff(weakest_closed[0]))
    print("strongest_margin=triple:" + str(strongest[1]) + ",pair:" + str(strongest[2][0])
          + ",upper:" + ff(strongest[2][1]) + ",margin:" + ff(strongest[0]))
    print("exact_masses_used_for_selection=no; physical_and_raw_masses_postselection_controls=yes")
    print("optimization_safe_checks=yes")
    print("checks=" + str(CHECKS))
    print("rows_sha256=" + semantic["rows_sha256"])
    print("semantic_sha256=" + digest)
    print("verdict=PASS")


if __name__ == "__main__":
    main()
