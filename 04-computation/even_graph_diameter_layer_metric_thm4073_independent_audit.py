#!/usr/bin/env python3
"""Independent quotient-level hostile audit for THM-4073.

This path scans every edge subset, recognizes Eulerian graphs and simple
cycles from degree/connectivity conditions, canonicalizes graphs by a direct
permutation scan, builds each simple orbit graph, and runs BFS there.  It does
not use a cycle-space basis or cyclic-order generation.

No floating-point arithmetic or Python assertions are used.
"""

from collections import deque
from itertools import combinations, permutations


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def complete_edges(n):
    return tuple(combinations(range(n), 2))


def degrees(mask, n):
    result = [0] * n
    for position, (left, right) in enumerate(complete_edges(n)):
        if (mask >> position) & 1:
            result[left] += 1
            result[right] += 1
    return result


def connected_nonisolated_support(mask, n):
    degree_list = degrees(mask, n)
    present = {vertex for vertex, degree in enumerate(degree_list) if degree}
    if not present:
        return False
    adjacency = [[] for _ in range(n)]
    for position, (left, right) in enumerate(complete_edges(n)):
        if (mask >> position) & 1:
            adjacency[left].append(right)
            adjacency[right].append(left)
    reached = {next(iter(present))}
    frontier = list(reached)
    while frontier:
        vertex = frontier.pop()
        for neighbor in adjacency[vertex]:
            if neighbor not in reached:
                reached.add(neighbor)
                frontier.append(neighbor)
    return reached == present


def recognize_universes(n):
    eulerian = []
    cycles = {}
    for mask in range(1 << len(complete_edges(n))):
        degree_list = degrees(mask, n)
        if all(degree % 2 == 0 for degree in degree_list):
            eulerian.append(mask)
        support_size = sum(degree != 0 for degree in degree_list)
        if support_size >= 3 \
                and all(degree in (0, 2) for degree in degree_list) \
                and connected_nonisolated_support(mask, n):
            cycles.setdefault(support_size, []).append(mask)
    return tuple(eulerian), {length: tuple(values)
                             for length, values in cycles.items()}


def relabel(mask, n, ordering):
    edges = complete_edges(n)
    edge_position = {edge: position for position, edge in enumerate(edges)}
    result = 0
    for position, (left, right) in enumerate(edges):
        if (mask >> position) & 1:
            image = tuple(sorted((ordering[left], ordering[right])))
            result |= 1 << edge_position[image]
    return result


def canonical_table(states, n):
    orderings = tuple(permutations(range(n)))
    return {
        state: min(relabel(state, n, ordering) for ordering in orderings)
        for state in states
    }


def orbit_graph(states, generators, canonical):
    classes = tuple(sorted(set(canonical.values())))
    class_number = {representative: index
                    for index, representative in enumerate(classes)}
    adjacency = [set() for _ in classes]
    for state in states:
        source = class_number[canonical[state]]
        for generator in generators:
            target = class_number[canonical[state ^ generator]]
            if source != target:
                adjacency[source].add(target)
                adjacency[target].add(source)
    return classes, adjacency


def distances_from(adjacency, source):
    distances = {source: 0}
    queue = deque([source])
    while queue:
        vertex = queue.popleft()
        for neighbor in adjacency[vertex]:
            if neighbor not in distances:
                distances[neighbor] = distances[vertex] + 1
                queue.append(neighbor)
    return distances


def graph_diameter(adjacency):
    return max(
        max(distances_from(adjacency, source).values())
        for source in range(len(adjacency))
    )


def main():
    print("THM-4073 independent quotient audit")
    print("universe=all edge subsets; cycles=connected 2-regular recognition")
    print("quotient=direct permutation canon; arithmetic=integer/GF(2)")
    print()

    metric_checks = 0
    subset_scans = 0
    transition_gates = 0
    relabel_gates = 0
    print("quotient rows: n D vertices edges deg_empty ecc_empty diameter cycle_distances")
    for n in range(3, 7):
        states, by_length = recognize_universes(n)
        subset_scans += 1 << len(complete_edges(n))
        expected_states = 1 << ((n - 1) * (n - 2) // 2)
        check(len(states) == expected_states, "Eulerian subset count failed")
        canonical = canonical_table(states, n)
        relabel_gates += len(states) * factorial(n)
        for diameter in range(2, n):
            generators = tuple(
                cycle
                for length in range(3, diameter + 2)
                for cycle in by_length[length]
            )
            transition_gates += len(states) * len(generators)
            classes, adjacency = orbit_graph(states, generators, canonical)
            empty = classes.index(canonical[0])
            distances = distances_from(adjacency, empty)
            check(len(distances) == len(classes), "orbit graph disconnected")
            profile = []
            for length in range(3, n + 1):
                target_classes = {canonical[cycle] for cycle in by_length[length]}
                check(len(target_classes) == 1, "cycle class split by recognition")
                target = classes.index(next(iter(target_classes)))
                value = distances[target]
                predicted = (length - 2 + diameter - 2) // (diameter - 1)
                check(value == predicted, "quotient distance formula failed")
                profile.append(value)
                metric_checks += 1
            check(len(adjacency[empty]) == diameter - 1,
                  "empty-class degree formula failed")
            edge_count = sum(map(len, adjacency)) // 2
            eccentricity = max(distances.values())
            diameter_value = graph_diameter(adjacency)
            check(eccentricity >= profile[-1], "eccentricity lower bound failed")
            check(diameter_value >= eccentricity, "diameter lower bound failed")
            print("  %d %d %d %d %d %d %d %s"
                  % (n, diameter, len(classes), edge_count,
                     len(adjacency[empty]), eccentricity, diameter_value,
                     tuple(profile)))

    check(metric_checks == 30, "wrong independent metric-check count")
    check(subset_scans == 33864, "wrong edge-subset scan count")
    check(transition_gates == 433754, "wrong quotient transition count")
    check(relabel_gates == 745164, "wrong direct relabel count")
    print()
    print("metric_checks=%d edge_subsets=%d transition_gates=%d relabel_gates=%d"
          % (metric_checks, subset_scans, transition_gates, relabel_gates))
    print("hostile_n6_D3_D4_equal_diameter=3")
    print("PASS: independent simple-orbit graphs reproduce THM-4073")


def factorial(value):
    result = 1
    for factor in range(2, value + 1):
        result *= factor
    return result


if __name__ == "__main__":
    main()
