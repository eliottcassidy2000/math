#!/usr/bin/env python3
"""Primary exact audit for THM-4073.

Cycles are generated from cyclic vertex orders.  The Eulerian universe is
generated independently from the anchored-triangle basis of the cycle space.
Breadth-first search then computes the word metric for every cycle-length
layer through n=7.  A small exact orbit computation also checks the weighted
commuting lift and the failure of commutation after Boolean support collapse.

No floating-point arithmetic or Python assertions are used.
"""

from collections import deque
from itertools import combinations, permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def complete_edges(n):
    return tuple(combinations(range(n), 2))


def cycle_mask(order, edge_position):
    result = 0
    for index, vertex in enumerate(order):
        edge = tuple(sorted((vertex, order[(index + 1) % len(order)])))
        result ^= 1 << edge_position[edge]
    return result


def cycles_by_length(n):
    """Generate each unoriented labelled simple cycle exactly once."""
    edge_position = {edge: index for index, edge in enumerate(complete_edges(n))}
    result = {}
    for length in range(3, n + 1):
        masks = set()
        for vertices in combinations(range(n), length):
            anchor = vertices[0]
            for tail in permutations(vertices[1:]):
                order = (anchor,) + tail
                if order[1] > order[-1]:
                    continue
                masks.add(cycle_mask(order, edge_position))
        expected = combination(n, length) * factorial(length - 1) // 2
        require(len(masks) == expected, "cycle enumeration count failed")
        result[length] = tuple(sorted(masks))
    return result


def factorial(value):
    result = 1
    for factor in range(2, value + 1):
        result *= factor
    return result


def combination(n, k):
    return factorial(n) // (factorial(k) * factorial(n - k))


def anchored_triangle_basis(n):
    edge_position = {edge: index for index, edge in enumerate(complete_edges(n))}
    return tuple(
        cycle_mask((0, left, right), edge_position)
        for left, right in combinations(range(1, n), 2)
    )


def span(generators):
    values = [0]
    for generator in generators:
        values += [value ^ generator for value in values]
    require(len(set(values)) == len(values), "triangle basis is dependent")
    return frozenset(values)


def word_distances(generators):
    distances = {0: 0}
    queue = deque([0])
    while queue:
        state = queue.popleft()
        for generator in generators:
            target = state ^ generator
            if target not in distances:
                distances[target] = distances[state] + 1
                queue.append(target)
    return distances


def permuted_mask(mask, n, ordering):
    edges = complete_edges(n)
    edge_position = {edge: index for index, edge in enumerate(edges)}
    result = 0
    for position, (left, right) in enumerate(edges):
        if (mask >> position) & 1:
            image = tuple(sorted((ordering[left], ordering[right])))
            result |= 1 << edge_position[image]
    return result


def canonicalizer(n):
    orderings = tuple(permutations(range(n)))
    cache = {}

    def canonical(mask):
        if mask not in cache:
            cache[mask] = min(permuted_mask(mask, n, ordering)
                              for ordering in orderings)
        return cache[mask]

    return canonical


def matrix_product(left, right):
    size = len(left)
    return [
        [sum(left[row][middle] * right[middle][column]
             for middle in range(size))
         for column in range(size)]
        for row in range(size)
    ]


def orbit_operator(canonical_classes, cycles, canonical):
    return [
        [sum(canonical(source ^ cycle) == target for cycle in cycles)
         for target in canonical_classes]
        for source in canonical_classes
    ]


def simple_off_diagonal_support(matrix):
    return [
        [int(row != column and matrix[row][column] != 0)
         for column in range(len(matrix))]
        for row in range(len(matrix))
    ]


def main():
    print("THM-4073 primary exact metric audit")
    print("universe=anchored-triangle span; cycles=cyclic-order generation")
    print("arithmetic=integer/GF(2); assertions=disabled-by-construction")
    print()

    total_metric_checks = 0
    total_transitions = 0
    print("cycle-distance profiles; each row lists d([0],[C_3]),...,d([0],[C_n])")
    for n in range(3, 8):
        universe = span(anchored_triangle_basis(n))
        expected_states = 1 << ((n - 1) * (n - 2) // 2)
        require(len(universe) == expected_states, "wrong Eulerian universe size")
        by_length = cycles_by_length(n)
        rows = []
        for diameter in range(2, n):
            generators = tuple(
                cycle
                for length in range(3, diameter + 2)
                for cycle in by_length[length]
            )
            distances = word_distances(generators)
            require(set(distances) == set(universe), "layer Cayley graph disconnected")
            total_transitions += len(universe) * len(generators)
            profile = []
            for length in range(3, n + 1):
                observed = {distances[cycle] for cycle in by_length[length]}
                require(len(observed) == 1, "cycle orbit has varying distance")
                value = next(iter(observed))
                predicted = (length - 2 + diameter - 2) // (diameter - 1)
                require(value == predicted, "exact cycle-distance formula failed")
                profile.append(value)
                total_metric_checks += 1
            rows.append((diameter, tuple(profile)))
        print("  n=%d states=%d cycles=%d rows=%s"
              % (n, len(universe), sum(map(len, by_length.values())), rows))

    require(total_metric_checks == 55, "wrong primary metric-check count")
    require(total_transitions == 84024922, "wrong transition census")
    print("metric_checks=%d transition_gates=%d" %
          (total_metric_checks, total_transitions))
    print()

    n = 4
    universe = span(anchored_triangle_basis(n))
    canonical = canonicalizer(n)
    classes = tuple(sorted({canonical(state) for state in universe}))
    by_length = cycles_by_length(n)
    matrix_three = orbit_operator(classes, by_length[3], canonical)
    matrix_four = orbit_operator(classes, by_length[4], canonical)
    require(matrix_product(matrix_three, matrix_four)
            == matrix_product(matrix_four, matrix_three),
            "weighted orbit operators failed to commute")
    support_three = simple_off_diagonal_support(matrix_three)
    support_four = simple_off_diagonal_support(matrix_four)
    require(matrix_product(support_three, support_four)
            != matrix_product(support_four, support_three),
            "Boolean support collapse unexpectedly commuted")
    require(matrix_three == [[0, 4, 0], [1, 0, 3], [0, 4, 0]],
            "unexpected triangle orbit operator")
    require(matrix_four == [[0, 0, 3], [0, 3, 0], [1, 0, 2]],
            "unexpected four-cycle orbit operator")
    print("n=4 weighted orbit operators")
    print("  M3=%s" % matrix_three)
    print("  M4=%s" % matrix_four)
    print("  weighted_commute=True boolean_support_commute=False")
    print()
    print("PASS: 55 exact metric cases and the typed commuting lift survived")


if __name__ == "__main__":
    main()
