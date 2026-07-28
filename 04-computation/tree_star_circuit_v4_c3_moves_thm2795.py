#!/usr/bin/env python3
"""Exact referee for THM-2795.

The script checks the binary star-circuit/dwell-time theorem on every vertex
order of every nonisomorphic tree through seven vertices, verifies the local
V4 diamond and C2/C3 -> S3 relations, and enumerates the graceful-order move
graphs through eight vertices.  Every truth-bearing check uses explicit
exceptions, so optimized execution performs the same audit.
"""

from __future__ import annotations

from functools import reduce
from itertools import permutations

import networkx as nx


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def star_vectors(vertex_count: int, edges):
    stars = [0] * vertex_count
    for edge_index, (u, v) in enumerate(edges):
        bit = 1 << edge_index
        stars[u] ^= bit
        stars[v] ^= bit
    return tuple(stars)


def labels_from_order(order):
    labels = [0] * len(order)
    for position, vertex in enumerate(order):
        labels[vertex] = position
    return tuple(labels)


def prefix_walk(order, stars):
    states = [0]
    state = 0
    for vertex in order:
        state ^= stars[vertex]
        states.append(state)
    return tuple(states)


def dwell_times(states, edge_count: int):
    # Internal thresholds are c_1,...,c_m.
    return tuple(
        sum((state >> edge) & 1 for state in states[1:-1])
        for edge in range(edge_count)
    )


def edge_differences(edges, labels):
    return tuple(abs(labels[u] - labels[v]) for u, v in edges)


def is_graceful(edges, order) -> bool:
    labels = labels_from_order(order)
    return sorted(edge_differences(edges, labels)) == list(
        range(1, len(edges) + 1)
    )


def circuit_and_dwell_audit(max_vertices: int = 7):
    shape_counts = {}
    order_counts = {}
    diamond_count = 0
    for vertex_count in range(2, max_vertices + 1):
        shapes = list(nx.generators.nonisomorphic_trees(vertex_count))
        shape_counts[vertex_count] = len(shapes)
        local_orders = 0
        for graph in shapes:
            edges = tuple(graph.edges())
            edge_count = len(edges)
            stars = star_vectors(vertex_count, edges)
            require(
                reduce(
                    lambda left, right: left ^ right, stars, 0
                )
                == 0,
                "vertex-star sum stopped vanishing",
            )
            for subset in range(1, (1 << vertex_count) - 1):
                value = 0
                for vertex in range(vertex_count):
                    if (subset >> vertex) & 1:
                        value ^= stars[vertex]
                require(
                    value != 0,
                    "vertex stars acquired a proper subcircuit",
                )
            if vertex_count >= 3:
                require(
                    len(set(stars)) == vertex_count
                    and all(star != 0 for star in stars),
                    "nontrivial tree lost distinct nonzero star generators",
                )

            for order in permutations(range(vertex_count)):
                states = prefix_walk(order, stars)
                require(
                    states[0] == states[-1] == 0,
                    "star-circuit walk stopped closing",
                )
                require(
                    len(set(states[:-1])) == vertex_count,
                    "a proper prefix block returned to an old cut state",
                )
                labels = labels_from_order(order)
                require(
                    dwell_times(states, edge_count)
                    == edge_differences(edges, labels),
                    "cut-cube dwell time stopped equalling edge difference",
                )

                if vertex_count >= 3:
                    for position in range(vertex_count - 1):
                        left = order[position]
                        right = order[position + 1]
                        base = states[position]
                        first = base ^ stars[left]
                        second = base ^ stars[right]
                        finish = base ^ stars[left] ^ stars[right]
                        require(
                            len({base, first, second, finish}) == 4,
                            "adjacent swap stopped spanning an affine V4",
                        )

                        swapped = list(order)
                        swapped[position], swapped[position + 1] = (
                            swapped[position + 1],
                            swapped[position],
                        )
                        swapped_states = prefix_walk(tuple(swapped), stars)
                        require(
                            states[: position + 1]
                            == swapped_states[: position + 1]
                            and states[position + 2 :]
                            == swapped_states[position + 2 :]
                            and swapped_states[position + 1] == second,
                            "adjacent swap changed more than one prefix state",
                        )

                        old_dwell = dwell_times(states, edge_count)
                        new_dwell = dwell_times(
                            swapped_states, edge_count
                        )
                        for edge in range(edge_count):
                            changed = old_dwell[edge] != new_dwell[edge]
                            expected = bool(
                                ((stars[left] ^ stars[right]) >> edge) & 1
                            )
                            require(
                                changed == expected
                                and abs(
                                    old_dwell[edge] - new_dwell[edge]
                                )
                                <= 1,
                                "V4 swap lost its signed unit dwell law",
                            )
                        diamond_count += 1
                local_orders += 1
        order_counts[vertex_count] = local_orders
    return shape_counts, order_counts, diamond_count


def rotate_left_triple(order, position):
    result = list(order)
    block = result[position : position + 3]
    result[position : position + 3] = block[1:] + block[:1]
    return tuple(result)


def rotate_right_triple(order, position):
    result = list(order)
    block = result[position : position + 3]
    result[position : position + 3] = block[-1:] + block[:-1]
    return tuple(result)


def swap_adjacent(order, position):
    result = list(order)
    result[position], result[position + 1] = (
        result[position + 1],
        result[position],
    )
    return tuple(result)


def local_s3_audit():
    base = (0, 1, 2)

    def s(order):
        return swap_adjacent(order, 0)

    def r(order):
        return rotate_left_triple(order, 0)

    value = base
    for _ in range(2):
        value = s(value)
    require(value == base, "local binary move lost order two")

    value = base
    for _ in range(3):
        value = r(value)
    require(value == base, "local ternary move lost order three")

    value = base
    for _ in range(2):
        value = s(r(value))
    require(
        value == base,
        "local C2/C3 action lost the S3 quotient relation (sr)^2",
    )

    orbit = set()
    for rotation in range(3):
        value = base
        for _ in range(rotation):
            value = r(value)
        orbit.add(value)
        orbit.add(s(value))
    require(len(orbit) == 6, "local C2/C3 orbit stopped being S3")
    return tuple(sorted(orbit))


def component_sizes(vertices, neighbor_function):
    vertex_set = set(vertices)
    seen = set()
    sizes = []
    edge_twice = 0
    for vertex in vertices:
        edge_twice += sum(
            neighbor in vertex_set
            for neighbor in set(neighbor_function(vertex))
        )
        if vertex in seen:
            continue
        stack = [vertex]
        seen.add(vertex)
        size = 0
        while stack:
            current = stack.pop()
            size += 1
            for neighbor in set(neighbor_function(current)):
                if neighbor in vertex_set and neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        sizes.append(size)
    return tuple(sorted(sizes)), edge_twice // 2


def graceful_move_census(max_vertices: int = 8):
    binary = {}
    binary_ternary = {}
    graceful_totals = {}
    for vertex_count in range(2, max_vertices + 1):
        total_graceful = 0
        binary_components = []
        mixed_components = []
        binary_edges = 0
        mixed_edges = 0
        for graph in nx.generators.nonisomorphic_trees(vertex_count):
            edges = tuple(graph.edges())
            graceful_orders = tuple(
                order
                for order in permutations(range(vertex_count))
                if is_graceful(edges, order)
            )
            total_graceful += len(graceful_orders)

            def binary_neighbors(order):
                return tuple(
                    swap_adjacent(order, position)
                    for position in range(vertex_count - 1)
                )

            def mixed_neighbors(order):
                result = list(binary_neighbors(order))
                for position in range(vertex_count - 2):
                    result.append(rotate_left_triple(order, position))
                    result.append(rotate_right_triple(order, position))
                return tuple(result)

            sizes, edges_binary = component_sizes(
                graceful_orders, binary_neighbors
            )
            sizes_mixed, edges_mixed = component_sizes(
                graceful_orders, mixed_neighbors
            )
            binary_components.extend(sizes)
            mixed_components.extend(sizes_mixed)
            binary_edges += edges_binary
            mixed_edges += edges_mixed

        graceful_totals[vertex_count] = total_graceful
        binary[vertex_count] = (
            len(binary_components),
            binary_edges,
            sum(size == 1 for size in binary_components),
            max(binary_components),
        )
        binary_ternary[vertex_count] = (
            len(mixed_components),
            mixed_edges,
            sum(size == 1 for size in mixed_components),
            max(mixed_components),
        )

    require(
        graceful_totals
        == {
            2: 2,
            3: 4,
            4: 16,
            5: 68,
            6: 376,
            7: 2184,
            8: 15096,
        },
        "graceful-order totals changed",
    )
    require(
        binary
        == {
            2: (1, 1, 0, 2),
            3: (2, 2, 0, 2),
            4: (6, 12, 4, 6),
            5: (18, 76, 12, 24),
            6: (94, 532, 72, 120),
            7: (402, 4052, 260, 720),
            8: (2304, 34270, 1554, 5040),
        },
        "binary V4-move census changed",
    )
    require(
        binary_ternary
        == {
            2: (1, 1, 0, 2),
            3: (1, 4, 0, 4),
            4: (6, 24, 4, 6),
            5: (14, 176, 8, 24),
            6: (62, 1316, 40, 120),
            7: (334, 10216, 202, 720),
            8: (2000, 88696, 1244, 5040),
        },
        f"binary/ternary move census changed: {binary_ternary}",
    )
    return graceful_totals, binary, binary_ternary


def p4_and_star_boundaries():
    path_edges = ((1, 0), (1, 2), (0, 3))
    isolated = (0, 2, 1, 3)
    require(is_graceful(path_edges, isolated), "P4 hostile is nongraceful")
    neighbors = [
        swap_adjacent(isolated, position) for position in range(3)
    ]
    neighbors += [
        rotate_left_triple(isolated, position)
        for position in range(2)
    ]
    neighbors += [
        rotate_right_triple(isolated, position)
        for position in range(2)
    ]
    require(
        all(not is_graceful(path_edges, order) for order in neighbors),
        "P4 hostile acquired a graceful C2/C3 move",
    )

    # For K_(1,m), gracefulness forces the centre to one of the two ends.
    # Local moves preserving gracefulness cannot pass it through the interior
    # once m>=3.
    for leaf_count in range(3, 8):
        vertex_count = leaf_count + 1
        centre = 0
        star_edges = tuple((centre, leaf) for leaf in range(1, vertex_count))
        graceful_orders = tuple(
            order
            for order in permutations(range(vertex_count))
            if is_graceful(star_edges, order)
        )
        require(
            len(graceful_orders) == 2 * factorial(leaf_count)
            and all(
                order[0] == centre or order[-1] == centre
                for order in graceful_orders
            ),
            "star extreme-centre boundary changed",
        )
        sizes, _ = component_sizes(
            graceful_orders,
            lambda order: tuple(
                [swap_adjacent(order, position)
                 for position in range(vertex_count - 1)]
                + [rotate_left_triple(order, position)
                   for position in range(vertex_count - 2)]
                + [rotate_right_triple(order, position)
                   for position in range(vertex_count - 2)]
            ),
        )
        require(
            sizes == (factorial(leaf_count), factorial(leaf_count)),
            "star graceful move graph stopped having two extreme components",
        )
    return isolated


def factorial(value: int) -> int:
    result = 1
    for factor in range(2, value + 1):
        result *= factor
    return result


def main() -> None:
    shapes, orders, diamonds = circuit_and_dwell_audit()
    local_orbit = local_s3_audit()
    graceful, binary, mixed = graceful_move_census()
    p4 = p4_and_star_boundaries()

    print("THM2795 TREE STAR-CIRCUIT / V4-C3 EXACT REFEREE")
    print("all-order shapes through n=7:", shapes)
    print("all-order counts through n=7:", orders)
    print("verified affine V4 diamonds:", diamonds)
    print("local C2/C3 orbit:", local_orbit, "group=S3 quotient")
    print("graceful-order totals:", graceful)
    print("binary move census (components,edges,isolated,max):", binary)
    print("binary+ternary census (components,edges,isolated,max):", mixed)
    print("P4 isolated graceful order:", p4)
    print("star boundary n=4..8: two components of (n-1)! each")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
