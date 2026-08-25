#!/usr/bin/env python3
"""Independent perfect-matching audit of THM-4116's cubic four-pole atlas.

This path never assigns edge colors recursively.  For a fixed boundary word,
the internal edges of each color must instead be a perfect matching on the
vertices not already covered by a semiedge of that color.  Compatible triples
of such matchings are counted exactly.  Standard-library integers only.
"""

from collections import Counter
from functools import lru_cache
from itertools import combinations


COLORS = range(3)
REPRESENTATIVES = (
    (0, 0, 0, 0),
    (0, 0, 1, 1),
    (0, 1, 0, 1),
    (0, 1, 1, 0),
)


def edge(u, v):
    return (u, v) if u < v else (v, u)


def flower_edges(k):
    assert k >= 5 and k % 2 == 1
    center = lambda i: 4 * (i % k)
    b_leaf = lambda i: 4 * (i % k) + 1
    c_leaf = lambda i: 4 * (i % k) + 2
    d_leaf = lambda i: 4 * (i % k) + 3
    edges = set()
    for i in range(k):
        for leaf in (b_leaf(i), c_leaf(i), d_leaf(i)):
            edges.add(edge(center(i), leaf))
        edges.add(edge(b_leaf(i), b_leaf(i + 1)))
    cycle = tuple(c_leaf(i) for i in range(k)) + tuple(
        d_leaf(i) for i in range(k)
    )
    for index, vertex in enumerate(cycle):
        edges.add(edge(vertex, cycle[(index + 1) % len(cycle)]))
    return tuple(sorted(edges))


def graph_definitions():
    k4 = tuple(combinations(range(4), 2))
    petersen = (
        (0, 1), (0, 4), (0, 5), (1, 2), (1, 6),
        (2, 3), (2, 7), (3, 4), (3, 8), (4, 9),
        (5, 7), (5, 8), (6, 8), (6, 9), (7, 9),
    )
    blanusa1 = tuple(sorted((
        (0, 5), (1, 17), (2, 14), (3, 8), (4, 17),
        (6, 11), (7, 17), (9, 13), (10, 15), (12, 16),
        *((i, i + 1) for i in range(16)),
        (0, 16),
    )))
    blanusa2 = (
        (0, 1), (0, 2), (0, 14), (1, 5), (1, 11),
        (2, 3), (2, 6), (3, 4), (3, 9), (4, 5),
        (4, 7), (5, 6), (6, 8), (7, 8), (7, 17),
        (8, 9), (9, 15), (10, 11), (10, 14), (10, 16),
        (11, 12), (12, 13), (12, 17), (13, 14), (13, 15),
        (15, 16), (16, 17),
    )
    return (
        ("K4", tuple(range(4)), tuple(k4)),
        ("Petersen", tuple(range(10)), tuple(petersen)),
        ("J5", tuple(range(20)), flower_edges(5)),
        ("BlanusaI", tuple(range(18)), tuple(blanusa1)),
        ("BlanusaII", tuple(range(18)), tuple(blanusa2)),
    )


def matching_extension_count(vertices, edges, removed_edge, word):
    """Count one boundary word via three compatible perfect matchings."""
    u, v = removed_edge
    adjacency = {vertex: set() for vertex in vertices}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    boundary_live = tuple(
        sorted(adjacency[u] - {v}) + sorted(adjacency[v] - {u})
    )
    live_vertices = frozenset(set(vertices) - {u, v})
    internal_edges = frozenset(
        e for e in edges if e[0] in live_vertices and e[1] in live_vertices
    )
    live_adjacency = {
        vertex: frozenset(adjacency[vertex] & live_vertices)
        for vertex in live_vertices
    }

    boundary_by_color = {color: [] for color in COLORS}
    for vertex, color in zip(boundary_live, word):
        boundary_by_color[color].append(vertex)
    if any(
        len(points) != len(set(points))
        for points in boundary_by_color.values()
    ):
        return 0

    @lru_cache(maxsize=None)
    def perfect_matchings(remaining):
        remaining = frozenset(remaining)
        if not remaining:
            return (frozenset(),)
        if len(remaining) % 2:
            return ()
        first = min(remaining)
        results = []
        for partner in sorted(live_adjacency[first] & remaining):
            smaller = remaining - {first, partner}
            for rest in perfect_matchings(smaller):
                results.append(rest | {edge(first, partner)})
        return tuple(results)

    color_matchings = []
    for color in COLORS:
        uncovered = live_vertices - set(boundary_by_color[color])
        color_matchings.append(perfect_matchings(frozenset(uncovered)))

    count = 0
    for matching0 in color_matchings[0]:
        for matching1 in color_matchings[1]:
            if matching0 & matching1:
                continue
            used = matching0 | matching1
            for matching2 in color_matchings[2]:
                if used & matching2:
                    continue
                if used | matching2 == internal_edges:
                    count += 1
    return count


def atlas(vertices, edges):
    edges = tuple(sorted(edge(*item) for item in edges))
    adjacency = {vertex: set() for vertex in vertices}
    for u, v in edges:
        adjacency[u].add(v)
        adjacency[v].add(u)
    assert all(len(adjacency[vertex]) == 3 for vertex in vertices)
    signatures = Counter()
    for removed_edge in edges:
        signature = tuple(
            matching_extension_count(
                vertices, edges, removed_edge, representative
            )
            for representative in REPRESENTATIVES
        )
        signatures[signature] += 1
    return signatures


def main():
    expected = {
        "K4": Counter({(0, 1, 0, 1): 6}),
        "Petersen": Counter({(2, 2, 0, 0): 15}),
        "J5": Counter({
            (4, 4, 0, 0): 5,
            (6, 6, 0, 0): 10,
            (10, 10, 0, 0): 10,
            (12, 12, 0, 0): 5,
        }),
        "BlanusaI": Counter({(4, 4, 0, 0): 23, (6, 6, 0, 0): 4}),
        "BlanusaII": Counter({(4, 4, 0, 0): 25, (6, 6, 0, 0): 2}),
    }
    observed = {}
    edge_sides = 0
    for name, vertices, edges in graph_definitions():
        observed[name] = atlas(vertices, edges)
        assert observed[name] == expected[name]
        edge_sides += len(edges)
    assert edge_sides == 105
    fields = "; ".join(
        f"{name}={tuple(sorted(observed[name].items()))}"
        for name in expected
    )
    print(f"Independent matching atlas: edge_sides={edge_sides}; {fields}")


if __name__ == "__main__":
    main()
