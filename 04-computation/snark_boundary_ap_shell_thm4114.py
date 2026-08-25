#!/usr/bin/env python3
"""Exact boundary-state and AP odd-shell synchronization audit for THM-4114.

The graph half checks use labelled edges and labelled colors.  The AP checks
work on the two-torsion quotient from THM-4110; no floating point arithmetic
is used.
"""

from collections import Counter
from itertools import combinations


COLORS = range(3)


def boundary_color_counts(vertices, edges, side, boundary_order):
    """Return the labelled boundary extension-count vector for one side."""
    vertices = set(vertices)
    side = set(side)
    edge_set = set(edges)
    boundary_order = tuple(boundary_order)
    assert set(boundary_order) <= edge_set
    assert all((u in side) != (v in side) for u, v in boundary_order)

    internal = sorted(e for e in edges if e[0] in side and e[1] in side)
    active = internal + list(boundary_order)
    assert len(active) == len(set(active))
    assert all(u in vertices and v in vertices for u, v in active)

    used = {v: set() for v in side}
    assigned = {}
    counts = Counter()

    def search(index):
        if index == len(active):
            word = tuple(assigned[e] for e in boundary_order)
            counts[word] += 1
            return

        edge = active[index]
        endpoints = tuple(v for v in edge if v in side)
        for color in COLORS:
            if any(color in used[v] for v in endpoints):
                continue
            assigned[edge] = color
            for v in endpoints:
                used[v].add(color)
            search(index + 1)
            for v in endpoints:
                used[v].remove(color)
            del assigned[edge]

    search(0)
    return counts


def full_coloring_count(vertices, edges):
    """Independently count proper labelled 3-edge-colorings."""
    vertices = tuple(vertices)
    edges = tuple(sorted(edges))
    used = {v: set() for v in vertices}
    count = 0

    def search(index):
        nonlocal count
        if index == len(edges):
            count += 1
            return
        u, v = edges[index]
        for color in COLORS:
            if color in used[u] or color in used[v]:
                continue
            used[u].add(color)
            used[v].add(color)
            search(index + 1)
            used[u].remove(color)
            used[v].remove(color)

    search(0)
    return count


def multiplicity_histogram(counts):
    return tuple(sorted(Counter(counts.values()).items()))


def graph_boundary_audit(name, vertices, edges, side, boundary_order):
    complement = set(vertices) - set(side)
    left = boundary_color_counts(vertices, edges, side, boundary_order)
    right = boundary_color_counts(vertices, edges, complement, boundary_order)
    dot = sum(left[word] * right[word] for word in set(left) | set(right))
    direct = full_coloring_count(vertices, edges)
    assert dot == direct
    print(
        f"{name}: cut={list(boundary_order)}; "
        f"support=({len(left)},{len(right)}); "
        f"mass=({sum(left.values())},{sum(right.values())}); "
        f"dot={dot}; "
        f"multiplicity_hist=({multiplicity_histogram(left)},"
        f"{multiplicity_histogram(right)})"
    )
    return left, right


def component_count(vertices, edges):
    parent = {v: v for v in vertices}

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    def union(u, v):
        u = find(u)
        v = find(v)
        if u != v:
            parent[v] = u

    for u, v in edges:
        union(u, v)
    return len({find(v) for v in vertices})


def survivor_count_by_bits(vertices, edges):
    """Brute-force eps_1=0 and equality along edges."""
    vertices = tuple(vertices)
    assert vertices[0] == 1
    survivors = 0
    for mask in range(1 << (len(vertices) - 1)):
        epsilon = {1: 0}
        for index, v in enumerate(vertices[1:]):
            epsilon[v] = (mask >> index) & 1
        if all(epsilon[u] == epsilon[v] for u, v in edges):
            survivors += 1
    return survivors


def exhaustive_component_law(max_q=6):
    graph_cases = 0
    bit_gates = 0
    for q in range(2, max_q + 1):
        odds = tuple(range(1, 2 * q, 2))
        possible = tuple(combinations(odds, 2))
        for mask in range(1 << len(possible)):
            edges = tuple(
                edge for index, edge in enumerate(possible) if (mask >> index) & 1
            )
            expected = 1 << (component_count(odds, edges) - 1)
            observed = survivor_count_by_bits(odds, edges)
            assert observed == expected
            graph_cases += 1
            bit_gates += 1 << (q - 1)
    return graph_cases, bit_gates


def ap13_six_edge_census():
    odds = tuple(range(1, 14, 2))
    possible = tuple(combinations(odds, 2))
    hist = Counter()
    for edges in combinations(possible, 6):
        hist[1 << (component_count(odds, edges) - 1)] += 1
    assert sum(hist.values()) == 54264
    assert hist == Counter({1: 16807, 2: 32417, 4: 5005, 8: 35})
    return odds, possible, hist


def ap13_short_edge_floor(odds, possible):
    gates = 0
    minima = {}
    for edge_count in range(6):
        minimum = None
        for edges in combinations(possible, edge_count):
            sheets = 1 << (component_count(odds, edges) - 1)
            minimum = sheets if minimum is None else min(minimum, sheets)
            assert sheets >= 1 << (6 - edge_count)
            gates += 1
        minima[edge_count] = minimum
    assert minima == {0: 64, 1: 32, 2: 16, 3: 8, 4: 4, 5: 2}
    return gates, minima


def ap_phase_residual_audit(q):
    """Verify the explicit THM-4110 half-toggle packets modulo 2."""
    speeds = tuple(range(1, 2 * q))
    odds = tuple(v for v in speeds if v % 2)
    unequal_v2 = tuple(
        (u, v)
        for u, v in combinations(speeds, 2)
        if two_adic_valuation(u) != two_adic_valuation(v)
    )
    gates = 0
    for mask in range(1 << (q - 1)):
        epsilon = {v: 0 for v in speeds}
        for index, odd in enumerate(odds[1:]):
            epsilon[odd] = (mask >> index) & 1
        for u, v in unequal_v2:
            # Twice the residual commutator is v*eps_u-u*eps_v.
            assert (v * epsilon[u] - u * epsilon[v]) % 2 == 0
            gates += 1
    return len(unequal_v2), 1 << (q - 1), gates


def two_adic_valuation(value):
    valuation = 0
    while value % 2 == 0:
        value //= 2
        valuation += 1
    return valuation


def main():
    k4_vertices = tuple(range(4))
    k4_edges = tuple(combinations(k4_vertices, 2))
    k4_cut = ((0, 2), (0, 3), (1, 2), (1, 3))
    k4_left, k4_right = graph_boundary_audit(
        "K4", k4_vertices, k4_edges, {0, 1}, k4_cut
    )
    assert len(set(k4_left) & set(k4_right)) == 6

    petersen_vertices = tuple(range(10))
    petersen_edges = (
        (0, 1), (0, 4), (0, 5), (1, 2), (1, 6),
        (2, 3), (2, 7), (3, 4), (3, 8), (4, 9),
        (5, 7), (5, 8), (6, 8), (6, 9), (7, 9),
    )
    petersen_cut = ((0, 4), (0, 5), (1, 2), (1, 6))
    pet_left, pet_right = graph_boundary_audit(
        "Petersen", petersen_vertices, petersen_edges, {0, 1}, petersen_cut
    )
    expected_left = {
        (a, b, c, d)
        for a in COLORS for b in COLORS for c in COLORS for d in COLORS
        if a != b and c != d and {a, b} == {c, d}
    }
    expected_right = {(a, a, b, b) for a in COLORS for b in COLORS}
    assert set(pet_left) == expected_left
    assert set(pet_right) == expected_right
    assert set(pet_left).isdisjoint(pet_right)

    graph_cases, bit_gates = exhaustive_component_law()
    odds, possible, sheet_hist = ap13_six_edge_census()
    short_gates, minima = ap13_short_edge_floor(odds, possible)
    edge_count, phase_sheets, phase_gates = ap_phase_residual_audit(7)

    star = tuple((1, odd) for odd in odds[1:])
    path = tuple(zip(odds[:-1], odds[1:]))
    hostile = ((1, 3), (1, 5), (1, 7), (1, 9), (1, 11), (3, 5))
    assert survivor_count_by_bits(odds, star) == 1
    assert survivor_count_by_bits(odds, path) == 1
    assert survivor_count_by_bits(odds, hostile) == 2

    print(
        "AP13: six_subsets=54264; "
        f"sheet_hist={tuple(sorted(sheet_hist.items()))}"
    )
    print(
        "AP component law: "
        f"q=2..6 graph_cases={graph_cases}; bit_gates={bit_gates}; "
        f"AP13_short_subset_gates={short_gates}; minima={tuple(minima.items())}"
    )
    print(
        "AP13 phase packets: "
        f"unequal_v2_edges={edge_count}; sheets={phase_sheets}; "
        f"commutator_gates={phase_gates}; "
        "star=1; path=1; cyclic_hostile=2"
    )


if __name__ == "__main__":
    main()
