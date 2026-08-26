#!/usr/bin/env python3
"""Independent literal-cycle referee for THM-4200.

This path does not use SymPy, inclusion--exclusion, or the linear-forest
contraction formula.  It enumerates cycles directly, obtains each quartic gap
from five values, and certifies positivity in the Newton binomial basis.
"""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from itertools import combinations, permutations
from math import comb


TYPES: tuple[tuple[str, tuple[tuple[int, int], ...], int], ...] = (
    ("K1,4", ((0, 1), (0, 2), (0, 3), (0, 4)), 6),
    ("paw", ((0, 1), (0, 2), (0, 3), (1, 2)), 6),
    ("fork", ((0, 1), (0, 2), (0, 3), (1, 4)), 6),
    ("K1,3+K2", ((0, 1), (0, 2), (0, 3), (4, 5)), 6),
    ("K3+K2", ((0, 1), (0, 2), (1, 2), (3, 4)), 6),
    ("C4", ((0, 1), (1, 2), (2, 3), (0, 3)), 6),
    ("P5", ((0, 1), (1, 2), (2, 3), (3, 4)), 6),
    ("P4+K2", ((0, 1), (1, 2), (2, 3), (4, 5)), 6),
    ("2P3", ((0, 1), (1, 2), (3, 4), (4, 5)), 6),
    ("P3+2K2", ((0, 1), (1, 2), (3, 4), (5, 6)), 7),
    ("4K2", ((0, 1), (2, 3), (4, 5), (6, 7)), 8),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def edge_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def structure_signature(
    edges: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, int, tuple[int, ...]], ...]:
    active = {vertex for edge in edges for vertex in edge}
    adjacency = {vertex: set() for vertex in active}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    pieces = []
    unseen = set(active)
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        stack = [root]
        vertices = {root}
        while stack:
            vertex = stack.pop()
            for neighbor in adjacency[vertex]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    vertices.add(neighbor)
                    stack.append(neighbor)
        degree = tuple(sorted((len(adjacency[v]) for v in vertices), reverse=True))
        edge_count = sum(degree) // 2
        pieces.append((len(vertices), edge_count, degree))
    return tuple(sorted(pieces))


@lru_cache(maxsize=None)
def cycles(n: int, k: int) -> tuple[frozenset[tuple[int, int]], ...]:
    result: set[frozenset[tuple[int, int]]] = set()
    for vertices in combinations(range(n), k):
        root = vertices[0]
        for tail in permutations(vertices[1:]):
            order = (root,) + tail
            if order > (root,) + tuple(reversed(tail)):
                continue
            result.add(
                frozenset(edge_key(order[i], order[(i + 1) % k]) for i in range(k))
            )
    return tuple(result)


def cumulative_odd_count(edges: tuple[tuple[int, int], ...], n: int) -> int:
    negative = frozenset(edges)
    return sum(
        len(negative.intersection(cycle)) % 2
        for k in range(3, 7)
        for cycle in cycles(n, k)
    )


def candidate(n: int) -> int:
    return (n - 2) * (n**3 - 11 * n**2 + 41 * n - 50)


def forward_coefficients(values: list[int]) -> tuple[int, ...]:
    row = values[:]
    coefficients = []
    while row:
        coefficients.append(row[0])
        row = [row[i + 1] - row[i] for i in range(len(row) - 1)]
    return tuple(coefficients)


def newton_value(coefficients: tuple[int, ...], t: int) -> int:
    return sum(coefficient * comb(t, degree) for degree, coefficient in enumerate(coefficients))


def main() -> None:
    expected_signatures = {structure_signature(edges) for _, edges, _ in TYPES}
    require(len(expected_signatures) == 11, "hard-coded type signatures collide")
    complete_edges = tuple(combinations(range(8), 2))
    observed_signatures = {
        structure_signature(tuple(complete_edges[i] for i in indices))
        for indices in combinations(range(28), 4)
    }
    require(observed_signatures == expected_signatures, "independent K8 type inventory mismatch")

    print("THM-4200 independent literal-cycle audit")
    print("type inventory=11 signatures across all 20,475 four-edge subsets of K8")
    print("basis=Newton binom(t,j); degree bound=4 because every odd cycle contains a fixed edge")
    print()

    equality_rows = []
    for name, edges, first_n in TYPES:
        gaps = [
            cumulative_odd_count(edges, order) - candidate(order)
            for order in range(first_n, first_n + 5)
        ]
        coefficients = forward_coefficients(gaps)
        require(len(coefficients) == 5, "quartic interpolation length mismatch")
        require(all(value >= 0 for value in coefficients), f"negative Newton coefficient: {name}")
        if coefficients[0] == 0:
            require(name == "K1,4" and all(value > 0 for value in coefficients[1:]), "bad equality row")
            equality_rows.append((name, first_n))
        else:
            require(all(value > 0 for value in coefficients), f"nonstrict Newton row: {name}")

        control_t = 5
        direct_control = cumulative_odd_count(edges, first_n + control_t) - candidate(first_n + control_t)
        require(
            direct_control == newton_value(coefficients, control_t),
            f"sixth-point interpolation mismatch: {name}",
        )
        print(
            f"{name:9s} first_n={first_n} gaps={tuple(gaps)} "
            f"newton={coefficients} sixth={direct_control}"
        )

    require(equality_rows == [("K1,4", 6)], "independent equality ledger mismatch")
    star = frozenset(TYPES[0][1])
    center_cut = frozenset((0, vertex) for vertex in range(1, 6))
    switched = star.symmetric_difference(center_cut)
    require(switched == frozenset({(0, 5)}), "star center switch is not a single edge")

    print()
    print("equality ledger=[K1,4 at n=6], center switch={(0,5)}")
    print("result=PASS; literal cycles independently exclude frustration index four")


if __name__ == "__main__":
    main()
