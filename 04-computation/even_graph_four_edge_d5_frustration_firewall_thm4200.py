#!/usr/bin/env python3
"""Primary exact audit for THM-4200.

Enumerate the eleven unlabeled four-edge simple-graph types, derive their
cycle-parity profiles by exact inclusion--exclusion, and certify the D=5 gap
by nonnegative coefficients after the first admissible order is shifted to
zero.  Literal simple-cycle counts provide finite hostile controls.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations, permutations

import sympy as sp


N = sp.symbols("n", integer=True, positive=True)
T = sp.symbols("t", integer=True, nonnegative=True)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def edge_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def active_vertices(edges: tuple[tuple[int, int], ...]) -> tuple[int, ...]:
    return tuple(sorted({v for edge in edges for v in edge}))


def components(
    edges: tuple[tuple[int, int], ...],
) -> list[tuple[set[int], tuple[tuple[int, int], ...]]]:
    vertices = set(active_vertices(edges))
    adjacency = {v: set() for v in vertices}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    result = []
    while vertices:
        root = min(vertices)
        stack = [root]
        block = {root}
        vertices.remove(root)
        while stack:
            v = stack.pop()
            for w in adjacency[v]:
                if w in vertices:
                    vertices.remove(w)
                    block.add(w)
                    stack.append(w)
        block_edges = tuple(edge for edge in edges if edge[0] in block)
        result.append((block, block_edges))
    return result


def structure_signature(
    edges: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, int, tuple[int, ...]], ...]:
    pieces = []
    for vertices, block_edges in components(edges):
        degree = defaultdict(int)
        for a, b in block_edges:
            degree[a] += 1
            degree[b] += 1
        pieces.append(
            (len(vertices), len(block_edges), tuple(sorted(degree.values(), reverse=True)))
        )
    return tuple(sorted(pieces))


SIGNATURE_NAMES = {
    ((5, 4, (4, 1, 1, 1, 1)),): "K1,4",
    ((4, 4, (3, 2, 2, 1)),): "paw",
    ((5, 4, (3, 2, 1, 1, 1)),): "fork",
    ((2, 1, (1, 1)), (4, 3, (3, 1, 1, 1))): "K1,3+K2",
    ((2, 1, (1, 1)), (3, 3, (2, 2, 2))): "K3+K2",
    ((4, 4, (2, 2, 2, 2)),): "C4",
    ((5, 4, (2, 2, 2, 1, 1)),): "P5",
    ((2, 1, (1, 1)), (4, 3, (2, 2, 1, 1))): "P4+K2",
    ((3, 2, (2, 1, 1)), (3, 2, (2, 1, 1))): "2P3",
    ((2, 1, (1, 1)), (2, 1, (1, 1)), (3, 2, (2, 1, 1))): "P3+2K2",
    (
        (2, 1, (1, 1)),
        (2, 1, (1, 1)),
        (2, 1, (1, 1)),
        (2, 1, (1, 1)),
    ): "4K2",
}


def four_edge_types() -> list[tuple[str, tuple[tuple[int, int], ...]]]:
    complete_edges = tuple(combinations(range(8), 2))
    representatives: dict[
        tuple[tuple[int, int, tuple[int, ...]], ...], tuple[tuple[int, int], ...]
    ] = {}
    multiplicities = defaultdict(int)
    for indices in combinations(range(len(complete_edges)), 4):
        edges = tuple(complete_edges[i] for i in indices)
        signature = structure_signature(edges)
        representatives.setdefault(signature, edges)
        multiplicities[signature] += 1
    require(set(representatives) == set(SIGNATURE_NAMES), "four-edge type inventory mismatch")
    require(sum(multiplicities.values()) == sp.binomial(28, 4), "K8 inventory count mismatch")
    return [
        (SIGNATURE_NAMES[signature], representatives[signature])
        for signature in SIGNATURE_NAMES
    ]


def containing_cycles(edges: tuple[tuple[int, int], ...], k: int) -> sp.Expr:
    """Number of unoriented simple k-cycles in K_n containing all given edges."""

    j = len(edges)
    pieces = components(edges)
    degree = defaultdict(int)
    for a, b in edges:
        degree[a] += 1
        degree[b] += 1
    if max(degree.values()) > 2:
        return sp.Integer(0)

    cyclic = [
        (vertices, block_edges)
        for vertices, block_edges in pieces
        if len(vertices) == len(block_edges)
    ]
    if cyclic:
        if len(cyclic) == 1 and len(pieces) == 1 and j == k:
            return sp.Integer(1)
        return sp.Integer(0)

    component_count = len(pieces)
    vertex_count = len(degree)
    if k < vertex_count:
        return sp.Integer(0)
    return (
        sp.Integer(2) ** (component_count - 1)
        * sp.factorial(k - j - 1)
        * sp.binomial(N - vertex_count, k - vertex_count)
    )


def negative_k_cycles(edges: tuple[tuple[int, int], ...], k: int) -> sp.Expr:
    total = sp.Integer(0)
    for size in range(1, len(edges) + 1):
        coefficient = (-2) ** (size - 1)
        for subset in combinations(edges, size):
            total += coefficient * containing_cycles(tuple(subset), k)
    return sp.expand_func(total).expand()


def simple_cycles(n: int, k: int) -> tuple[frozenset[tuple[int, int]], ...]:
    cycles: set[frozenset[tuple[int, int]]] = set()
    for chosen in combinations(range(n), k):
        root = chosen[0]
        rest = chosen[1:]
        for order in permutations(rest):
            tour = (root,) + order
            reverse = (root,) + tuple(reversed(order))
            if tour > reverse:
                continue
            cycles.add(
                frozenset(edge_key(tour[i], tour[(i + 1) % k]) for i in range(k))
            )
    return tuple(cycles)


def direct_negative_count(
    edges: tuple[tuple[int, int], ...], cycles: tuple[frozenset[tuple[int, int]], ...]
) -> int:
    edge_set = frozenset(edges)
    return sum(len(edge_set.intersection(cycle)) % 2 for cycle in cycles)


def main() -> None:
    types = four_edge_types()
    candidate = (N - 2) * (N**3 - 11 * N**2 + 41 * N - 50)
    print("THM-4200 primary exact audit")
    print(f"four-edge isomorphism types={len(types)}; K8 labelled sets={sp.binomial(28, 4)}")
    print("method=parity inclusion-exclusion + exact linear-forest contraction")
    print()

    equality_rows = []
    for name, edges in types:
        profile = {k: sp.factor(negative_k_cycles(edges, k)) for k in range(3, 7)}
        cumulative = sp.factor(sum(profile.values()))
        gap = sp.factor(cumulative - candidate)
        first_n = max(6, len(active_vertices(edges)))
        shifted = sp.Poly(sp.expand(gap.subs(N, T + first_n)), T)
        coefficients = shifted.all_coeffs()
        require(all(coefficient >= 0 for coefficient in coefficients), f"negative shift coefficient: {name}")
        zeros = [value for value in range(first_n, first_n + 9) if gap.subs(N, value) == 0]
        if shifted.eval(0) == 0:
            require(name == "K1,4" and first_n == 6, "unexpected zero boundary")
            require(all(coefficient > 0 for coefficient in coefficients[:-1]), "nonunique star zero")
            equality_rows.append((name, first_n))
        else:
            require(shifted.eval(0) > 0, f"nonstrict boundary: {name}")
            require(not zeros, f"unexpected finite zero: {name}")

        for order in range(first_n, min(10, first_n + 2) + 1):
            cycle_bank = {k: simple_cycles(order, k) for k in range(3, 7)}
            for k in range(3, 7):
                direct = direct_negative_count(edges, cycle_bank[k])
                symbolic = int(profile[k].subs(N, order))
                require(direct == symbolic, f"direct mismatch {name} n={order} k={k}")

        print(
            f"{name:9s} v={len(active_vertices(edges))} first_n={first_n} "
            f"shift_gap={shifted.as_expr()}"
        )
        print("  " + " ".join(f"c{k}={profile[k]}" for k in range(3, 7)))

    require(equality_rows == [("K1,4", 6)], "equality ledger mismatch")
    star = dict(types)["K1,4"]
    center = 0
    star_at_six = frozenset(star)
    center_cut = frozenset(edge_key(center, vertex) for vertex in range(6) if vertex != center)
    switched = star_at_six.symmetric_difference(center_cut)
    require(len(switched) == 1, "K1,4 boundary does not switch to one edge")

    print()
    print("equality ledger=[K1,4 at n=6]")
    print(f"center switch leaves={tuple(sorted(switched))}")
    print("direct controls=all c3..c6 rows at the first three admissible orders")
    print("result=PASS; every frustration-four class is strictly above A_n")


if __name__ == "__main__":
    main()
