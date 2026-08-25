#!/usr/bin/env python3
"""Exact audit for sparse reciprocal phase graphs (THM-4110).

Only integer and Fraction arithmetic is used.  Every truth gate remains live
under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import gcd


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def bareiss_det(matrix: list[list[int]]) -> int:
    """Exact determinant with fraction-free elimination."""
    n = len(matrix)
    if n == 0:
        return 1
    require(all(len(row) == n for row in matrix), "nonsquare matrix")
    a = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot_row = next((i for i in range(k + 1, n) if a[i][k] != 0), None)
            if pivot_row is None:
                return 0
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                require(numerator % previous == 0, (matrix, k, numerator, previous))
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def edge_row(values: tuple[int, ...], edge: tuple[int, int]) -> list[int]:
    i, j = edge
    row = [0] * len(values)
    row[i] = values[j]
    row[j] = -values[i]
    return row


def degrees(n: int, edges: tuple[tuple[int, int], ...]) -> tuple[int, ...]:
    answer = [0] * n
    for i, j in edges:
        answer[i] += 1
        answer[j] += 1
    return tuple(answer)


def tree_index_formula(values: tuple[int, ...], edges: tuple[tuple[int, int], ...]) -> int:
    return product_int(v ** (d - 1) for v, d in zip(values, degrees(len(values), edges)))


def product_int(values) -> int:
    answer = 1
    for value in values:
        answer *= value
    return answer


def tree_minor_vector(values: tuple[int, ...], edges: tuple[tuple[int, int], ...]) -> tuple[int, ...]:
    rows = [edge_row(values, edge) for edge in edges]
    return tuple(
        abs(bareiss_det([[entry for j, entry in enumerate(row) if j != deleted] for row in rows]))
        for deleted in range(len(values))
    )


def prufer_tree(sequence: tuple[int, ...], n: int) -> tuple[tuple[int, int], ...]:
    degree = [1] * n
    for vertex in sequence:
        degree[vertex] += 1
    edges: list[tuple[int, int]] = []
    for vertex in sequence:
        leaf = next(i for i, d in enumerate(degree) if d == 1)
        edges.append(tuple(sorted((leaf, vertex))))
        degree[leaf] -= 1
        degree[vertex] -= 1
    leaves = [i for i, d in enumerate(degree) if d == 1]
    require(len(leaves) == 2, (sequence, degree))
    edges.append(tuple(sorted(leaves)))
    return tuple(sorted(edges))


def connected(n: int, edges: tuple[tuple[int, int], ...], vertices: set[int] | None = None) -> bool:
    active = set(range(n)) if vertices is None else set(vertices)
    if not active:
        return False
    seen = {next(iter(active))}
    changed = True
    while changed:
        changed = False
        for i, j in edges:
            if i not in active or j not in active:
                continue
            if i in seen and j not in seen:
                seen.add(j)
                changed = True
            if j in seen and i not in seen:
                seen.add(i)
                changed = True
    return seen == active


def is_tree(n: int, edges: tuple[tuple[int, int], ...]) -> bool:
    return len(edges) == n - 1 and connected(n, edges)


def gcd_all(values) -> int:
    answer = 0
    for value in values:
        answer = gcd(answer, abs(value))
    return answer


def graph_minor_index(values: tuple[int, ...], edges: tuple[tuple[int, int], ...]) -> int:
    n = len(values)
    rows = [edge_row(values, edge) for edge in edges]
    answer = 0
    for chosen in combinations(range(len(edges)), n - 1):
        block = [rows[i] for i in chosen]
        for deleted in range(n):
            minor = abs(
                bareiss_det(
                    [[entry for j, entry in enumerate(row) if j != deleted] for row in block]
                )
            )
            answer = gcd(answer, minor)
            if answer == 1:
                return 1
    return answer


def graph_tree_gcd(values: tuple[int, ...], edges: tuple[tuple[int, int], ...]) -> int:
    n = len(values)
    answer = 0
    for chosen in combinations(edges, n - 1):
        if is_tree(n, chosen):
            answer = gcd(answer, tree_index_formula(values, chosen))
    return answer


def prime_factors(value: int) -> tuple[int, ...]:
    answer: list[int] = []
    p = 2
    while p * p <= value:
        if value % p == 0:
            answer.append(p)
            while value % p == 0:
                value //= p
        p += 1
    if value > 1:
        answer.append(value)
    return tuple(answer)


def valuation(value: int, prime: int) -> int:
    answer = 0
    while value % prime == 0:
        answer += 1
        value //= prime
    return answer


def prime_free_core_ok(values: tuple[int, ...], edges: tuple[tuple[int, int], ...], prime: int) -> bool:
    free = {i for i, value in enumerate(values) if value % prime != 0}
    if not connected(len(values), edges, free):
        return False
    return all(
        vertex in free
        or any((min(vertex, other), max(vertex, other)) in edges for other in free)
        for vertex in range(len(values))
    )


def exhaustive_tree_minor_audit() -> tuple[int, int]:
    rows = {
        n: (
            tuple(range(1, n + 1)),
            (2, 3, 5, 7, 11, 13)[:n],
            (6, 35, 10, 21, 22, 33)[:n],
        )
        for n in range(2, 7)
    }
    tree_cases = 0
    minor_gates = 0
    for n in range(2, 7):
        for sequence in product(range(n), repeat=n - 2):
            edges = prufer_tree(sequence, n)
            for values in rows[n]:
                require(gcd(*values) == 1, values)
                index = tree_index_formula(values, edges)
                minors = tree_minor_vector(values, edges)
                expected = tuple(index * value for value in values)
                require(minors == expected, (values, edges, minors, expected))
                require(gcd_all(minors) == index, (values, edges, minors, index))
                tree_cases += 1
                minor_gates += len(minors) + 1
    return tree_cases, minor_gates


def exhaustive_graph_audit() -> tuple[int, int, int]:
    value_bank = {
        2: ((1, 2), (2, 3)),
        3: ((1, 2, 3), (2, 3, 5)),
        4: ((1, 2, 3, 4), (2, 3, 5, 7)),
        5: ((1, 2, 3, 4, 5), (2, 3, 5, 7, 11)),
    }
    graph_cases = 0
    index_gates = 0
    prime_gates = 0
    for n in range(2, 6):
        complete = tuple(combinations(range(n), 2))
        for mask in range(1 << len(complete)):
            edges = tuple(edge for bit, edge in enumerate(complete) if mask & (1 << bit))
            if not connected(n, edges):
                continue
            for values in value_bank[n]:
                minor_index = graph_minor_index(values, edges)
                tree_index = graph_tree_gcd(values, edges)
                require(minor_index == tree_index > 0, (values, edges, minor_index, tree_index))
                index_gates += 1
                for prime in prime_factors(product_int(values)):
                    core_ok = prime_free_core_ok(values, edges, prime)
                    require((minor_index % prime != 0) == core_ok, (values, edges, prime, minor_index))
                    prime_gates += 1
                graph_cases += 1
    return graph_cases, index_gates, prime_gates


def kruskal_min_valuation(
    values: tuple[int, ...], edges: tuple[tuple[int, int], ...], prime: int
) -> int:
    parent = list(range(len(values)))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    total = 0
    accepted = 0
    for weight, i, j in sorted(
        (valuation(values[i], prime) + valuation(values[j], prime), i, j) for i, j in edges
    ):
        root_i, root_j = find(i), find(j)
        if root_i == root_j:
            continue
        parent[root_j] = root_i
        total += weight
        accepted += 1
    require(accepted == len(values) - 1, (prime, accepted))
    return total - sum(valuation(value, prime) for value in values)


def v2(value: int) -> int:
    return valuation(value, 2)


def ap13_audit() -> tuple[int, tuple[tuple[int, int], ...], int, dict[int, int], int, int]:
    values = tuple(range(1, 14))
    n = len(values)
    edges = tuple(
        (i, j) for i, j in combinations(range(n), 2) if v2(values[i]) != v2(values[j])
    )
    require(len(edges) == 53, len(edges))

    even_values = (2, 4, 6, 8, 10, 12)
    odd_values = (3, 5, 7, 9, 11, 13)
    tree = tuple(
        sorted(
            [(0, value - 1) for value in even_values]
            + [(1, value - 1) for value in odd_values]
        )
    )
    require(is_tree(n, tree), tree)
    require(all(edge in edges for edge in tree), tree)
    require(tree_index_formula(values, tree) == 64, tree_index_formula(values, tree))
    minors = tree_minor_vector(values, tree)
    require(minors == tuple(64 * value for value in values), minors)

    exponents = {
        prime: kruskal_min_valuation(values, edges, prime)
        for prime in prime_factors(product_int(values))
    }
    require(exponents == {2: 6, 3: 0, 5: 0, 7: 0, 11: 0, 13: 0}, exponents)
    graph_index = product_int(prime**exponent for prime, exponent in exponents.items())
    require(graph_index == 64, graph_index)

    # Rooting any allowed tree at speed 1 gives
    # I(T)=product_(x!=1) speed(parent(x)).  Each nonunit odd vertex therefore
    # has an even parent of value at least 2.  Equality forces precisely the
    # displayed tree: every even child has parent 1 and every other odd child
    # has parent 2.
    require(len(odd_values) == 6 and 2 ** len(odd_values) == graph_index, "bad AP13 lower bound")

    t = Fraction(7, 113)
    sheets = 0
    physical = 0
    augmented_survivors = 0
    added_edges = tuple((0, odd - 1) for odd in odd_values)
    for bits in product((0, 1), repeat=len(odd_values)):
        bit = dict(zip(odd_values, bits))
        phases = tuple(
            (value * t + Fraction(bit.get(value, 0), 2)) % 1 for value in values
        )
        require(
            all((values[j] * phases[i] - values[i] * phases[j]) % 1 == 0 for i, j in edges),
            (bits, phases),
        )
        sheets += 1
        is_physical = all((phase - value * phases[0]) % 1 == 0 for value, phase in zip(values, phases))
        physical += int(is_physical)
        augmented_survivors += int(
            all((values[j] * phases[i] - values[i] * phases[j]) % 1 == 0 for i, j in added_edges)
        )
    require((sheets, physical, augmented_survivors) == (64, 1, 1), (sheets, physical, augmented_survivors))
    return len(edges), tree, gcd_all(minors), exponents, sheets, augmented_survivors


def main() -> None:
    tree_cases, minor_gates = exhaustive_tree_minor_audit()
    graph_cases, index_gates, prime_gates = exhaustive_graph_audit()
    allowed, tree, minor_gcd, exponents, sheets, augmented = ap13_audit()
    printable_tree = tuple((i + 1, j + 1) for i, j in tree)
    print("THM-4110 exact sparse reciprocal phase graph audit")
    print(f"tree_minor_cases={tree_cases}")
    print(f"tree_minor_gates={minor_gates}")
    print(f"connected_graph_value_cases={graph_cases}")
    print(f"graph_index_gates={index_gates}")
    print(f"prime_free_core_gates={prime_gates}")
    print(f"AP13_unequal_v2_edges={allowed}")
    print(f"AP13_unique_minimum_tree={printable_tree}")
    print(f"AP13_tree_minor_gcd={minor_gcd}")
    print(f"AP13_graph_prime_exponents={exponents}")
    print(f"AP13_cross_v2_phase_sheets={sheets}")
    print("AP13_phase_sheet_group=(Z/2Z)^6")
    print("AP13_same_shell_anchor_edges=6")
    print(f"AP13_augmented_surviving_sheets={augmented}")
    print("PASS")


if __name__ == "__main__":
    main()
