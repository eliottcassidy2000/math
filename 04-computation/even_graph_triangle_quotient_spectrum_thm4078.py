#!/usr/bin/env python3
"""Primal orbit-graph audit for THM-4078.

This path constructs the binary cycle space from anchored triangles,
canonicalizes every Eulerian graph under all vertex permutations through
n=6, and builds the multiplicity-weighted triangle and four-cycle operators
directly on isomorphism classes.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from math import factorial

import sympy


N_MAX = 6
EXPECTED_ORBIT_COUNTS = {3: 2, 4: 3, 5: 7, 6: 16}
EXPECTED_CHARPOLY = {
    3: "(lambda - 1)*(lambda + 1)",
    4: "lambda*(lambda - 4)*(lambda + 4)",
    5: "lambda*(lambda - 10)*(lambda - 4)*(lambda - 2)*(lambda + 2)*(lambda + 4)*(lambda + 10)",
    6: "lambda**4*(lambda - 20)*(lambda - 12)*(lambda - 8)*(lambda - 4)**3*(lambda + 4)**3*(lambda + 8)*(lambda + 12)*(lambda + 20)",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def edge_system(n: int) -> tuple[list[tuple[int, int]], dict[tuple[int, int], int]]:
    edges = list(itertools.combinations(range(n), 2))
    return edges, {edge: index for index, edge in enumerate(edges)}


def edge_mask(edge_index: dict[tuple[int, int], int], pairs: list[tuple[int, int]]) -> int:
    mask = 0
    for a, b in pairs:
        edge = (a, b) if a < b else (b, a)
        mask |= 1 << edge_index[edge]
    return mask


def simple_cycles(n: int, length: int, edge_index: dict[tuple[int, int], int]) -> list[int]:
    cycles: list[int] = []
    for vertices in itertools.combinations(range(n), length):
        first = vertices[0]
        for tail in itertools.permutations(vertices[1:]):
            if tail[0] > tail[-1]:
                continue
            order = (first,) + tail
            pairs = [
                (order[index], order[(index + 1) % length])
                for index in range(length)
            ]
            cycles.append(edge_mask(edge_index, pairs))
    expected = factorial(n) // (2 * length * factorial(n - length))
    require(len(cycles) == expected and len(set(cycles)) == expected, f"cycle generation failed n={n},k={length}")
    return cycles


def permutation_edge_maps(n: int, edges: list[tuple[int, int]], edge_index: dict[tuple[int, int], int]) -> list[tuple[int, ...]]:
    maps: list[tuple[int, ...]] = []
    for permutation in itertools.permutations(range(n)):
        image: list[int] = []
        for a, b in edges:
            x, y = permutation[a], permutation[b]
            image.append(edge_index[(x, y) if x < y else (y, x)])
        maps.append(tuple(image))
    return maps


def permute_mask(mask: int, edge_map: tuple[int, ...]) -> int:
    image = 0
    remaining = mask
    while remaining:
        low = remaining & -remaining
        source = low.bit_length() - 1
        image |= 1 << edge_map[source]
        remaining ^= low
    return image


def cycle_space(n: int, edge_index: dict[tuple[int, int], int]) -> list[int]:
    basis = [
        edge_mask(edge_index, [(0, i), (i, j), (j, 0)])
        for i in range(1, n)
        for j in range(i + 1, n)
    ]
    states = [0]
    for vector in basis:
        states += [state ^ vector for state in states]
    require(len(states) == 2 ** ((n - 1) * (n - 2) // 2), f"cycle-space size failed n={n}")
    return states


def quotient_data(n: int) -> tuple[list[int], dict[int, int], list[int], int]:
    edges, edge_index = edge_system(n)
    states = cycle_space(n, edge_index)
    maps = permutation_edge_maps(n, edges, edge_index)
    canonical: dict[int, int] = {}
    relabel_gates = 0
    for state in states:
        best = None
        for edge_map in maps:
            image = permute_mask(state, edge_map)
            best = image if best is None or image < best else best
            relabel_gates += 1
        require(best is not None, "empty permutation family")
        canonical[state] = best
    representatives = sorted(set(canonical.values()))
    orbit_index = {representative: index for index, representative in enumerate(representatives)}
    state_to_orbit = {state: orbit_index[canonical[state]] for state in states}
    orbit_sizes = [0] * len(representatives)
    for state in states:
        orbit_sizes[state_to_orbit[state]] += 1
    return representatives, state_to_orbit, orbit_sizes, relabel_gates


def weighted_operator(representatives: list[int], state_to_orbit: dict[int, int], generators: list[int]) -> list[list[int]]:
    size = len(representatives)
    matrix = [[0] * size for _ in range(size)]
    for row, representative in enumerate(representatives):
        for generator in generators:
            matrix[row][state_to_orbit[representative ^ generator]] += 1
    return matrix


def boolean_support(matrix: list[list[int]]) -> list[list[int]]:
    return [
        [1 if row != column and matrix[row][column] else 0 for column in range(len(matrix))]
        for row in range(len(matrix))
    ]


def matrix_product_entry(left: list[list[int]], right: list[list[int]], row: int, column: int) -> int:
    return sum(left[row][middle] * right[middle][column] for middle in range(len(left)))


def main() -> None:
    total_relabel_gates = 0
    rows: list[dict[str, object]] = []

    for n in range(3, N_MAX + 1):
        _, edge_index = edge_system(n)
        representatives, state_to_orbit, orbit_sizes, relabel_gates = quotient_data(n)
        total_relabel_gates += relabel_gates
        require(len(representatives) == EXPECTED_ORBIT_COUNTS[n], f"orbit count changed n={n}")

        triangles = simple_cycles(n, 3, edge_index)
        m3 = weighted_operator(representatives, state_to_orbit, triangles)
        degree = len(triangles)
        require(all(sum(row) == degree for row in m3), f"row sum failed n={n}")
        for i in range(len(m3)):
            for j in range(len(m3)):
                require(orbit_sizes[i] * m3[i][j] == orbit_sizes[j] * m3[j][i], f"reversibility failed n={n},i={i},j={j}")

        symbolic = sympy.Matrix(m3)
        polynomial = sympy.factor(symbolic.charpoly().as_expr())
        polynomial_text = str(polynomial)
        require(polynomial_text == EXPECTED_CHARPOLY[n], f"characteristic polynomial changed n={n}: {polynomial_text}")
        eigenvalues = symbolic.eigenvals()
        second = degree - 2 * (n - 2)
        require(eigenvalues.get(degree) == 1, f"top eigenvalue failed n={n}")
        require(eigenvalues.get(second) == 1, f"simple quotient gap failed n={n}")
        require(eigenvalues.get(-degree) == 1, f"bipartite eigenvalue failed n={n}")

        commutator_pair: tuple[int, int] | None = None
        if n >= 4:
            four_cycles = simple_cycles(n, 4, edge_index)
            m4 = weighted_operator(representatives, state_to_orbit, four_cycles)
            b3 = boolean_support(m3)
            b4 = boolean_support(m4)
            empty_index = state_to_orbit[0]
            triangle_index = state_to_orbit[triangles[0]]
            left = matrix_product_entry(b3, b4, empty_index, triangle_index)
            right = matrix_product_entry(b4, b3, empty_index, triangle_index)
            require((left, right) == (0, 1), f"Boolean hostile failed n={n}")
            commutator_pair = (left, right)

        rows.append(
            {
                "n": n,
                "states": len(state_to_orbit),
                "orbits": len(representatives),
                "triangle_degree": degree,
                "charpoly": polynomial_text,
                "commutator_0_C3": commutator_pair,
            }
        )

    require(total_relabel_gates == 745164, "relabel aggregate changed")
    digest = hashlib.sha256(
        json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    print("THM-4078 primary primal-orbit audit")
    print(f"n_range=3..{N_MAX} direct_relabel_gates={total_relabel_gates}")
    print(f"rows={json.dumps(rows, sort_keys=True, separators=(',', ':'))}")
    print(f"semantic_sha256={digest}")
    print("PASS: quotient matrices, exact spectra, reversibility, and Boolean hostile")


if __name__ == "__main__":
    main()
