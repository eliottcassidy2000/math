#!/usr/bin/env python3
"""Exact finite companion for THM-2708.

The general theorem is proved representation-theoretically.  This script
exhausts all small C3-invariant binary circulant blocks, checks the F4
Hermitian standard-block formula in the presence of fixed vertices, verifies
gain switching on every labelled forest through five vertices, and audits the
sharp balanced 3 K3 / unbalanced C9 voltage-lift pair.  No floating point
arithmetic or truth-bearing Python ``assert`` is used.
"""

from __future__ import annotations

from itertools import product

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# F4 = F2[w]/(w^2+w+1), encoded as a+b*w in the two bits of an int.
ONE = 1
OMEGA = 2


def f4_add(a: int, b: int) -> int:
    return a ^ b


def f4_mul(a: int, b: int) -> int:
    a0, a1 = a & 1, (a >> 1) & 1
    b0, b1 = b & 1, (b >> 1) & 1
    constant = (a0 * b0) ^ (a1 * b1)
    omega = (a0 * b1) ^ (a1 * b0) ^ (a1 * b1)
    return constant | (omega << 1)


def f4_pow(a: int, exponent: int) -> int:
    result = ONE
    base = a
    exponent = exponent % 3 if a != 0 else exponent
    while exponent:
        if exponent & 1:
            result = f4_mul(result, base)
        base = f4_mul(base, base)
        exponent >>= 1
    return result


def f4_conjugate(a: int) -> int:
    return f4_mul(a, a)


def f4_inverse(a: int) -> int:
    require(a != 0, "F4 inverse of zero")
    return f4_conjugate(a)


def f4_rank(matrix: list[list[int]]) -> int:
    if not matrix:
        return 0
    work = [row[:] for row in matrix]
    rows, columns = len(work), len(work[0])
    rank = 0
    for column in range(columns):
        pivot = next((i for i in range(rank, rows) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = f4_inverse(work[rank][column])
        work[rank] = [f4_mul(inverse, value) for value in work[rank]]
        for i in range(rows):
            if i == rank or work[i][column] == 0:
                continue
            factor = work[i][column]
            work[i] = [
                f4_add(left, f4_mul(factor, right))
                for left, right in zip(work[i], work[rank])
            ]
        rank += 1
    return rank


def f4_det(matrix: list[list[int]]) -> int:
    n = len(matrix)
    work = [row[:] for row in matrix]
    determinant = ONE
    for column in range(n):
        pivot = next((i for i in range(column, n) if work[i][column]), None)
        if pivot is None:
            return 0
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        determinant = f4_mul(determinant, pivot_value)
        inverse = f4_inverse(pivot_value)
        work[column] = [f4_mul(inverse, value) for value in work[column]]
        for i in range(column + 1, n):
            factor = work[i][column]
            if factor:
                work[i] = [
                    f4_add(left, f4_mul(factor, right))
                    for left, right in zip(work[i], work[column])
                ]
    return determinant


def gf2_rank(rows: list[int], columns: int) -> int:
    work = rows[:]
    rank = 0
    for column in range(columns):
        pivot = next((i for i in range(rank, len(work))
                      if (work[i] >> column) & 1), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for i in range(len(work)):
            if i != rank and ((work[i] >> column) & 1):
                work[i] ^= work[rank]
        rank += 1
    return rank


def gf2_rows(matrix: list[list[int]]) -> list[int]:
    return [sum((entry & 1) << j for j, entry in enumerate(row))
            for row in matrix]


def evaluate_circulant(coefficients: tuple[int, int, int]) -> int:
    value = 0
    for exponent, coefficient in enumerate(coefficients):
        if coefficient:
            value = f4_add(value, f4_pow(OMEGA, exponent))
    return value


def build_full_matrix(
    fixed_block: tuple[tuple[int, ...], ...],
    couplings: tuple[tuple[int, ...], ...],
    self_coefficients: tuple[tuple[int, int, int], ...],
    off_coefficients: dict[tuple[int, int], tuple[int, int, int]],
) -> list[list[int]]:
    fixed_count = len(fixed_block)
    free_count = len(self_coefficients)
    size = fixed_count + 3 * free_count
    matrix = [[0] * size for _ in range(size)]

    for i in range(fixed_count):
        for j in range(fixed_count):
            matrix[i][j] = fixed_block[i][j] & 1

    def node(orbit: int, sheet: int) -> int:
        return fixed_count + 3 * orbit + sheet

    for f in range(fixed_count):
        for orbit in range(free_count):
            for sheet in range(3):
                matrix[f][node(orbit, sheet)] = couplings[f][orbit] & 1
                matrix[node(orbit, sheet)][f] = couplings[f][orbit] & 1

    for orbit, coefficients in enumerate(self_coefficients):
        require(coefficients[1] == coefficients[2], "symmetric self block")
        for source in range(3):
            for target in range(3):
                matrix[node(orbit, source)][node(orbit, target)] = (
                    coefficients[(target - source) % 3]
                )

    for (left, right), coefficients in off_coefficients.items():
        require(left < right, "oriented free-orbit pair")
        for source in range(3):
            for target in range(3):
                value = coefficients[(target - source) % 3]
                matrix[node(left, source)][node(right, target)] = value
                matrix[node(right, target)][node(left, source)] = value

    for i in range(size):
        for j in range(size):
            require(matrix[i][j] == matrix[j][i], "full matrix symmetry")
    return matrix


def hermitian_block(
    self_coefficients: tuple[tuple[int, int, int], ...],
    off_coefficients: dict[tuple[int, int], tuple[int, int, int]],
) -> list[list[int]]:
    count = len(self_coefficients)
    block = [[0] * count for _ in range(count)]
    for i, coefficients in enumerate(self_coefficients):
        block[i][i] = evaluate_circulant(coefficients)
    for (left, right), coefficients in off_coefficients.items():
        value = evaluate_circulant(coefficients)
        block[left][right] = value
        block[right][left] = f4_conjugate(value)
    for i in range(count):
        for j in range(count):
            require(block[j][i] == f4_conjugate(block[i][j]),
                    "F4 Hermitian symmetry")
    return block


def standard_kernel_dimension(
    full_matrix: list[list[int]], fixed_count: int, free_count: int
) -> int:
    size = fixed_count + 3 * free_count
    norm_rows: list[int] = []
    for fixed in range(fixed_count):
        norm_rows.append(1 << fixed)
    for orbit in range(free_count):
        mask = sum(1 << (fixed_count + 3 * orbit + sheet)
                   for sheet in range(3))
        norm_rows.extend([mask] * 3)
    stacked = gf2_rows(full_matrix) + norm_rows
    return size - gf2_rank(stacked, size)


def verify_case(
    fixed_block: tuple[tuple[int, ...], ...],
    couplings: tuple[tuple[int, ...], ...],
    self_coefficients: tuple[tuple[int, int, int], ...],
    off_coefficients: dict[tuple[int, int], tuple[int, int, int]],
) -> int:
    full = build_full_matrix(
        fixed_block, couplings, self_coefficients, off_coefficients
    )
    block = hermitian_block(self_coefficients, off_coefficients)
    observed = standard_kernel_dimension(
        full, len(fixed_block), len(self_coefficients)
    )
    predicted = 2 * (len(block) - f4_rank(block))
    require(observed == predicted, "Hermitian standard-block formula")
    determinant = f4_det(block)
    require(determinant == f4_conjugate(determinant),
            "Hermitian determinant belongs to F2")
    return observed


def symmetric_binary_matrices(size: int):
    positions = [(i, j) for i in range(size) for j in range(i, size)]
    for values in product((0, 1), repeat=len(positions)):
        matrix = [[0] * size for _ in range(size)]
        for (i, j), value in zip(positions, values):
            matrix[i][j] = matrix[j][i] = value
        yield tuple(tuple(row) for row in matrix)


def prufer_tree(sequence: tuple[int, ...], size: int) -> tuple[tuple[int, int], ...]:
    if size == 1:
        return ()
    degree = [1] * size
    for vertex in sequence:
        degree[vertex] += 1
    edges: list[tuple[int, int]] = []
    for vertex in sequence:
        leaf = next(i for i in range(size) if degree[i] == 1)
        edges.append((min(leaf, vertex), max(leaf, vertex)))
        degree[leaf] -= 1
        degree[vertex] -= 1
    remaining = [i for i in range(size) if degree[i] == 1]
    require(len(remaining) == 2, "Pruefer terminal pair")
    edges.append((min(remaining), max(remaining)))
    return tuple(sorted(edges))


def gain_matrix(size: int, edge_gains: dict[tuple[int, int], int]) -> list[list[int]]:
    matrix = [[0] * size for _ in range(size)]
    for (left, right), exponent in edge_gains.items():
        value = f4_pow(OMEGA, exponent)
        matrix[left][right] = value
        matrix[right][left] = f4_conjugate(value)
    return matrix


def gauge_matrix(matrix: list[list[int]], heights: tuple[int, ...]) -> list[list[int]]:
    size = len(matrix)
    result = [[0] * size for _ in range(size)]
    for i in range(size):
        for j in range(size):
            left = f4_pow(OMEGA, heights[i])
            right = f4_pow(OMEGA, -heights[j])
            result[i][j] = f4_mul(left, f4_mul(matrix[i][j], right))
    return result


def connected_components(adjacency: list[list[int]]) -> list[tuple[int, ...]]:
    unseen = set(range(len(adjacency)))
    components: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        stack = [start]
        unseen.remove(start)
        component = []
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbor, value in enumerate(adjacency[vertex]):
                if value and neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return components


def triangle_lift(gains: tuple[int, int, int]) -> list[list[int]]:
    adjacency = [[0] * 9 for _ in range(9)]

    def node(orbit: int, sheet: int) -> int:
        return 3 * orbit + (sheet % 3)

    for left, right, gain in ((0, 1, gains[0]),
                              (1, 2, gains[1]),
                              (2, 0, gains[2])):
        for sheet in range(3):
            u, v = node(left, sheet), node(right, sheet + gain)
            adjacency[u][v] = adjacency[v][u] = 1
    return adjacency


def smith_diagonal(matrix: list[list[int]]) -> tuple[int, ...]:
    normal = smith_normal_form(sp.Matrix(matrix), domain=sp.ZZ)
    return tuple(sorted(abs(int(normal[i, i])) for i in range(normal.rows)))


def main() -> None:
    require(f4_mul(OMEGA, OMEGA) == 3, "omega square")
    require(f4_mul(OMEGA, 3) == 1, "omega cube")
    require(f4_conjugate(OMEGA) == 3, "F4 conjugation")

    self_options = (
        (0, 0, 0),
        (1, 0, 0),
        (0, 1, 1),
        (1, 1, 1),
    )
    off_options = tuple(product((0, 1), repeat=3))

    free_case_count = 0
    free_positive_count = 0
    for free_count in range(1, 4):
        pairs = tuple((i, j) for i in range(free_count)
                      for j in range(i + 1, free_count))
        for diagonal in product(self_options, repeat=free_count):
            for off_values in product(off_options, repeat=len(pairs)):
                off = dict(zip(pairs, off_values))
                observed = verify_case((), (), diagonal, off)
                free_case_count += 1
                free_positive_count += observed > 0

    fixed_case_count = 0
    for fixed_block in symmetric_binary_matrices(2):
        for coupling_values in product((0, 1), repeat=4):
            couplings = (
                coupling_values[:2],
                coupling_values[2:],
            )
            for diagonal in product(self_options, repeat=2):
                for coefficients in off_options:
                    verify_case(
                        fixed_block,
                        couplings,
                        diagonal,
                        {(0, 1): coefficients},
                    )
                    fixed_case_count += 1

    forest_case_count = 0
    for size in range(1, 6):
        sequences = ((),) if size <= 2 else product(range(size), repeat=size - 2)
        for sequence in sequences:
            edges = prufer_tree(tuple(sequence), size)
            for gains in product(range(3), repeat=len(edges)):
                edge_gains = dict(zip(edges, gains))
                matrix = gain_matrix(size, edge_gains)
                heights = [None] * size
                heights[0] = 0
                stack = [0]
                while stack:
                    parent = stack.pop()
                    for left, right in edges:
                        if parent not in (left, right):
                            continue
                        child = right if parent == left else left
                        if heights[child] is not None:
                            continue
                        gain = edge_gains[(left, right)]
                        if parent == left:
                            heights[child] = (heights[parent] + gain) % 3
                        else:
                            heights[child] = (heights[parent] - gain) % 3
                        stack.append(child)
                require(all(value is not None for value in heights),
                        "tree gauge potentials")
                switched = gauge_matrix(matrix, tuple(int(x) for x in heights))
                for i in range(size):
                    for j in range(size):
                        expected = 1 if (min(i, j), max(i, j)) in edge_gains else 0
                        require(switched[i][j] == expected,
                                "forest gains switch to one")
                forest_case_count += 1

    balanced_count = 0
    unbalanced_count = 0
    balanced_matrix = None
    unbalanced_matrix = None
    for gains in product(range(3), repeat=3):
        g01, g12, g20 = (f4_pow(OMEGA, exponent) for exponent in gains)
        block = [
            [0, g01, f4_conjugate(g20)],
            [f4_conjugate(g01), 0, g12],
            [g20, f4_conjugate(g12), 0],
        ]
        holonomy = sum(gains) % 3
        determinant = f4_det(block)
        require(determinant == (0 if holonomy == 0 else 1),
                "triangle trace-holonomy determinant")
        for heights in product(range(3), repeat=3):
            require(f4_det(gauge_matrix(block, heights)) == determinant,
                    "gain switching preserves determinant")

        lift = triangle_lift(gains)
        components = connected_components(lift)
        degree_word = tuple(sum(row) for row in lift)
        require(degree_word == (2,) * 9, "triangle lift is two-regular")
        determinant_integer = abs(int(sp.Matrix(lift).det()))
        standard_dimension = standard_kernel_dimension(lift, 0, 3)
        if holonomy == 0:
            require(tuple(sorted(map(len, components))) == (3, 3, 3),
                    "balanced lift is three triangles")
            require(determinant_integer == 8, "balanced lift determinant")
            require(standard_dimension == 2, "balanced lift one standard plane")
            balanced_count += 1
            balanced_matrix = lift
        else:
            require(tuple(sorted(map(len, components))) == (9,),
                    "unbalanced lift is one nine-cycle")
            require(determinant_integer == 2, "unbalanced lift determinant")
            require(standard_dimension == 0, "unbalanced lift no standard plane")
            unbalanced_count += 1
            unbalanced_matrix = lift

    require(balanced_count == 9 and unbalanced_count == 18,
            "triangle holonomy census")
    require(balanced_matrix is not None and unbalanced_matrix is not None,
            "sharp controls exist")
    balanced_smith = smith_diagonal(balanced_matrix)
    unbalanced_smith = smith_diagonal(unbalanced_matrix)
    require(balanced_smith == (1, 1, 1, 1, 1, 1, 2, 2, 2),
            "three-triangle Smith form")
    require(unbalanced_smith == (1, 1, 1, 1, 1, 1, 1, 1, 2),
            "nine-cycle Smith form")

    print("THM-2708 C3 HERMITIAN GAIN-HOLONOMY AUDIT")
    print(f"free_orbit_matrices_checked={free_case_count}")
    print(f"free_matrices_with_standard_kernel={free_positive_count}")
    print(f"fixed_coupling_matrices_checked={fixed_case_count}")
    print(f"forest_gain_assignments_checked={forest_case_count}")
    print("triangle_gain_assignments=27 balanced=9 unbalanced=18")
    print("balanced_lift=3K3 det=8 smith=(1^6,2^3) standard_multiplicity=1")
    print("unbalanced_lift=C9 det=2 smith=(1^8,2) standard_multiplicity=0")
    print("standard_multiplicity=F4_nullity_of_Hermitian_block")
    print("representative_switching=diagonal_F4_gauge")
    print("tree_case=balanced_gain_specialization")
    print("surface_gate_requires_full_rank_boundary_lattice")
    print("gain_holonomy_is_not_physical_LRC_homology")
    print("PASS")


if __name__ == "__main__":
    main()
