#!/usr/bin/env python3
"""Exact finite companion for THM-2711.

The theorem is a lattice/discriminant-form argument.  This script exhausts
small C3-invariant matrices modulo four, computes the induced order-two
pairing, and verifies the D4 obstruction, the mod-two-identical/mod-four-
different hostile, and a sharp stable-metabolizer positive overlattice.
There is no floating point arithmetic and no truth-bearing Python assert.
"""

from __future__ import annotations

from itertools import product

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def determinant_bareiss(matrix: list[list[int]]) -> int:
    n = len(matrix)
    if n == 0:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for column in range(n - 1):
        pivot = next((i for i in range(column, n) if work[i][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        pivot_value = work[column][column]
        for i in range(column + 1, n):
            for j in range(column + 1, n):
                numerator = (
                    work[i][j] * pivot_value
                    - work[i][column] * work[column][j]
                )
                require(numerator % previous == 0, "Bareiss exact division")
                work[i][j] = numerator // previous
        previous = pivot_value
        for i in range(column + 1, n):
            work[i][column] = 0
    return sign * work[n - 1][n - 1]


def smith_diagonal(matrix: list[list[int]]) -> tuple[int, ...]:
    normal = smith_normal_form(sp.Matrix(matrix), domain=sp.ZZ)
    return tuple(sorted(abs(int(normal[i, i])) for i in range(normal.rows)))


def apply_permutation(vector: int, permutation: tuple[int, ...]) -> int:
    result = 0
    for source, target in enumerate(permutation):
        if (vector >> source) & 1:
            result |= 1 << target
    return result


def invariant_matrix(matrix: list[list[int]], permutation: tuple[int, ...]) -> bool:
    size = len(matrix)
    return all(
        matrix[permutation[i]][permutation[j]] == matrix[i][j]
        for i in range(size)
        for j in range(size)
    )


def kernel_mod_two(matrix: list[list[int]]) -> tuple[int, ...]:
    size = len(matrix)
    kernel = []
    for vector in range(1 << size):
        if all(
            sum(matrix[i][j] * ((vector >> j) & 1) for j in range(size)) % 2 == 0
            for i in range(size)
        ):
            kernel.append(vector)
    return tuple(kernel)


def pairing_bit(matrix: list[list[int]], left: int, right: int) -> int:
    size = len(matrix)
    numerator = sum(
        ((left >> i) & 1) * matrix[i][j] * ((right >> j) & 1)
        for i in range(size)
        for j in range(size)
    )
    require(numerator % 2 == 0, "order-two pairing numerator is even")
    return (numerator // 2) % 2


def standard_planes(
    kernel: tuple[int, ...], permutation: tuple[int, ...]
) -> tuple[frozenset[int], ...]:
    kernel_set = set(kernel)
    planes: set[frozenset[int]] = set()
    for vector in kernel:
        if vector == 0:
            continue
        image = apply_permutation(vector, permutation)
        image2 = apply_permutation(image, permutation)
        if image == vector or vector ^ image ^ image2:
            continue
        plane = frozenset((0, vector, image, vector ^ image))
        require(plane <= kernel_set, "standard plane remains in kernel")
        planes.add(plane)
    return tuple(sorted(planes, key=lambda plane: tuple(sorted(plane))))


def plane_is_isotropic(matrix: list[list[int]], plane: frozenset[int]) -> bool:
    return all(pairing_bit(matrix, left, right) == 0
               for left in plane for right in plane)


def span(generators: tuple[int, ...]) -> frozenset[int]:
    values = {0}
    for generator in generators:
        values |= {value ^ generator for value in tuple(values)}
    return frozenset(values)


def orthogonal(
    matrix: list[list[int]], kernel: tuple[int, ...], subgroup: frozenset[int]
) -> frozenset[int]:
    return frozenset(
        vector
        for vector in kernel
        if all(pairing_bit(matrix, vector, member) == 0 for member in subgroup)
    )


def two_orbit_matrix(parameters: tuple[int, ...]) -> list[list[int]]:
    a0, b0, a1, b1, c0, c1, c2 = parameters
    matrix = [[0] * 6 for _ in range(6)]
    for orbit, (diagonal, off_diagonal) in enumerate(((a0, b0), (a1, b1))):
        for source in range(3):
            for target in range(3):
                matrix[3 * orbit + source][3 * orbit + target] = (
                    diagonal if source == target else off_diagonal
                )
    coefficients = (c0, c1, c2)
    for source in range(3):
        for target in range(3):
            value = coefficients[(target - source) % 3]
            matrix[source][3 + target] = value
            matrix[3 + target][source] = value
    return matrix


def pairing_matrix_on_plane(
    matrix: list[list[int]], plane: frozenset[int]
) -> tuple[tuple[int, int], tuple[int, int]]:
    nonzero = sorted(plane - {0})
    left = nonzero[0]
    right = next(value for value in nonzero[1:] if value != left)
    return (
        (pairing_bit(matrix, left, left), pairing_bit(matrix, left, right)),
        (pairing_bit(matrix, right, left), pairing_bit(matrix, right, right)),
    )


def main() -> None:
    two_orbit_permutation = (1, 2, 0, 4, 5, 3)
    matrices_checked = 0
    nonsingular_count = 0
    standard_count = 0
    isotropic_standard_count = 0
    standard_plane_count = 0
    isotropic_plane_count = 0

    for parameters in product(range(4), repeat=7):
        matrix = two_orbit_matrix(parameters)
        require(invariant_matrix(matrix, two_orbit_permutation),
                "two-orbit C3 invariance")
        matrices_checked += 1
        if determinant_bareiss(matrix) == 0:
            continue
        nonsingular_count += 1
        kernel = kernel_mod_two(matrix)
        planes = standard_planes(kernel, two_orbit_permutation)
        isotropic_planes = tuple(
            plane for plane in planes if plane_is_isotropic(matrix, plane)
        )
        standard_plane_count += len(planes)
        isotropic_plane_count += len(isotropic_planes)
        standard_count += bool(planes)
        isotropic_standard_count += bool(isotropic_planes)
        for plane in planes:
            for left in plane:
                for right in plane:
                    require(
                        pairing_bit(
                            matrix,
                            apply_permutation(left, two_orbit_permutation),
                            apply_permutation(right, two_orbit_permutation),
                        )
                        == pairing_bit(matrix, left, right),
                        "C3 invariance of mod-four pairing",
                    )

    # Same mod-two standard plane, different mod-four isotropy.
    one_orbit_permutation = (1, 2, 0)
    symplectic = [[-2 if i == j else 0 for j in range(3)] for i in range(3)]
    isotropic = [[-2 if i == j else 2 for j in range(3)] for i in range(3)]
    require(
        [[value % 2 for value in row] for row in symplectic]
        == [[value % 2 for value in row] for row in isotropic],
        "mod-two-identical hostile",
    )
    symplectic_plane = standard_planes(
        kernel_mod_two(symplectic), one_orbit_permutation
    )
    isotropic_plane = standard_planes(
        kernel_mod_two(isotropic), one_orbit_permutation
    )
    require(len(symplectic_plane) == len(isotropic_plane) == 1,
            "one standard plane in hostile pair")
    require(not plane_is_isotropic(symplectic, symplectic_plane[0]),
            "symplectic half of mod-four hostile")
    require(plane_is_isotropic(isotropic, isotropic_plane[0]),
            "isotropic half of mod-four hostile")
    require(abs(determinant_bareiss(symplectic)) == 8, "symplectic determinant")
    require(abs(determinant_bareiss(isotropic)) == 32, "isotropic determinant")

    # D4 triality: the unique W is alternating nondegenerate, so it cannot be
    # contained in a C3-stable metabolizer.
    d4 = [
        [-2, 1, 1, 1],
        [1, -2, 0, 0],
        [1, 0, -2, 0],
        [1, 0, 0, -2],
    ]
    d4_permutation = (0, 2, 3, 1)
    require(invariant_matrix(d4, d4_permutation), "D4 triality")
    d4_kernel = kernel_mod_two(d4)
    d4_planes = standard_planes(d4_kernel, d4_permutation)
    require(len(d4_planes) == 1, "D4 unique standard plane")
    require(not plane_is_isotropic(d4, d4_planes[0]),
            "D4 standard plane is symplectic")
    d4_pairing = pairing_matrix_on_plane(d4, d4_planes[0])
    require(d4_pairing == ((0, 1), (1, 0)), "D4 symplectic matrix")
    d4_fixed_nonzero = tuple(
        vector for vector in d4_kernel
        if vector and apply_permutation(vector, d4_permutation) == vector
    )
    require(not d4_fixed_nonzero, "D4 has no stable metabolizer line")
    require(abs(determinant_bareiss(d4)) == 4, "D4 determinant")
    require(smith_diagonal(d4) == (1, 1, 2, 2), "D4 Smith form")

    # Positive control: a fixed +1 direction and two C3 triples of -2 roots.
    # The diagonal subgroup across the triples is a stable metabolizer, and
    # adjoining its half-vectors gives the odd unimodular lattice I_(1,6).
    positive = sp.diag(1, *([-2] * 6))
    positive_list = [[int(positive[i, j]) for j in range(7)] for i in range(7)]
    positive_permutation = (0, 2, 3, 1, 5, 6, 4)
    require(invariant_matrix(positive_list, positive_permutation),
            "positive-control C3 invariance")
    generators = tuple((1 << (1 + i)) | (1 << (4 + i)) for i in range(3))
    metabolizer = span(generators)
    positive_kernel = kernel_mod_two(positive_list)
    require(metabolizer <= set(positive_kernel), "metabolizer lies in A[2]")
    require(all(
        pairing_bit(positive_list, left, right) == 0
        for left in metabolizer for right in metabolizer
    ), "positive metabolizer isotropic")
    require(orthogonal(positive_list, positive_kernel, metabolizer) == metabolizer,
            "positive metabolizer is self-orthogonal")
    require(
        {apply_permutation(value, positive_permutation) for value in metabolizer}
        == set(metabolizer),
        "positive metabolizer C3 stable",
    )
    positive_w = span((generators[0] ^ generators[1],
                       generators[0] ^ generators[2]))
    require(positive_w in set(standard_planes(positive_kernel, positive_permutation)),
            "positive metabolizer contains standard plane")

    transform = sp.zeros(7)
    transform[0, 0] = 1
    column = 1
    for i in range(3):
        e, f = 1 + i, 4 + i
        transform[e, column] = sp.Rational(1, 2)
        transform[f, column] = sp.Rational(1, 2)
        column += 1
        transform[e, column] = sp.Rational(1, 2)
        transform[f, column] = sp.Rational(-1, 2)
        column += 1
    overlattice_gram = sp.simplify(transform.T * positive * transform)
    require(overlattice_gram == sp.diag(1, *([-1] * 6)),
            "positive unimodular overlattice Gram")
    require(abs(transform.det()) == sp.Rational(1, 8),
            "positive overlattice index eight")
    require(abs(determinant_bareiss(positive_list)) == 64,
            "positive boundary determinant")

    print("THM-2711 C3 DISCRIMINANT-METABOLIZER MOD-FOUR AUDIT")
    print(f"two_orbit_mod4_matrices_checked={matrices_checked}")
    print(f"nonsingular_matrices={nonsingular_count}")
    print(f"matrices_with_standard_plane={standard_count}")
    print(f"matrices_with_isotropic_standard_plane={isotropic_standard_count}")
    print(f"standard_planes_total={standard_plane_count}")
    print(f"isotropic_standard_planes_total={isotropic_plane_count}")
    print("mod2_identical_pair_dets=8,32 standard_isotropy=NO,YES")
    print("d4_det=4 smith=(1,1,2,2) lambda_W=((0,1),(1,0)) stable_metabolizer=NO")
    print("positive_L_det=64 H_order=8 H_equals_Hperp=YES contains_W=YES")
    print("positive_overlattice_gram=diag(1,-1^6) index=8")
    print("order_two_pairing=lambda(x,y)=x^T*M*y/2_mod2")
    print("boundary_realization_and_reflection_not_tested")
    print("PASS")


if __name__ == "__main__":
    main()
