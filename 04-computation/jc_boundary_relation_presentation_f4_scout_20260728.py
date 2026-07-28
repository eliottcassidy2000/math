#!/usr/bin/env python3
"""Exact scout for the C3 boundary-presentation Kummer gate.

The script checks the F4 identity

    nullity(A^dagger G A)
      = nullity(A) + dim(rad(im(A)))

for every matrix A of size at most 2 by 2 and every nonsingular Hermitian G
of the corresponding target size.  It also checks four integral presentation
controls which separate boundary relations, saturation, an isotropic-Gram
false positive, and a nonsingular exclusion.  There is no floating point
arithmetic and no truth-bearing Python ``assert`` statement.
"""

from __future__ import annotations

from itertools import product

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# F4 = F2[w]/(w^2+w+1), encoded as a0 + a1*w in the two bits of an int.
def f4_add(left: int, right: int) -> int:
    return left ^ right


def f4_mul(left: int, right: int) -> int:
    left_0, left_1 = left & 1, (left >> 1) & 1
    right_0, right_1 = right & 1, (right >> 1) & 1
    constant = (left_0 * right_0) ^ (left_1 * right_1)
    omega = (left_0 * right_1) ^ (left_1 * right_0) ^ (left_1 * right_1)
    return constant | (omega << 1)


def f4_conjugate(value: int) -> int:
    return f4_mul(value, value)


def f4_inverse(value: int) -> int:
    require(value != 0, "zero has no inverse")
    for candidate in range(1, 4):
        if f4_mul(value, candidate) == 1:
            return candidate
    raise RuntimeError("missing F4 inverse")


def matrix_shape(matrix: list[list[int]]) -> tuple[int, int]:
    return len(matrix), len(matrix[0]) if matrix else 0


def matrix_dagger(matrix: list[list[int]]) -> list[list[int]]:
    rows, columns = matrix_shape(matrix)
    return [
        [f4_conjugate(matrix[row][column]) for row in range(rows)]
        for column in range(columns)
    ]


def matrix_multiply(
    left: list[list[int]], right: list[list[int]]
) -> list[list[int]]:
    left_rows, shared = matrix_shape(left)
    right_rows, right_columns = matrix_shape(right)
    require(shared == right_rows, "matrix product shape")
    answer = [[0] * right_columns for _ in range(left_rows)]
    for row in range(left_rows):
        for column in range(right_columns):
            value = 0
            for middle in range(shared):
                value = f4_add(
                    value, f4_mul(left[row][middle], right[middle][column])
                )
            answer[row][column] = value
    return answer


def matrix_vector(
    matrix: list[list[int]], vector: tuple[int, ...]
) -> tuple[int, ...]:
    rows, columns = matrix_shape(matrix)
    require(columns == len(vector), "matrix-vector shape")
    answer = []
    for row in range(rows):
        value = 0
        for column in range(columns):
            value = f4_add(value, f4_mul(matrix[row][column], vector[column]))
        answer.append(value)
    return tuple(answer)


def matrix_rank(matrix: list[list[int]]) -> int:
    rows, columns = matrix_shape(matrix)
    work = [row[:] for row in matrix]
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = f4_inverse(work[pivot_row][column])
        work[pivot_row] = [f4_mul(scale, value) for value in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or work[row][column] == 0:
                continue
            coefficient = work[row][column]
            work[row] = [
                f4_add(work[row][entry], f4_mul(coefficient, work[pivot_row][entry]))
                for entry in range(columns)
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def matrix_nullity(matrix: list[list[int]]) -> int:
    return matrix_shape(matrix)[1] - matrix_rank(matrix)


def all_vectors(dimension: int) -> tuple[tuple[int, ...], ...]:
    return tuple(product(range(4), repeat=dimension))


def hermitian_value(
    left: tuple[int, ...],
    gram: list[list[int]],
    right: tuple[int, ...],
) -> int:
    gram_right = matrix_vector(gram, right)
    value = 0
    for left_entry, right_entry in zip(left, gram_right):
        value = f4_add(value, f4_mul(f4_conjugate(left_entry), right_entry))
    return value


def image_vectors(matrix: list[list[int]]) -> frozenset[tuple[int, ...]]:
    return frozenset(
        matrix_vector(matrix, vector)
        for vector in all_vectors(matrix_shape(matrix)[1])
    )


def power_of_four_dimension(cardinality: int) -> int:
    dimension = 0
    while cardinality > 1:
        require(cardinality % 4 == 0, "F4-space cardinality")
        cardinality //= 4
        dimension += 1
    return dimension


def image_radical_dimension(
    matrix: list[list[int]], gram: list[list[int]]
) -> int:
    image = image_vectors(matrix)
    radical = tuple(
        vector
        for vector in image
        if all(hermitian_value(vector, gram, other) == 0 for other in image)
    )
    return power_of_four_dimension(len(radical))


def gram_shadow(
    matrix: list[list[int]], gram: list[list[int]]
) -> list[list[int]]:
    return matrix_multiply(matrix_multiply(matrix_dagger(matrix), gram), matrix)


def nonsingular_hermitian_matrices(size: int) -> tuple[list[list[int]], ...]:
    matrices = []
    for entries in product(range(4), repeat=size * size):
        matrix = [list(entries[row * size : (row + 1) * size]) for row in range(size)]
        if matrix != matrix_dagger(matrix):
            continue
        if matrix_rank(matrix) != size:
            continue
        matrices.append(matrix)
    return tuple(matrices)


def smith_data(matrix: list[list[int]]) -> tuple[int, int, int]:
    integer_matrix = sp.Matrix(matrix)
    rational_rank = int(integer_matrix.rank())
    normal = smith_normal_form(integer_matrix, domain=sp.ZZ)
    diagonal = [
        abs(int(normal[index, index]))
        for index in range(min(normal.rows, normal.cols))
        if normal[index, index] != 0
    ]
    relation_rank = integer_matrix.cols - rational_rank
    coker_two_rank = sum(value % 2 == 0 for value in diagonal)
    binary_matrix = [[value % 2 for value in row] for row in matrix]
    mod_two_kernel_rank = integer_matrix.cols - matrix_rank(binary_matrix)
    require(
        mod_two_kernel_rank == relation_rank + coker_two_rank,
        "universal-coefficient/Tor dimension identity",
    )
    return relation_rank, coker_two_rank, mod_two_kernel_rank


def integer_gram(
    boundary_map: list[list[int]], ambient_gram: list[list[int]]
) -> list[list[int]]:
    delta = sp.Matrix(boundary_map)
    gram = delta.T * sp.Matrix(ambient_gram) * delta
    return [[int(gram[row, column]) for column in range(gram.cols)] for row in range(gram.rows)]


def mod_four_pairing_bit(
    matrix: list[list[int]], left: tuple[int, ...], right: tuple[int, ...]
) -> int:
    numerator = sum(
        left[row] * matrix[row][column] * right[column]
        for row in range(len(matrix))
        for column in range(len(matrix))
    )
    require(numerator % 2 == 0, "mod-four pairing numerator")
    return (numerator // 2) % 2


def main() -> None:
    exhaustive_cases = 0
    hermitian_counts: dict[int, int] = {}
    for target_size in (1, 2):
        grams = nonsingular_hermitian_matrices(target_size)
        hermitian_counts[target_size] = len(grams)
        for source_size in (1, 2):
            for entries in product(range(4), repeat=target_size * source_size):
                matrix = [
                    list(entries[row * source_size : (row + 1) * source_size])
                    for row in range(target_size)
                ]
                for gram in grams:
                    shadow = gram_shadow(matrix, gram)
                    nullity_a = matrix_nullity(matrix)
                    nullity_b = matrix_nullity(shadow)
                    radical = image_radical_dimension(matrix, gram)
                    require(
                        nullity_b == nullity_a + radical,
                        "Hermitian shadow nullity identity",
                    )
                    exhaustive_cases += 1

    controls = {
        "relation": ([[1, 1]], [[1]], (1, 1, 0)),
        "saturation": ([[0]], [[1]], (1, 1, 0)),
        "isotropic_radical_hostile": ([[1], [1]], [[1, 0], [0, 1]], (0, 1, 1)),
        "nonsingular": ([[1]], [[1]], (0, 0, 0)),
    }
    control_results: dict[str, tuple[int, int, int]] = {}
    for name, (matrix, gram, expected) in controls.items():
        observed = (
            matrix_nullity(matrix),
            matrix_nullity(gram_shadow(matrix, gram)),
            image_radical_dimension(matrix, gram),
        )
        require(observed == expected, f"{name} control")
        control_results[name] = observed

    identity = [[1 if row == column else 0 for column in range(3)] for row in range(3)]
    relation_map = [identity[row] + identity[row] for row in range(3)]
    saturation_map = [[2 if row == column else 0 for column in range(3)] for row in range(3)]
    hostile_map = (
        [[0, 0, 0]]
        + identity
        + identity
    )
    nonsingular_map = identity
    presentation_results = {
        "relation": smith_data(relation_map),
        "saturation": smith_data(saturation_map),
        "isotropic_radical_hostile": smith_data(hostile_map),
        "nonsingular": smith_data(nonsingular_map),
    }
    require(presentation_results["relation"] == (3, 0, 3), "relation branch")
    require(presentation_results["saturation"] == (0, 3, 3), "saturation branch")
    require(
        presentation_results["isotropic_radical_hostile"] == (0, 0, 0),
        "hostile has no Kummer kernel",
    )
    require(presentation_results["nonsingular"] == (0, 0, 0), "nonsingular branch")

    ambient_hostile = sp.diag(1, *([-1] * 6))
    hostile_gram = integer_gram(
        hostile_map,
        [[int(ambient_hostile[row, column]) for column in range(7)] for row in range(7)],
    )
    require(
        hostile_gram == [[-2 if row == column else 0 for column in range(3)] for row in range(3)],
        "integral hostile Gram",
    )
    standard_left = (1, 1, 0)
    standard_right = (0, 1, 1)
    require(
        mod_four_pairing_bit(hostile_gram, standard_left, standard_left) == 0,
        "hostile alternating diagonal",
    )
    require(
        mod_four_pairing_bit(hostile_gram, standard_left, standard_right) == 1,
        "hostile symplectic cross-pairing",
    )

    print("JC BOUNDARY-RELATION PRESENTATION/F4 SCOUT")
    print(f"nonsingular_hermitian_counts={hermitian_counts}")
    print(f"exhaustive_A_G_cases={exhaustive_cases}")
    for name in ("relation", "saturation", "isotropic_radical_hostile", "nonsingular"):
        print(f"f4_{name}_nullA_nullB_rad={control_results[name]}")
    for name in ("relation", "saturation", "isotropic_radical_hostile", "nonsingular"):
        print(f"Z_{name}_Krank_Tor2rank_mod2kernel={presentation_results[name]}")
    print("hostile_integral_gram=diag(-2,-2,-2)")
    print("hostile_standard_mod4_pairing=((0,1),(1,0))")
    print("geometric_boundary_realization_and_S3_reflection_not_tested")
    print("PASS")


if __name__ == "__main__":
    main()
