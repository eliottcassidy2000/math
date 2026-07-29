#!/usr/bin/env python3
"""Exact companion for THM-2953.

The script realizes the rank-one Reynolds hostile on the real
augmentation representation of C_7.  All truth-bearing arithmetic uses
fractions and explicit exceptions; optimized execution retains every gate.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def identity(size: int) -> list[list[Fraction]]:
    return [
        [Fraction(row == column) for column in range(size)]
        for row in range(size)
    ]


def transpose(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    return [list(column) for column in zip(*matrix)]


def multiply(
    left: list[list[Fraction]],
    right: list[list[Fraction]],
) -> list[list[Fraction]]:
    require(
        len(left[0]) == len(right),
        "matrix multiplication shape mismatch",
    )
    right_t = transpose(right)
    return [
        [sum(a * b for a, b in zip(row, column)) for column in right_t]
        for row in left
    ]


def subtract(
    left: list[list[Fraction]],
    right: list[list[Fraction]],
) -> list[list[Fraction]]:
    require(
        len(left) == len(right) and len(left[0]) == len(right[0]),
        "matrix subtraction shape mismatch",
    )
    return [
        [a - b for a, b in zip(left_row, right_row)]
        for left_row, right_row in zip(left, right)
    ]


def matrix_rank(matrix: list[list[Fraction]]) -> int:
    work = [row[:] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [entry / scale for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [
                entry - scale * pivot_entry
                for entry, pivot_entry in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    size = len(matrix)
    require(all(len(row) == size for row in matrix), "determinant not square")
    work = [row[:] for row in matrix]
    answer = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        pivot_value = work[column][column]
        answer *= pivot_value
        for row in range(column + 1, size):
            if not work[row][column]:
                continue
            scale = work[row][column] / pivot_value
            for entry in range(column + 1, size):
                work[row][entry] -= scale * work[column][entry]
    return answer


def deleted_minor(
    matrix: list[list[Fraction]],
    deleted_row: int,
    deleted_column: int,
) -> Fraction:
    return determinant(
        [
            [
                entry
                for column, entry in enumerate(row)
                if column != deleted_column
            ]
            for row_index, row in enumerate(matrix)
            if row_index != deleted_row
        ]
    )


def shift_matrix(size: int) -> list[list[Fraction]]:
    matrix = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for column in range(size):
        matrix[(column + 1) % size][column] = Fraction(1)
    return matrix


def augmentation_basis(size: int) -> list[list[Fraction]]:
    basis = [[Fraction(0) for _ in range(size - 1)] for _ in range(size)]
    for column in range(size - 1):
        basis[column][column] = Fraction(1)
        basis[size - 1][column] = Fraction(-1)
    return basis


def augmentation_coordinates(
    operator: list[list[Fraction]],
) -> list[list[Fraction]]:
    size = len(operator)
    basis = augmentation_basis(size)
    image = multiply(operator, basis)
    coordinates = [row[:] for row in image[: size - 1]]
    require(
        multiply(basis, coordinates) == image,
        "operator did not preserve augmentation",
    )
    return coordinates


def reynolds_laplacian(size: int) -> list[list[Fraction]]:
    shift = shift_matrix(size)
    shift_inverse = transpose(shift)
    answer = identity(size)
    for row in range(size):
        for column in range(size):
            answer[row][column] = (
                2 * answer[row][column]
                - shift[row][column]
                - shift_inverse[row][column]
            ) / size
    return answer


def fraction_text(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def main() -> None:
    prime = 7
    ambient_shift = shift_matrix(prime)
    averaged = reynolds_laplacian(prime)
    charged = augmentation_coordinates(averaged)
    charged_shift = augmentation_coordinates(ambient_shift)

    require(
        multiply(charged, charged_shift)
        == multiply(charged_shift, charged),
        "Reynolds average lost cyclic equivariance",
    )
    require(matrix_rank(charged) == 6, "charged average lost full rank")
    charged_determinant = determinant(charged)
    require(
        charged_determinant == Fraction(1, 2401),
        "charged determinant changed",
    )

    cofactors = [
        deleted_minor(charged, row, column)
        for row in range(6)
        for column in range(6)
    ]
    nonzero_cofactors = sum(value != 0 for value in cofactors)
    fifth_energy = sum(value * value for value in cofactors)
    canonical_cofactor = deleted_minor(charged, 0, 0)
    require(nonzero_cofactors == 30, "fifth-minor census changed")
    require(
        canonical_cofactor == Fraction(3, 2401),
        "canonical fifth minor changed",
    )
    require(
        fifth_energy == Fraction(36, 823543),
        "fifth-minor squared energy changed",
    )

    # x=e0+e1-2e2 lies in ker(<e0-e1,.>), but its cyclic shift does not.
    kernel_vector = [Fraction(1), Fraction(1), Fraction(-2)] + [
        Fraction(0)
    ] * 4
    shifted_vector = [
        kernel_vector[(index - 1) % prime] for index in range(prime)
    ]
    response = lambda vector: vector[0] - vector[1]
    require(sum(kernel_vector) == 0, "hostile vector left augmentation")
    require(response(kernel_vector) == 0, "hostile vector left the kernel")
    require(
        response(shifted_vector) != 0,
        "rank-one kernel accidentally became shift-stable",
    )

    general_controls = {}
    for odd_prime in (3, 5, 7, 11):
        control = augmentation_coordinates(reynolds_laplacian(odd_prime))
        value = determinant(control)
        expected = Fraction(odd_prime ** 3, odd_prime**odd_prime)
        require(value == expected, "general prime determinant law failed")
        general_controls[odd_prime] = value

    print("THM2953 CYCLIC CHARGED CORANK / REYNOLDS BOUNDARY")
    print("p=7 augmentation_dim=6 real_character_pairs=3")
    print("stable_kernel_nullity=EVEN rank_gate=one_nonzero_5minor_implies_rank6")
    print("Reynolds_kernel=intersection_of_seven_rotated_kernels")
    print("hostile_original_rank=1 kernel_shift_stable=NO")
    print(
        "averaged_rank=6 determinant="
        f"{fraction_text(charged_determinant)} equivariant=YES"
    )
    print(
        f"fifth_minors_nonzero={nonzero_cofactors}/36 "
        f"canonical={fraction_text(canonical_cofactor)} "
        f"energy={fraction_text(fifth_energy)}"
    )
    print(
        "general_det="
        + ",".join(
            f"p{odd_prime}:{fraction_text(value)}"
            for odd_prime, value in general_controls.items()
        )
    )
    print("scope=joint_rotated_injectivity_not_owner_or_carrier_survival")
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
