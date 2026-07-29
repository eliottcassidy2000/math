#!/usr/bin/env python3
"""Exact companion for THM-2951.

Only integer and Fraction arithmetic is truth-bearing.  In particular, every
gate uses ``require`` so optimized mode checks the same statements.
"""

from fractions import Fraction
from itertools import combinations, permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def parity_of_permutation(values):
    inversions = sum(
        values[i] > values[j]
        for i in range(len(values))
        for j in range(i + 1, len(values))
    )
    return -1 if inversions % 2 else 1


def matrix_rank(rows):
    matrix = [[Fraction(value) for value in row] for row in rows]
    if not matrix:
        return 0
    row_count = len(matrix)
    column_count = len(matrix[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, row_count)
                if matrix[row][column] != 0
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [entry / scale for entry in matrix[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or matrix[row][column] == 0:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                matrix[row][index] - scale * matrix[pivot_row][index]
                for index in range(column_count)
            ]
        pivot_row += 1
    return pivot_row


def determinant(matrix):
    size = len(matrix)
    require(all(len(row) == size for row in matrix), "determinant not square")
    work = [[Fraction(value) for value in row] for row in matrix]
    value = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            value = -value
        pivot_value = work[column][column]
        value *= pivot_value
        for row in range(column + 1, size):
            if work[row][column] == 0:
                continue
            scale = work[row][column] / pivot_value
            for index in range(column, size):
                work[row][index] -= scale * work[column][index]
    return value


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def matrix_inverse(matrix):
    size = len(matrix)
    require(all(len(row) == size for row in matrix), "inverse not square")
    work = []
    for row, row_values in enumerate(matrix):
        work.append(
            [Fraction(value) for value in row_values]
            + [Fraction(int(row == column)) for column in range(size)]
        )
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        require(pivot is not None, "matrix unexpectedly singular")
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [entry / scale for entry in work[column]]
        for row in range(size):
            if row == column or work[row][column] == 0:
                continue
            scale = work[row][column]
            work[row] = [
                work[row][index] - scale * work[column][index]
                for index in range(2 * size)
            ]
    return [row[size:] for row in work]


def scalar_matrix_product(scalar, matrix):
    return [[Fraction(scalar) * entry for entry in row] for row in matrix]


def matrix_product(left, right):
    require(
        left and right and len(left[0]) == len(right),
        "matrix product dimension mismatch",
    )
    return [
        [
            sum(
                Fraction(left[row][middle]) * Fraction(right[middle][column])
                for middle in range(len(right))
            )
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def compound_matrix(matrix, degree):
    row_sets = list(combinations(range(len(matrix)), degree))
    column_sets = list(combinations(range(len(matrix[0])), degree))
    return [
        [
            determinant(
                [[matrix[row][column] for column in columns] for row in rows]
            )
            for columns in column_sets
        ]
        for rows in row_sets
    ]


def contraction_sign(form_pair, five_tuple):
    first, second = form_pair
    if first not in five_tuple or second not in five_tuple:
        return 0, None
    positions = [
        five_tuple.index(first),
        five_tuple.index(second),
    ] + [
        index
        for index, value in enumerate(five_tuple)
        if value not in form_pair
    ]
    output = tuple(value for value in five_tuple if value not in form_pair)
    return parity_of_permutation(positions), output


def gaussian_mul(left, right):
    a, b = left
    c, d = right
    return a * c - b * d, a * d + b * c


def gaussian_conjugate(value):
    return value[0], -value[1]


def gaussian_div(left, right):
    denominator = right[0] * right[0] + right[1] * right[1]
    require(denominator != 0, "Gaussian division by zero")
    numerator = gaussian_mul(left, gaussian_conjugate(right))
    return numerator[0] / denominator, numerator[1] / denominator


def gaussian_product(values):
    result = (Fraction(1), Fraction(0))
    for value in values:
        result = gaussian_mul(result, value)
    return result


def half_trace_vector(values):
    traces = []
    for second, third in product((0, 1), repeat=2):
        choices = (
            values[0],
            values[1] if second == 0 else gaussian_conjugate(values[1]),
            values[2] if third == 0 else gaussian_conjugate(values[2]),
        )
        traces.append(2 * gaussian_product(choices)[0])
    return tuple(traces)


def format_integral_tuple(values):
    require(
        all(Fraction(value).denominator == 1 for value in values),
        "nonintegral tuple in transcript",
    )
    return "(" + ",".join(str(Fraction(value).numerator) for value in values) + ")"


# Signed-pair character calculation.
pair_permutations = list(permutations(range(3)))
character_numerator = 0
signed_pair_order = 0
for flips in product((0, 1), repeat=3):
    for pair_permutation in pair_permutations:
        signed_pair_order += 1
        fixed_v = sum(
            pair_permutation[pair] == pair
            and (entry ^ flips[pair]) == entry
            for pair in range(3)
            for entry in (0, 1)
        )
        determinant_character = (-1) ** sum(flips)
        character_fifth = determinant_character * fixed_v

        fixed_balanced = 0
        for epsilon in product((0, 1), repeat=3):
            image = [None, None, None]
            for pair in range(3):
                image[pair_permutation[pair]] = epsilon[pair] ^ flips[pair]
            fixed_balanced += tuple(image) == epsilon
        character_balanced = (
            parity_of_permutation(pair_permutation) * fixed_balanced
        )
        character_numerator += character_fifth * character_balanced

require(signed_pair_order == 48, "wrong signed-pair group order")
require(character_numerator == 0, "signed-pair Hom character not zero")


# Pure and projected two-form contractions.
fifth_basis = list(combinations(range(6), 5))
third_basis = list(combinations(range(6), 3))
two_form_basis = list(combinations(range(6), 2))
balanced_basis = [
    basis
    for basis in third_basis
    if sorted(coordinate // 2 for coordinate in basis) == [0, 1, 2]
]
unbalanced_basis = [basis for basis in third_basis if basis not in balanced_basis]

unbalanced_constraints = []
for fifth in fifth_basis:
    for third in unbalanced_basis:
        row = []
        for form_pair in two_form_basis:
            sign, output = contraction_sign(form_pair, fifth)
            row.append(sign if output == third else 0)
        unbalanced_constraints.append(row)

projected_columns = []
for form_pair in two_form_basis:
    column = []
    for fifth in fifth_basis:
        sign, output = contraction_sign(form_pair, fifth)
        for balanced in balanced_basis:
            column.append(sign if output == balanced else 0)
    projected_columns.append(column)
projected_matrix = [list(row) for row in zip(*projected_columns)]

unbalanced_rank = matrix_rank(unbalanced_constraints)
projected_rank = matrix_rank(projected_matrix)
within_pair_forms = [
    form_pair
    for form_pair in two_form_basis
    if form_pair[0] // 2 == form_pair[1] // 2
]
require(len(fifth_basis) == 6, "wrong fifth basis dimension")
require(len(third_basis) == 20, "wrong third basis dimension")
require(len(balanced_basis) == 8, "wrong balanced basis dimension")
require(len(two_form_basis) == 15, "wrong two-form dimension")
require(unbalanced_rank == 15, "pure-contraction kernel is nonzero")
require(projected_rank == 12, "wrong projected-contraction rank")
require(len(within_pair_forms) == 3, "wrong within-pair kernel count")
for form_pair in within_pair_forms:
    column = projected_columns[two_form_basis.index(form_pair)]
    require(not any(column), "within-pair form survives balanced projection")


# Full fifth-compound reconstruction and scalar hostile.
identity = [
    [Fraction(int(row == column)) for column in range(6)]
    for row in range(6)
]
rotation = [
    [0, -1, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1],
]
reconstruction_controls = 0
compound_controls = []
raw_fifth_basis = list(combinations(range(6), 5))
kappa = [[Fraction(0) for _ in range(6)] for _ in range(6)]
for wedge_index, wedge_basis in enumerate(raw_fifth_basis):
    missing = next(index for index in range(6) if index not in wedge_basis)
    kappa[missing][wedge_index] = Fraction((-1) ** missing)
kappa_inverse = matrix_inverse(kappa)
for matrix in (identity, rotation):
    norm = determinant(matrix)
    raw_fifth = compound_matrix(matrix, 5)
    fifth = matrix_product(matrix_product(kappa, raw_fifth), kappa_inverse)
    require(determinant(fifth) == norm**5, "det(C5)=N^5 failed")
    reconstructed = transpose(scalar_matrix_product(norm, matrix_inverse(fifth)))
    require(reconstructed == matrix, "fifth-compound reconstruction failed")
    compound_controls.append(fifth)
    reconstruction_controls += 1

fifth_identity, fifth_rotation = compound_controls
energy_identity = sum(entry * entry for row in fifth_identity for entry in row)
energy_rotation = sum(entry * entry for row in fifth_rotation for entry in row)
require(energy_identity == energy_rotation == 6, "fifth energy hostile failed")
require(
    fifth_identity[2][2] == fifth_rotation[2][2] == 1,
    "selected fifth-cofactor hostile failed",
)

trace_identity = half_trace_vector(
    ((Fraction(1), Fraction(0)),) * 3
)
trace_rotation = half_trace_vector(
    (
        (Fraction(0), Fraction(1)),
        (Fraction(1), Fraction(0)),
        (Fraction(1), Fraction(0)),
    )
)
require(trace_identity == (2, 2, 2, 2), "identity trace frame failed")
require(trace_rotation == (0, 0, 0, 0), "rotation trace frame failed")


# Eight exact controls of h=N^3/product(c_i).
gaussian_values = (
    (Fraction(1), Fraction(1)),
    (Fraction(2), Fraction(1)),
    (Fraction(3), Fraction(2)),
)
pair_norms = [
    value[0] * value[0] + value[1] * value[1]
    for value in gaussian_values
]
total_norm = pair_norms[0] * pair_norms[1] * pair_norms[2]
fifth_weights = tuple(
    gaussian_div((total_norm, Fraction(0)), value)
    for value in gaussian_values
)
balanced_weight_controls = 0
for epsilon in product((0, 1), repeat=3):
    chosen_values = tuple(
        gaussian_values[index]
        if epsilon[index] == 0
        else gaussian_conjugate(gaussian_values[index])
        for index in range(3)
    )
    chosen_fifth = tuple(
        fifth_weights[index]
        if epsilon[index] == 0
        else gaussian_conjugate(fifth_weights[index])
        for index in range(3)
    )
    left = gaussian_product(chosen_values)
    right = gaussian_div(
        (Fraction(total_norm**3), Fraction(0)),
        gaussian_product(chosen_fifth),
    )
    require(left == right, "balanced fifth-weight reconstruction failed")
    balanced_weight_controls += 1


# Singular boundaries.
rank_four = [
    [int(row == column and row >= 2) for column in range(6)]
    for row in range(6)
]
rank_five = [
    [int(row == column and row >= 1) for column in range(6)]
    for row in range(6)
]
rank_four_compound = compound_matrix(rank_four, 5)
rank_five_compound = compound_matrix(rank_five, 5)
require(matrix_rank(rank_four_compound) == 0, "rank-four C5 is nonzero")
require(matrix_rank(rank_five_compound) == 1, "rank-five C5 rank is not one")


print("THM2951 FIFTH-COMPOUND / V4 PHASE AUDIT")
print(
    "basis_dims fifth={} third={} balanced={} two_forms={}".format(
        len(fifth_basis),
        len(third_basis),
        len(balanced_basis),
        len(two_form_basis),
    )
)
print(
    "W_order={} character_inner_numerator={} Hom_W=0".format(
        signed_pair_order,
        character_numerator,
    )
)
print(
    "contraction_unbalanced_rank={} pure_kernel={}".format(
        unbalanced_rank,
        len(two_form_basis) - unbalanced_rank,
    )
)
print(
    "projected_contraction_rank={} kernel={}".format(
        projected_rank,
        len(two_form_basis) - projected_rank,
    )
)
print(
    "reconstruction_controls={} detC_equals_N5=PASS".format(
        reconstruction_controls
    )
)
print(
    "balanced_weight_controls={} spectral_formula=PASS".format(
        balanced_weight_controls
    )
)
print(
    "scalar_hostile N=1 E5={} selected_cofactor={}".format(
        energy_identity,
        fifth_identity[2][2],
    )
)
print("trace_identity={}".format(format_integral_tuple(trace_identity)))
print("trace_rotation={}".format(format_integral_tuple(trace_rotation)))
print(
    "singular_C5_ranks rank4={} rank5={}".format(
        matrix_rank(rank_four_compound),
        matrix_rank(rank_five_compound),
    )
)
print("ALL CHECKS PASS")
