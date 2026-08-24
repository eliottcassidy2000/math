#!/usr/bin/env python3
"""Exact controls for the S^6 octonion structure and arXiv:2608.20875.

The script has three independent jobs, all over integers or rationals:

* exhibit nonzero Nijenhuis tensor for the standard octonionic almost-complex
  structure on S^6;
* replay the exact order-seven Lyapunov separator from Kressner--Vandereycken;
* test whether its skew witness is in either summand of
  Lambda^2(R^7) = Lambda^2_7 + Lambda^2_14 under the standard G2 form.

It uses only the Python standard library.  Coordinates are zero based in the
code and one based in the printed mathematical output.
"""

from fractions import Fraction


def require(condition, message):
    """Keep every certificate gate active under ``python -O``."""
    if not condition:
        raise ArithmeticError(message)


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def multiply(left, right):
    return [
        [
            sum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def add(left, right):
    return [
        [x + y for x, y in zip(left_row, right_row)]
        for left_row, right_row in zip(left, right)
    ]


def subtract(left, right):
    return [
        [x - y for x, y in zip(left_row, right_row)]
        for left_row, right_row in zip(left, right)
    ]


def matrix_rank(matrix):
    work = [[Fraction(entry) for entry in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [entry / scale for entry in work[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [
                x - scale * y for x, y in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
    return pivot_row


def parity_basis(size, symmetric):
    basis = []
    if symmetric:
        for i in range(size):
            unit = [[0] * size for _ in range(size)]
            unit[i][i] = 1
            basis.append(unit)
    for i in range(size):
        for j in range(i + 1, size):
            unit = [[0] * size for _ in range(size)]
            unit[i][j] = 1
            unit[j][i] = 1 if symmetric else -1
            basis.append(unit)
    return basis


def lyapunov_restriction(matrix, symmetric):
    """Coordinates of X -> AX + XA^T in the parity-adapted basis."""
    size = len(matrix)
    basis = parity_basis(size, symmetric)
    restriction = [[0] * len(basis) for _ in basis]
    for column, source in enumerate(basis):
        target = add(
            multiply(matrix, source),
            multiply(source, transpose(matrix)),
        )
        coordinates = []
        if symmetric:
            coordinates.extend(target[i][i] for i in range(size))
        coordinates.extend(
            target[i][j]
            for i in range(size)
            for j in range(i + 1, size)
        )
        for row, entry in enumerate(coordinates):
            restriction[row][column] = entry
    return restriction


def weighted_gram(matrix, weights):
    return [
        [
            sum(
                matrix[k][i] * weights[k] * matrix[k][j]
                for k in range(len(matrix))
            )
            for j in range(len(matrix))
        ]
        for i in range(len(matrix))
    ]


def ldlt(matrix):
    """Exact no-pivot LDL^T; positive pivots certify positive definiteness."""
    size = len(matrix)
    lower = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    pivots = []
    for i in range(size):
        lower[i][i] = Fraction(1)
        pivot = Fraction(matrix[i][i]) - sum(
            lower[i][k] ** 2 * pivots[k] for k in range(i)
        )
        require(pivot, f"zero LDL pivot at index {i}")
        pivots.append(pivot)
        for j in range(i + 1, size):
            lower[j][i] = (
                Fraction(matrix[j][i])
                - sum(
                    lower[j][k] * lower[i][k] * pivots[k]
                    for k in range(i)
                )
            ) / pivot
    return lower, pivots


def quadratic_form(vector, matrix):
    return sum(
        vector[i] * matrix[i][j] * vector[j]
        for i in range(len(vector))
        for j in range(len(vector))
    )


# Kressner--Vandereycken's integer counterexample and upper-triangle K data.
A = [
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [6, -5, -13, 0, 0, 0, 0],
    [-14, -18, 0, 0, 0, 0, 0],
    [12, -12, 11, 0, 0, 0, 0],
    [0, 0, 0, -6, 0, -14, 0],
    [0, 0, 0, 6, -18, 0, 0],
]
K_UPPER = [
    4, -9, 41, 71, 4, 6,
    -11, -59, 105, -12, 3,
    59, -10, 36, -8,
    9, 0, 2,
    -15, 0,
    4,
]


def audit_lyapunov_certificate():
    separator = 1196
    symmetric_restriction = lyapunov_restriction(A, symmetric=True)
    skew_restriction = lyapunov_restriction(A, symmetric=False)
    symmetric_weights = [1] * 7 + [2] * 21
    skew_weights = [2] * 21
    symmetric_gram = weighted_gram(
        symmetric_restriction, symmetric_weights
    )
    separator_matrix = [
        [
            (separator * symmetric_weights[i] if i == j else 0)
            - symmetric_gram[i][j]
            for j in range(28)
        ]
        for i in range(28)
    ]
    lower, pivots = ldlt(separator_matrix)
    diagonal = [
        [pivots[i] if i == j else Fraction(0) for j in range(28)]
        for i in range(28)
    ]
    require(
        multiply(multiply(lower, diagonal), transpose(lower))
        == separator_matrix,
        "exact LDL reconstruction failed",
    )
    require(min(pivots) >= 82, "symmetric separator is not positive definite")

    skew_gram = weighted_gram(skew_restriction, skew_weights)
    denominator = sum(2 * entry * entry for entry in K_UPPER)
    numerator = quadratic_form(K_UPPER, skew_gram)
    surplus = numerator - separator * denominator
    require(
        (denominator, numerator, surplus) == (53836, 64387950, 94),
        "skew Lyapunov witness arithmetic changed",
    )

    commutator = subtract(
        multiply(A, transpose(A)), multiply(transpose(A), A)
    )
    commutator_rank = matrix_rank(commutator)
    commutator_norm_squared = sum(
        entry * entry for row in commutator for entry in row
    )
    require(
        (commutator_rank, commutator_norm_squared) == (7, 819518),
        "nonnormality hostile changed",
    )
    print(
        "LYAPUNOV_SEPARATOR=PASS"
        f"|symmetric_pivot_min={min(pivots)}"
        f"|skew_fraction={numerator}/{denominator}"
        f"|surplus_over_{separator}={surplus}"
    )
    print(
        "A_IS_NORMAL=NO"
        f"|commutator_rank={commutator_rank}"
        f"|commutator_frobenius_squared={commutator_norm_squared}"
    )


# phi=e123+e145+e167+e246-e257-e347-e356.
PHI_SORTED = {
    (0, 1, 2): 1,
    (0, 3, 4): 1,
    (0, 5, 6): 1,
    (1, 3, 5): 1,
    (1, 4, 6): -1,
    (2, 3, 6): -1,
    (2, 4, 5): -1,
}


def permutation_sign(indices):
    inversions = sum(
        indices[i] > indices[j]
        for i in range(len(indices))
        for j in range(i + 1, len(indices))
    )
    return -1 if inversions % 2 else 1


def phi(i, j, k):
    if len({i, j, k}) < 3:
        return 0
    ordered = (i, j, k)
    sorted_indices = tuple(sorted(ordered))
    permutation = [sorted_indices.index(index) for index in ordered]
    return PHI_SORTED.get(sorted_indices, 0) * permutation_sign(permutation)


def vector_add(*vectors):
    return [sum(entries) for entries in zip(*vectors)]


def vector_scale(scalar, vector):
    return [scalar * entry for entry in vector]


def dot(left, right):
    return sum(x * y for x, y in zip(left, right))


def cross(left, right):
    return [
        sum(
            left[i] * right[j] * phi(i, j, k)
            for i in range(7)
            for j in range(7)
        )
        for k in range(7)
    ]


def unit(index):
    return [int(i == index) for i in range(7)]


def audit_s6_nijenhuis():
    point = unit(6)

    def almost_complex(vector):
        return cross(point, vector)

    # At x=e7, tangent vectors are e1,...,e6.
    for index in range(6):
        tangent = unit(index)
        require(
            almost_complex(almost_complex(tangent)) == vector_scale(-1, tangent),
            f"J squared is not -I on tangent basis vector {index + 1}",
        )

    # On the unit sphere, (nabla_U J)V is the tangent projection of U x V.
    def covariant_j(U, V):
        product = cross(U, V)
        return vector_add(product, vector_scale(-dot(product, point), point))

    def nijenhuis(U, V):
        # Torsion-free connection form of N_J, with N_J=[J,J].
        return vector_add(
            covariant_j(almost_complex(U), V),
            vector_scale(-1, covariant_j(almost_complex(V), U)),
            almost_complex(covariant_j(V, U)),
            vector_scale(-1, almost_complex(covariant_j(U, V))),
        )

    witness = nijenhuis(unit(0), unit(1))
    require(witness == vector_scale(4, unit(3)), "Nijenhuis witness changed")
    print("OCTONIONIC_S6=PASS|J_squared=-I|N_J(e1,e2)=4e4")


def audit_g2_decomposition():
    pairs = [(i, j) for i in range(7) for j in range(i + 1, 7)]
    # Lambda^2_7 is {i_v phi}; these seven rows have Gram matrix 3 I.
    contraction_basis = [
        [phi(v, i, j) for i, j in pairs] for v in range(7)
    ]
    gram = [
        [dot(left, right) for right in contraction_basis]
        for left in contraction_basis
    ]
    require(
        gram == [[3 if i == j else 0 for j in range(7)] for i in range(7)],
        "contraction basis does not have Gram matrix 3I",
    )
    coefficients = [
        Fraction(dot(row, K_UPPER), 3) for row in contraction_basis
    ]
    expected = [
        Fraction(2, 3), Fraction(3), Fraction(17, 3),
        Fraction(-67, 3), Fraction(80, 3), Fraction(-55, 3),
        Fraction(-160, 3),
    ]
    require(coefficients == expected, "Lambda^2_7 projection changed")
    total_norm_squared = dot(K_UPPER, K_UPPER)
    projection_7_norm_squared = 3 * dot(coefficients, coefficients)
    projection_14_norm_squared = total_norm_squared - projection_7_norm_squared
    require(
        (
            total_norm_squared,
            projection_7_norm_squared,
            projection_14_norm_squared,
        ) == (26918, 13296, 13622),
        "G2 skew-witness norm split changed",
    )

    # Decompose the paper's A itself as End(R7)=1+27+7+14.  Matrix Frobenius
    # norm gives a factor two on the upper-triangle skew coordinates.
    symmetric_part = [
        [Fraction(A[i][j] + A[j][i], 2) for j in range(7)]
        for i in range(7)
    ]
    skew_upper = [
        Fraction(A[i][j] - A[j][i], 2) for i, j in pairs
    ]
    trace_average = Fraction(sum(A[i][i] for i in range(7)), 7)
    scalar_norm_squared = 7 * trace_average ** 2
    symmetric_norm_squared = sum(
        entry * entry for row in symmetric_part for entry in row
    )
    traceless_norm_squared = symmetric_norm_squared - scalar_norm_squared
    skew_coefficients = [
        Fraction(dot(row, skew_upper), 3) for row in contraction_basis
    ]
    skew_7_norm_squared = 6 * dot(skew_coefficients, skew_coefficients)
    skew_norm_squared = 2 * dot(skew_upper, skew_upper)
    skew_14_norm_squared = skew_norm_squared - skew_7_norm_squared
    require(
        (
            scalar_norm_squared,
            traceless_norm_squared,
            skew_7_norm_squared,
            skew_14_norm_squared,
        ) == (Fraction(729, 7), Fraction(6677, 7), 237, 456),
        "G2 channel split of A changed",
    )
    print(
        "G2_TWO_FORM_SPLIT=PASS"
        f"|upper_norm_squared={total_norm_squared}"
        f"|projection_7={projection_7_norm_squared}"
        f"|projection_14={projection_14_norm_squared}"
        "|witness_is_pure_7=NO|witness_is_pure_14=NO"
    )
    print(
        "A_G2_CHANNELS=PASS"
        f"|one={scalar_norm_squared}|twenty_seven={traceless_norm_squared}"
        f"|seven={skew_7_norm_squared}|fourteen={skew_14_norm_squared}"
        "|all_nonzero=YES"
    )


def main():
    audit_s6_nijenhuis()
    audit_lyapunov_certificate()
    audit_g2_decomposition()


if __name__ == "__main__":
    main()
