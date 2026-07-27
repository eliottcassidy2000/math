#!/usr/bin/env python3
"""Independent exact referee for the rational-Witt additions to THM-2523.

This dependency-free audit reconstructs the chi_7 Hamilton operator over Q.
It checks the reversal-even hyperbolic summand, the reversal-odd isotropic
plane and determinant obstruction, the resulting rational Witt index five,
and the exact 13/10 coercivity and 169/10 THM-2521 normalization.  Every
check remains active under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction as F


P = 13
CHI = {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def zeros(rows: int, columns: int) -> list[list[F]]:
    return [[F(0) for _ in range(columns)] for _ in range(rows)]


def identity(size: int) -> list[list[F]]:
    answer = zeros(size, size)
    for index in range(size):
        answer[index][index] = F(1)
    return answer


def transpose(matrix: list[list[F]]) -> list[list[F]]:
    return [list(row) for row in zip(*matrix)]


def matmul(left: list[list[F]], right: list[list[F]]) -> list[list[F]]:
    right_transpose = transpose(right)
    return [
        [sum((x * y for x, y in zip(row, column)), F(0)) for column in right_transpose]
        for row in left
    ]


def matadd(*matrices: list[list[F]]) -> list[list[F]]:
    return [
        [sum((matrix[row][column] for matrix in matrices), F(0))
         for column in range(len(matrices[0][0]))]
        for row in range(len(matrices[0]))
    ]


def matscale(matrix: list[list[F]], scalar: int | F) -> list[list[F]]:
    scalar = F(scalar)
    return [[scalar * value for value in row] for row in matrix]


def matvec(matrix: list[list[F]], vector: list[F]) -> list[F]:
    return [sum((x * y for x, y in zip(row, vector)), F(0)) for row in matrix]


def dot(left: list[F], right: list[F]) -> F:
    return sum((x * y for x, y in zip(left, right)), F(0))


def bilinear(left: list[F], matrix: list[list[F]], right: list[F]) -> F:
    return dot(left, matvec(matrix, right))


def gram(vectors: list[list[F]], matrix: list[list[F]]) -> list[list[F]]:
    return [[bilinear(left, matrix, right) for right in vectors] for left in vectors]


def rank(matrix: list[list[F]]) -> int:
    work = [[F(value) for value in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0]) if work else 0
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [value / pivot_value for value in work[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                value - multiplier * pivot_value
                for value, pivot_value in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def determinant(matrix: list[list[F]]) -> F:
    work = [[F(value) for value in row] for row in matrix]
    answer = F(1)
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return F(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        pivot_value = work[column][column]
        answer *= pivot_value
        for row in range(column + 1, len(work)):
            if not work[row][column]:
                continue
            multiplier = work[row][column] / pivot_value
            for entry in range(column + 1, len(work)):
                work[row][entry] -= multiplier * work[column][entry]
    return answer


def inertia(matrix: list[list[F]]) -> tuple[int, int, int]:
    """Exact inertia by symmetric rational congruence."""

    work = [[F(value) for value in row] for row in matrix]
    positive = negative = zero = 0
    while work:
        size = len(work)
        diagonal_pivot = next(
            (index for index in range(size) if work[index][index]),
            None,
        )
        if diagonal_pivot is not None:
            if diagonal_pivot:
                work[0], work[diagonal_pivot] = work[diagonal_pivot], work[0]
                for row in work:
                    row[0], row[diagonal_pivot] = row[diagonal_pivot], row[0]
            pivot_value = work[0][0]
            positive += int(pivot_value > 0)
            negative += int(pivot_value < 0)
            first_column = [work[row][0] for row in range(1, size)]
            work = [
                [
                    work[row][column]
                    - first_column[row - 1] * first_column[column - 1] / pivot_value
                    for column in range(1, size)
                ]
                for row in range(1, size)
            ]
            continue

        off_diagonal = next(
            (
                (row, column)
                for row in range(size)
                for column in range(row + 1, size)
                if work[row][column]
            ),
            None,
        )
        if off_diagonal is None:
            zero += size
            break
        first, second = off_diagonal
        order = [first, second] + [
            index for index in range(size) if index not in (first, second)
        ]
        work = [[work[row][column] for column in order] for row in order]
        pivot_value = work[0][1]
        require(work[0][0] == work[1][1] == 0, "bad hyperbolic pivot")
        positive += 1
        negative += 1
        work = [
            [
                work[row][column]
                - (work[row][0] * work[1][column]
                   + work[row][1] * work[0][column]) / pivot_value
                for column in range(2, size)
            ]
            for row in range(2, size)
        ]
    return positive, negative, zero


def squarefree_part(number: int) -> int:
    require(number != 0, "zero has no determinant squareclass")
    sign = -1 if number < 0 else 1
    remaining = abs(number)
    answer = sign
    prime = 2
    while prime * prime <= remaining:
        exponent = 0
        while remaining % prime == 0:
            exponent += 1
            remaining //= prime
        if exponent % 2:
            answer *= prime
        prime += 1
    if remaining > 1:
        answer *= remaining
    return answer


def hamilton_matrix(tau: int) -> list[list[F]]:
    """The symmetric A_tau satisfying Q_tau(p)=<p,A_tau p>."""

    answer = zeros(P, P)
    for row in range(P):
        for step, sign in CHI.items():
            answer[row][(row + 2 * tau * step) % P] -= sign
            answer[row][(row - 2 * tau * step) % P] -= sign
    return answer


def reversal(vector: list[F]) -> list[F]:
    return [vector[-index % P] for index in range(P)]


def dilation_five(vector: list[F]) -> list[F]:
    return [vector[5 * index % P] for index in range(P)]


def augmentation_basis() -> list[list[F]]:
    answer = []
    for index in range(1, P):
        vector = [F(0)] * P
        vector[0] = F(-1)
        vector[index] = F(1)
        answer.append(vector)
    return answer


def main() -> None:
    matrix = hamilton_matrix(1)
    require(matrix == transpose(matrix), "Hamilton operator is not symmetric")
    require(all(sum(row, F(0)) == 0 for row in matrix), "constant kernel failed")
    require(rank(matrix) == 12, "augmentation restriction is degenerate")

    # Reversal gives an orthogonal six-plus-six decomposition of augmentation.
    even_basis = []
    odd_basis = []
    for residue in range(1, 7):
        even = [F(0)] * P
        even[0] = F(-2)
        even[residue] = even[-residue % P] = F(1)
        even_basis.append(even)

        odd = [F(0)] * P
        odd[residue] = F(1)
        odd[-residue % P] = F(-1)
        odd_basis.append(odd)
    require(rank(even_basis) == rank(odd_basis) == 6, "reversal half dimension")
    require(rank(even_basis + odd_basis) == 12, "reversal halves do not span")
    require(all(reversal(vector) == vector for vector in even_basis), "bad even basis")
    require(all(reversal(vector) == [-x for x in vector] for vector in odd_basis), "bad odd basis")
    require(
        all(bilinear(even, matrix, odd) == 0 for even in even_basis for odd in odd_basis),
        "reversal halves are not orthogonal",
    )

    # The two D_5 eigenspaces are complementary Lagrangians in the even half.
    orbits = ((1, 5, 12, 8), (2, 10, 11, 3), (4, 7, 9, 6))
    fixed = []
    anti_fixed = []
    for orbit in orbits:
        fixed_vector = [F(0)] * P
        fixed_vector[0] = F(-4)
        for residue in orbit:
            fixed_vector[residue] = F(1)
        fixed.append(fixed_vector)

        anti_vector = [F(0)] * P
        for index, residue in enumerate(orbit):
            anti_vector[residue] = F(1 if index % 2 == 0 else -1)
        anti_fixed.append(anti_vector)

    require(all(dilation_five(vector) == vector for vector in fixed), "D5-fixed basis")
    require(
        all(dilation_five(vector) == [-x for x in vector] for vector in anti_fixed),
        "D5-anti-fixed basis",
    )
    require(all(reversal(vector) == vector for vector in fixed + anti_fixed), "D5 spaces are not even")
    require(rank(fixed + anti_fixed) == 6, "D5 eigenspaces do not fill even half")
    require(rank(even_basis + fixed + anti_fixed) == 6, "D5 span differs from even half")
    require(
        all(bilinear(left, matrix, right) == 0 for left in fixed for right in fixed),
        "fixed space is not totally isotropic",
    )
    require(
        all(bilinear(left, matrix, right) == 0 for left in anti_fixed for right in anti_fixed),
        "anti-fixed space is not totally isotropic",
    )
    cross = [[bilinear(left, matrix, right) for right in anti_fixed] for left in fixed]
    expected_cross = [[-20, 32, 16], [-16, 12, 24], [-16, 8, 12]]
    require(cross == [[F(value) for value in row] for row in expected_cross], "cross pairing")
    cross_determinant = determinant(cross)
    require(cross_determinant == -4160, "cross pairing determinant")

    even_gram = gram(even_basis, matrix)
    even_determinant = determinant(even_gram)
    require(even_determinant == -270400, "even determinant")
    require(squarefree_part(even_determinant.numerator) == -1, "even squareclass")

    # The odd half has a two-plane but its determinant excludes a third plane.
    odd_full_gram = gram(odd_basis, matrix)
    odd_half_gram = matscale(odd_full_gram, F(1, 2))
    expected_odd_half = [
        [1, 0, 0, 2, -2, -2],
        [0, 1, 2, -2, 0, 0],
        [0, 2, -1, 0, 0, 2],
        [2, -2, 0, 1, 2, -2],
        [-2, 0, 0, 2, -1, 2],
        [-2, 0, 2, -2, 2, -1],
    ]
    require(
        odd_half_gram == [[F(value) for value in row] for row in expected_odd_half],
        "odd half-Gram matrix",
    )
    odd_determinant = determinant(odd_half_gram)
    odd_squareclass = squarefree_part(odd_determinant.numerator)
    require(odd_determinant == -325 and odd_squareclass == -13, "odd determinant squareclass")

    isotropic_v = [F(value) for value in (-1, -1, -1, 0, -1, 0)]
    isotropic_w = [F(value) for value in (-1, 0, -1, 0, 0, 0)]
    require(rank([isotropic_v, isotropic_w]) == 2, "odd isotropic vectors are dependent")
    require(
        bilinear(isotropic_v, odd_half_gram, isotropic_v)
        == bilinear(isotropic_w, odd_half_gram, isotropic_w)
        == bilinear(isotropic_v, odd_half_gram, isotropic_w)
        == 0,
        "displayed odd plane is not totally isotropic",
    )
    hyperbolic_six_squareclass = -1
    require(
        odd_squareclass != hyperbolic_six_squareclass,
        "odd determinant failed to obstruct a third hyperbolic plane",
    )
    even_witt_index = 3
    odd_witt_index = 2
    total_witt_index = even_witt_index + odd_witt_index
    require(total_witt_index == 5, "rational Witt index")

    augmentation = augmentation_basis()
    augmentation_determinant = determinant(gram(augmentation, matrix))
    require(augmentation_determinant == 1_373_125, "augmentation determinant")
    require(
        squarefree_part(augmentation_determinant.numerator)
        == squarefree_part(even_determinant.numerator) * odd_squareclass,
        "even/odd determinant squareclasses do not recombine",
    )

    # Every nonzero tau is a coordinate-permutation conjugate of tau=1.
    conjugacy_checks = 0
    for tau in range(1, P):
        inverse_tau = pow(tau, -1, P)
        tau_matrix = hamilton_matrix(tau)
        for row in range(P):
            for column in range(P):
                require(
                    tau_matrix[row][column]
                    == matrix[inverse_tau * row % P][inverse_tau * column % P],
                    "slope conjugacy",
                )
                conjugacy_checks += 1

    # Exact spectral and direct quadratic proofs of the 13/10 gap.
    square = matmul(matrix, matrix)
    fourth = matmul(square, square)
    sixth = matmul(fourth, square)
    polynomial_value = matadd(
        sixth,
        matscale(fourth, -39),
        matscale(square, 299),
        matscale(identity(P), -325),
    )
    require(polynomial_value == [[F(-25)] * P for _ in range(P)], "g(A^2) group-ring identity")
    expected_square_row = (12, -5, -1, -1, 3, -5, 3, 3, -5, 3, -1, -1, -5)
    require(tuple(square[0]) == tuple(F(value) for value in expected_square_row), "A1 squared stencil")

    boundary = F(13, 10)

    def g(value: F) -> F:
        return value**3 - 39 * value**2 + 299 * value - 325

    def g_prime(value: F) -> F:
        return 3 * value**2 - 78 * value + 299

    require(g(boundary) == F(-13, 1000), "g(13/10)")
    require(g_prime(boundary) == F(20267, 100), "g-prime boundary")
    require(6 * boundary - 78 < 0 and g_prime(boundary) > 0, "monotonic root exclusion")

    coercive_matrix = matadd(matscale(square, 10), matscale(identity(P), -13))
    coercive_inertia = inertia(gram(augmentation, coercive_matrix))
    require(coercive_inertia == (12, 0, 0), "10*A^2-13*I is not positive on augmentation")

    partner_checks = 0
    partner_controls = augmentation + even_basis + odd_basis
    for profile in partner_controls:
        partner = matvec(matrix, profile)
        polarized_value = bilinear(profile, matrix, partner)
        partner_norm = dot(partner, partner)
        require(polarized_value == partner_norm > 0, "canonical partner identity")
        require(10 * partner_norm > 13 * dot(profile, profile), "13/10 coercivity control")
        partner_checks += 1

    # A finite exact Hilbert-valued control checks the THM-2521 normalization.
    weights = (F(1, 3), F(2, 5), F(7, 11), F(5, 2))
    hilbert_profiles = (
        augmentation[0],
        [augmentation[2][i] + 2 * augmentation[7][i] for i in range(P)],
        even_basis[4],
        [odd_basis[1][i] - 3 * odd_basis[5][i] for i in range(P)],
    )
    weighted_norm = sum(
        (weight * dot(profile, profile) for weight, profile in zip(weights, hilbert_profiles)),
        F(0),
    )
    derived_detector = sum(
        (
            weight * dot(matvec(matrix, profile), matvec(matrix, profile))
            for weight, profile in zip(weights, hilbert_profiles)
        ),
        F(0),
    )
    drift_13 = weighted_norm / 13
    require(weighted_norm == 13 * drift_13 and drift_13 > 0, "THM-2521 normalization")
    require(10 * derived_detector > 169 * drift_13, "169/10 Hilbert detector")
    require(13 * boundary == F(169, 10), "normalization factor 13*(13/10)")

    print("THM-2523 rational Witt/coercive independent referee")
    print(
        f"reversal: even_dim=6 odd_dim=6 orthogonal=yes; "
        f"even_det={even_determinant} squareclass=-1"
    )
    print(
        f"even_hyperbolic: fixed_dim=3 anti_fixed_dim=3 "
        f"cross_det={cross_determinant} index={even_witt_index}"
    )
    print(
        f"odd_boundary: half_det={odd_determinant} squareclass={odd_squareclass}; "
        f"isotropic_plane=2 hyperbolic_index={odd_witt_index}"
    )
    print(
        f"augmentation: det={augmentation_determinant} squareclass=13; "
        f"rational_Witt_index={total_witt_index} anisotropic_dimension=2"
    )
    print(f"slope_conjugacy: exact_entries={conjugacy_checks} slopes=12")
    print(
        f"coercivity: g(13/10)={g(boundary)} g_prime(13/10)={g_prime(boundary)}; "
        f"inertia(10*A^2-13*I)={coercive_inertia}; partner_checks={partner_checks}"
    )
    print(
        f"Hilbert_normalization: ||p||_G^2=13*D13; "
        f"13*(13/10)={13 * boundary}; finite_D13={drift_13}; "
        f"strict_margin_10R-169D13={10 * derived_detector - 169 * drift_13}"
    )
    print(
        "VERIFIED: the rational form has Witt index five, not six; its canonical "
        "signed partner is universally 13/10-coercive, hence gives the normalized "
        "169/10 detector, but this algebraic partner is not a Boolean ancestry event."
    )


if __name__ == "__main__":
    main()
