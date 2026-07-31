#!/usr/bin/env python3
"""Independent exact hostile audit of the THM-2775 generator frame.

This script derives the D3 matrices from the four root variables instead of
copying the candidate matrices.  It also records the diagonal-translation
kernel of the half-Hadamard coordinates and independently checks the local
110 inertia identification from the THM-2769 polynomial.
"""

from itertools import permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def transpose(matrix):
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0])))


def matmul(left, right):
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(len(right)))
              for j in range(len(right[0])))
        for i in range(len(left))
    )


def matscale_exact(matrix, denominator):
    require(all(entry % denominator == 0
                for row in matrix for entry in row),
            "derived root-action matrix stopped being integral")
    return tuple(tuple(entry // denominator for entry in row)
                 for row in matrix)


def identity(size):
    return tuple(tuple(int(i == j) for j in range(size))
                 for i in range(size))


def matpow(matrix, exponent):
    out = identity(len(matrix))
    for _ in range(exponent):
        out = matmul(matrix, out)
    return out


def determinant_three(matrix):
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def permutation_matrix(substitution):
    """Matrix P with new root coordinate r'_i=r_{substitution[i]}."""

    return tuple(
        tuple(int(column == substitution[row]) for column in range(4))
        for row in range(4)
    )


def derive_d_action(substitution, d_matrix):
    """Derive M from D P = M D, using D D^T=4I."""

    numerator = matmul(matmul(d_matrix, permutation_matrix(substitution)),
                       transpose(d_matrix))
    return matscale_exact(numerator, 4)


def generated_group(generators):
    unit = identity(len(generators[0]))
    group = {unit}
    frontier = [unit]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            for candidate in (matmul(generator, current),
                              matmul(current, generator)):
                if candidate not in group:
                    group.add(candidate)
                    frontier.append(candidate)
    return frozenset(group)


def all_even_signed_permutations():
    matrices = set()
    for permutation in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            if signs[0] * signs[1] * signs[2] != 1:
                continue
            rows = []
            for row_index in range(3):
                row = [0, 0, 0]
                row[permutation[row_index]] = signs[row_index]
                rows.append(tuple(row))
            matrices.add(tuple(rows))
    return frozenset(matrices)


def row_times_matrix(row, matrix):
    return tuple(sum(row[k] * matrix[k][j] for k in range(len(row)))
                 for j in range(len(row)))


def induced_state_action(matrix, state_rows):
    action = []
    for row in state_rows:
        image = row_times_matrix(row, matrix)
        require(image in state_rows,
                "centered-root state left the Hadamard state set")
        action.append(state_rows.index(image))
    return tuple(action)


def permutation_order(permutation):
    unit = tuple(range(len(permutation)))
    current = unit
    for exponent in range(1, 25):
        current = tuple(permutation[current[index]]
                        for index in range(len(permutation)))
        if current == unit:
            return exponent
    raise RuntimeError("small permutation order escaped audit bound")


def main():
    # d_1=r1+r2-r3-r4, d_2=r1+r3-r2-r4,
    # d_3=r1+r4-r2-r3.
    d_matrix = (
        (1, 1, -1, -1),
        (1, -1, 1, -1),
        (1, -1, -1, 1),
    )
    require(matmul(d_matrix, transpose(d_matrix)) == (
        (4, 0, 0), (0, 4, 0), (0, 0, 4)
    ), "opposite-edge rows stopped being orthogonal")

    # Substitutions are written zero-based and mean r'_i=r_{p(i)}.
    sheet_x_s = (1, 0, 2, 3)       # (12)
    sheet_y_s = (0, 2, 3, 1)       # (234)
    sheet_x_a = (1, 0, 3, 2)       # (12)(34)
    sheet_y_a = (1, 2, 0, 3)       # (123)

    x_s = derive_d_action(sheet_x_s, d_matrix)
    y_s = derive_d_action(sheet_y_s, d_matrix)
    x_a = derive_d_action(sheet_x_a, d_matrix)
    y_a = derive_d_action(sheet_y_a, d_matrix)

    require(x_s == (
        (1, 0, 0), (0, 0, -1), (0, -1, 0)
    ), "derived X_S disagrees with THM-2775")
    require(y_s == (
        (0, 1, 0), (0, 0, 1), (1, 0, 0)
    ), "derived Y_S disagrees with THM-2775")
    require(x_a == (
        (1, 0, 0), (0, -1, 0), (0, 0, -1)
    ), "derived X_A disagrees with THM-2775")
    require(y_a == (
        (0, 0, -1), (1, 0, 0), (0, -1, 0)
    ), "derived Y_A disagrees with THM-2775")

    unit_three = identity(3)
    face_s = matmul(x_s, y_s)
    face_a = matmul(x_a, y_a)
    require(matpow(x_s, 2) == unit_three
            and matpow(y_s, 3) == unit_three
            and matpow(face_s, 4) == unit_three
            and matpow(face_s, 2) != unit_three,
            "derived S-frame triangle orders are not (2,3,4)")
    require(matpow(x_a, 2) == unit_three
            and matpow(y_a, 3) == unit_three
            and matpow(face_a, 3) == unit_three
            and face_a != unit_three,
            "derived A-frame triangle orders are not (2,3,3)")
    face_square = matpow(face_s, 2)
    require(face_square == (
        (-1, 0, 0), (0, -1, 0), (0, 0, 1)
    ), "order-four face square stopped being diag(-1,-1,1)")

    group_s = generated_group((x_s, y_s))
    group_a = generated_group((x_a, y_a))
    require(group_s == all_even_signed_permutations()
            and len(group_s) == 24,
            "derived S-frame image stopped being W(D3)")
    require(len(group_a) == 12
            and all(determinant_three(element) == 1
                    for element in group_a),
            "derived A-frame image stopped being the determinant-one half")

    state_rows = (
        (1, 1, 1),
        (1, -1, -1),
        (-1, 1, -1),
        (-1, -1, 1),
    )
    require(induced_state_action(x_s, state_rows) == sheet_x_s,
            "X_S stopped inducing the literal sheet transposition")
    require(induced_state_action(y_s, state_rows) == sheet_y_s,
            "Y_S stopped inducing the literal sheet three-cycle")
    require(induced_state_action(x_a, state_rows) == sheet_x_a,
            "X_A stopped inducing the literal double transposition")
    require(induced_state_action(y_a, state_rows) == sheet_y_a,
            "Y_A stopped inducing the literal sheet three-cycle")
    require({permutation_order(induced_state_action(matrix, state_rows))
             for matrix in (x_s, x_a)} == {2}
            and {permutation_order(induced_state_action(matrix, state_rows))
                 for matrix in (y_s, y_a)} == {3},
            "induced sheet orders changed")

    # H D=4I-J.  Thus the h_i are exactly centered roots 4r_i-sum(r_j);
    # the absolute common translation is deliberately absent.
    centered_matrix = matmul(state_rows, d_matrix)
    four_i_minus_j = tuple(
        tuple(3 if i == j else -1 for j in range(4))
        for i in range(4)
    )
    require(centered_matrix == four_i_minus_j,
            "half-Hadamard forms stopped being centered root states")
    diagonal = ((1,), (1,), (1,), (1,))
    require(matmul(d_matrix, diagonal) == ((0,), (0,), (0,)),
            "opposite-edge coordinates stopped killing common translation")
    require(matmul(centered_matrix, diagonal)
            == ((0,), (0,), (0,), (0,)),
            "centered states unexpectedly retained common translation")

    # Independent Newton/local audit of
    # U^3-4U^2+16tU-64t^2 at t=0.
    newton_points = ((0, 2), (1, 1), (2, 0), (3, 0))
    require(newton_points[1][1] - newton_points[0][1] == -1
            and newton_points[2][1] - newton_points[1][1] == -1
            and newton_points[3][1] - newton_points[2][1] == 0,
            "THM-2769 Newton slopes changed")
    root_valuations = (1, 1, 0)
    residual_quadratic = (1, -4, 16)
    residual_discriminant = (
        residual_quadratic[1] ** 2
        - 4 * residual_quadratic[0] * residual_quadratic[2]
    )
    require(residual_discriminant == -48,
            "two valuation-one roots stopped being residually distinct")
    require(3 * 4 ** 2 - 8 * 4 == 16,
            "valuation-zero root U=4 stopped being simple")
    parity_row = tuple(value % 2 for value in root_valuations)
    require(parity_row == (1, 1, 0),
            "THM-2769 local parity row stopped being 110")
    require(face_square == tuple(
        tuple((-1 if parity_row[i] else 1) if i == j else 0
              for j in range(3))
        for i in range(3)
    ), "face square stopped matching the local 110 inertia signs")

    even_code = {
        word for word in product((0, 1), repeat=3)
        if sum(word) % 2 == 0
    }
    require(even_code == {
        (0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1)
    }, "even valuation code changed")

    print("THM-2775 INDEPENDENT HOSTILE AUDIT")
    print("root_derived_XS=(d1,-d3,-d2) XA=(d1,-d2,-d3)")
    print("root_derived_YS=(d2,d3,d1) YA=(-d3,d1,-d2)")
    print("triangle_orders=S:(2,3,4) A:(2,3,3)")
    print("images=S:W(D3)=S4_order24 A:det_plus_order12")
    print("centered_state_actions=XS:(12),YS:(234),XA:(12)(34),YA:(123)")
    print("half_Hadamard_identity=H*D=4I-J")
    print("translation_hostile=diagonal_kernel_dimension1")
    print("absolute_root_recovery_requires_common_sum_e1")
    print("S_face_square=diag(-1,-1,1)")
    print("THM2769_Newton_valuations=(1,1,0)")
    print("THM2769_residual_discriminant=-48 simple_U4_derivative=16")
    print("face_square_equals_110_local_Kummer_inertia")
    print("boundary_code={000,110,101,011}")
    print("VERDICT=ACCEPT_WITH_CENTERED_ROOT_WORDING_GUARD")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
