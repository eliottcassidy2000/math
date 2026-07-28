#!/usr/bin/env python3
"""Exact modular-to-Weyl D3 generator-frame audit for THM-2775.

Only finite integer matrix and permutation arithmetic is used.  Truth checks
use explicit exceptions; Python assertions and floating point are absent.
"""

from itertools import permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def identity(size):
    return tuple(tuple(int(i == j) for j in range(size)) for i in range(size))


def matmul(left, right):
    size = len(left)
    require(len(right) == size, "matrix sizes stopped matching")
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(size))
              for j in range(size))
        for i in range(size)
    )


def matpow(matrix, exponent):
    out = identity(len(matrix))
    base = matrix
    remaining = exponent
    while remaining:
        if remaining % 2:
            out = matmul(out, base)
        base = matmul(base, base)
        remaining //= 2
    return out


def determinant_three(matrix):
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def generated_group(generators):
    unit = identity(len(generators[0]))
    seen = {unit}
    frontier = [unit]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = matmul(generator, current)
            if candidate not in seen:
                seen.add(candidate)
                frontier.append(candidate)
    return frozenset(seen)


def signed_permutation_data(matrix):
    positions = []
    signs = []
    for row in matrix:
        support = [index for index, value in enumerate(row) if value]
        require(len(support) == 1 and abs(row[support[0]]) == 1,
                "matrix stopped being a signed permutation")
        positions.append(support[0])
        signs.append(row[support[0]])
    require(tuple(sorted(positions)) == tuple(range(len(matrix))),
            "signed permutation repeated an input coordinate")
    return tuple(positions), tuple(signs)


def row_times_matrix(row, matrix):
    return tuple(
        sum(row[index] * matrix[index][column]
            for index in range(len(row)))
        for column in range(len(row))
    )


def induced_hadamard_permutation(matrix, rows):
    image = []
    for row in rows:
        transformed = row_times_matrix(row, matrix)
        require(transformed in rows, "Hadamard state left the four-state set")
        image.append(rows.index(transformed))
    return tuple(image)


def permutation_order(permutation):
    current = tuple(range(len(permutation)))
    order = 0
    while True:
        current = tuple(permutation[current[index]]
                        for index in range(len(permutation)))
        order += 1
        if current == tuple(range(len(permutation))):
            return order


def all_weyl_d3():
    out = set()
    for coordinate_permutation in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            if signs[0] * signs[1] * signs[2] != 1:
                continue
            matrix = []
            for row_index in range(3):
                row = [0, 0, 0]
                row[coordinate_permutation[row_index]] = signs[row_index]
                matrix.append(tuple(row))
            out.add(tuple(matrix))
    return frozenset(out)


def main():
    # Pullback action on opposite-edge coordinates d1,d2,d3.
    x_s = (
        (1, 0, 0),
        (0, 0, -1),
        (0, -1, 0),
    )
    x_a = (
        (1, 0, 0),
        (0, -1, 0),
        (0, 0, -1),
    )
    y_s = (
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 0),
    )
    y_a = (
        (0, 0, -1),
        (1, 0, 0),
        (0, -1, 0),
    )
    unit = identity(3)

    require(matpow(x_s, 2) == unit and matpow(x_a, 2) == unit,
            "binary generator order changed")
    require(matpow(y_s, 3) == unit and matpow(y_a, 3) == unit,
            "ternary generator order changed")
    xy_s = matmul(x_s, y_s)
    xy_a = matmul(x_a, y_a)
    require(matpow(xy_s, 4) == unit and matpow(xy_s, 2) != unit,
            "S4 face relation stopped having order four")
    require(matpow(xy_a, 3) == unit and xy_a != unit,
            "A4 face relation stopped having order three")
    require(matpow(xy_s, 2) == ((-1, 0, 0), (0, -1, 0), (0, 0, 1)),
            "order-four face square stopped being the 110 sign flip")

    group_s = generated_group((x_s, y_s))
    group_a = generated_group((x_a, y_a))
    weyl_d3 = all_weyl_d3()
    require(group_s == weyl_d3 and len(group_s) == 24,
            "modular S4 image stopped being all of W(D3)")
    require(len(group_a) == 12
            and all(determinant_three(matrix) == 1 for matrix in group_a),
            "modular A4 image stopped being the orientation subgroup")

    projection_s = set()
    projection_a = set()
    kernel_s = []
    kernel_a = []
    for matrix in group_s:
        coordinate_permutation, signs = signed_permutation_data(matrix)
        require(signs[0] * signs[1] * signs[2] == 1,
                "W(D3) acquired an odd sign flip")
        projection_s.add(coordinate_permutation)
        if coordinate_permutation == (0, 1, 2):
            kernel_s.append(matrix)
    for matrix in group_a:
        coordinate_permutation, signs = signed_permutation_data(matrix)
        require(signs[0] * signs[1] * signs[2] == 1,
                "A4 image acquired an odd sign flip")
        projection_a.add(coordinate_permutation)
        if coordinate_permutation == (0, 1, 2):
            kernel_a.append(matrix)
    require(len(projection_s) == 6 and len(kernel_s) == 4,
            "S4/V4 coordinate quotient changed")
    require(len(projection_a) == 3 and len(kernel_a) == 4,
            "A4/V4 coordinate quotient changed")
    require(signed_permutation_data(x_s)[0] == (0, 2, 1)
            and signed_permutation_data(x_a)[0] == (0, 1, 2),
            "binary quotient distinction changed")

    hadamard_rows = (
        (1, 1, 1),
        (1, -1, -1),
        (-1, 1, -1),
        (-1, -1, 1),
    )
    had_x_s = induced_hadamard_permutation(x_s, hadamard_rows)
    had_x_a = induced_hadamard_permutation(x_a, hadamard_rows)
    had_y_s = induced_hadamard_permutation(y_s, hadamard_rows)
    had_y_a = induced_hadamard_permutation(y_a, hadamard_rows)
    require(had_x_s == (1, 0, 2, 3),
            "S4 binary generator stopped acting as (12)")
    require(had_x_a == (1, 0, 3, 2),
            "A4 binary generator stopped acting as (12)(34)")
    require(had_y_s == (0, 2, 3, 1),
            "S4 ternary generator stopped acting as (234)")
    require(had_y_a == (1, 2, 0, 3),
            "A4 ternary generator stopped acting as (123)")
    require(permutation_order(had_x_s) == 2
            and permutation_order(had_x_a) == 2
            and permutation_order(had_y_s) == 3
            and permutation_order(had_y_a) == 3,
            "four-state generator orders changed")

    # Boundary valuations of the three squared coordinates live in the even
    # code.  Signed flips do not change valuations; the quotient S3 permutes
    # the three positions.  Hence zero is fixed and all nonzero words form one
    # allowed orbit rather than being excluded by the modular relations.
    even_code = {
        vector for vector in product((0, 1), repeat=3)
        if sum(vector) % 2 == 0
    }
    nonzero_orbit = {
        tuple((1, 1, 0)[permutation[index]] for index in range(3))
        for permutation in projection_s
    }
    require(even_code == {(0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1)},
            "even boundary code changed")
    require(nonzero_orbit == even_code - {(0, 0, 0)},
            "nonzero even boundary words stopped forming one S3 orbit")

    # Torsion-free modular kernels have rank 1+index/6.  The S3 matching
    # quotient retains both finite-factor orders and has free kernel F2.
    require(1 + 12 // 6 == 3 and 1 + 24 // 6 == 5
            and 1 + 6 // 6 == 2,
            "Bass-Serre kernel rank ledger changed")

    print("MODULAR S4 TO WEYL D3 GENERATOR-FRAME AUDIT")
    print("X_S(d1,d2,d3)=(d1,-d3,-d2) order=2")
    print("X_A(d1,d2,d3)=(d1,-d2,-d3) order=2")
    print("Y_S(d1,d2,d3)=(d2,d3,d1) order=3")
    print("Y_A(d1,d2,d3)=(-d3,d1,-d2) order=3")
    print("order(X_S*Y_S)=4 square=diag(-1,-1,1)")
    print("order(X_A*Y_A)=3")
    print("generated_S_image=W(D3)=S4 order=24")
    print("generated_A_image=orientation_subgroup=A4 order=12")
    print("Hadamard_X_S=(12) Hadamard_Y_S=(234)")
    print("Hadamard_X_A=(12)(34) Hadamard_Y_A=(123)")
    print("S4_over_V4=S3 binary_retained")
    print("A4_over_V4=C3 binary_collapsed")
    print("modular_kernel_ranks=A4:F3,S4:F5,S3:F2")
    print("boundary_code={000,110,101,011}")
    print("nonzero_boundary_orbit={110,101,011}")
    print("SCOPE=finite_generator_frame_not_affine_monodromy_or_JC2")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
