#!/usr/bin/env python3
"""Exact standard-library audit for THM-3412.

The program constructs the full numerator differential

    nabla_q = g(f' d/dt - d/dx) + q g'

on B/(g,t)^q over Q.  Exact Fraction row reduction checks the finite kernel
and cokernel, the graded four-term complex, transition persistence, local
primitive, Jordan arms, collision controls, and nonsplit rational controls.
No floating point and no Python ``assert`` statement is used.
"""

from fractions import Fraction
from hashlib import sha256
from math import comb
import json


EXPECTED_SEMANTIC_SHA256 = "88110e08603773bebcfd88e2fe224badfa014bf3950d12800e5f5fd65007c3ae"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def trim(poly):
    values = [Fraction(value) for value in poly]
    while values and not values[-1]:
        values.pop()
    return tuple(values)


ZERO = ()
ONE = (Fraction(1),)
X = (Fraction(0), Fraction(1))


def poly_add(left, right):
    size = max(len(left), len(right))
    return trim([
        (left[i] if i < len(left) else 0)
        + (right[i] if i < len(right) else 0)
        for i in range(size)
    ])


def poly_scale(poly, scalar):
    scalar = Fraction(scalar)
    return trim([scalar * value for value in poly])


def poly_sub(left, right):
    return poly_add(left, poly_scale(right, -1))


def poly_mul(left, right):
    if not left or not right:
        return ZERO
    answer = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            answer[i + j] += left_value * right_value
    return trim(answer)


def poly_pow(poly, exponent):
    require(exponent >= 0, ("negative polynomial exponent", exponent))
    answer = ONE
    base = trim(poly)
    current = exponent
    while current:
        if current & 1:
            answer = poly_mul(answer, base)
        base = poly_mul(base, base)
        current //= 2
    return answer


def poly_deriv(poly):
    return trim([i * poly[i] for i in range(1, len(poly))])


def poly_divmod(dividend, divisor):
    divisor = trim(divisor)
    require(divisor, "division by zero polynomial")
    remainder = list(trim(dividend))
    quotient = [Fraction(0)] * max(1, len(remainder) - len(divisor) + 1)
    while remainder and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] += coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] -= coefficient * value
        remainder = list(trim(remainder))
    return trim(quotient), trim(remainder)


def poly_mod(poly, modulus):
    return poly_divmod(poly, modulus)[1]


def monomial(degree, coefficient=1):
    require(degree >= 0, ("negative monomial degree", degree))
    return trim([0] * degree + [Fraction(coefficient)])


def matrix_zero(rows, columns):
    return [[Fraction(0) for _ in range(columns)] for _ in range(rows)]


def matrix_identity(size):
    return [[Fraction(row == column) for column in range(size)]
            for row in range(size)]


def matrix_add(left, right):
    require(len(left) == len(right), "matrix row mismatch")
    if not left:
        return []
    require(len(left[0]) == len(right[0]), "matrix column mismatch")
    return [[left[row][column] + right[row][column]
             for column in range(len(left[0]))]
            for row in range(len(left))]


def matrix_scale(matrix, scalar):
    scalar = Fraction(scalar)
    return [[scalar * value for value in row] for row in matrix]


def matrix_mul(left, right):
    if not left:
        return []
    require(right, "right matrix has no rows")
    require(len(left[0]) == len(right), "matrix product mismatch")
    return [[
        sum(left[row][inner] * right[inner][column]
            for inner in range(len(right)))
        for column in range(len(right[0]))
    ] for row in range(len(left))]


def matrix_power(matrix, exponent):
    require(exponent >= 0, ("negative matrix exponent", exponent))
    require(matrix and len(matrix) == len(matrix[0]), "power needs square matrix")
    answer = matrix_identity(len(matrix))
    base = matrix
    current = exponent
    while current:
        if current & 1:
            answer = matrix_mul(answer, base)
        base = matrix_mul(base, base)
        current //= 2
    return answer


def matrix_is_zero(matrix):
    return all(not value for row in matrix for value in row)


def rref(matrix):
    rows = [list(map(Fraction, row)) for row in matrix]
    if not rows:
        return rows, []
    column_count = len(rows[0])
    pivot_columns = []
    pivot_row = 0
    for column in range(column_count):
        pivot = next((row for row in range(pivot_row, len(rows))
                      if rows[row][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        factor = rows[pivot_row][column]
        rows[pivot_row] = [value / factor for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                rows[row][entry] - factor * rows[pivot_row][entry]
                for entry in range(column_count)
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return rows, pivot_columns


def matrix_rank(matrix):
    return len(rref(matrix)[1])


def matrix_nullity(matrix):
    require(matrix, "nullity needs explicit rows")
    return len(matrix[0]) - matrix_rank(matrix)


def nullspace_columns(matrix):
    reduced, pivots = rref(matrix)
    column_count = len(matrix[0])
    free_columns = [column for column in range(column_count)
                    if column not in pivots]
    answer = []
    for free in free_columns:
        vector = [Fraction(0)] * column_count
        vector[free] = Fraction(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row][free]
        answer.append(vector)
    return answer


def matrix_from_columns(columns, row_count=None):
    if not columns:
        require(row_count is not None, "empty columns need row count")
        return matrix_zero(row_count, 0)
    rows = len(columns[0])
    require(all(len(column) == rows for column in columns), "ragged columns")
    return [[columns[column][row] for column in range(len(columns))]
            for row in range(rows)]


def matrix_columns(matrix):
    if not matrix:
        return []
    return [[matrix[row][column] for row in range(len(matrix))]
            for column in range(len(matrix[0]))]


def solve_in_column_basis(columns, target):
    require(columns, "empty coordinate basis")
    row_count = len(target)
    variable_count = len(columns)
    rows = [[columns[column][row] for column in range(variable_count)]
            + [target[row]] for row in range(row_count)]
    reduced, pivots = rref(rows)
    coefficient_pivots = [pivot for pivot in pivots if pivot < variable_count]
    require(len(coefficient_pivots) == variable_count,
            ("dependent coordinate basis", coefficient_pivots, variable_count))
    for row in reduced:
        if all(not row[column] for column in range(variable_count)):
            require(not row[-1], ("target outside basis", target))
    answer = [Fraction(0)] * variable_count
    for row, pivot in enumerate(pivots):
        if pivot < variable_count:
            answer[pivot] = reduced[row][-1]
    return answer


def restrict_operator(operator, invariant_columns):
    basis_matrix = matrix_from_columns(invariant_columns)
    images = matrix_columns(matrix_mul(operator, basis_matrix))
    coordinates = [solve_in_column_basis(invariant_columns, image)
                   for image in images]
    return matrix_from_columns(coordinates)


def polynomial_of_matrix(coefficients, matrix):
    require(matrix and len(matrix) == len(matrix[0]), "square matrix required")
    size = len(matrix)
    answer = matrix_zero(size, size)
    power = matrix_identity(size)
    for coefficient in coefficients:
        answer = matrix_add(answer, matrix_scale(power, coefficient))
        power = matrix_mul(power, matrix)
    return answer


def shifted(matrix, scalar):
    return matrix_add(matrix, matrix_scale(matrix_identity(len(matrix)), -scalar))


def check_nilpotent_blocks(matrix, scalar, blocks, label):
    require(sum(blocks) == len(matrix), (label, "block dimension", blocks, len(matrix)))
    nilpotent = shifted(matrix, scalar)
    maximum = max(blocks, default=0)
    profile = []
    for exponent in range(1, maximum + 1):
        power = matrix_power(nilpotent, exponent)
        actual = matrix_nullity(power)
        expected = sum(min(exponent, block) for block in blocks)
        require(actual == expected,
                (label, "Jordan nullity", exponent, actual, expected, blocks))
        profile.append(actual)
    require(matrix_is_zero(matrix_power(nilpotent, maximum)),
            (label, "top power nonzero", maximum))
    return profile


def check_primary_blocks(matrix, scalar, blocks, label):
    """Check the blocks at one eigenvalue while other eigenvalues remain."""
    require(sum(blocks) <= len(matrix),
            (label, "primary dimension", blocks, len(matrix)))
    primary = shifted(matrix, scalar)
    maximum = max(blocks, default=0)
    profile = []
    for exponent in range(1, maximum + 1):
        actual = matrix_nullity(matrix_power(primary, exponent))
        expected = sum(min(exponent, block) for block in blocks)
        require(actual == expected,
                (label, "primary nullity", exponent, actual, expected, blocks))
        profile.append(actual)
    return profile


def quotient_basis(g, q):
    degree = len(g) - 1
    require(degree >= 1 and g[-1] == 1, ("g must be monic", g))
    return [(j, a) for j in range(q) for a in range(degree * (q - j))]


def quotient_index(g, q):
    basis = quotient_basis(g, q)
    return basis, {item: index for index, item in enumerate(basis)}


def add_slice(vector, index, g, q, t_degree, polynomial, scalar=1):
    if t_degree < 0 or t_degree >= q:
        return
    modulus = poly_pow(g, q - t_degree)
    remainder = poly_mod(poly_scale(polynomial, scalar), modulus)
    for x_degree, value in enumerate(remainder):
        if value:
            vector[index[(t_degree, x_degree)]] += value


def nabla_matrix(f, g, q):
    basis, index = quotient_index(g, q)
    dimension = len(basis)
    matrix = matrix_zero(dimension, dimension)
    f_prime = poly_deriv(f)
    g_prime = poly_deriv(g)
    gf_prime = poly_mul(g, f_prime)
    for column, (j, a) in enumerate(basis):
        vector = [Fraction(0)] * dimension
        xa = monomial(a)
        if j:
            add_slice(vector, index, g, q, j - 1,
                      poly_mul(gf_prime, xa), j)
        if a:
            add_slice(vector, index, g, q, j,
                      poly_mul(g, monomial(a - 1)), -a)
        add_slice(vector, index, g, q, j,
                  poly_mul(g_prime, xa), q)
        for row, value in enumerate(vector):
            matrix[row][column] = value
    return matrix


def p_matrix(f, g, q):
    basis, index = quotient_index(g, q)
    dimension = len(basis)
    matrix = matrix_zero(dimension, dimension)
    for column, (j, a) in enumerate(basis):
        vector = [Fraction(0)] * dimension
        add_slice(vector, index, g, q, j, poly_mul(f, monomial(a)))
        if j + 1 < q:
            add_slice(vector, index, g, q, j + 1, monomial(a))
        for row, value in enumerate(vector):
            matrix[row][column] = value
    return matrix


def transition_matrix(g, q):
    source, _ = quotient_index(g, q)
    _, target_index = quotient_index(g, q + 1)
    target_dimension = len(quotient_basis(g, q + 1))
    matrix = matrix_zero(target_dimension, len(source))
    for column, (j, a) in enumerate(source):
        vector = [Fraction(0)] * target_dimension
        add_slice(vector, target_index, g, q + 1, j,
                  poly_mul(g, monomial(a)))
        for row, value in enumerate(vector):
            matrix[row][column] = value
    return matrix


def rectangle_basis(c, q):
    degree = len(c) - 1
    return [(j, a) for j in range(q) for a in range(degree)]


def reduction_matrix(g, c, q):
    source, _ = quotient_index(g, q)
    target = rectangle_basis(c, q)
    target_index = {item: index for index, item in enumerate(target)}
    matrix = matrix_zero(len(target), len(source))
    for column, (j, a) in enumerate(source):
        remainder = poly_mod(monomial(a), c)
        for degree, value in enumerate(remainder):
            if value:
                matrix[target_index[(j, degree)]][column] += value
    return matrix


def rectangle_p_matrix(f, c, q):
    basis = rectangle_basis(c, q)
    index = {item: position for position, item in enumerate(basis)}
    matrix = matrix_zero(len(basis), len(basis))
    for column, (j, a) in enumerate(basis):
        remainder = poly_mod(poly_mul(f, monomial(a)), c)
        for degree, value in enumerate(remainder):
            if value:
                matrix[index[(j, degree)]][column] += value
        if j + 1 < q:
            matrix[index[(j + 1, a)]][column] += 1
    return matrix


def verify_finite_exact(f, g, c, q, label):
    differential = nabla_matrix(f, g, q)
    reduction = reduction_matrix(g, c, q)
    dimension = len(differential)
    target_dimension = q * (len(c) - 1)
    require(matrix_is_zero(matrix_mul(reduction, differential)),
            (label, q, "rho*nabla nonzero"))
    require(matrix_rank(reduction) == target_dimension,
            (label, q, "rho not onto"))
    require(matrix_rank(differential) == dimension - target_dimension,
            (label, q, "wrong nabla rank", matrix_rank(differential),
             dimension, target_dimension))
    return differential, nullspace_columns(differential)


def sparse_clean(poly):
    return {key: Fraction(value) for key, value in poly.items() if value}


def sparse_add(left, right):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = answer.get(key, Fraction(0)) + value
    return sparse_clean(answer)


def sparse_scale(poly, scalar):
    scalar = Fraction(scalar)
    return sparse_clean({key: scalar * value for key, value in poly.items()})


def sparse_mul(left, right):
    answer = {}
    for (left_y, left_t), left_value in left.items():
        for (right_y, right_t), right_value in right.items():
            key = (left_y + right_y, left_t + right_t)
            answer[key] = answer.get(key, Fraction(0)) + left_value * right_value
    return sparse_clean(answer)


def sparse_pow(poly, exponent):
    answer = {(0, 0): Fraction(1)}
    base = sparse_clean(poly)
    current = exponent
    while current:
        if current & 1:
            answer = sparse_mul(answer, base)
        base = sparse_mul(base, base)
        current //= 2
    return answer


def sparse_dy(poly):
    return sparse_clean({(y_degree - 1, t_degree): y_degree * value
                         for (y_degree, t_degree), value in poly.items()
                         if y_degree})


def sparse_dt(poly):
    return sparse_clean({(y_degree, t_degree - 1): t_degree * value
                         for (y_degree, t_degree), value in poly.items()
                         if t_degree})


Y_SPARSE = {(1, 0): Fraction(1)}
T_SPARSE = {(0, 1): Fraction(1)}
U_SPARSE = sparse_add(Y_SPARSE, T_SPARSE)


def primitive_numerator(e, q):
    answer = {}
    for j in range(q):
        coefficient = Fraction(((-1) ** j) * comb(q - 1, j), e * q - j - 1)
        term = sparse_mul(
            sparse_pow(U_SPARSE, q - 1 - j),
            {(j + 1, 0): coefficient},
        )
        answer = sparse_add(answer, term)
    return sparse_clean(answer)


def monomial_nabla_sparse(poly, e, q):
    first = sparse_scale(
        sparse_mul({(e, 0): Fraction(1)},
                   sparse_add(sparse_dt(poly), sparse_scale(sparse_dy(poly), -1))),
        1,
    )
    correction = sparse_scale(
        sparse_mul({(e - 1, 0): Fraction(1)}, poly), q * e)
    return sparse_add(first, correction)


def sparse_to_quotient_vector(poly, g, q):
    _, index = quotient_index(g, q)
    vector = [Fraction(0)] * len(index)
    by_t = {}
    for (x_degree, t_degree), value in poly.items():
        coefficients = by_t.setdefault(t_degree, [])
        if len(coefficients) <= x_degree:
            coefficients.extend([Fraction(0)] * (x_degree + 1 - len(coefficients)))
        coefficients[x_degree] += value
    for t_degree, coefficients in by_t.items():
        add_slice(vector, index, g, q, t_degree, trim(coefficients))
    return vector


def graded_monomial_complex(e, q):
    # G_q=Q[y,t]/(y^e,t^q), differential qe*y^(e-1).
    g_basis = [(j, a) for j in range(q) for a in range(e)]
    g_index = {item: index for index, item in enumerate(g_basis)}
    c_basis = [(j, a) for j in range(q) for a in range(e - 1)]
    differential = matrix_zero(len(g_basis), len(g_basis))
    injection = matrix_zero(len(g_basis), len(c_basis))
    reduction = matrix_zero(len(c_basis), len(g_basis))
    for column, (j, a) in enumerate(g_basis):
        if a == 0:
            differential[g_index[(j, e - 1)]][column] = q * e
        if a < e - 1:
            reduction[c_basis.index((j, a))][column] = 1
    for column, (j, a) in enumerate(c_basis):
        injection[g_index[(j, a + 1)]][column] = 1
    require(matrix_rank(injection) == len(c_basis), (e, q, "graded injection"))
    require(matrix_is_zero(matrix_mul(differential, injection)),
            (e, q, "graded d*i"))
    require(matrix_rank(differential) == q, (e, q, "graded rank"))
    require(matrix_is_zero(matrix_mul(reduction, differential)),
            (e, q, "graded rho*d"))
    require(matrix_rank(reduction) == len(c_basis), (e, q, "graded rho"))
    return len(c_basis)


def multiplication_modulus_matrix(f, modulus):
    degree = len(modulus) - 1
    require(degree >= 1 and modulus[-1] == 1, "monic modulus required")
    matrix = matrix_zero(degree, degree)
    for column in range(degree):
        remainder = poly_mod(poly_mul(f, monomial(column)), modulus)
        for row, value in enumerate(remainder):
            matrix[row][column] = value
    return matrix


def require_annihilator(matrix, polynomial, proper_polynomial, label):
    require(matrix_is_zero(polynomial_of_matrix(polynomial, matrix)),
            (label, "annihilator fails", polynomial))
    require(not matrix_is_zero(polynomial_of_matrix(proper_polynomial, matrix)),
            (label, "proper candidate already kills", proper_polynomial))


def audit_one_root():
    cells = []
    exact_cells = 0
    graded_cells = 0
    primitive_cells = 0
    persistence_cells = 0
    ephemeral_cells = 0
    hostile = None
    for e in range(2, 6):
        g = monomial(e)
        c = monomial(e - 1)
        f = X
        previous_kernel = None
        previous_q = None
        for q in range(1, 5):
            differential, kernel = verify_finite_exact(f, g, c, q, ("one-root", e))
            exact_cells += 1
            p_ambient = p_matrix(f, g, q)
            p_kernel = restrict_operator(p_ambient, kernel)
            length = q * (e - 1)
            profile = check_nilpotent_blocks(
                p_kernel, 0, [length], ("one-root-kernel", e, q))

            coker_p = rectangle_p_matrix(f, c, q)
            coker_blocks = [q + e - 2 - 2 * j
                            for j in range(min(q, e - 1))]
            coker_profile = check_nilpotent_blocks(
                coker_p, 0, coker_blocks, ("one-root-coker", e, q))

            graded_dimension = graded_monomial_complex(e, q)
            require(graded_dimension == length, (e, q, "graded dimension"))
            graded_cells += 1

            numerator = primitive_numerator(e, q)
            expected_derivative = {(e, q - 1): Fraction(1)}
            require(monomial_nabla_sparse(numerator, e, q) == expected_derivative,
                    (e, q, "primitive derivative", numerator,
                     monomial_nabla_sparse(numerator, e, q)))
            vector = sparse_to_quotient_vector(numerator, g, q)
            require(matrix_is_zero(matrix_mul(
                differential, matrix_from_columns([vector]))),
                (e, q, "primitive not a cycle"))
            cyclic_columns = []
            current = matrix_from_columns([vector])
            for _ in range(length):
                cyclic_columns.append(matrix_columns(current)[0])
                current = matrix_mul(p_ambient, current)
            require(matrix_rank(matrix_from_columns(cyclic_columns)) == length,
                    (e, q, "primitive not cyclic", length))
            require(matrix_is_zero(current), (e, q, "primitive chain too long"))
            primitive_cells += 1

            if previous_kernel is not None:
                transition = transition_matrix(g, previous_q)
                image_ambient = matrix_mul(
                    transition, matrix_from_columns(previous_kernel))
                image_columns = matrix_columns(image_ambient)
                require(matrix_rank(image_ambient) == (q - 1) * (e - 1),
                        (e, q, "kernel transition rank"))
                image_coordinates = [solve_in_column_basis(kernel, column)
                                     for column in image_columns]
                image_in_kernel = matrix_from_columns(image_coordinates)
                require(matrix_rank(image_in_kernel) == (q - 1) * (e - 1),
                        (e, q, "kernel coordinates rank"))
                # A cyclic target has one submodule of this dimension; the
                # quotient must be the constant arm of length e-1.
                for exponent in range(1, e):
                    power = matrix_power(p_kernel, exponent)
                    combined = [power[row] + image_in_kernel[row]
                                for row in range(len(power))]
                    quotient_nullity = len(p_kernel) - matrix_rank(combined)
                    require(quotient_nullity == min(exponent, e - 1),
                            (e, q, "persistence quotient", exponent,
                             quotient_nullity))
                persistence_cells += 1

                reduction_next = reduction_matrix(g, c, q)
                transition_previous = transition_matrix(g, q - 1)
                require(matrix_is_zero(matrix_mul(reduction_next,
                                                  transition_previous)),
                        (e, q, "coker transition not zero"))
                ephemeral_cells += 1

            if e == 3 and q == 2:
                hostile = {
                    "e": e,
                    "q": q,
                    "kernel_blocks": [length],
                    "graded_kernel_blocks": coker_blocks,
                    "kernel_dimension": length,
                    "liftable_symbol_dimension": e - 1,
                    "graded_kernel_dimension": graded_dimension,
                    "deficit": graded_dimension - (e - 1),
                }

            cells.append({
                "e": e,
                "q": q,
                "kernel_block": length,
                "kernel_profile": profile,
                "coker_blocks": coker_blocks,
                "coker_profile": coker_profile,
                "deficit": (q - 1) * (e - 1),
            })
            previous_kernel = kernel
            previous_q = q
    require(hostile == {
        "e": 3,
        "q": 2,
        "kernel_blocks": [4],
        "graded_kernel_blocks": [3, 1],
        "kernel_dimension": 4,
        "liftable_symbol_dimension": 2,
        "graded_kernel_dimension": 4,
        "deficit": 2,
    }, ("hostile mismatch", hostile))
    return {
        "cells": cells,
        "exact_cells": exact_cells,
        "graded_cells": graded_cells,
        "primitive_cells": primitive_cells,
        "persistence_cells": persistence_cells,
        "ephemeral_cells": ephemeral_cells,
        "hostile": hostile,
    }


def product_polynomials(factors):
    answer = ONE
    for factor in factors:
        answer = poly_mul(answer, factor)
    return answer


def audit_control(label, f, g, c, block_checker, eta_ann,
                  eta_proper, q_values=(1, 2, 3)):
    records = []
    for q in q_values:
        _, kernel = verify_finite_exact(f, g, c, q, label)
        p_kernel = restrict_operator(p_matrix(f, g, q), kernel)
        block_summary = block_checker(p_kernel, q)
        require(len(kernel) == q * (len(c) - 1),
                (label, q, "kernel dimension", len(kernel)))
        records.append({"q": q, "dimension": len(kernel),
                        "blocks": block_summary})

    eta_matrix = multiplication_modulus_matrix(f, g)
    require_annihilator(eta_matrix, eta_ann, eta_proper, (label, "eta"))
    return records


def audit_named_controls():
    xm1 = poly_sub(X, ONE)
    xp1 = poly_add(X, ONE)
    x2m1 = poly_sub(poly_pow(X, 2), ONE)
    x2p1 = poly_add(poly_pow(X, 2), ONE)

    distinct_g = poly_mul(poly_pow(X, 2), poly_pow(xm1, 3))
    distinct_c = poly_mul(X, poly_pow(xm1, 2))

    def distinct_checker(matrix, q):
        zero = check_primary_blocks(matrix, 0, [q], ("distinct-zero", q))
        one = check_primary_blocks(matrix, 1, [2 * q], ("distinct-one", q))
        return {"at_0": [q], "at_1": [2 * q],
                "profiles": [zero, one]}

    distinct_eta = poly_mul(poly_pow(X, 2), poly_pow(xm1, 3))
    distinct_eta_proper = poly_mul(poly_pow(X, 2), poly_pow(xm1, 2))

    equal_g = poly_pow(x2m1, 2)
    equal_c = x2m1

    def equal_checker(matrix, q):
        profile = check_nilpotent_blocks(matrix, 1, [q, q], ("equal", q))
        return {"at_1": [q, q], "profile": profile}

    equal_eta = poly_pow(xm1, 2)
    equal_eta_proper = xm1

    unequal_g = poly_mul(poly_pow(xm1, 2), poly_pow(xp1, 3))
    unequal_c = poly_mul(xm1, poly_pow(xp1, 2))

    def unequal_checker(matrix, q):
        profile = check_nilpotent_blocks(matrix, 1, [q, 2 * q],
                                         ("unequal", q))
        return {"at_1": [q, 2 * q], "profile": profile}

    unequal_eta = poly_pow(xm1, 3)
    unequal_eta_proper = poly_pow(xm1, 2)

    nonsplit_g = poly_pow(x2p1, 2)
    nonsplit_c = x2p1

    def nonsplit_checker(matrix, q):
        irreducible = matrix_add(matrix_mul(matrix, matrix),
                                 matrix_identity(len(matrix)))
        profile = []
        for exponent in range(1, q + 1):
            actual = matrix_nullity(matrix_power(irreducible, exponent))
            expected = 2 * min(exponent, q)
            require(actual == expected,
                    ("nonsplit", q, exponent, actual, expected))
            profile.append(actual)
        require(matrix_is_zero(matrix_power(irreducible, q)),
                ("nonsplit top", q))
        return {"at_T2_plus_1": [q], "profile": profile}

    nonsplit_eta = poly_pow(x2p1, 2)
    nonsplit_eta_proper = x2p1

    def nonsplit_collision_checker(matrix, q):
        profile = check_nilpotent_blocks(matrix, -1, [q, q],
                                         ("nonsplit-collision", q))
        return {"at_minus_1": [q, q], "profile": profile}

    nonsplit_collision_eta = poly_pow(xp1, 2)
    nonsplit_collision_eta_proper = xp1

    return {
        "distinct": audit_control(
            "distinct", X, distinct_g, distinct_c, distinct_checker,
            distinct_eta, distinct_eta_proper),
        "equal_collision": audit_control(
            "equal-collision", poly_pow(X, 2), equal_g, equal_c,
            equal_checker, equal_eta, equal_eta_proper),
        "unequal_collision": audit_control(
            "unequal-collision", poly_pow(X, 2), unequal_g, unequal_c,
            unequal_checker, unequal_eta, unequal_eta_proper),
        "nonsplit": audit_control(
            "nonsplit", X, nonsplit_g, nonsplit_c, nonsplit_checker,
            nonsplit_eta, nonsplit_eta_proper),
        "nonsplit_collision": audit_control(
            "nonsplit-collision", poly_pow(X, 2), nonsplit_g, nonsplit_c,
            nonsplit_collision_checker,
            nonsplit_collision_eta, nonsplit_collision_eta_proper),
    }


def semantic_payload():
    one_root = audit_one_root()
    controls = audit_named_controls()
    return {
        "theorem": "THM-3412",
        "field": "Q exact Fraction arithmetic",
        "operator": "g(f' dt-dx)+qg'",
        "one_root": one_root,
        "controls": controls,
        "control_interpretation": {
            "distinct": "blocks q at 0 and 2q at 1",
            "equal_collision": "two q-blocks at 1; eta exponent 2",
            "unequal_collision": "q and 2q blocks at 1; eta exponent 3",
            "nonsplit": "one (T^2+1)^q rational arm",
            "nonsplit_collision": "two (T+1)^q arms over Q",
        },
        "scope": "linear-z response only; no JC(2) implication",
    }


def main():
    payload = semantic_payload()
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    semantic = sha256(encoded).hexdigest()
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest mismatch", semantic,
             EXPECTED_SEMANTIC_SHA256))

    one_root = payload["one_root"]
    hostile = one_root["hostile"]
    print("THM-3412 HAMILTONIAN PRINCIPAL-PART / PRUFER AUDIT")
    print("arithmetic: exact Fraction row reduction; assertion-independent")
    print("operator: nabla_q = g(f' partial_t-partial_x)+qg'")
    print(f"one-root finite exact cells: {one_root['exact_cells']}")
    print(f"graded four-term cells: {one_root['graded_cells']}")
    print(f"localized primitive/cyclic-generator cells: {one_root['primitive_cells']}")
    print(f"persistent kernel-transition cells: {one_root['persistence_cells']}")
    print(f"zero finite-cokernel-transition cells: {one_root['ephemeral_cells']}")
    print("sharp e=3,q=2 hostile: kernel [4], graded kernel/coker [3,1], "
          f"lift deficit {hostile['deficit']}")
    print("distinct values: blocks [q] at 0 and [2q] at 1")
    print("equal collision: two [q] arms at 1; diagonal eta exponent 2")
    print("unequal collision: [q,2q] arms at 1; diagonal eta exponent 3")
    print("nonsplit Q control: one (T^2+1)^q arm")
    print("nonsplit collided Q control: two (T+1)^q arms")
    print("scope: linear-z response only; no JC(2) implication")
    print(f"semantic_sha256={semantic}")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
