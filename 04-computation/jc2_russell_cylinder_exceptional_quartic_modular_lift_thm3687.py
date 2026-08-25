"""Finite-field controls for actual lifts at the THM-3683 exceptional quartic.

This is a research probe, not a characteristic-zero theorem.  It specializes
the four exceptional parameters at split good primes, solves the actual-ring
J0 equation in the raw target packet, and (optionally) tests the coupled
J1/J2 system.  Feasibility modulo a prime is positive evidence only; failure
at one prime is not a characteristic-zero obstruction.
"""

import argparse
import gc

import sympy as sp
from flint import nmod_mat, nmod_poly


x, r = sp.symbols("x r")
F_QUARTIC = (
    72783360 * r**4
    - 77822208 * r**3
    - 28419741 * r**2
    + 7849770 * r
    - 1276420
)

P = x**2 * (x**2 - 1) ** 2
R1 = P * (1 - x**2)
R2 = P * (4 - 9 * x)
Q1 = x**5 + sp.Rational(9, 2) * x**4 - 2 * x**3 - sp.Rational(27, 4) * x**2 + x - sp.Rational(3, 4)
Q6 = Q1 - sp.Rational(259, 36) * P
P_OF_R = sp.Rational(520, 9) * r**2 - sp.Rational(1688, 81) * r - sp.Rational(5717, 729)
Q_R = sp.Poly(sp.expand(Q6 + P_OF_R * R1 + r * R2), x, domain=sp.QQ.poly_ring(r))


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def inv(value, prime):
    return pow(int(value) % prime, -1, prime)


def rational_mod(value, prime):
    value = sp.cancel(value)
    numerator, denominator = sp.fraction(value)
    return int(numerator) % prime * inv(int(denominator), prime) % prime


def specialize_q(prime, root):
    coefficients = [0] * (Q_R.degree() + 1)
    for (degree,), coefficient in Q_R.terms():
        coefficients[degree] = rational_mod(coefficient.as_expr().subs(r, root), prime)
    return nmod_poly(coefficients, prime)


def compiler(prime, root):
    one = nmod_poly([1], prime)
    X = nmod_poly([0, 1], prime)
    Q = specialize_q(prime, root)
    D = one + X * X * Q
    B = (D - one) * (D + 2) * (D + 2)
    C = X * D * (D + 2)
    E = Q * (D + 3)
    dB = 3 * X * X * D * (D + 2)
    dC = 2 * X * X * X * (D + 1)
    dE = 2 * (D + 1)
    require(C * C * E == B * (B + 4), "compiler relation")
    require((B.degree(), C.degree(), E.degree()) == (30, 21, 18), "generic restriction degrees")
    return B, C, E, dB, dC, dE


def packet(generators, deltas, weights, cutoff, prime):
    powers = []
    for generator, weight in zip(generators, weights):
        row = [nmod_poly([1], prime)]
        for _ in range(cutoff // weight):
            row.append(row[-1] * generator)
        powers.append(row)
    answer = []
    for i in range(cutoff // weights[0] + 1):
        for j in range(cutoff // weights[1] + 1):
            for k in range(cutoff // weights[2] + 1):
                metadata = (i, j, k)
                if sum(a * b for a, b in zip(metadata, weights)) > cutoff:
                    continue
                restriction = powers[0][i] * powers[1][j] * powers[2][k]
                delta = nmod_poly([], prime)
                for coordinate, exponent in enumerate(metadata):
                    if not exponent:
                        continue
                    term = nmod_poly([exponent], prime)
                    for index, coordinate_exponent in enumerate(metadata):
                        term *= powers[index][coordinate_exponent - int(index == coordinate)]
                    delta += term * deltas[coordinate]
                answer.append((metadata, restriction, delta))
    return answer


def coefficient(polynomial, degree):
    return int(polynomial[degree]) if degree <= polynomial.degree() else 0


def matrix_from_columns(columns, rhs, prime):
    row_count = max(rhs.degree(), *(column.degree() for column in columns)) + 1
    operator_data = []
    augmented_data = []
    for degree in range(row_count):
        row = [coefficient(column, degree) for column in columns]
        operator_data.extend(row)
        augmented_data.extend([*row, -coefficient(rhs, degree) % prime])
    operator = nmod_mat(row_count, len(columns), operator_data, prime)
    augmented = nmod_mat(row_count, len(columns) + 1, augmented_data, prime)
    return operator, augmented


def solve_j0(prime, root, cutoff):
    B, C, E, dB, dC, dE = compiler(prime, root)
    monomials = packet((B, C, E), (dB, dC, dE), (30, 21, 18), cutoff, prime)
    restrictions = [item[1] for item in monomials]
    columns = [C.derivative() * polynomial for polynomial in restrictions]
    columns += [-E.derivative() * polynomial for polynomial in restrictions]
    operator, augmented = matrix_from_columns(columns, -nmod_poly([1], prime), prime)
    operator_rank = operator.rank()
    reduced, augmented_rank = augmented.rref()
    feasible = operator_rank == augmented_rank
    solution = [0] * len(columns)
    pivot_columns = []
    pivot_rows = []
    if feasible:
        for row in range(augmented_rank):
            pivot = next((column for column in range(len(columns) + 1) if reduced[row, column]), None)
            require(pivot is not None and pivot < len(columns), "J0 pivot")
            pivot_columns.append(pivot)
            solution[pivot] = int(reduced[row, len(columns)])
        selected_data = []
        for row in range(operator.nrows()):
            selected_data.extend(int(operator[row, column]) for column in pivot_columns)
        selected = nmod_mat(operator.nrows(), len(pivot_columns), selected_data, prime)
        transpose_reduced, transpose_rank = selected.transpose().rref()
        require(transpose_rank == augmented_rank, "J0 selected-column rank")
        for row in range(transpose_rank):
            pivot_rows.append(
                next(
                    column
                    for column in range(operator.nrows())
                    if transpose_reduced[row, column]
                )
            )
    size = len(monomials)
    G_values, F_values = solution[:size], solution[size:]
    F1 = nmod_poly([], prime)
    G1 = nmod_poly([], prime)
    delta_F1 = nmod_poly([], prime)
    delta_G1 = nmod_poly([], prime)
    for value, (_metadata, restriction, delta) in zip(F_values, monomials):
        F1 += value * restriction
        delta_F1 += value * delta
    for value, (_metadata, restriction, delta) in zip(G_values, monomials):
        G1 += value * restriction
        delta_G1 += value * delta
    if feasible:
        require(C.derivative() * G1 - F1 * E.derivative() == nmod_poly([1], prime), "J0 residual")
    return {
        "compiler": (B, C, E, dB, dC, dE),
        "monomials": monomials,
        "rank": operator_rank,
        "augmented": augmented_rank,
        "rows": operator.nrows(),
        "columns": operator.ncols(),
        "F1": F1,
        "G1": G1,
        "delta_F1": delta_F1,
        "delta_G1": delta_G1,
        "pivot_columns": pivot_columns,
        "pivot_rows": pivot_rows,
    }


def coupled_ranks(prime, root, predecessor, cutoff):
    B, C, E, dB, dC, dE = predecessor["compiler"]
    F1, G1 = predecessor["F1"], predecessor["G1"]
    delta_F1, delta_G1 = predecessor["delta_F1"], predecessor["delta_G1"]
    monomials = packet((B, C, E), (dB, dC, dE), (30, 21, 18), cutoff, prime)
    restrictions = [item[1] for item in monomials]
    zero = nmod_poly([], prime)
    known_j1 = 2 * C.derivative() * dE + F1.derivative() * G1 - F1 * G1.derivative() - 2 * dC * E.derivative()
    known_j2 = (
        3 * C.derivative() * delta_G1
        + 2 * F1.derivative() * dE
        + dC.derivative() * G1
        - F1 * dE.derivative()
        - 2 * dC * G1.derivative()
        - 3 * delta_F1 * E.derivative()
    )
    columns = []
    for polynomial in restrictions:
        columns.append((-2 * polynomial * E.derivative(), polynomial.derivative() * G1 - 2 * polynomial * G1.derivative()))
    for polynomial in restrictions:
        columns.append((2 * C.derivative() * polynomial, 2 * F1.derivative() * polynomial - F1 * polynomial.derivative()))
    for polynomial in restrictions:
        columns.append((zero, -3 * polynomial * E.derivative()))
    for polynomial in restrictions:
        columns.append((zero, 3 * C.derivative() * polynomial))
    rows_j1 = max(known_j1.degree(), *(pair[0].degree() for pair in columns)) + 1
    rows_j2 = max(known_j2.degree(), *(pair[1].degree() for pair in columns)) + 1
    operator_data = []
    augmented_data = []
    for degree in range(rows_j1):
        row = [coefficient(pair[0], degree) for pair in columns]
        operator_data.extend(row)
        augmented_data.extend([*row, -coefficient(known_j1, degree) % prime])
    for degree in range(rows_j2):
        row = [coefficient(pair[1], degree) for pair in columns]
        operator_data.extend(row)
        augmented_data.extend([*row, -coefficient(known_j2, degree) % prime])
    operator = nmod_mat(rows_j1 + rows_j2, len(columns), operator_data, prime)
    augmented = nmod_mat(rows_j1 + rows_j2, len(columns) + 1, augmented_data, prime)
    ranks = (operator.rank(), augmented.rank())
    size = len(restrictions)
    j1_data = []
    for row in range(rows_j1):
        j1_data.extend(int(operator[row, column]) for column in range(2 * size))
    j1_operator = nmod_mat(rows_j1, 2 * size, j1_data, prime)
    low_reach = max(
        pair[1].degree()
        for pair in columns[2 * size :]
    ) + 1
    tail_rows = [*range(rows_j1), *range(rows_j1 + low_reach, rows_j1 + rows_j2)]
    tail_data = []
    tail_augmented_data = []
    for row in tail_rows:
        values = [int(operator[row, column]) for column in range(2 * size)]
        tail_data.extend(values)
        tail_augmented_data.extend([*values, int(augmented[row, len(columns)])])
    tail_operator = nmod_mat(len(tail_rows), 2 * size, tail_data, prime)
    tail_augmented = nmod_mat(len(tail_rows), 2 * size + 1, tail_augmented_data, prime)
    structure = {
        "j1_rank": j1_operator.rank(),
        "low_reach": low_reach,
        "tail_rows": len(tail_rows),
        "tail_rank": tail_operator.rank(),
        "tail_augmented": tail_augmented.rank(),
    }
    return len(monomials), rows_j1, rows_j2, len(columns), ranks, structure


def structured_coupled_solution(prime, root, predecessor, cutoff=375):
    """Solve J1,J2 by a 395 + 174 + 395 selected-direction factorization.

    The first/third stable blocks use the same Bezout operator.  We solve J1,
    use 177 unreachable high J2 rows plus the one omitted low-row residual to
    form a 178-coordinate response of rank 174, select 174 kernel directions,
    and reuse the Bezout square for low J2.
    """
    B, C, E, dB, dC, dE = predecessor["compiler"]
    F1, G1 = predecessor["F1"], predecessor["G1"]
    delta_F1, delta_G1 = predecessor["delta_F1"], predecessor["delta_G1"]
    monomials = packet((B, C, E), (dB, dC, dE), (30, 21, 18), cutoff, prime)
    restrictions = [item[1] for item in monomials]
    size = len(restrictions)
    known_j1 = 2 * C.derivative() * dE + F1.derivative() * G1 - F1 * G1.derivative() - 2 * dC * E.derivative()
    known_j2 = (
        3 * C.derivative() * delta_G1
        + 2 * F1.derivative() * dE
        + dC.derivative() * G1
        - F1 * dE.derivative()
        - 2 * dC * G1.derivative()
        - 3 * delta_F1 * E.derivative()
    )
    l0_columns = [-(polynomial * E.derivative()) for polynomial in restrictions]
    l0_columns += [C.derivative() * polynomial for polynomial in restrictions]
    m_columns = [polynomial.derivative() * G1 - 2 * polynomial * G1.derivative() for polynomial in restrictions]
    m_columns += [2 * F1.derivative() * polynomial - F1 * polynomial.derivative() for polynomial in restrictions]
    low_rows = max(column.degree() for column in l0_columns) + 1
    rows_j1 = max(known_j1.degree(), *(column.degree() for column in l0_columns)) + 1
    rows_j2 = max(known_j2.degree(), *(column.degree() for column in m_columns), low_rows - 1) + 1
    require(rows_j1 == low_rows, "structured common low row count")

    l0_data = []
    for degree in range(low_rows):
        l0_data.extend(coefficient(column, degree) for column in l0_columns)
    l0 = nmod_mat(low_rows, 2 * size, l0_data, prime)
    l0_reduced, l0_rank = l0.rref()
    require(l0_rank == low_rows - 1, "structured Bezout rank")
    pivot_columns = []
    for row in range(l0_rank):
        pivot_columns.append(
            next(column for column in range(2 * size) if l0_reduced[row, column])
        )
    selected_data = []
    for row in range(low_rows):
        selected_data.extend(int(l0[row, column]) for column in pivot_columns)
    selected = nmod_mat(low_rows, l0_rank, selected_data, prime)
    transpose_reduced, transpose_rank = selected.transpose().rref()
    require(transpose_rank == l0_rank, "structured Bezout selected rank")
    pivot_rows = []
    for row in range(transpose_rank):
        pivot_rows.append(next(column for column in range(low_rows) if transpose_reduced[row, column]))
    omitted_rows = [row for row in range(low_rows) if row not in set(pivot_rows)]
    require(len(omitted_rows) == 1, "structured one omitted low row")
    omitted_row = omitted_rows[0]

    square_data = []
    for row in pivot_rows:
        square_data.extend(int(l0[row, column]) for column in pivot_columns)
    square = nmod_mat(l0_rank, l0_rank, square_data, prime)
    free_columns = [column for column in range(2 * size) if column not in set(pivot_columns)]
    rhs_data = []
    inv2 = inv(2, prime)
    for row in pivot_rows:
        rhs_data.append(-coefficient(known_j1, row) * inv2 % prime)
        rhs_data.extend(-int(l0[row, column]) % prime for column in free_columns)
    rhs = nmod_mat(l0_rank, 1 + len(free_columns), rhs_data, prime)
    solved = square.solve(rhs)

    def vector_from_solution(rhs_column, free_column=None):
        answer = [0] * (2 * size)
        for row, column in enumerate(pivot_columns):
            answer[column] = int(solved[row, rhs_column])
        if free_column is not None:
            answer[free_column] = 1
        return answer

    base_x = vector_from_solution(0)
    directions = [
        vector_from_solution(index + 1, column)
        for index, column in enumerate(free_columns)
    ]

    # Evaluate every affine J1 solution on J2 in two dense finite-field matrix
    # products.  This avoids thousands of separate Python-level dot products.
    scenario_count = 1 + len(free_columns)
    x_matrix_data = []
    free_position = {column: index for index, column in enumerate(free_columns)}
    pivot_position = {column: index for index, column in enumerate(pivot_columns)}
    for variable in range(2 * size):
        if variable in pivot_position:
            solved_row = pivot_position[variable]
            x_matrix_data.extend(int(solved[solved_row, column]) for column in range(scenario_count))
        else:
            position = free_position[variable]
            x_matrix_data.extend(int(column == position + 1) for column in range(scenario_count))
    x_matrix = nmod_mat(2 * size, scenario_count, x_matrix_data, prime)
    m_data = []
    for row in range(rows_j2):
        m_data.extend(coefficient(column, row) for column in m_columns)
    m_matrix = nmod_mat(rows_j2, 2 * size, m_data, prime)
    mx = m_matrix * x_matrix

    inv3 = inv(3, prime)
    low_rhs_data = []
    for row in pivot_rows:
        for scenario in range(scenario_count):
            total = int(mx[row, scenario])
            if scenario == 0:
                total += coefficient(known_j2, row)
            low_rhs_data.append(-total * inv3 % prime)
    low_rhs = nmod_mat(l0_rank, scenario_count, low_rhs_data, prime)
    z_solved = square.solve(low_rhs)
    base_z_pivot = [int(z_solved[row, 0]) for row in range(l0_rank)]
    direction_z_pivots = [
        [int(z_solved[row, scenario + 1]) for row in range(l0_rank)]
        for scenario in range(len(free_columns))
    ]

    high_rows = list(range(low_rows, rows_j2))
    responses = []
    l0_omitted = [int(l0[omitted_row, column]) for column in pivot_columns]
    for scenario in range(scenario_count):
        answer = []
        for row in high_rows:
            total = int(mx[row, scenario])
            if scenario == 0:
                total += coefficient(known_j2, row)
            answer.append(total % prime)
        omitted = int(mx[omitted_row, scenario])
        omitted += 3 * sum(
            int(z_solved[row, scenario]) * value
            for row, value in enumerate(l0_omitted)
        )
        if scenario == 0:
            omitted += coefficient(known_j2, omitted_row)
        answer.append(omitted % prime)
        responses.append(answer)
    base_response = responses[0]
    direction_responses = responses[1:]
    require(len(base_response) == 178, "structured response dimension")
    response_data = []
    for row in range(len(base_response)):
        response_data.extend(direction[row] for direction in direction_responses)
    response_matrix = nmod_mat(len(base_response), len(directions), response_data, prime)
    response_reduced, response_rank = response_matrix.rref()
    require(response_rank == 174, "structured response rank")
    selected_direction_indices = []
    for row in range(response_rank):
        selected_direction_indices.append(
            next(column for column in range(len(directions)) if response_reduced[row, column])
        )
    response_selected_data = []
    for row in range(len(base_response)):
        response_selected_data.extend(
            direction_responses[column][row] for column in selected_direction_indices
        )
    response_selected = nmod_mat(
        len(base_response), response_rank, response_selected_data, prime
    )
    response_transpose_reduced, response_transpose_rank = response_selected.transpose().rref()
    require(response_transpose_rank == response_rank, "structured response selected rank")
    response_rows = []
    for row in range(response_rank):
        response_rows.append(
            next(
                column
                for column in range(len(base_response))
                if response_transpose_reduced[row, column]
            )
        )
    response_square_data = []
    for row in response_rows:
        response_square_data.extend(
            direction_responses[column][row] for column in selected_direction_indices
        )
    response_square = nmod_mat(response_rank, response_rank, response_square_data, prime)
    response_rhs = nmod_mat(
        response_rank, 1, [-base_response[row] % prime for row in response_rows], prime
    )
    parameters = response_square.solve(response_rhs)

    x_values = list(base_x)
    z_pivot_values = list(base_z_pivot)
    for position, direction_index in enumerate(selected_direction_indices):
        scalar = int(parameters[position, 0])
        x_values = [(a + scalar * b) % prime for a, b in zip(x_values, directions[direction_index])]
        z_pivot_values = [
            (a + scalar * b) % prime
            for a, b in zip(z_pivot_values, direction_z_pivots[direction_index])
        ]
    z_values = [0] * (2 * size)
    for value, column in zip(z_pivot_values, pivot_columns):
        z_values[column] = value

    j1_residual = known_j1 + 2 * sum(
        (value * column for value, column in zip(x_values, l0_columns)),
        nmod_poly([], prime),
    )
    j2_residual = known_j2
    j2_residual += sum(
        (value * column for value, column in zip(x_values, m_columns)),
        nmod_poly([], prime),
    )
    j2_residual += 3 * sum(
        (value * column for value, column in zip(z_values, l0_columns)),
        nmod_poly([], prime),
    )
    require(not j1_residual and not j2_residual, "structured full residuals")
    return {
        "monomials": monomials,
        "pivot_columns": pivot_columns,
        "pivot_rows": pivot_rows,
        "omitted_row": omitted_row,
        "free_columns": free_columns,
        "selected_direction_indices": selected_direction_indices,
        "response_rows": response_rows,
        "x_values": x_values,
        "z_values": z_values,
        "low_rows": low_rows,
        "rows_j2": rows_j2,
        "response_rank": response_rank,
    }


def roots_mod_prime(prime):
    return [value for value in range(prime) if int(F_QUARTIC.subs(r, value)) % prime == 0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=137)
    parser.add_argument("--j0-cutoff", type=int, default=195)
    parser.add_argument("--coupled-cutoff", type=int)
    parser.add_argument("--root", type=int, action="append")
    args = parser.parse_args()
    prime = args.prime
    roots = args.root if args.root is not None else roots_mod_prime(prime)
    require(roots, "quartic has no selected roots modulo prime")
    print(f"prime={prime};roots={','.join(map(str, roots))}")
    for root in roots:
        require(int(F_QUARTIC.subs(r, root)) % prime == 0, "selected quartic root")
        predecessor = solve_j0(prime, root, args.j0_cutoff)
        print(
            f"root={root};J0_cutoff={args.j0_cutoff};monomials={len(predecessor['monomials'])};"
            f"shape={predecessor['rows']}x{predecessor['columns']};rank={predecessor['rank']};"
            f"augmented={predecessor['augmented']};feasible={predecessor['rank']==predecessor['augmented']}"
        )
        if args.coupled_cutoff is not None and predecessor["rank"] == predecessor["augmented"]:
            packet_size, rows_j1, rows_j2, columns, ranks, structure = coupled_ranks(
                prime, root, predecessor, args.coupled_cutoff
            )
            print(
                f"root={root};J12_cutoff={args.coupled_cutoff};monomials={packet_size};"
                f"shape={rows_j1}+{rows_j2}x{columns};rank={ranks[0]};augmented={ranks[1]};"
                f"feasible={ranks[0]==ranks[1]}"
            )
            print(
                f"root={root};J12_structure;j1_rank={structure['j1_rank']};"
                f"low_reach={structure['low_reach']};tail_rows={structure['tail_rows']};"
                f"tail_rank={structure['tail_rank']};tail_augmented={structure['tail_augmented']}"
            )
        del predecessor
        gc.collect()


if __name__ == "__main__":
    main()
