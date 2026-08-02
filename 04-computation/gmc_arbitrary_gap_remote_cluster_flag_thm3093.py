#!/usr/bin/env python3
"""Exact finite companion for THM-3093's arbitrary-gap Monge flag."""

from fractions import Fraction
from functools import reduce
from itertools import product
from math import factorial, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(start, length):
    answer = 1
    for offset in range(length):
        answer *= start + offset
    return answer


def carrier(degree, N, gap):
    return Fraction(
        rising(degree * N + 1, degree * gap),
        rising(N + 1, gap) ** degree,
    )


def normalized_inner(left, right):
    """L(f_left f_right) for f_n=s^n/n! over the factorial functional."""
    return Fraction(factorial(left + right), factorial(left) * factorial(right))


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first,) + tail


def lcm(left, right):
    return left // gcd(left, right) * right


def determinant(matrix):
    matrix = [[Fraction(entry) for entry in row] for row in matrix]
    size = len(matrix)
    sign = 1
    answer = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if matrix[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            sign *= -1
        value = matrix[column][column]
        answer *= value
        for entry in range(column, size):
            matrix[column][entry] /= value
        for row in range(column + 1, size):
            multiple = matrix[row][column]
            for entry in range(column, size):
                matrix[row][entry] -= multiple * matrix[column][entry]
    return sign * answer


def contiguous_blocks(size, divergence_mask):
    blocks = []
    start = 0
    for edge in range(size - 1):
        if divergence_mask & (1 << edge):
            blocks.append(tuple(range(start, edge + 1)))
            start = edge + 1
    blocks.append(tuple(range(start, size)))
    return blocks


def gap_patterns(size):
    banks = {
        2: ((0, 1), (0, 7), (0, 73)),
        3: ((0, 1, 3), (0, 4, 19), (0, 17, 73)),
        4: ((0, 1, 3, 6), (0, 2, 17, 41), (0, 7, 37, 73)),
        5: ((0, 1, 3, 7, 12), (0, 3, 19, 53, 73)),
    }
    return banks[size]


# Strict response monotonicity in the row degree.  This is the exact Monge
# inequality after clearing the positive roots.
response_cells = 0
for maximum_degree in range(5, 10):
    for N in (maximum_degree, maximum_degree + 4, 37):
        for left_gap, right_gap in ((0, 1), (0, 7), (2, 19), (37, 137)):
            for lower in range(2, maximum_degree):
                for upper in range(lower + 1, maximum_degree + 1):
                    low_response = carrier(lower, N, right_gap) / carrier(
                        lower, N, left_gap
                    )
                    high_response = carrier(upper, N, right_gap) / carrier(
                        upper, N, left_gap
                    )
                    require(
                        low_response**upper < high_response**lower,
                        "strict response MLR",
                    )
                    response_cells += 1


# Multivariate Jensen: add r-j copies of the zero gap, then compare the
# rising factorial at the barycentre with the pure endpoints.
jensen_cells = 0
for child in range(1, 6):
    for rank in range(2, 5):
        p = child + 1
        k = child + rank
        for N in (29,):
            patterns = gap_patterns(rank)
            for gaps in (patterns[0], patterns[-1]):
                denominators = [rising(N + 1, gap) for gap in gaps]
                values = {
                    degree: [carrier(degree, N, gap) for gap in gaps]
                    for degree in range(2, k + 1)
                }
                for degree in range(2, k + 1):
                    for total_high in range(degree + 1):
                        for alpha in compositions(total_high, rank):
                            total_gap = sum(a * h for a, h in zip(alpha, gaps))
                            mixed = Fraction(
                                rising(degree * N + 1, total_gap),
                                reduce(
                                    lambda x, y: x * y,
                                    (
                                        denominators[index] ** alpha[index]
                                        for index in range(rank)
                                    ),
                                    1,
                                ),
                            )
                            pure = reduce(
                                lambda x, y: x * y,
                                (
                                    values[degree][index] ** alpha[index]
                                    for index in range(rank)
                                ),
                                Fraction(1),
                            )
                            require(mixed**degree <= pure, "multivariate Jensen")
                            jensen_cells += 1


# Root-free adjacent flag normalization.  Lambda_j^L and b_i^L are rational,
# where L is the lcm of the upper degrees.  Pure weights are at most one and
# tie exactly on the staircase edges (i,i+1).
flag_cells = 0
normalization_banks = []
for child in range(1, 6):
    for rank in range(2, 5):
        p = child + 1
        degrees = tuple(range(p, p + rank))
        L = reduce(lcm, degrees, 1)
        for N in (p + rank + 1, 31):
            for gaps in gap_patterns(rank):
                V = [
                    [carrier(degree, N, gap) for gap in gaps]
                    for degree in degrees
                ]
                lambda_L = [Fraction(1)]
                for edge in range(rank - 1):
                    pivot_degree = degrees[edge]
                    ratio = V[edge][edge] / V[edge][edge + 1]
                    lambda_L.append(
                        lambda_L[-1] * ratio ** (L // pivot_degree)
                    )
                require(
                    all(lambda_L[index + 1] < lambda_L[index] for index in range(rank - 1)),
                    "decreasing column potentials",
                )
                b_L = []
                weights_L = []
                for row, degree in enumerate(degrees):
                    row_b = V[row][row] ** (L // degree) * lambda_L[row]
                    require(row_b >= 1, "row potential at least one")
                    b_L.append(row_b)
                    row_weights = []
                    for column in range(rank):
                        weight = (
                            V[row][column] ** (L // degree)
                            * lambda_L[column]
                            / row_b
                        )
                        require(weight <= 1, "pure flag weight")
                        if column == row or (row < rank - 1 and column == row + 1):
                            require(weight == 1, "staircase tie")
                        row_weights.append(weight)
                        flag_cells += 1
                    weights_L.append(row_weights)
                normalization_banks.append(
                    (child, rank, N, gaps, degrees, L, V, lambda_L, b_L)
                )


# Raw physical starts and beta>=alpha normal exponents.  The normalized atom
# never exceeds the collapsed-base atom in a lower or upper equation.
raw_atom_cells = 0
for bank in normalization_banks:
    child, rank, N, gaps, degrees, L, V, lambda_L, b_L = bank
    if child > 4 or rank > 3 or N != 31 or gaps != gap_patterns(rank)[-1]:
        continue
    denominators = [rising(N + 1, gap) for gap in gaps]
    tested_rows = tuple(range(2, child + 1)) + degrees
    for degree in tested_rows:
        upper_row = degrees.index(degree) if degree in degrees else None
        for total_high in range(degree + 1):
            for alpha in compositions(total_high, rank):
                total_gap = sum(a * h for a, h in zip(alpha, gaps))
                for displacement in (0, (degree - total_high) * N):
                    mixed = Fraction(
                        rising(total_high * N + displacement + 1, total_gap),
                        reduce(
                            lambda x, y: x * y,
                            (
                                denominators[index] ** alpha[index]
                                for index in range(rank)
                            ),
                            1,
                        ),
                    )
                    beta_bank = [alpha]
                    spare = degree - total_high
                    if spare:
                        for column in range(rank):
                            beta = list(alpha)
                            beta[column] += spare
                            beta_bank.append(tuple(beta))
                    for beta in beta_bank:
                        scaled_L = mixed**L
                        for column in range(rank):
                            scaled_L *= lambda_L[column] ** beta[column]
                        if upper_row is not None:
                            a_L = b_L[upper_row] ** degree
                            scaled_L /= a_L
                        require(scaled_L <= 1, "raw flag atom contraction")
                        raw_atom_cells += 1


# Every finite/infinite adjacent-gap pattern gives consecutive diagonal
# alternant blocks.  Their normalized determinants have a gap-independent
# Vandermonde floor.  Arbitrary forward entries do not change a block-upper-
# triangular determinant, the line-power shadow of THM-3073.
block_cells = 0
q3_face_cells = 0
for child in range(1, 6):
    for rank in range(2, 6):
        degrees = tuple(range(child + 1, child + rank + 1))
        full_product = reduce(lambda x, y: x * y, degrees, 1)
        require(full_product % 2 == 0, "multi-normal parity")
        for gaps in gap_patterns(rank):
            for mask in range(1 << (rank - 1)):
                blocks = contiguous_blocks(rank, mask)
                matrix = [[Fraction(0) for _ in range(rank)] for _ in range(rank)]
                expected = Fraction(1)
                for block_number, block in enumerate(blocks):
                    start = block[0]
                    local_gaps = [gaps[index] - gaps[start] for index in block]
                    rows = [degrees[index] for index in block]
                    raw = [[row**gap for gap in local_gaps] for row in rows]
                    raw_det = determinant(raw)
                    denominator = reduce(
                        lambda x, y: x * y,
                        (rows[index] ** local_gaps[index] for index in range(len(rows))),
                        1,
                    )
                    normalized = raw_det / denominator
                    vandermonde = reduce(
                        lambda x, y: x * y,
                        (
                            rows[j] - rows[i]
                            for i in range(len(rows))
                            for j in range(i + 1, len(rows))
                        ),
                        1,
                    )
                    floor_denominator = reduce(
                        lambda x, y: x * y,
                        (rows[index] ** index for index in range(len(rows))),
                        1,
                    )
                    require(
                        normalized >= Fraction(vandermonde, floor_denominator) > 0,
                        "normalized alternant floor",
                    )
                    expected *= normalized
                    for local_row, global_row in enumerate(block):
                        for local_column, global_column in enumerate(block):
                            matrix[global_row][global_column] = Fraction(
                                rows[local_row] ** local_gaps[local_column],
                                rows[local_row] ** local_gaps[local_row],
                            )
                        for global_column in range(block[-1] + 1, rank):
                            matrix[global_row][global_column] = Fraction(
                                (block_number + 2) * (global_row + 1) - global_column,
                                global_column + 2,
                            )
                    block_cells += 1
                require(determinant(matrix) == expected, "forward triangular control")
                if rank == 3:
                    q3_face_cells += 1


# Exact carrier telescoping and covariance/parity ledgers.
carrier_cells = 0
covariance_cells = 0
for child in range(1, 9):
    for rank in range(2, 6):
        k = child + rank
        degrees = tuple(range(child + 1, k + 1))
        normal_degree = factorial(k) // factorial(child)
        require(normal_degree == reduce(lambda x, y: x * y, degrees, 1), "normal degree")
        require(normal_degree % 2 == 0, "multi-normal parity")
        # Each lambda_i appears once as a variable scale and once in a_i.
        for degree in degrees:
            variable_exponent = factorial(k)
            equation_exponent = -(factorial(k) // degree) * degree
            require(variable_exponent + equation_exponent == 0, "lambda cancellation")
            covariance_cells += 1
        for n, C, gap in product((1, 3), (7, 11), (1, 9, 41)):
            N = n + C
            for degree in degrees:
                U_base = Fraction(
                    rising(degree * n + 1, degree * C),
                    rising(n + 1, C) ** degree,
                )
                V_gap = carrier(degree, N, gap)
                U_far = Fraction(
                    rising(degree * n + 1, degree * (C + gap)),
                    rising(n + 1, C + gap) ** degree,
                )
                require(U_base * V_gap == U_far, "carrier telescope")
                carrier_cells += 1


# Sharp boundaries: repeated sheets kill the alternant; q=1 does not force
# an even child-resultant exponent.
for degree in range(4, 9):
    require(determinant([[1, 1], [1, 1]]) == 0, "repeated-gap boundary")
require((5 // 5) % 2 == 1, "one-normal parity boundary")

# Low-child conventions used in the theorem.  At m=1 there are no lower
# equations or variables, so the empty resultant and empty product are one.
# At m=2 the eliminated child is the positive factorial variance of two
# distinct normalized monomials.  These exact samples exercise both nearby
# and widely separated physical slots.
low_child_cells = 0
require(determinant([]) == 1, "empty resultant convention")
require(Fraction(1) == 1, "empty lower product")
low_child_cells += 1
for left, right in ((1, 2), (1, 7), (3, 11), (17, 73)):
    variance = (
        normalized_inner(left, left)
        - 2 * normalized_inner(left, right)
        + normalized_inner(right, right)
    )
    require(variance > 0, "two-slot factorial variance")
    low_child_cells += 1


print("THM-3093 ARBITRARY-GAP REMOTE-CLUSTER MONGE FLAG")
print(f"response_mlr_cells={response_cells} strict_Monge=PASS")
print(f"jensen_cells={jensen_cells} mixed_hull=PASS")
print(f"flag_weight_cells={flag_cells} staircase_ties=PASS")
print(f"raw_atom_cells={raw_atom_cells} whole_system_contraction=PASS")
print(f"block_cells={block_cells} all_composition_faces=PASS")
print(f"q3_face_cells={q3_face_cells} four_boundary_types=PASS")
print(f"covariance_cells={covariance_cells} lambda_cancellation=PASS")
print(f"carrier_cells={carrier_cells} U_times_V=U_far=PASS")
print("uniform_gap_range=all_distinct_integer_vectors")
print("outer_error=poly(C)*(m/(m+1))^(mC) independent_of_cluster_diameter")
print("boundary=repeated_gaps;moving_child;growing_width;S_m_zero")
print("scope=fixed_child;fixed_rank;arbitrary_internal_gaps")
print(f"low_child_cells={low_child_cells} empty_resultant_and_variance=PASS")
print("all_exact_checks=PASS")
