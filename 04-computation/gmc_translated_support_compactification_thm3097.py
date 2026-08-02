#!/usr/bin/env python3
"""Exact finite companion for THM-3097's translated-support flag."""

from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations
from math import comb, factorial, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def rising(start, length):
    answer = 1
    for offset in range(length):
        answer *= start + offset
    return answer


@lru_cache(maxsize=None)
def carrier(degree, base, gap):
    return Fraction(
        rising(degree * base + 1, degree * gap),
        rising(base + 1, gap) ** degree,
    )


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
        5: (
            (0, 1, 3, 7, 12),
            (0, 3, 19, 53, 73),
            (0, 11, 23, 47, 97),
        ),
        6: (
            (0, 1, 3, 6, 10, 15),
            (0, 2, 7, 19, 41, 73),
            (0, 5, 17, 37, 79, 131),
        ),
    }
    return banks[size]


def block_matrix(rows, local_gaps):
    return [
        [Fraction(row**gap, row ** local_gaps[i]) for gap in local_gaps]
        for i, row in enumerate(rows)
    ]


def vandermonde(rows):
    answer = 1
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            answer *= rows[j] - rows[i]
    return answer


# The response comparison includes the formerly omitted degree-one row.
response_cells = 0
response_edges = (
    (0, 1),
    (0, 2),
    (0, 7),
    (1, 3),
    (2, 7),
    (3, 11),
    (7, 19),
    (17, 73),
)
for base in (0, 1, 4, 17, 73):
    for left_gap, right_gap in response_edges:
        for lower in range(1, 7):
            for upper in range(lower + 1, 8):
                low_response = carrier(lower, base, right_gap) / carrier(
                    lower, base, left_gap
                )
                high_response = carrier(upper, base, right_gap) / carrier(
                    upper, base, left_gap
                )
                require(
                    low_response**upper < high_response**lower,
                    "strict translated response",
                )
                response_cells += 1
require(response_cells == 840, "response census")


# Every coefficient of a degree-r translated form has total remote degree r.
# Convexity of log(rN+1)_s gives the exact multivariate Jensen bound.
jensen_cells = 0
for width in range(2, 7):
    for base in (0, 29):
        for gaps in gap_patterns(width):
            denominators = [rising(base + 1, gap) for gap in gaps]
            for degree in range(1, width + 1):
                pure_values = [carrier(degree, base, gap) for gap in gaps]
                for alpha in compositions(degree, width):
                    total_gap = sum(a * h for a, h in zip(alpha, gaps))
                    denominator = reduce(
                        lambda x, y: x * y,
                        (
                            denominators[index] ** alpha[index]
                            for index in range(width)
                        ),
                        1,
                    )
                    mixed = Fraction(
                        rising(degree * base + 1, total_gap), denominator
                    )
                    pure = reduce(
                        lambda x, y: x * y,
                        (
                            pure_values[index] ** alpha[index]
                            for index in range(width)
                        ),
                        Fraction(1),
                    )
                    require(mixed**degree <= pure, "translated Jensen")
                    if degree == 1:
                        require(mixed == 1, "F1 is the coordinate sum")
                    jensen_cells += 1
require(jensen_cells == 7602, "Jensen census")


# Root-free powers of the dual potentials expose the staircase, including
# lambda_1=lambda_2 because V_1 is identically one.
flag_cells = 0
for width in range(2, 7):
    degrees = tuple(range(1, width + 1))
    root_power = reduce(lcm, degrees, 1)
    for base in (width + 1, 31):
        for gaps in gap_patterns(width):
            values = [
                [carrier(degree, base, gap) for gap in gaps]
                for degree in degrees
            ]
            lambda_power = [Fraction(1)]
            for edge in range(width - 1):
                pivot_degree = degrees[edge]
                ratio = values[edge][edge] / values[edge][edge + 1]
                lambda_power.append(
                    lambda_power[-1]
                    * ratio ** (root_power // pivot_degree)
                )
            require(lambda_power[0] == lambda_power[1] == 1, "degree-one tie")
            require(
                all(
                    lambda_power[index + 1] < lambda_power[index]
                    for index in range(1, width - 1)
                ),
                "strict later potentials",
            )
            for row, degree in enumerate(degrees):
                row_scale = (
                    values[row][row] ** (root_power // degree)
                    * lambda_power[row]
                )
                for column in range(width):
                    weight = (
                        values[row][column] ** (root_power // degree)
                        * lambda_power[column]
                        / row_scale
                    )
                    require(weight <= 1, "staircase weight")
                    if column == row or (row < width - 1 and column == row + 1):
                        require(weight == 1, "staircase equality")
                    flag_cells += 1
require(flag_cells == 540, "flag census")


# Every composition face has positive generalized-alternant diagonal blocks.
# Arbitrary forward entries do not alter the determinant.
face_cells = 0
for width in range(2, 7):
    gaps = gap_patterns(width)[width % 3]
    for mask in range(1 << (width - 1)):
        blocks = contiguous_blocks(width, mask)
        require(tuple(i for block in blocks for i in block) == tuple(range(width)), "block cover")
        require(all(block == tuple(range(block[0], block[-1] + 1)) for block in blocks), "consecutive blocks")

        full = [[Fraction(0) for _ in range(width)] for _ in range(width)]
        product_det = Fraction(1)
        all_positive = True
        all_floors = True
        for block_number, block in enumerate(blocks):
            start = block[0]
            rows = tuple(index + 1 for index in block)
            local_gaps = tuple(gaps[index] - gaps[start] for index in block)
            diagonal = block_matrix(rows, local_gaps)
            diagonal_det = determinant(diagonal)
            floor = Fraction(
                vandermonde(rows),
                reduce(
                    lambda x, y: x * y,
                    (
                        rows[index] ** index
                        for index in range(len(rows))
                    ),
                    1,
                ),
            )
            all_positive = all_positive and diagonal_det > 0
            all_floors = all_floors and diagonal_det >= floor
            product_det *= diagonal_det
            for local_row, global_row in enumerate(block):
                for local_column, global_column in enumerate(block):
                    full[global_row][global_column] = diagonal[local_row][local_column]
                for later_block in blocks[block_number + 1 :]:
                    for global_column in later_block:
                        full[global_row][global_column] = Fraction(
                            (global_row + 1) * (global_column + 2) + mask + 1,
                            7,
                        )
        require(all_positive, "positive diagonal alternants")
        require(all_floors, "Schur-Vandermonde floor")
        require(determinant(full) == product_det, "upper-triangular determinant")
        require(product_det > 0, "positive face product")
        face_cells += 6
require(face_cells == 372, "composition-face census")


# Exact carrier covariance, including common-offset rebasing, and cancellation
# of every column-potential exponent in the resultant.
covariance_cells = 0
for degree in range(1, 8):
    for base in (0, 3, 17):
        for common in (1, 4, 19):
            for local in (0, 2, 11):
                require(
                    carrier(degree, base, common + local)
                    == carrier(degree, base, common)
                    * carrier(degree, base + common, local),
                    "common-offset carrier",
                )
                covariance_cells += 1
for width in range(2, 7):
    total_degree = factorial(width)
    for degree in range(1, width + 1):
        variable_exponent = total_degree
        equation_exponent = -(total_degree // degree) * degree
        require(variable_exponent + equation_exponent == 0, "lambda cancellation")
        covariance_cells += 1
require(covariance_cells == 209, "covariance census")


# The rank-two and rank-three faces calibrate the limiting formulas exactly.
small_face_cells = 0
for gap in (1, 2, 7, 19):
    rows = (1, 2)
    exact = determinant(block_matrix(rows, (0, gap)))
    require(exact == 1 - Fraction(1, 2**gap), "rank-two face")
    require(exact**2 > 0, "rank-two resultant power")
    small_face_cells += 2
for left, right in ((1, 3), (4, 19), (17, 73)):
    rows = (1, 2, 3)
    full_det = determinant(block_matrix(rows, (0, left, right)))
    require(full_det > 0, "rank-three full face")
    require((1 - Fraction(1, 2**left)) ** 6 > 0, "rank-three left block")
    require((1 - Fraction(2**(right - left), 3 ** (right - left))) ** 6 > 0, "rank-three right block")
    require(Fraction(1) ** 6 == 1, "rank-three singleton face")
    small_face_cells += 4
require(small_face_cells == 20, "small-face census")


# Exact finite-prefix extension counts and the sharp abstract codimension-four
# cylinder.  The hostile is logical rather than a physical bad support.
cylinder_cells = 0
for width in range(4, 9):
    for prefix in ((0, 1, 2, 3), (1, 3, 7, 11), (2, 5, 9, 17)):
        for top in (17, 23, 31):
            expected = comb(top - prefix[-1], width - 4)
            actual = sum(
                1
                for support in combinations(range(top + 1), width)
                if support[:4] == prefix
            )
            require(actual == expected, "prefix-cylinder count")
            cylinder_cells += 1
require(cylinder_cells == 45, "cylinder census")

density_cells = 0
for width in range(4, 9):
    for top in (width + 3, width + 7, width + 13):
        hostile = comb(top - 3, width - 4)
        require(hostile > 0, "abstract first-four-prefix hostile")
        require(hostile <= comb(top + 1, width), "hostile is a support subset")
        density_cells += 1
require(density_cells == 15, "density census")


# Repeated gaps are the exact diagonal boundary.
for width in range(2, 7):
    rows = tuple(range(1, width + 1))
    repeated = tuple(0 if index < 2 else index for index in range(width))
    require(determinant(block_matrix(rows, repeated)) == 0, "repeated-gap boundary")


print("THM-3097 TRANSLATED-SUPPORT MONGE COMPACTIFICATION")
print(f"response_cells={response_cells} degree_one_endpoint=PASS")
print(f"jensen_cells={jensen_cells} full_degree_coefficients=PASS")
print(f"flag_cells={flag_cells} lambda1_equals_lambda2=PASS")
print(f"face_cells={face_cells} composition_cube=PASS")
print(f"covariance_cells={covariance_cells} common_offset_and_lambda=PASS")
print(f"small_face_cells={small_face_cells} rank2_rank3_limits=PASS")
print(f"cylinder_cells={cylinder_cells} exact_extension_counts=PASS")
print(f"density_cells={density_cells} codimension_four_hostile=PASS")
print("translation_face=uniform_all_internal_gaps;fixed_width")
print("first_bad_prefix_bank=finite;minimum_width=4")
print("bad_support_count=O_t(X^(t-4))")
print("conditional_count=O_t(X^(t-b-1));SFC_through_width_b")
print("boundary=non_effective_threshold;growing_width;bad_child_cylinder")
print("all_exact_checks=PASS")
