#!/usr/bin/env python3
"""Exact companion for THM-3085's multi-normal fixed-gap cluster tail."""

from fractions import Fraction
from itertools import combinations
from math import factorial, prod


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def bareiss_determinant(matrix):
    matrix = [list(map(int, row)) for row in matrix]
    size = len(matrix)
    if size == 1:
        return matrix[0][0]
    sign = 1
    denominator = 1
    for pivot_index in range(size - 1):
        if matrix[pivot_index][pivot_index] == 0:
            swap = next(
                (row for row in range(pivot_index + 1, size) if matrix[row][pivot_index]),
                None,
            )
            require(swap is not None, "singular determinant control")
            matrix[pivot_index], matrix[swap] = matrix[swap], matrix[pivot_index]
            sign *= -1
        pivot = matrix[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = matrix[row][column] * pivot - matrix[row][pivot_index] * matrix[pivot_index][column]
                require(numerator % denominator == 0, "Bareiss exact division")
                matrix[row][column] = numerator // denominator
        denominator = pivot
        for row in range(pivot_index + 1, size):
            matrix[row][pivot_index] = 0
    return sign * matrix[-1][-1]


def power_determinant(rows, gaps):
    return bareiss_determinant([[row**gap for gap in gaps] for row in rows])


def normal_determinant(rows, gaps):
    last = gaps[-1]
    return bareiss_determinant(
        [[-(row**last)] + [row**gap - row**last for gap in gaps[:-1]] for row in rows]
    )


# The complete physical layer bank through width fourteen.  The 5,296 count
# excludes exactly the declared grouped survivors.
layer_cells = 0
survivor_cells = 0
equality_cells = 0
for child_width in range(3, 11):
    first_upper = child_width + 1
    rho = Fraction(child_width, first_upper) ** child_width
    for degree in range(2, 15):
        for normal_degree in range(degree + 1):
            for high_factors in range(normal_degree + 1):
                lower_survivor = degree <= child_width and normal_degree == high_factors == 0
                upper_child_survivor = (
                    degree == first_upper and normal_degree == high_factors == 0
                )
                upper_normal_survivor = (
                    degree >= first_upper and normal_degree == high_factors == degree
                )
                if lower_survivor or upper_child_survivor or upper_normal_survivor:
                    survivor_cells += 1
                    continue
                if degree <= child_width:
                    numerator = 1 if high_factors == 0 else high_factors**high_factors
                    base = Fraction(numerator, first_upper**normal_degree)
                else:
                    numerator = 1 if high_factors == 0 else high_factors**high_factors
                    base = Fraction(
                        numerator * first_upper ** (degree - normal_degree),
                        degree**degree,
                    )
                require(base <= rho, "layer base exceeds universal gap")
                equality_cells += int(base == rho)
                layer_cells += 1

require(layer_cells == 5296, "non-survivor layer count")
require(survivor_cells == 112, "survivor layer count")
require(equality_cells == 16, "sharp-gap equality count")


# Generalized Vandermonde signs for every base-three cluster with two through
# six normal directions and gaps chosen from {0,...,9} with zero anchored.
vandermonde_cells = 0
for normal_rank in range(2, 7):
    rows = tuple(range(4, 4 + normal_rank))
    for tail in combinations(range(1, 10), normal_rank - 1):
        gaps = (0,) + tail
        power_det = power_determinant(rows, gaps)
        normal_det = normal_determinant(rows, gaps)
        require(power_det > 0, "generalized Vandermonde positivity")
        require(normal_det == (-1) ** normal_rank * power_det, "normal determinant sign")
        vandermonde_cells += 1
require(vandermonde_cells == 381, "Vandermonde control count")


# Conditional general-node covariance and final all-width exponent ledgers.
covariance_cells = 0
for child_width in range(3, 11):
    for final_width in range(child_width + 1, 15):
        normal_rank = final_width - child_width
        first_upper = child_width + 1
        full_degree = factorial(final_width)
        variable_power = -Fraction(normal_rank * full_degree, first_upper)
        equation_power = sum(
            Fraction(degree, first_upper) * Fraction(full_degree, degree)
            for degree in range(first_upper, final_width + 1)
        )
        require(variable_power + equation_power == 0, "U_p covariance cancellation")
        require(
            prod(range(first_upper, final_width + 1))
            == full_degree // factorial(child_width),
            "normal degree product",
        )
        require(factorial(child_width) % 2 == 0, "final normal sign exponent")
        covariance_cells += 1

width_cells = 0
for final_width in range(4, 15):
    full_degree = factorial(final_width)
    require(full_degree // 6 == prod(range(4, final_width + 1)), "base exponent")
    require(6 * (full_degree // 6) == full_degree, "gap-symbol exponent")
    for degree in range(4, final_width + 1):
        require(full_degree % degree == 0, "U multidegree")
    width_cells += 1


# The clean six-slot consecutive-gap control.
rows_456 = (4, 5, 6)
gaps_012 = (0, 1, 2)
det_456 = power_determinant(rows_456, gaps_012)
require(det_456 == 2, "K6 Vandermonde")
require(prod(rows_456) == 120, "K6 exponential base")
k6_power = Fraction(factorial(6), 2) * sum(
    Fraction(degree - 1, degree) for degree in rows_456
)
require(k6_power == 858, "K6 inverse-power exponent")


print("THM-3085 MULTI-NORMAL FIXED-GAP CLUSTER")
print(f"layer_non_survivors={layer_cells} survivors={survivor_cells} sharp_equalities={equality_cells}")
print("general_gap=rho_m=(m/(m+1))^m;base3_rho=27/64")
print(f"generalized_vandermonde_cells={vandermonde_cells} det_power_positive=PASS")
print(f"covariance_cells={covariance_cells} widths=4..14:{width_cells}")
print("carrier=Sm^(K!/m!)*E^(m!)*prod_r U_r^(K!/r)")
print("base3= S3^(K!/6)*E_4..K^6*prod_r U_r^(K!/r)")
print("normal_symbol=det[r^h]^(K!);distinct_gaps=positive")
print("K6_gaps_012=det:2;base:120;power:858;symbol:2^720")
print("boundaries=repeated_gaps;zero_child;moving_gaps;arbitrary_support")
print("all_exact_checks=PASS")
