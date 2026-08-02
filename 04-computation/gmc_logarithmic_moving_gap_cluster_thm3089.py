#!/usr/bin/env python3
"""Exact companion for THM-3089's square-root moving-gap cone."""

from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations
from math import factorial, prod


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(start, length):
    answer = 1
    for offset in range(length):
        answer *= start + offset
    return answer


def weak_compositions(total, length):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, length - 1):
            yield (first,) + tail


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
            require(swap is not None, "singular determinant")
            matrix[pivot_index], matrix[swap] = matrix[swap], matrix[pivot_index]
            sign *= -1
        pivot = matrix[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    matrix[row][column] * pivot
                    - matrix[row][pivot_index] * matrix[pivot_index][column]
                )
                require(numerator % denominator == 0, "Bareiss exact division")
                matrix[row][column] = numerator // denominator
        denominator = pivot
        for row in range(pivot_index + 1, size):
            matrix[row][pivot_index] = 0
    return sign * matrix[-1][-1]


def permanent(matrix):
    size = len(matrix)
    states = {0: 1}
    for row in range(size):
        updated = {}
        for mask, value in states.items():
            for column in range(size):
                if not (mask >> column) & 1:
                    new_mask = mask | (1 << column)
                    updated[new_mask] = updated.get(new_mask, 0) + value * matrix[row][column]
        states = updated
    return states[(1 << size) - 1]


def find_positive_permutation(matrix):
    size = len(matrix)

    def extend(row, used, answer):
        if row == size:
            return tuple(answer)
        for column in range(size):
            if column not in used and matrix[row][column] > 0:
                found = extend(row + 1, used | {column}, answer + [column])
                if found is not None:
                    return found
        return None

    return extend(0, set(), [])


def decompose_regular(matrix, degree):
    remaining = [row[:] for row in matrix]
    answer = []
    for _ in range(degree):
        matching = find_positive_permutation(remaining)
        require(matching is not None, "regular bipartite matching")
        answer.append(matching)
        for row, column in enumerate(matching):
            remaining[row][column] -= 1
    require(all(entry == 0 for row in remaining for entry in row), "Birkhoff decomposition")
    return answer


def add(left, right):
    return left[0] + right[0], left[1] + right[1]


def scale(value, interval):
    if value >= 0:
        return value * interval[0], value * interval[1]
    return value * interval[1], value * interval[0]


def log_bounds(value, terms=45):
    require(value > 0, "log domain")
    if value < 1:
        lower, upper = log_bounds(1 / value, terms)
        return -upper, -lower
    power_two = 0
    reduced = value
    while reduced >= 2:
        reduced /= 2
        power_two += 1
    z = (reduced - 1) / (reduced + 1)
    partial = Fraction(0)
    zpower = z
    for index in range(terms):
        partial += 2 * zpower / (2 * index + 1)
        zpower *= z * z
    tail = 2 * zpower / ((2 * terms + 1) * (1 - z * z))
    core = partial, partial + tail
    if not power_two:
        return core
    z2 = Fraction(1, 3)
    partial2 = Fraction(0)
    zpower2 = z2
    for index in range(terms):
        partial2 += 2 * zpower2 / (2 * index + 1)
        zpower2 *= z2 * z2
    tail2 = 2 * zpower2 / ((2 * terms + 1) * (1 - z2 * z2))
    return add(core, scale(power_two, (partial2, partial2 + tail2)))


# Uniform exact factorial-ratio controls.  The normalized all-high
# coefficient for alpha is compared to its line-power model r^S.
ratio_cells = 0
for child in range(3, 6):
    for cluster_size in range(2, 4):
        parent = child + cluster_size
        for gap_tail in combinations(range(1, 6), cluster_size - 1):
            gaps = (0,) + gap_tail
            height = gaps[-1]
            for C in (50, 113):
                N = C + 2
                for degree in sorted({child + 1, parent}):
                    for alpha in weak_compositions(degree, cluster_size):
                        shift = sum(a * h for a, h in zip(alpha, gaps))
                        quotient = Fraction(rising(degree * N + 1, shift))
                        for multiplicity, gap in zip(alpha, gaps):
                            quotient /= rising(N + 1, gap) ** multiplicity
                        quotient /= degree**shift
                        lower, upper = log_bounds(quotient)
                        bound = Fraction(degree * height * (height + 1), N)
                        require(lower >= -bound and upper <= bound, "factorial distortion")
                        ratio_cells += 1


# Uniform total-positivity controls.  The permanent/determinant ratio stays
# bounded as the integer gaps move.  In this exact bank its maximum occurs at
# consecutive gaps for ranks two through five.
permanent_cells = 0
permanent_maxima = []
expected_maxima = (Fraction(9), Fraction(375), Fraction(168268, 3), Fraction(241529660, 9))
for rank, expected_maximum in zip(range(2, 6), expected_maxima):
    rows = tuple(range(4, rank + 4))
    maximum = Fraction(0)
    maximizing_gaps = None
    for gap_tail in combinations(range(1, 31), rank - 1):
        gaps = (0,) + gap_tail
        matrix = [[row**gap for gap in gaps] for row in rows]
        determinant = bareiss_determinant(matrix)
        require(determinant >= 1, "integer alternant lower bound")
        ratio = Fraction(permanent(matrix), determinant)
        vandermonde = prod(
            rows[right] - rows[left]
            for left in range(rank)
            for right in range(left + 1, rank)
        )
        explicit_schur_bound = Fraction(
            factorial(rank)
            * prod(row**index for index, row in enumerate(rows)),
            vandermonde,
        )
        require(ratio <= explicit_schur_bound, "explicit Schur condition bound")
        if ratio > maximum:
            maximum = ratio
            maximizing_gaps = gaps
        permanent_cells += 1
    require(maximum == expected_maximum, "permanent/determinant maximum")
    require(maximizing_gaps == tuple(range(rank)), "consecutive-gap maximum")
    permanent_maxima.append(maximum)


# Every nonnegative integer exponent matrix with constant row and column sum
# decomposes into permutation matrices.  Such a monomial is bounded by the
# corresponding power of the permanent.
birkhoff_cells = 0
for rank in range(2, 6):
    rows = tuple(range(4, rank + 4))
    gaps = tuple(index * (index + 1) // 2 for index in range(rank))
    matrix = [[row**gap for gap in gaps] for row in rows]
    matrix_permanent = permanent(matrix)
    permutation_bank = tuple(permutations(range(rank)))[: 3 * rank]
    for degree in range(1, 4):
        for chosen in combinations_with_replacement(range(len(permutation_bank)), degree):
            exponent_matrix = [[0 for _ in range(rank)] for _ in range(rank)]
            for index in chosen:
                for row, column in enumerate(permutation_bank[index]):
                    exponent_matrix[row][column] += 1
            require(
                all(sum(row) == degree for row in exponent_matrix),
                "regular row degree",
            )
            require(
                all(sum(exponent_matrix[row][column] for row in range(rank)) == degree for column in range(rank)),
                "regular column degree",
            )
            decomposition = decompose_regular(exponent_matrix, degree)
            require(len(decomposition) == degree, "decomposition length")
            monomial = 1
            for row in range(rank):
                for column in range(rank):
                    monomial *= matrix[row][column] ** exponent_matrix[row][column]
            require(monomial <= matrix_permanent**degree, "permanent monomial bound")
            birkhoff_cells += 1


# The auxiliary positive direct basis v and the fixed-pivot normal basis
# u=(w,z_1,...,z_(q-1)) are related by one gap-independent unimodular map:
# v_j=z_j for j<q and v_q=-w-sum z_j.
normal_transform_cells = 0
for rank in range(2, 7):
    transform = [[0 for _ in range(rank)] for _ in range(rank)]
    for index in range(rank - 1):
        transform[index][index + 1] = 1
    transform[-1][0] = -1
    for column in range(1, rank):
        transform[-1][column] = -1
    require(
        bareiss_determinant(transform) == (-1) ** rank,
        "direct-to-normal determinant",
    )
    gap_patterns = (
        tuple(range(rank)),
        tuple(index * index for index in range(rank)),
        tuple(index * (index + 1) // 2 for index in range(rank)),
    )
    for gaps in gap_patterns:
        for degree in range(4, rank + 4):
            direct_row = [degree**gap for gap in gaps]
            normal_row = [
                sum(direct_row[row] * transform[row][column] for row in range(rank))
                for column in range(rank)
            ]
            expected = [-degree ** gaps[-1]] + [
                degree**gaps[index] - degree ** gaps[-1]
                for index in range(rank - 1)
            ]
            require(normal_row == expected, "direct-to-normal row")
            normal_transform_cells += 1


# Resultant coefficient degree and the logarithmic conditioning invoice.
# Appending a degree multiplies the previous invoice, so the endpoint cone
# automatically controls every internal cluster prefix.
degree_cells = 0
for child in range(3, 9):
    previous_A_over_log = 0
    for cluster_size in range(2, 6):
        parent = child + cluster_size
        normal_degree = factorial(parent) // factorial(child)
        total_coefficient_degree = sum(normal_degree // degree for degree in range(child + 1, parent + 1))
        A_over_log = parent * total_coefficient_degree
        require(A_over_log > previous_A_over_log, "prefix conditioning monotonicity")
        previous_A_over_log = A_over_log
        degree_cells += 1


# At H comparable to sqrt(C), the exact first correction is a Gaussian
# variance multiplier.  These exact identities verify the new boundary face
# and its strict negativity whenever two distinct moving slopes are active.
sqrt_boundary_cells = 0
sqrt_boundary_negative = 0
for rank in range(2, 5):
    slopes = tuple(index * index for index in range(rank))
    for degree in range(rank, rank + 4):
        for alpha in weak_compositions(degree, rank):
            weighted_sum = sum(a * slope for a, slope in zip(alpha, slopes))
            quadratic = Fraction(weighted_sum * weighted_sum, 2 * degree) - Fraction(
                sum(a * slope * slope for a, slope in zip(alpha, slopes)), 2
            )
            variance = -Fraction(
                sum(
                    alpha[left]
                    * alpha[right]
                    * (slopes[left] - slopes[right]) ** 2
                    for left in range(rank)
                    for right in range(left + 1, rank)
                ),
                2 * degree,
            )
            require(quadratic == variance, "square-root Gaussian identity")
            active = sum(value > 0 for value in alpha)
            if active >= 2:
                require(quadratic < 0, "nontrivial Gaussian damping")
                sqrt_boundary_negative += 1
            else:
                require(quadratic == 0, "one-sheet boundary equality")
            sqrt_boundary_cells += 1


k12_child = 3
k12_parent = 12
k12_D = factorial(k12_parent) // factorial(k12_child)
k12_T = sum(k12_D // degree for degree in range(k12_child + 1, k12_parent + 1))
k12_A_over_log = k12_parent * k12_T

print("THM-3089 SQUARE-ROOT MOVING-GAP CLUSTER CONE")
print(f"factorial_ratio_cells={ratio_cells} exact_log_bounds=PASS")
print(f"permanent_ratio_cells={permanent_cells} maxima=9,375,168268/3,241529660/9")
print("explicit_schur_condition_bound=PASS")
print(f"birkhoff_condition_cells={birkhoff_cells} monomial_bound=PASS")
print(f"normal_transform_cells={normal_transform_cells} determinant_and_rows=PASS")
print(f"conditioning_degree_cells={degree_cells} prefix_monotonicity=PASS")
print("uniform_relative_error=O(H^2/C);physical_cone=H=o(sqrt(C))")
print(f"K12_m3_Tz={k12_T} A_over_logK={k12_A_over_log}")
print(
    f"sqrt_boundary_identities={sqrt_boundary_cells} "
    f"strict_gaussian_cells={sqrt_boundary_negative}"
)
print("composition_extension=same_entropy_chamber;moving_internal_gaps")
print("scope=fixed_width;distinct_gaps;no_H_comparable_to_sqrtC_claim")
print("all_exact_checks=PASS")
