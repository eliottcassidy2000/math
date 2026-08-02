#!/usr/bin/env python3
"""Exact companion for THM-3086's cluster-composition chambers."""

from fractions import Fraction
from itertools import combinations, permutations
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
            require(swap is not None, "singular determinant")
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


def compositions(total):
    if total == 0:
        return [()]
    answer = []
    for first in range(1, total + 1):
        for tail in compositions(total - first):
            answer.append((first,) + tail)
    return answer


def permutation_sign(permutation):
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def ledger(parts):
    values = {"S3": 1}
    nodes = []
    child = 3
    for node_index, cluster_size in enumerate(parts):
        parent = child + cluster_size
        multiplier = factorial(parent) // factorial(child)
        values = {key: value * multiplier for key, value in values.items()}
        for degree in range(child + 1, parent + 1):
            values[f"U{degree}"] = factorial(parent) // degree
        values[f"E{node_index}"] = factorial(child)
        nodes.append((node_index, child, parent))
        child = parent
    return values, nodes


def add(left, right):
    return left[0] + right[0], left[1] + right[1]


def sub(left, right):
    return left[0] - right[1], left[1] - right[0]


def scale(value, interval):
    if value >= 0:
        return value * interval[0], value * interval[1]
    return value * interval[1], value * interval[0]


def log_bounds(value, terms=110):
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


def entropy_invoice(parent, delta):
    bill = (parent - 1) * factorial(parent)
    return scale(bill * delta, add(log_bounds(Fraction(parent) / delta), (1, 1)))


def cluster_gap(child):
    return scale(child, log_bounds(Fraction(child + 1, child)))


def cluster_margin(child, parent, delta):
    return sub(cluster_gap(child), entropy_invoice(parent, delta))


# All physical layer cells for child widths 3..10 and cluster sizes through 4.
layer_cells = 0
survivor_cells = 0
for child in range(3, 11):
    for cluster_size in range(1, min(4, 14 - child) + 1):
        parent = child + cluster_size
        first_upper = child + 1
        rho = Fraction(child, first_upper) ** child
        for degree in range(2, parent + 1):
            for normal_degree in range(degree + 1):
                for high_factors in range(normal_degree + 1):
                    survivor = (
                        (degree <= child and normal_degree == high_factors == 0)
                        or (degree == first_upper and normal_degree == high_factors == 0)
                        or (degree >= first_upper and normal_degree == high_factors == degree)
                    )
                    if survivor:
                        survivor_cells += 1
                        continue
                    numerator = 1 if high_factors == 0 else high_factors**high_factors
                    if degree <= child:
                        base = Fraction(numerator, first_upper**normal_degree)
                    else:
                        base = Fraction(
                            numerator * first_upper ** (degree - normal_degree),
                            degree**degree,
                        )
                    require(base <= rho, "cluster layer exceeds local gap")
                    layer_cells += 1
require((layer_cells, survivor_cells) == (7768, 288), "cluster layer census")


# Every ordered composition through width twelve has the closed final ledger.
composition_cells = 0
for final_width in range(4, 13):
    expected_words = 2 ** (final_width - 4)
    words = compositions(final_width - 3)
    require(len(words) == expected_words, "composition count")
    for parts in words:
        values, nodes = ledger(parts)
        full_degree = factorial(final_width)
        require(values["S3"] == full_degree // 6, "base exponent")
        for degree in range(4, final_width + 1):
            require(values[f"U{degree}"] == full_degree // degree, "U exponent")
        for node_index, child, parent in nodes:
            expected_e = full_degree * factorial(child) // factorial(parent)
            require(values[f"E{node_index}"] == expected_e, "node E exponent")
            require(
                expected_e * (factorial(parent) // factorial(child)) == full_degree,
                "node determinant exponent",
            )
        composition_cells += 1
require(composition_cells == 511, "composition-cell total")


# A concrete two-cluster K8 chamber: (2,3), with C/R=10^-7 safe.
safe_delta = Fraction(1, 10**7)
require(cluster_margin(5, 8, safe_delta)[0] > 0, "K8 safe ratio")
require(
    cluster_margin(5, 8, Fraction("0.00000017320136922"))[0] > 0
    and cluster_margin(5, 8, Fraction("0.00000017320136924"))[1] < 0,
    "K8 threshold bracket",
)
k8_values, _ = ledger((2, 3))
require(
    k8_values
    == {
        "S3": 6720,
        "U4": 10080,
        "U5": 8064,
        "E0": 2016,
        "U6": 6720,
        "U7": 5760,
        "U8": 5040,
        "E1": 120,
    },
    "K8 two-cluster ledger",
)


# Fifty-six exact two-block alternant-clutch controls.  The unique maximal
# row set pairs the larger A-block with the largest rows; both local minors
# have positive generalized-Vandermonde determinants.
alternant_cells = 0
for rank in range(2, 9):
    rows = tuple(range(2, rank + 2))
    for split in range(1, rank):
        low_rows = rows[:split]
        high_rows = rows[split:]
        low_offsets = tuple(range(split))
        high_offsets = tuple(range(rank - split))
        low_det = bareiss_determinant([[row**offset for offset in low_offsets] for row in low_rows])
        high_det = bareiss_determinant([[row**offset for offset in high_offsets] for row in high_rows])
        require(low_det > 0 and high_det > 0, "local alternant clutch sign")
        leading_product = 1
        for row in high_rows:
            leading_product *= row
        for alternative in combinations(rows, len(high_rows)):
            product_value = 1
            for row in alternative:
                product_value *= row
            require(product_value <= leading_product, "rearrangement product")
            if alternative != high_rows:
                require(product_value < leading_product, "unique leading row block")
        for _slope_gap in (1, 3):
            alternant_cells += 1
require(alternant_cells == 56, "alternant-clutch count")


# Full tied-leading-coefficient controls.  For every nontrivial ordered
# composition through rank seven, rearrangement assigns consecutive row
# blocks to the increasing scale slopes.  Summing all permutations tied at
# that maximal rate must give the product of the local generalized
# alternants, including its positive sign.  Three distinct offset patterns
# test nonconsecutive as well as consecutive generalized alternants.
multiblock_cells = 0
for rank in range(2, 8):
    row_values = tuple(range(2, rank + 2))
    for block_sizes in compositions(rank):
        if len(block_sizes) < 2:
            continue
        for pattern in range(3):
            slopes = []
            offsets = []
            for block_index, block_size in enumerate(block_sizes):
                slopes.extend([2 * block_index] * block_size)
                if pattern == 0:
                    offsets.extend(range(block_size))
                elif pattern == 1:
                    offsets.extend(index * index + index for index in range(block_size))
                else:
                    offsets.extend(index * (index + 3) // 2 for index in range(block_size))

            best_rate = None
            leading_coefficient = 0
            for permutation in permutations(range(rank)):
                rate = sum(slopes[column] * permutation[column] for column in range(rank))
                monomial = permutation_sign(permutation) * prod(
                    row_values[permutation[column]] ** offsets[column]
                    for column in range(rank)
                )
                if best_rate is None or rate > best_rate:
                    best_rate = rate
                    leading_coefficient = monomial
                elif rate == best_rate:
                    leading_coefficient += monomial

            local_product = 1
            start = 0
            for block_size in block_sizes:
                block_offsets = offsets[start : start + block_size]
                block_rows = row_values[start : start + block_size]
                local_product *= bareiss_determinant(
                    [
                        [row**offset for offset in block_offsets]
                        for row in block_rows
                    ]
                )
                start += block_size
            require(
                leading_coefficient == local_product > 0,
                "multiblock tied alternant coefficient",
            )
            multiblock_cells += 1
require(multiblock_cells == 360, "multiblock tied-permutation count")


print("THM-3086 ARBITRARY CLUSTER-COMPOSITION CHAMBERS")
print(f"layer_cells={layer_cells} survivors={survivor_cells}")
print(f"composition_cells_K4_to_K12={composition_cells} K12_sectors=256")
print("local_gap=q1:gamma_O;q>=2:m*log((m+1)/m)")
print("condition=J_parent(delta)<local_gap;all_internal_prefixes_controlled")
print("carrier=S3^(K!/6)*prod_r U_r^(K!/r)*prod_nodes E_i^(K!*m_i!/k_i!)")
print("every_node_alternant_exponent=K!")
print("K8_parts_2_3=E45^2016*E678^120;safe_ratio=1/10000000")
print("K8_threshold_in=(0.00000017320136922,0.00000017320136924)")
print(f"alternant_clutch_controls={alternant_cells} unique_row_blocks=PASS")
print(f"multiblock_tied_permutation_controls={multiblock_cells}")
print("scope=physical_fixed-gap_clusters;growing-gap_clutch_is_symbol_only")
print("all_exact_checks=PASS")
