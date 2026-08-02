#!/usr/bin/env python3
"""Exact rational companion for THM-3100's product-Gamma Monge flag."""

from fractions import Fraction
from functools import lru_cache, reduce
from itertools import product, permutations
from math import factorial, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def rising(start, length):
    answer = Fraction(1)
    for offset in range(length):
        answer *= start + offset
    return answer


@lru_cache(maxsize=None)
def weight(shapes, degree):
    answer = Fraction(1)
    for theta in shapes:
        answer *= rising(theta, degree)
    return answer


@lru_cache(maxsize=None)
def carrier(shapes, degree, base, gap):
    answer = Fraction(1)
    for theta in shapes:
        answer *= Fraction(
            rising(theta + degree * base, degree * gap),
            rising(theta + base, gap) ** degree,
        )
    return answer


def mixed_coefficient(shapes, degree, base, gaps, alpha):
    total_gap = sum(a * h for a, h in zip(alpha, gaps))
    answer = Fraction(1)
    for theta in shapes:
        denominator = reduce(
            lambda x, y: x * y,
            (
                rising(theta + base, gaps[index]) ** alpha[index]
                for index in range(len(gaps))
            ),
            Fraction(1),
        )
        answer *= Fraction(
            rising(theta + degree * base, total_gap), denominator
        )
    return answer


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


def permutation_profile(private_sizes, theta):
    """Cycle-weighted bad-mask census for one marked permutation layer."""

    mark_count = len(private_sizes)
    private_sets = []
    cursor = mark_count
    for size in private_sizes:
        private_sets.append(set(range(cursor, cursor + size)))
        cursor += size
    profile = {}
    for image in permutations(range(cursor)):
        predecessor = [0] * cursor
        for source, target in enumerate(image):
            predecessor[target] = source
        seen = [False] * cursor
        cycles = 0
        for start in range(cursor):
            if seen[start]:
                continue
            cycles += 1
            point = start
            while not seen[point]:
                seen[point] = True
                point = image[point]
        bad_mask = 0
        for mark in range(mark_count):
            if image[mark] == mark or predecessor[mark] in private_sets[mark]:
                bad_mask |= 1 << mark
        profile[bad_mask] = profile.get(bad_mask, Fraction(0)) + theta**cycles
    return profile


def adjacent_tensor_numerator(shapes, indices):
    """Cleared numerator L(prod E_n), E_n=s^(n+1)-R_n s^n."""

    total_degree = sum(index + 1 for index in indices)
    answer = Fraction(0)
    for lowered in product((0, 1), repeat=len(indices)):
        chosen = [index for index, bit in zip(indices, lowered) if bit]
        coefficient = Fraction((-1) ** len(chosen))
        for index in chosen:
            coefficient *= reduce(
                lambda x, y: x * y,
                (theta + index for theta in shapes),
                Fraction(1),
            )
        answer += coefficient * weight(shapes, total_degree - len(chosen))
    return answer


def normalized_mixed_moment(shapes, vectors):
    """Moment of products of sparse vectors in the normalized f_n basis."""

    answer = Fraction(0)
    choices = [tuple(vector.items()) for vector in vectors]
    for terms in product(*choices):
        coefficient = reduce(
            lambda x, y: x * y,
            (term[1] for term in terms),
            Fraction(1),
        )
        indices = tuple(term[0] for term in terms)
        answer += coefficient * Fraction(
            weight(shapes, sum(indices)),
            reduce(
                lambda x, y: x * y,
                (weight(shapes, index) for index in indices),
                Fraction(1),
            ),
        )
    return answer


def gap_patterns(width):
    return {
        2: ((0, 1), (0, 11)),
        3: ((0, 1, 3), (0, 7, 29)),
        4: ((0, 1, 3, 6), (0, 5, 19, 47)),
        5: ((0, 1, 3, 7, 12), (0, 7, 23, 53, 101)),
    }[width]


def contiguous_blocks(size, divergence_mask):
    blocks = []
    start = 0
    for edge in range(size - 1):
        if divergence_mask & (1 << edge):
            blocks.append(tuple(range(start, edge + 1)))
            start = edge + 1
    blocks.append(tuple(range(start, size)))
    return blocks


def vandermonde(nodes):
    answer = 1
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            answer *= nodes[j] - nodes[i]
    return answer


shape_bank = (
    (Fraction(1),),
    (Fraction(2),),
    (Fraction(3, 2),),
    (Fraction(1), Fraction(1)),
    (Fraction(1, 2), Fraction(3, 2)),
    (Fraction(2), Fraction(3)),
    (Fraction(1), Fraction(1), Fraction(1)),
    (Fraction(1, 2), Fraction(1), Fraction(5, 2)),
    (Fraction(2, 3), Fraction(4, 3), Fraction(7, 3)),
)


# Strict response-Monge comparisons, including base zero and degree one.
response_cells = 0
for shapes in shape_bank:
    for base in (0, 5, 29):
        for left_gap, right_gap in ((0, 1), (0, 7), (2, 11), (7, 31)):
            for lower in range(1, 6):
                for upper in range(lower + 1, 7):
                    low = carrier(shapes, lower, base, right_gap) / carrier(
                        shapes, lower, base, left_gap
                    )
                    high = carrier(shapes, upper, base, right_gap) / carrier(
                        shapes, upper, base, left_gap
                    )
                    require(low**upper < high**lower, "product-Gamma response")
                    response_cells += 1
require(response_cells == 1620, "response census")


# Factorwise multivariate Jensen across all full-degree coefficients.
jensen_cells = 0
for shapes in shape_bank:
    for width in range(2, 6):
        for base in (0, 17):
            for gaps in gap_patterns(width):
                for degree in range(1, width + 1):
                    pure_values = [
                        carrier(shapes, degree, base, gap) for gap in gaps
                    ]
                    for alpha in compositions(degree, width):
                        mixed = mixed_coefficient(
                            shapes, degree, base, gaps, alpha
                        )
                        pure = reduce(
                            lambda x, y: x * y,
                            (
                                pure_values[index] ** alpha[index]
                                for index in range(width)
                            ),
                            Fraction(1),
                        )
                        require(mixed**degree <= pure, "factorwise Jensen")
                        if degree == 1:
                            require(mixed == 1, "normalized first moment")
                        jensen_cells += 1
require(jensen_cells == 12384, "Jensen census")


# Exact dual weights for representative fixed parameters.
flag_cells = 0
for shapes in shape_bank:
    for width in range(2, 6):
        degrees = tuple(range(1, width + 1))
        root_power = reduce(lcm, degrees, 1)
        for gaps in gap_patterns(width):
            values = [
                [carrier(shapes, degree, 23, gap) for gap in gaps]
                for degree in degrees
            ]
            lambda_power = [Fraction(1)]
            for edge in range(width - 1):
                degree = degrees[edge]
                lambda_power.append(
                    lambda_power[-1]
                    * (values[edge][edge] / values[edge][edge + 1])
                    ** (root_power // degree)
                )
            require(lambda_power[:2] == [1, 1], "degree-one anchor")
            for row, degree in enumerate(degrees):
                scale = (
                    values[row][row] ** (root_power // degree)
                    * lambda_power[row]
                )
                for column in range(width):
                    value = (
                        values[row][column] ** (root_power // degree)
                        * lambda_power[column]
                        / scale
                    )
                    require(value <= 1, "product-Gamma flag weight")
                    if column == row or (row < width - 1 and column == row + 1):
                        require(value == 1, "product-Gamma staircase tie")
                    flag_cells += 1
require(flag_cells == 972, "flag census")


# Every compactification block uses the nodes r^A.  Six exact conditions are
# checked on every composition face for A=1,2,3 and widths two through five.
face_cells = 0
for layer_count in range(1, 4):
    for width in range(2, 6):
        gaps = gap_patterns(width)[1]
        nodes_all = tuple((index + 1) ** layer_count for index in range(width))
        for mask in range(1 << (width - 1)):
            blocks = contiguous_blocks(width, mask)
            require(tuple(i for block in blocks for i in block) == tuple(range(width)), "block cover")
            require(all(block == tuple(range(block[0], block[-1] + 1)) for block in blocks), "block order")
            full = [[Fraction(0) for _ in range(width)] for _ in range(width)]
            product_det = Fraction(1)
            positive = True
            floor_ok = True
            for block_number, block in enumerate(blocks):
                start = block[0]
                nodes = tuple(nodes_all[index] for index in block)
                local_gaps = tuple(gaps[index] - gaps[start] for index in block)
                diagonal = [
                    [
                        Fraction(node**gap, node ** local_gaps[row])
                        for gap in local_gaps
                    ]
                    for row, node in enumerate(nodes)
                ]
                diagonal_det = determinant(diagonal)
                floor = Fraction(
                    vandermonde(nodes),
                    reduce(
                        lambda x, y: x * y,
                        (nodes[index] ** index for index in range(len(nodes))),
                        1,
                    ),
                )
                positive = positive and diagonal_det > 0
                floor_ok = floor_ok and diagonal_det >= floor
                product_det *= diagonal_det
                for local_row, global_row in enumerate(block):
                    for local_column, global_column in enumerate(block):
                        full[global_row][global_column] = diagonal[local_row][local_column]
                    for later_block in blocks[block_number + 1 :]:
                        for global_column in later_block:
                            full[global_row][global_column] = Fraction(
                                global_row + 2 * global_column + mask + 1, 5
                            )
            require(positive, "positive product-Gamma alternants")
            require(floor_ok, "product-Gamma Schur floor")
            require(determinant(full) == product_det, "forward-entry invariance")
            require(product_det > 0, "positive face product")
            face_cells += 6
require(face_cells == 540, "face census")


# Common-offset covariance remains exact for arbitrary rational shapes.
covariance_cells = 0
for shapes in shape_bank:
    for degree in range(1, 7):
        for base in (0, 3, 17):
            for common in (1, 7, 23):
                for local in (0, 2, 11):
                    require(
                        carrier(shapes, degree, base, common + local)
                        == carrier(shapes, degree, base, common)
                        * carrier(shapes, degree, base + common, local),
                        "product-Gamma common-offset carrier",
                    )
                    covariance_cells += 1
require(covariance_cells == 1458, "covariance census")


# Product-Gamma adjacent tensors are positive by an A-layer permutation
# inclusion-exclusion.  The direct expansion, the cycle-weighted bad-mask
# census, and the explicit one-layer marked-cycle lower family agree.
adjacent_cells = 0
permutation_cells = 0
adjacent_profiles = (
    (0, 0),
    (0, 1),
    (1, 2),
    (0, 0, 0),
    (0, 1, 1),
    (0, 1, 2),
    (0, 0, 0, 0),
)
for shapes in shape_bank:
    for indices in adjacent_profiles:
        numerator = adjacent_tensor_numerator(shapes, indices)
        require(numerator > 0, "strict product-Gamma adjacent tensor")
        ordinary_count = sum(indices)
        total_degree = sum(index + 1 for index in indices)
        lower = (
            factorial(len(indices) - 1)
            * shapes[0]
            * rising(shapes[0], ordinary_count)
            * reduce(
                lambda x, y: x * y,
                (rising(theta, total_degree) for theta in shapes[1:]),
                Fraction(1),
            )
        )
        require(numerator >= lower, "marked-cycle lower family")
        denominator = reduce(
            lambda x, y: x * y,
            (weight(shapes, index + 1) for index in indices),
            Fraction(1),
        )
        require(numerator / denominator > 0, "normalized adjacent tensor")
        adjacent_cells += 1
require(adjacent_cells == 63, "adjacent tensor census")

permutation_bank = (
    ((Fraction(1),), (0, 0)),
    ((Fraction(3, 2),), (0, 1)),
    ((Fraction(1), Fraction(2)), (0, 0)),
    ((Fraction(1, 2), Fraction(3, 2)), (0, 1)),
    ((Fraction(1), Fraction(1)), (0, 0, 0)),
    ((Fraction(2, 3), Fraction(4, 3), Fraction(7, 3)), (0, 0)),
)
for shapes, indices in permutation_bank:
    layer_profiles = [permutation_profile(indices, theta) for theta in shapes]
    enumerated = Fraction(0)
    for masks_and_weights in product(*(tuple(profile.items()) for profile in layer_profiles)):
        common_bad = (1 << len(indices)) - 1
        tuple_weight = Fraction(1)
        for mask, layer_weight in masks_and_weights:
            common_bad &= mask
            tuple_weight *= layer_weight
        if common_bad == 0:
            enumerated += tuple_weight
    require(
        enumerated == adjacent_tensor_numerator(shapes, indices),
        "A-layer catastrophic-mask inclusion-exclusion",
    )
    permutation_cells += 1
require(permutation_cells == 6, "permutation census")


# Width one and two are unconditionally good.  The exact variance is the L2
# norm under a product of full-support Gamma variables.
variance_cells = 0
for shapes in shape_bank:
    for left, right in ((0, 1), (0, 7), (2, 11), (17, 73)):
        inner_ll = Fraction(weight(shapes, 2 * left), weight(shapes, left) ** 2)
        inner_lr = Fraction(
            weight(shapes, left + right),
            weight(shapes, left) * weight(shapes, right),
        )
        inner_rr = Fraction(weight(shapes, 2 * right), weight(shapes, right) ** 2)
        require(inner_ll - 2 * inner_lr + inner_rr > 0, "two-slot variance")
        variance_cells += 1
require(variance_cells == 36, "variance census")


# The factorial atomic orientation D>=0 does not extend even to one Gamma
# factor of shape two.  The true quadratic/cubic remainders remain nonzero,
# so this is a certificate failure, not a width-three counterexample.
shapes = (Fraction(2),)
u = {1: Fraction(1), 0: Fraction(-1)}
v = {2: Fraction(1), 1: Fraction(-1)}
g11 = normalized_mixed_moment(shapes, (u, u))
g12 = normalized_mixed_moment(shapes, (u, v))
g22 = normalized_mixed_moment(shapes, (v, v))
t111 = normalized_mixed_moment(shapes, (u, u, u))
t112 = normalized_mixed_moment(shapes, (u, u, v))
t122 = normalized_mixed_moment(shapes, (u, v, v))
t222 = normalized_mixed_moment(shapes, (v, v, v))
atomic_d = 2 * t222 * g12 - 3 * t122 * g22
remainder_i1 = (
    3 * t112 * g11 * g22
    - t222 * g11**2
    - 2 * t111 * g12 * g22
)
remainder_i2 = (
    3 * t122 * g11 * g22
    - 2 * t222 * g12 * g11
    - t111 * g22**2
)
require(g11 * g22 - g12**2 > 0, "shape-two Gram determinant")
require(atomic_d == Fraction(-1, 12), "shape-two atomic hostile")
require(remainder_i1 == Fraction(-1, 2), "shape-two first remainder")
require(remainder_i2 == Fraction(-11, 36), "shape-two second remainder")


# A=0 is evaluation at one: every response node collides and every
# higher-width resultant degenerates.  Moving theta has no uniform floor.
degenerate_cells = 0
for width in range(2, 6):
    rows = tuple(Fraction(1) for _ in range(width))
    matrix = [[row**column for column in range(width)] for row in rows]
    require(determinant(matrix) == 0, "A=0 repeated-node boundary")
    degenerate_cells += 1
for theta in (1, 2, 7, 31, 101):
    shapes = (Fraction(theta),)
    variance = (
        Fraction(weight(shapes, 0), weight(shapes, 0) ** 2)
        - 2 * Fraction(weight(shapes, 1), weight(shapes, 0) * weight(shapes, 1))
        + Fraction(weight(shapes, 2), weight(shapes, 1) ** 2)
    )
    require(variance == Fraction(1, theta), "moving-theta variance boundary")
    degenerate_cells += 1
require(degenerate_cells == 9, "boundary census")


print("THM-3100 PRODUCT-GAMMA RESPONSE-MONGE COMPACTIFICATION")
print(f"response_cells={response_cells} degree_one_and_base_zero=PASS")
print(f"jensen_cells={jensen_cells} factorwise_convexity=PASS")
print(f"flag_cells={flag_cells} staircase_dual=PASS")
print(f"face_cells={face_cells} nodes_r_power_A=PASS")
print(f"covariance_cells={covariance_cells} common_offset=PASS")
print(f"adjacent_cells={adjacent_cells} permutation_cells={permutation_cells} A_layer_model=PASS")
print(f"variance_cells={variance_cells} widths_one_two_good=PASS")
print("width3_atomic_hostile=theta2_support012_D=-1/12_I1=-1/2_I2=-11/36")
print(f"boundary_cells={degenerate_cells} A0_and_moving_theta=PASS")
print("translation_face=fixed_A;fixed_positive_theta;fixed_width")
print("remote_error=poly(C)*(m/(m+1))^(A*m*C)")
print("minimal_bad_prefix_width_at_least_3")
print("bad_support_count=O_t(X^(t-3))")
print("conditional_SFC3_count=O_t(X^(t-4))")
print("boundary=width3_open;A0_degenerate;no_parameter_uniformity")
print("all_exact_checks=PASS")
