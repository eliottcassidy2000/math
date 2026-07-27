#!/usr/bin/env python3
"""Dependency-free exact companion for THM-2590.

The input is the thirteen septimal slice rows printed in proved THM-2585.
All arithmetic is over F_13.  The THM-2586 part is deliberately only the
finite coefficient-incidence splice: the script never identifies its
sigma={b} physical carrier with THM-2585's sigma={a} carrier.
"""

from collections import Counter, defaultdict
from itertools import combinations
from math import comb


P = 13
CHECKS = 0


def check(condition, message="exact check failed"):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


# Rows are indexed by target section q; columns are ell=0,...,6.
Y_ROWS = (
    (0, 9, 9, 0, 0, 4, 4),
    (5, 9, 7, 11, 11, 4, 9),
    (3, 11, 7, 10, 10, 11, 1),
    (9, 11, 4, 2, 10, 5, 10),
    (11, 5, 4, 10, 12, 11, 8),
    (11, 1, 11, 2, 8, 1, 1),
    (6, 9, 12, 8, 4, 4, 7),
    (7, 6, 9, 9, 5, 1, 4),
    (2, 12, 12, 5, 11, 2, 12),
    (2, 5, 2, 1, 3, 9, 8),
    (4, 3, 8, 3, 11, 9, 2),
    (10, 12, 2, 3, 3, 6, 2),
    (8, 4, 9, 2, 2, 6, 4),
)


def reduce_phi7(coefficients):
    """Reduce degree <=6 modulo Phi_7=1+z+...+z^6."""
    check(len(coefficients) == 7)
    return tuple((coefficients[i] - coefficients[6]) % P for i in range(6))


def add_vectors(left, right):
    return tuple((x + y) % P for x, y in zip(left, right))


def determinant_mod(matrix):
    matrix = [[entry % P for entry in row] for row in matrix]
    size = len(matrix)
    check(all(len(row) == size for row in matrix))
    determinant = 1
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if matrix[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            determinant = -determinant
        value = matrix[column][column] % P
        determinant = determinant * value % P
        inverse = pow(value, -1, P)
        for row in range(column + 1, size):
            factor = matrix[row][column] * inverse % P
            if factor:
                for j in range(column, size):
                    matrix[row][j] = (
                        matrix[row][j] - factor * matrix[column][j]
                    ) % P
    return determinant % P


def rref_mod(matrix):
    matrix = [[entry % P for entry in row] for row in matrix]
    row = 0
    pivots = []
    for column in range(len(matrix[0])):
        pivot = next(
            (candidate for candidate in range(row, len(matrix))
             if matrix[candidate][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        inverse = pow(matrix[row][column], -1, P)
        matrix[row] = [(entry * inverse) % P for entry in matrix[row]]
        for other in range(len(matrix)):
            if other == row:
                continue
            factor = matrix[other][column]
            if factor:
                matrix[other] = [
                    (matrix[other][j] - factor * matrix[row][j]) % P
                    for j in range(len(matrix[0]))
                ]
        pivots.append(column)
        row += 1
        if row == len(matrix):
            break
    return tuple(tuple(row) for row in matrix), tuple(pivots)


def multiply_phi7(left, right):
    """Multiply in F_13[z]/(Phi_7), using the six-term power basis."""
    coefficients = [0] * 11
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            coefficients[i + j] = (coefficients[i + j] + x * y) % P
    for degree in range(10, 5, -1):
        leading = coefficients[degree] % P
        coefficients[degree] = 0
        for shift in range(1, 7):
            coefficients[degree - shift] = (
                coefficients[degree - shift] - leading
            ) % P
    return tuple(coefficients[:6])


def polynomial_product(left, right):
    coefficients = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            coefficients[i + j] = (coefficients[i + j] + x * y) % P
    return tuple(coefficients)


def multiplication_norm(poly):
    columns = []
    for column in range(6):
        basis = tuple(int(row == column) for row in range(6))
        columns.append(multiply_phi7(poly, basis))
    matrix = [[columns[column][row] for column in range(6)] for row in range(6)]
    return determinant_mod(matrix)


def divisible_by_quadratic(poly, middle):
    """Test divisibility by z^2+middle*z+1 using exact Horner reduction."""
    constant = 0
    linear = 0
    for coefficient in reversed(poly):
        # (constant+linear*z)z = -linear+(constant-middle*linear)z.
        constant, linear = (
            (coefficient - linear) % P,
            (constant - middle * linear) % P,
        )
    return constant == 0 and linear == 0


def factor_profile(poly):
    return tuple(
        middle for middle in (3, 5, 6)
        if divisible_by_quadratic(poly, middle)
    )


def zeta_power(exponent):
    exponent %= 7
    if exponent < 6:
        return tuple(int(index == exponent) for index in range(6))
    return (P - 1,) * 6


def owner_substitution(poly, kappa):
    """Apply the R_7 automorphism z -> z^kappa."""
    output = (0,) * 6
    for exponent, coefficient in enumerate(poly):
        term = zeta_power(kappa * exponent)
        output = tuple(
            (output[index] + coefficient * term[index]) % P
            for index in range(6)
        )
    return output


def support(mask):
    return tuple(index for index in range(13) if (mask >> index) & 1)


def subset_sum(mask, columns):
    total = (0,) * 6
    for q, column in enumerate(columns):
        if (mask >> q) & 1:
            total = add_vectors(total, column)
    return total


# The three factors are checked before and after quotienting.
full_factor_3 = (1, 3, 1)
full_factor_5 = (1, 5, 1)
full_factor_6 = (1, 6, 1)
check(
    polynomial_product(polynomial_product(full_factor_3, full_factor_5), full_factor_6)
    == (1,) * 7,
    "quadratic factorization of Phi_7 failed",
)
for middle in (3, 5, 6):
    check(
        all((root * root + middle * root + 1) % P for root in range(P)),
        "quadratic factor is not irreducible",
    )
factor_3 = reduce_phi7((1, 3, 1, 0, 0, 0, 0))
factor_5 = reduce_phi7((1, 5, 1, 0, 0, 0, 0))
factor_6 = reduce_phi7((1, 6, 1, 0, 0, 0, 0))
phi_product = multiply_phi7(multiply_phi7(factor_3, factor_5), factor_6)
check(phi_product == (0,) * 6, "quadratics do not multiply to Phi_7")

# In F_13[u]/(u^12), zeta_13=1+u and Omega=sum_m m*zeta_13^m=-u^11.
omega_u = tuple(
    sum(m * comb(m, degree) for m in range(max(1, degree), 13)) % P
    for degree in range(12)
)
check(omega_u == (0,) * 11 + (P - 1,), "Omega is not the socle generator")


# Linear section map L:F_13^13 -> R_7.
COLUMNS = tuple(reduce_phi7(row) for row in Y_ROWS)
for kappa in range(1, 7):
    inverse = pow(kappa, -1, 7)
    for coordinate in range(6):
        basis = tuple(int(index == coordinate) for index in range(6))
        check(owner_substitution(owner_substitution(basis, kappa), inverse) == basis)
LINEAR_MATRIX = tuple(
    tuple(COLUMNS[q][coordinate] for q in range(13))
    for coordinate in range(6)
)
RREF, PIVOTS = rref_mod(LINEAR_MATRIX)
EXPECTED_RREF = (
    (1, 0, 0, 0, 0, 0, 2, 0, 11, 3, 2, 1, 9),
    (0, 1, 0, 0, 0, 0, 9, 5, 2, 0, 6, 7, 4),
    (0, 0, 1, 0, 0, 0, 0, 12, 12, 2, 0, 7, 2),
    (0, 0, 0, 1, 0, 0, 11, 2, 12, 11, 6, 11, 6),
    (0, 0, 0, 0, 1, 0, 2, 1, 6, 12, 5, 8, 2),
    (0, 0, 0, 0, 0, 1, 10, 5, 9, 12, 9, 0, 0),
)
check(PIVOTS == (0, 1, 2, 3, 4, 5))
check(RREF == EXPECTED_RREF)

BASIS_DETERMINANT = determinant_mod([COLUMNS[q] for q in range(6)])
check(BASIS_DETERMINANT == 11)

KERNEL_BASIS = (
    (11, 4, 0, 2, 11, 3, 1, 0, 0, 0, 0, 0, 0),
    (0, 8, 1, 11, 12, 8, 0, 1, 0, 0, 0, 0, 0),
    (2, 11, 1, 1, 7, 4, 0, 0, 1, 0, 0, 0, 0),
    (10, 0, 11, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0),
    (11, 7, 0, 7, 8, 4, 0, 0, 0, 0, 1, 0, 0),
    (12, 6, 6, 2, 5, 0, 0, 0, 0, 0, 0, 1, 0),
    (4, 9, 11, 7, 11, 0, 0, 0, 0, 0, 0, 0, 1),
)
for vector in KERNEL_BASIS:
    image = tuple(
        sum(LINEAR_MATRIX[row][q] * vector[q] for q in range(13)) % P
        for row in range(6)
    )
    check(image == (0,) * 6)

minor_histogram = Counter(
    determinant_mod([COLUMNS[q] for q in indices])
    for indices in combinations(range(13), 6)
)
check(sum(minor_histogram.values()) == 1716)
check(minor_histogram[0] == 135)
check(sum(count for determinant, count in minor_histogram.items() if determinant) == 1581)


# Restrict L to the 13-dimensional Boolean cube.
fibres = defaultdict(list)
boolean_norm_histogram = Counter()
boolean_factor_histogram = Counter()
for mask in range(1 << 13):
    image = subset_sum(mask, COLUMNS)
    fibres[image].append(mask)
    norm = multiplication_norm(image)
    profile = factor_profile(image)
    check((norm != 0) == (profile == ()), "unit tests disagree")
    boolean_norm_histogram[norm] += 1
    boolean_factor_histogram[profile] += 1
    if mask:
        for kappa in range(1, 7):
            check(owner_substitution(image, kappa) != (0,) * 6)

check(len(fibres) == 8184)
check(max(len(fibre) for fibre in fibres.values()) == 2)
double_fibres = [fibre for fibre in fibres.values() if len(fibre) == 2]
check(len(double_fibres) == 8)
check(fibres[(0,) * 6] == [0])
check(boolean_norm_histogram == Counter({
    0: 181,
    1: 682,
    2: 635,
    3: 676,
    4: 682,
    5: 683,
    6: 627,
    7: 677,
    8: 651,
    9: 660,
    10: 682,
    11: 688,
    12: 668,
}))
check(boolean_factor_histogram == Counter({
    (): 8011,
    (3,): 59,
    (5,): 55,
    (6,): 63,
    (3, 5): 1,
    (3, 6): 1,
    (5, 6): 1,
    (3, 5, 6): 1,
}))

COLLISION_LEFT = frozenset((3, 4, 6, 7, 9, 10))
COLLISION_RIGHT = frozenset((1, 5, 8, 12))
COLLISION_FREE = frozenset((0, 2, 11))
collision_common_sets = set()
for left_mask, right_mask in double_fibres:
    left = set(support(left_mask))
    right = set(support(right_mask))
    common = frozenset(left & right)
    left_only = frozenset(left - common)
    right_only = frozenset(right - common)
    if (left_only, right_only) != (COLLISION_LEFT, COLLISION_RIGHT):
        left_only, right_only = right_only, left_only
    check((left_only, right_only) == (COLLISION_LEFT, COLLISION_RIGHT))
    check(common <= COLLISION_FREE)
    collision_common_sets.add(common)
EXPECTED_COMMONS = {
    frozenset(),
    frozenset((0,)),
    frozenset((2,)),
    frozenset((11,)),
    frozenset((0, 2)),
    frozenset((0, 11)),
    frozenset((2, 11)),
    frozenset((0, 2, 11)),
}
check(collision_common_sets == EXPECTED_COMMONS)


# THM-2586 admissible theta-zero selectors, inserted only as coefficient rows.
selector_records = []
for displacement in range(1, 13):
    for mask in range(1 << 7):
        selector = tuple(6 if (mask >> ell) & 1 else 0 for ell in range(7))
        if displacement == 6 and any(selector[ell] != 0 for ell in (4, 5, 6)):
            continue
        if displacement == 7 and any(selector[ell] != 6 for ell in (4, 5, 6)):
            continue
        full = tuple(Y_ROWS[selector[ell]][ell] for ell in range(7))
        image = reduce_phi7(full)
        check(
            (multiplication_norm(image) != 0) == (factor_profile(image) == ()),
            "selector unit tests disagree",
        )
        selector_records.append((displacement, selector, image))
        for kappa in range(1, 7):
            check(owner_substitution(image, kappa) != (0,) * 6)

check(len(selector_records) == 1312)
selector_images = {record[2] for record in selector_records}
check(len(selector_images) == 32)
selector_norm_histogram = Counter(multiplication_norm(image) for _, _, image in selector_records)
check(selector_norm_histogram == Counter({
    0: 40,
    1: 80,
    2: 80,
    3: 82,
    4: 82,
    5: 330,
    7: 208,
    9: 122,
    10: 42,
    11: 82,
    12: 164,
}))
check(sum(image == (0,) * 6 for _, _, image in selector_records) == 0)
check(sum(multiplication_norm(image) != 0 for _, _, image in selector_records) == 1272)

distinct_selector_profiles = Counter(factor_profile(image) for image in selector_images)
labeled_selector_profiles = Counter(
    factor_profile(image) for _, _, image in selector_records
)
check(distinct_selector_profiles == Counter({(): 31, (3,): 1}))
check(labeled_selector_profiles == Counter({(): 1272, (3,): 40}))

NONUNIT_SELECTOR = (9, 5, 8, 4, 0, 0)
check(
    tuple(
        sum(entry != 0 for entry in owner_substitution(NONUNIT_SELECTOR, kappa))
        for kappa in range(1, 7)
    ) == (4, 6, 6, 4, 4, 6)
)
for displacement, selector, image in selector_records:
    is_expected_nonunit = (
        displacement not in (6, 7)
        and selector[0] == 0
        and selector[2:5] == (6, 6, 6)
        and selector[6] == 0
    )
    check((image == NONUNIT_SELECTOR) == is_expected_nonunit)
    if image == NONUNIT_SELECTOR:
        check(factor_profile(image) == (3,))


# The actual priority selector of THM-2586: only (s,ell)=(7,4..6) uses q=6.
priority_norms = []
for displacement in range(1, 13):
    selector = tuple(
        6 if displacement == 7 and ell in (4, 5, 6) else 0
        for ell in range(7)
    )
    full = tuple(Y_ROWS[selector[ell]][ell] for ell in range(7))
    image = reduce_phi7(full)
    priority_norms.append(multiplication_norm(image))
check(Counter(priority_norms) == Counter({7: 11, 10: 1}))
check(all(norm != 0 for norm in priority_norms))


print("THM-2590 exact coefficient-incidence audit")
print("field=F13 quotient=F13[z]/Phi7")
print("phi7_factors=(z^2+3z+1)(z^2+5z+1)(z^2+6z+1)")
print("omega=-u^11 in F13[u]/(u^12); socle_embedding=1")
print(f"linear_rank={len(PIVOTS)} kernel_dimension={13-len(PIVOTS)}")
print(f"basis_q0_to_q5_determinant={BASIS_DETERMINANT}")
print("six_minors_total=1716 nonzero=1581 zero=135")
print(f"six_minor_histogram={sorted(minor_histogram.items())}")
print(f"kernel_basis={KERNEL_BASIS}")
print("boolean_subsets=8192 distinct_images=8184 max_fibre=2")
print("boolean_zero_fibre=empty_only nonempty_images_nonzero=8191")
print("boolean_double_fibres=8")
print(f"boolean_collision_left={tuple(sorted(COLLISION_LEFT))}")
print(f"boolean_collision_right={tuple(sorted(COLLISION_RIGHT))}")
print("boolean_collision_common=every_subset_of_(0,2,11)")
print("boolean_units=8011 nonzero_zero_divisors=180")
print("boolean_nonempty_owner_colours=49146/49146")
print(f"boolean_factor_histogram={sorted(boolean_factor_histogram.items())}")
print("selector_patterns=1312 distinct_images=32 zero=0 units=1272 nonunits=40")
print("selector_owner_colours=7872/7872")
print("selector_distinct_units=31 distinct_nonunits=1")
print(f"selector_nonunit_class={NONUNIT_SELECTOR} factor=z^2+3z+1")
print("priority_selector_norms=7:11,10:1 all_unit=1")
print(f"checks={CHECKS}")
print("PASS")
