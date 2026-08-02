#!/usr/bin/env python3
"""Exact companion for candidate THM-3075.

The theorem is analytic/algebraic.  This dependency-free companion checks its
three finite interfaces independently over Q:

1. the Laguerre/falling-factorial coefficient identity;
2. total nonnegativity of finite Toeplitz coefficient banks;
3. the row-transform and Cauchy--Binet identities, followed by direct
   arbitrary-offset Hankel determinants.

The final strict-prefix bank is explicitly finite evidence for the larger
open conjecture, not part of the all-order theorem.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(value, length):
    out = Fraction(1)
    for step in range(length):
        out *= value + step
    return out


def falling(value, length):
    out = Fraction(1)
    for step in range(length):
        out *= value - step
    return out


def determinant(matrix):
    """Bareiss determinant over Q."""
    a = [[Fraction(entry) for entry in row] for row in matrix]
    size = len(a)
    if size == 0:
        return Fraction(1)
    previous = Fraction(1)
    orientation = 1
    for column in range(size - 1):
        pivot = next((row for row in range(column, size) if a[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
            orientation *= -1
        value = a[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                a[row][target] = (
                    a[row][target] * value
                    - a[row][column] * a[column][target]
                ) / previous
        previous = value
    return orientation * a[-1][-1]


def sign(value):
    return int(value > 0) - int(value < 0)


def laguerre_coefficients(gap):
    return tuple(
        factorial(gap) * comb(gap, index) // factorial(index)
        for index in range(gap + 1)
    )


def coefficient(coefficients, index):
    if 0 <= index < len(coefficients):
        return coefficients[index]
    return 0


def reciprocal_gamma(base, index):
    """Gamma(base)/Gamma(base+index) on the audited integer strip."""
    if index >= 0:
        return Fraction(1, 1) / rising(base, index)
    return rising(base + index, -index)


def strict_gap_moment(base, gap, index):
    return rising(base + gap, index) / rising(base, index) ** 2


def transform_moment(base, gap, index):
    coefficients = laguerre_coefficients(gap)
    return sum(
        Fraction(coefficients[shift]) * reciprocal_gamma(base, index - shift)
        for shift in range(gap + 1)
    ) / rising(base, gap)


def feed(digest, *items):
    digest.update("|".join(str(item) for item in items).encode("ascii"))
    digest.update(b"\n")


print("THM-3075 LAGUERRE PF STRICT RECIPROCAL BETA GAP")
digest = sha256()

# (x+1)_m = sum_j m! binom(m,j)/j! x^(underline j).
laguerre_cells = 0
for gap in range(0, 11):
    coefficients = laguerre_coefficients(gap)
    require(coefficients[0] == factorial(gap), "Laguerre constant changed")
    require(coefficients[-1] == 1, "Laguerre leading coefficient changed")
    for value in tuple(Fraction(n, 3) for n in range(-8, 25)):
        left = rising(value + 1, gap)
        right = sum(
            Fraction(entry) * falling(value, index)
            for index, entry in enumerate(coefficients)
        )
        require(left == right, "falling-factorial expansion failed")
        feed(digest, "L", gap, value, left)
        laguerre_cells += 1
print(f"laguerre_falling_cells={laguerre_cells} gaps=0..10")

# Finite exact controls for the PF-infinity orientation.  Rows and columns are
# increasing, and T_(u,z)=c_(u-z), so this is the lower Toeplitz convention
# used in the proof.  All sampled minors must be nonnegative.
toeplitz_cells = 0
toeplitz_positive = 0
for gap in range(1, 9):
    coefficients = laguerre_coefficients(gap)
    row_sets = {
        order: tuple(combinations(range(7), order))
        for order in range(1, 5)
    }
    column_sets = {
        order: tuple(combinations(range(-gap, 7), order))
        for order in range(1, 5)
    }
    for order in range(1, 5):
        for rows in row_sets[order]:
            for columns in column_sets[order]:
                value = determinant(
                    [
                        [coefficient(coefficients, row - column) for column in columns]
                        for row in rows
                    ]
                )
                require(value >= 0, "Laguerre Toeplitz minor became negative")
                toeplitz_positive += int(value > 0)
                toeplitz_cells += 1
print(
    f"laguerre_toeplitz_minor_cells={toeplitz_cells} "
    f"positive={toeplitz_positive} negative=0"
)

# The exact row-transform identity.  The parameter restriction base>gap is
# the theorem's positive extended-Gamma strip, not merely a test filter.
transform_cells = 0
parameters = tuple(
    (gap, Fraction(gap) + tail)
    for gap in range(1, 7)
    for tail in (Fraction(1, 3), Fraction(1), Fraction(7, 2))
)
for gap, base in parameters:
    for index in range(0, 19):
        direct = strict_gap_moment(base, gap, index)
        transformed = transform_moment(base, gap, index)
        require(direct == transformed, "PF row-transform identity failed")
        feed(digest, "R", gap, base, index, direct)
        transform_cells += 1
print(f"row_transform_cells={transform_cells} parameters={len(parameters)}")

# Literal matrix multiplication and Cauchy--Binet on nonconsecutive controls.
# Reversing the columns is built into K.  Each summand is checked nonnegative,
# and the selected z_i=u_i term is checked strictly positive.
cauchy_cells = 0
cauchy_terms = 0
for gap, base in parameters[:12]:
    coefficients = laguerre_coefficients(gap)
    for rows, columns in (
        ((0, 2, 5), (0, 3, 7)),
        ((1, 4, 7), (1, 2, 8)),
        ((0, 1, 4, 8), (0, 2, 5, 9)),
    ):
        reverse_columns = tuple(reversed(columns))
        nodes = tuple(sorted({row - shift for row in rows for shift in range(gap + 1)}))
        a_matrix = [
            [coefficient(coefficients, row - node) for node in nodes]
            for row in rows
        ]
        k_matrix = [
            [reciprocal_gamma(base, node + column) for column in reverse_columns]
            for node in nodes
        ]
        product = [
            [
                sum(a_matrix[i][t] * k_matrix[t][j] for t in range(len(nodes)))
                / rising(base, gap)
                for j in range(len(rows))
            ]
            for i in range(len(rows))
        ]
        literal = [
            [strict_gap_moment(base, gap, row + column) for column in reverse_columns]
            for row in rows
        ]
        require(product == literal, "literal PF matrix product failed")
        direct = determinant(literal)
        expansion = Fraction(0)
        selected_term = None
        for chosen_positions in combinations(range(len(nodes)), len(rows)):
            left = determinant(
                [[a_matrix[i][position] for position in chosen_positions] for i in range(len(rows))]
            )
            right = determinant(
                [[k_matrix[position][j] for j in range(len(rows))] for position in chosen_positions]
            )
            require(left >= 0, "Cauchy--Binet Toeplitz term became negative")
            require(right > 0, "reversed reciprocal-Gamma term lost strict positivity")
            term = left * right / rising(base, gap) ** len(rows)
            require(term >= 0, "Cauchy--Binet term became negative")
            expansion += term
            chosen_nodes = tuple(nodes[position] for position in chosen_positions)
            if chosen_nodes == rows:
                selected_term = term
            cauchy_terms += 1
        require(expansion == direct, "Cauchy--Binet determinant identity failed")
        require(selected_term is not None and selected_term > 0,
                "distinguished strict Cauchy--Binet term failed")
        require(direct > 0, "reversed-column determinant is not positive")
        feed(digest, "C", gap, base, rows, columns, direct)
        cauchy_cells += 1
print(f"cauchy_binet_cells={cauchy_cells} terms={cauchy_terms}")

# Direct arbitrary-offset determinant path.  This does not call the transform.
minor_cells = 0
for gap, base in parameters:
    values = [strict_gap_moment(base, gap, index) for index in range(25)]
    for order in range(2, 7):
        expected = (-1) ** (order * (order - 1) // 2)
        offset_sets = tuple(combinations(range(8), order))
        controls = [(item, item) for item in offset_sets]
        controls += [
            (offset_sets[position], offset_sets[-1 - position])
            for position in range(min(24, len(offset_sets)))
        ]
        for rows, columns in controls:
            value = determinant(
                [[values[row + column] for column in columns] for row in rows]
            )
            require(sign(value) == expected,
                    "direct arbitrary-offset checkerboard sign failed")
            feed(digest, "H", gap, base, rows, columns, value)
            minor_cells += 1
print(f"direct_generalized_minor_cells={minor_cells} orders=2..6")

# Independent small controls for the abstract PF row-transform criterion, not
# tied to Laguerre coefficients.  Each coefficient polynomial is a product of
# positive two-term factors and therefore has a manifest bidiagonal Toeplitz
# factorization.
generic_pf_cells = 0
generic_factors = (
    (Fraction(1, 3), Fraction(2), Fraction(7, 2)),
    (Fraction(1, 5), Fraction(4, 3), Fraction(3), Fraction(8)),
    (Fraction(2, 7), Fraction(5, 4), Fraction(11, 3), Fraction(9, 2), Fraction(12)),
)
for factors in generic_factors:
    coefficients = [Fraction(1)]
    for factor in factors:
        updated = [Fraction(0)] * (len(coefficients) + 1)
        for index, entry in enumerate(coefficients):
            updated[index] += entry
            updated[index + 1] += factor * entry
        coefficients = updated
    width = len(coefficients) - 1
    base = Fraction(width + 2)
    values = [
        sum(
            coefficients[shift] * reciprocal_gamma(base, index - shift)
            for shift in range(width + 1)
        )
        for index in range(25)
    ]
    for order in range(2, 7):
        expected = (-1) ** (order * (order - 1) // 2)
        for rows, columns in (
            (tuple(range(order)), tuple(range(order))),
            (tuple(2 * index for index in range(order)),
             tuple(3 * index - int(index > 0) for index in range(order))),
        ):
            value = determinant(
                [[values[row + column] for column in columns] for row in rows]
            )
            require(sign(value) == expected, "generic PF row transform failed")
            generic_pf_cells += 1
print(f"generic_pf_transform_minor_cells={generic_pf_cells} families={len(generic_factors)}")

# Finite hostile scout for the larger OPEN strict-prefix conjecture.  The
# cumulative deficits S_j are positive integers; e_0=-S_0 and
# e_j=S_(j-1)-S_j.  This bank deliberately includes oscillating inventories
# and four widely separated rational meshes.
strict_prefix_cells = 0
meshes = (
    (Fraction(1, 31), Fraction(2, 7), Fraction(9, 2)),
    (Fraction(1, 7), Fraction(2, 3), Fraction(20)),
    (Fraction(3, 5), Fraction(11, 6), Fraction(13, 2), Fraction(40)),
    (Fraction(1, 101), Fraction(1, 3), Fraction(5, 4), Fraction(8), Fraction(100)),
)
def inventory_values(shapes, deficits, maximum):
    cumulative = tuple(-entry for entry in deficits)
    exponents = (cumulative[0],) + tuple(
        cumulative[index] - cumulative[index - 1]
        for index in range(1, len(cumulative))
    )
    require(all(sum(exponents[: index + 1]) < 0 for index in range(len(exponents))),
            "hostile scout admitted a non-strict prefix")
    values = [Fraction(1)]
    for index in range(maximum):
        ratio = Fraction(1)
        for shape, exponent in zip(shapes, exponents):
            factor = shape + index
            ratio = ratio * factor**exponent if exponent >= 0 else ratio / factor ** (-exponent)
        values.append(values[-1] * ratio)
    return values

for mesh_index, shapes in enumerate(meshes):
    deficit_bank = []
    for seed in range(1, 25):
        deficits = tuple(1 + ((seed * (2 * index + 1) + index * index + mesh_index) % 9)
                         for index in range(len(shapes)))
        if any(deficits[index] < deficits[index - 1] for index in range(1, len(deficits))):
            deficit_bank.append(deficits)
    for deficits in tuple(dict.fromkeys(deficit_bank)):
        values = inventory_values(shapes, deficits, 23)
        for order, rows, columns in (
            (3, (0, 1, 4), (0, 3, 8)),
            (4, (0, 1, 3, 8), (0, 2, 5, 9)),
            (5, (0, 1, 2, 6, 10), (0, 2, 4, 7, 11)),
            (6, (0, 1, 3, 5, 8, 11), (0, 2, 4, 6, 9, 12)),
        ):
            value = determinant(
                [[values[row + column] for column in columns] for row in rows]
            )
            expected = (-1) ** (order * (order - 1) // 2)
            require(sign(value) == expected, "finite strict-prefix hostile found")
            strict_prefix_cells += 1
print(
    f"open_strict_prefix_hostile_cells={strict_prefix_cells} "
    "orders=3..6 finite_only"
)

print(f"exact_value_digest={digest.hexdigest()}")
print("scope=integer_gap_two_shape_all_orders_for_a_gt_gap;general_strict_prefix_open")
print("all_exact_checks=PASS")
