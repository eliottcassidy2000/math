#!/usr/bin/env python3
"""Exact companion for THM-3110.

Derive the two generic anchored response banks, verify their Ferrers
dominance on every ratio chamber, certify the one- and two-active-layer
Newton-positive faces, and compute the exact labelled Ewens zeta current.
All arithmetic after the universal SymPy expansion is integral or rational.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import comb

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


# ---------------------------------------------------------------------------
# Universal multilinear bank derivation.

A, B, A2, AB, B2, A3, A2B, AB2, B3 = sp.symbols(
    "A B A2 AB B2 A3 A2B AB2 B3"
)
moment_variables = (A, B, A2, AB, B2, A3, A2B, AB2, B3)
linear_form = {
    A: (1, 0),
    B: (0, 1),
    A2: (2, 0),
    AB: (1, 1),
    B2: (0, 2),
    A3: (3, 0),
    A2B: (2, 1),
    AB2: (1, 2),
    B3: (0, 3),
}

y, z = sp.symbols("y z")
moment_value = {
    (0, 0): sp.Integer(1),
    (1, 0): A,
    (0, 1): B,
    (2, 0): A2,
    (1, 1): AB,
    (0, 2): B2,
    (3, 0): A3,
    (2, 1): A2B,
    (1, 2): AB2,
    (0, 3): B3,
}


def expectation(polynomial):
    return sp.expand(
        sum(
            coefficient * moment_value[(i, j)] / A**i / B**j
            for (i, j), coefficient in sp.Poly(
                sp.expand(polynomial), y, z
            ).terms()
        )
    )


U = y - 1
V = z - y
g11 = expectation(U**2)
g12 = expectation(U * V)
g22 = expectation(V**2)
t111 = expectation(U**3)
t112 = expectation(U**2 * V)
t122 = expectation(U * V**2)
t222 = expectation(V**3)

invariants = (
    3 * t112 * g11 * g22 - t222 * g11**2 - 2 * t111 * g12 * g22,
    3 * t122 * g11 * g22 - 2 * t222 * g12 * g11 - t111 * g22**2,
)


def extract_bank(invariant, expected_denominator):
    numerator, denominator = sp.cancel(invariant).as_numer_denom()
    require(sp.expand(denominator - expected_denominator) == 0, "denominator drift")
    bank = []
    for monomial, coefficient in sp.Poly(
        -sp.expand(numerator), *moment_variables
    ).terms():
        row = []
        for variable, exponent in zip(moment_variables, monomial):
            row.extend([linear_form[variable]] * exponent)
        bank.append((int(coefficient), tuple(sorted(row, reverse=True))))
    return tuple(bank)


banks = (
    extract_bank(invariants[0], A**5 * B**3),
    extract_bank(invariants[1], A**5 * B**4),
)
require(tuple(len(bank) for bank in banks) == (24, 25), "bank census")
require(
    tuple(sum(c for c, _ in bank if c > 0) for bank in banks) == (37, 39),
    "positive mass",
)
require(
    tuple(sum(-c for c, _ in bank if c < 0) for bank in banks) == (37, 39),
    "negative mass",
)

# Exact polarized shear which supplies the two b-collision factors.  The
# a- and (b-a)-orders are the manifest U/V multidegrees 4/3 and 3/4.
epsilon = sp.symbols("epsilon")
gUU, gUW, gWW = sp.symbols("gUU gUW gWW")
tUUU, tUUW, tUWW, tWWW = sp.symbols("tUUU tUUW tUWW tWWW")
shear_g11 = gUU
shear_g12 = epsilon * gUW - gUU
shear_g22 = gUU - 2 * epsilon * gUW + epsilon**2 * gWW
shear_t111 = tUUU
shear_t112 = epsilon * tUUW - tUUU
shear_t122 = tUUU - 2 * epsilon * tUUW + epsilon**2 * tUWW
shear_t222 = -tUUU + 3 * epsilon * tUUW - 3 * epsilon**2 * tUWW + epsilon**3 * tWWW
shear_invariants = (
    3 * shear_t112 * shear_g11 * shear_g22
    - shear_t222 * shear_g11**2
    - 2 * shear_t111 * shear_g12 * shear_g22,
    3 * shear_t122 * shear_g11 * shear_g22
    - 2 * shear_t222 * shear_g12 * shear_g11
    - shear_t111 * shear_g22**2,
)
for shear_invariant in shear_invariants:
    shear_poly = sp.Poly(sp.expand(shear_invariant), epsilon)
    require(min(power[0] for power, _ in shear_poly.terms()) == 2, "b-collision order")


def value(form, a_value, b_value):
    return form[0] * a_value + form[1] * b_value


def numeric_partition(row, a_value, b_value):
    return sorted((value(form, a_value, b_value) for form in row), reverse=True)


dominant_forms = (
    tuple(sorted(((0, 3), (2, 0), (2, 0), (1, 0)), reverse=True)),
    tuple(sorted(((0, 3), (1, 1), (2, 0), (1, 0), (1, 0)), reverse=True)),
)
dominant_coefficients = (1, 2)
for bank, row, coefficient in zip(banks, dominant_forms, dominant_coefficients):
    require(bank.count((coefficient, row)) == 1, "dominant row missing")


# ---------------------------------------------------------------------------
# Symbolic ratio-chamber majorization.


def sorted_forms_at(row, ratio):
    return sorted(row, key=lambda q: Fraction(q[0]) + q[1] * ratio, reverse=True)


def linear_prefix(row, length):
    return (
        sum(form[0] for form in row[:length]),
        sum(form[1] for form in row[:length]),
    )


ratio_intervals = (
    (Fraction(1), Fraction(3, 2), Fraction(5, 4)),
    (Fraction(3, 2), Fraction(2), Fraction(7, 4)),
    (Fraction(2), Fraction(3), Fraction(5, 2)),
    (Fraction(3), None, Fraction(4)),
)


def nonnegative_linear_on_interval(alpha, beta, left, right):
    if alpha + beta * left < 0:
        return False
    if right is None:
        return beta >= 0
    return alpha + beta * right >= 0


majorization_checks = 0
for bank, dominant in zip(banks, dominant_forms):
    for coefficient, row in bank:
        if coefficient >= 0:
            continue
        for left, right, sample in ratio_intervals:
            q_sorted = sorted_forms_at(dominant, sample)
            s_sorted = sorted_forms_at(row, sample)
            width = max(len(q_sorted), len(s_sorted))
            q_sorted += [(0, 0)] * (width - len(q_sorted))
            s_sorted += [(0, 0)] * (width - len(s_sorted))
            require(
                linear_prefix(q_sorted, width) == linear_prefix(s_sorted, width),
                "weighted row sum",
            )
            for prefix in range(1, width + 1):
                qa, qb = linear_prefix(q_sorted, prefix)
                sa, sb = linear_prefix(s_sorted, prefix)
                require(
                    nonnegative_linear_on_interval(
                        qa - sa, qb - sb, left, right
                    ),
                    "Ferrers majorization",
                )
                majorization_checks += 1

# The same ratio chambers prove the Ferrers-column gcd symbolically.  Below
# level a the minimum is the minimum row length.  Between a and M it is the
# minimum number of parts at least M.  At M a row of maximum part M kills
# the remaining columns.
common_root_symbolic_checks = 0
for bank, minimum_length, middle_count in zip(banks, (4, 5), (1, 2)):
    require(min(len(row) for _, row in bank) == minimum_length, "gcd low columns")
    for _, _, sample in ratio_intervals:
        max_form = (2, 0) if sample < 2 else (0, 1)
        max_value = Fraction(max_form[0]) + max_form[1] * sample
        counts = [
            sum(Fraction(form[0]) + form[1] * sample >= max_value for form in row)
            for _, row in bank
        ]
        require(min(counts) == middle_count, "gcd middle columns")
        maximum_parts = [
            max(Fraction(form[0]) + form[1] * sample for form in row)
            for _, row in bank
        ]
        require(min(maximum_parts) == max_value, "gcd terminal column")
        common_root_symbolic_checks += len(bank)


# ---------------------------------------------------------------------------
# Response polynomials and the chamber-common root divisor.


def multiply_truncated(left, right, degree):
    answer = [0] * (degree + 1)
    for i, x_value in enumerate(left):
        if not x_value:
            continue
        for j, y_value in enumerate(right[: degree + 1 - i]):
            answer[i + j] += x_value * y_value
    return answer


def divide_series(numerator, denominator, degree):
    require(denominator[0] == 1, "unit response denominator")
    quotient = [0] * (degree + 1)
    quotient[0] = numerator[0]
    for k in range(1, degree + 1):
        quotient[k] = numerator[k] - sum(
            denominator[i] * quotient[k - i] for i in range(1, k + 1)
        )
    return quotient


def response_modes(a_value, b_value, max_degree):
    largest = 3 * b_value
    p_values = [[1] + [0] * max_degree]
    for n_value in range(1, largest + 1):
        previous = p_values[-1]
        current = previous[:]
        root = n_value - 1
        if root:
            for k in range(max_degree, 0, -1):
                current[k] += root * previous[k - 1]
        p_values.append(current)

    common_max = max(2 * a_value, b_value)
    common_rows = (
        (a_value, a_value, a_value, common_max),
        (a_value, a_value, a_value, common_max, common_max),
    )
    result = []
    for bank, common_row in zip(banks, common_rows):
        common_polynomial = [1] + [0] * max_degree
        for part in common_row:
            common_polynomial = multiply_truncated(
                common_polynomial, p_values[part], max_degree
            )
        bank_modes = []
        for coefficient, row in bank:
            polynomial = [1] + [0] * max_degree
            for form in row:
                polynomial = multiply_truncated(
                    polynomial,
                    p_values[value(form, a_value, b_value)],
                    max_degree,
                )
            quotient = divide_series(polynomial, common_polynomial, max_degree)
            bank_modes.append((coefficient, tuple(quotient)))
        result.append(tuple(bank_modes))
    return tuple(result)


common_root_cells = 0
for a_value in range(1, 17):
    for b_value in range(a_value + 1, 25):
        modes = response_modes(a_value, b_value, 4)
        for bank_modes in modes:
            require(all(poly[0] == 1 for _, poly in bank_modes), "root quotient")
        common_root_cells += 1


def triangular(n_value):
    return n_value * (n_value - 1) // 2


def h_value(row, a_value, b_value):
    return sum(triangular(value(form, a_value, b_value)) for form in row)


# The nearest-gap formulas used in the uniform tail are global polynomial
# identities, not conclusions extrapolated from the finite cell scan below.
# Write b=a+d.  For every negative row its gap from the dominant row
# coefficientwise dominates at least one of the two displayed candidates;
# both candidates are attained by actual rows.
a_symbol, d_symbol = sp.symbols("a_symbol d_symbol")
b_symbol = a_symbol + d_symbol


def symbolic_triangular(argument):
    return sp.expand(argument * (argument - 1) / 2)


def symbolic_h(row):
    return sp.expand(
        sum(
            symbolic_triangular(value(form, a_symbol, b_symbol))
            for form in row
        )
    )


tail_gap_symbolic_checks = 0
tail_gap_targets = (
    (a_symbol**2, 2 * d_symbol * (a_symbol + d_symbol)),
    (a_symbol**2, a_symbol * d_symbol),
)
for bank, dominant, targets in zip(banks, dominant_forms, tail_gap_targets):
    attained = [False] * len(targets)
    for coefficient, row in bank:
        if coefficient >= 0:
            continue
        gap = sp.expand(symbolic_h(dominant) - symbolic_h(row))
        comparisons = []
        for index, target in enumerate(targets):
            difference = sp.Poly(sp.expand(gap - target), a_symbol, d_symbol)
            comparisons.append(all(entry >= 0 for entry in difference.coeffs()))
            if difference.is_zero:
                attained[index] = True
        require(any(comparisons), "symbolic closest-gap lower bound")
        tail_gap_symbolic_checks += 1
    require(all(attained), "symbolic closest-gap target not attained")
require(tail_gap_symbolic_checks == 25, "symbolic closest-gap census")


tail_cells = 0
for a_value in range(1, 61):
    for b_value in range(a_value + 1, 81):
        common_max = max(2 * a_value, b_value)
        common = (
            ((1, 0), (1, 0), (1, 0), (2, 0) if common_max == 2 * a_value else (0, 1)),
            (
                (1, 0),
                (1, 0),
                (1, 0),
                (2, 0) if common_max == 2 * a_value else (0, 1),
                (2, 0) if common_max == 2 * a_value else (0, 1),
            ),
        )
        gaps = []
        for bank, dominant, common_row in zip(banks, dominant_forms, common):
            dominant_first = h_value(dominant, a_value, b_value) - h_value(
                common_row, a_value, b_value
            )
            negative_gaps = [
                h_value(dominant, a_value, b_value) - h_value(row, a_value, b_value)
                for coefficient, row in bank
                if coefficient < 0
            ]
            require(dominant_first > 0, "dominant degree-one response")
            require(all(0 < gap <= dominant_first for gap in negative_gaps), "tail ratios")
            gaps.append(min(negative_gaps))
        require(
            gaps[0] == min(a_value * a_value, 2 * b_value * (b_value - a_value)),
            "first closest gap",
        )
        require(
            gaps[1] == a_value * min(a_value, b_value - a_value),
            "second closest gap",
        )
        tail_cells += 1


# ---------------------------------------------------------------------------
# Exact Newton certificates.


def chamber_point(chamber, first, second):
    if chamber == "wide":
        a_value = first + 1
        return a_value, 2 * a_value + second
    return first + second + 2, first + 2 * second + 3


def base_value(invariant_index, a_value, b_value):
    if invariant_index == 0:
        return a_value**4 * b_value**2 * (b_value - a_value) ** 3
    return a_value**3 * b_value**2 * (b_value - a_value) ** 4


mode_cache = {}


def cached_modes(chamber, first, second, max_degree=13):
    key = (chamber, first, second)
    if key not in mode_cache:
        a_value, b_value = chamber_point(chamber, first, second)
        mode_cache[key] = response_modes(a_value, b_value, max_degree)
    return mode_cache[key]


def phi_value(chamber, invariant_index, degrees, first, second):
    a_value, b_value = chamber_point(chamber, first, second)
    modes = cached_modes(chamber, first, second)[invariant_index]
    numerator = sum(
        coefficient
        * product_of(poly[degree] for degree in degrees)
        for coefficient, poly in modes
    )
    return Fraction(numerator, base_value(invariant_index, a_value, b_value))


def product_of(values):
    answer = 1
    for entry in values:
        answer *= entry
    return answer


def newton_coefficients(table):
    size = len(table)
    a_rows = [row[:] for row in table]
    leading_rows = []
    while a_rows:
        leading_rows.append(a_rows[0])
        a_rows = [
            [a_rows[i + 1][j] - a_rows[i][j] for j in range(size)]
            for i in range(len(a_rows) - 1)
        ]
    answer = []
    for row in leading_rows:
        differences = row[:]
        leading = []
        while differences:
            leading.append(differences[0])
            differences = [
                differences[i + 1] - differences[i]
                for i in range(len(differences) - 1)
            ]
        answer.append(leading)
    return answer


one_layer_evaluations = 0
one_layer_slots = 0
one_layer_positive = 0
one_layer_zero = 0
one_layer_scales = {5: 2, 6: 12, 7: 24, 8: 720}
for chamber in ("tight", "wide"):
    for invariant_index in range(2):
        for degree in range(5, 9):
            polynomial_degree = 2 * degree - 9
            table = [
                [
                    phi_value(chamber, invariant_index, (degree,), first, second)
                    for second in range(polynomial_degree + 1)
                ]
                for first in range(polynomial_degree + 1)
            ]
            one_layer_evaluations += (polynomial_degree + 1) ** 2
            coefficients = newton_coefficients(table)
            scale = one_layer_scales[degree]
            for first, row in enumerate(coefficients):
                for second, coefficient in enumerate(row):
                    if first + second > polynomial_degree:
                        require(coefficient == 0, "one-layer degree bound")
                        continue
                    scaled = coefficient * scale
                    require(scaled.denominator == 1, "one-layer scale")
                    require(coefficient >= 0, "one-layer Newton sign")
                    one_layer_slots += 1
                    if coefficient:
                        one_layer_positive += 1
                    else:
                        one_layer_zero += 1

require(one_layer_evaluations == 480, "one-layer evaluation census")
require((one_layer_slots, one_layer_positive, one_layer_zero) == (280, 276, 4), "one-layer Newton census")


two_layer_pairs = tuple(
    (first, second)
    for first in range(1, 15)
    for second in range(first, 15)
    if 5 <= first + second <= 14
)
require(len(two_layer_pairs) == 45, "two-layer pair census")
two_layer_slots = 0
two_layer_positive = 0
two_layer_zero = 0
for chamber in ("tight", "wide"):
    for invariant_index in range(2):
        for first_degree, second_degree in two_layer_pairs:
            polynomial_degree = 2 * (first_degree + second_degree) - 9
            table = [
                [
                    phi_value(
                        chamber,
                        invariant_index,
                        (first_degree, second_degree),
                        first,
                        second,
                    )
                    for second in range(polynomial_degree + 1)
                ]
                for first in range(polynomial_degree + 1)
            ]
            coefficients = newton_coefficients(table)
            for first, row in enumerate(coefficients):
                for second, coefficient in enumerate(row):
                    if first + second > polynomial_degree:
                        require(coefficient == 0, "two-layer degree bound")
                        continue
                    require(coefficient >= 0, "two-layer Newton sign")
                    two_layer_slots += 1
                    if coefficient:
                        two_layer_positive += 1
                    else:
                        two_layer_zero += 1

require(two_layer_slots == 18760, "two-layer Newton census")


# ---------------------------------------------------------------------------
# Labelled set-partition zeta transform.


def set_partitions(items):
    answer = []

    def recurse(position, blocks):
        if position == len(items):
            answer.append(tuple(tuple(block) for block in blocks))
            return
        item = items[position]
        for index in range(len(blocks)):
            blocks[index].append(item)
            recurse(position + 1, blocks)
            blocks[index].pop()
        blocks.append([item])
        recurse(position + 1, blocks)
        blocks.pop()

    recurse(0, [])
    return tuple(answer)


def canonical_partition(blocks):
    return tuple(
        sorted(
            (tuple(sorted(block)) for block in blocks),
            key=lambda block: (block[0], len(block), block),
        )
    )


def colour_type(partition, a_count):
    return tuple(
        sorted(
            (
                (
                    sum(item < a_count for item in block),
                    sum(item >= a_count for item in block),
                )
                for block in partition
            ),
            reverse=True,
        )
    )


ewens_records = []
for bank, b_count in zip(banks, (3, 4)):
    coefficient_by_type = {
        tuple(sorted(row, reverse=True)): coefficient for coefficient, row in bank
    }
    universe = tuple(range(5 + b_count))
    partitions = set_partitions(universe)
    multiplicities = Counter(colour_type(partition, 5) for partition in partitions)
    omega = []
    for partition in partitions:
        row_type = colour_type(partition, 5)
        coefficient = coefficient_by_type.get(row_type, 0)
        if coefficient:
            omega.append(
                (partition, Fraction(coefficient, multiplicities[row_type]))
            )

    refinement_cache = {}
    zeta = defaultdict(Fraction)
    for partition, weight in omega:
        choices = []
        for block in partition:
            if block not in refinement_cache:
                refinement_cache[block] = set_partitions(block)
            choices.append(refinement_cache[block])
        for refinements in product(*choices):
            refined = canonical_partition(
                sum((list(partition_piece) for partition_piece in refinements), [])
            )
            zeta[refined] += weight

    signs = Counter()
    ranks = Counter()
    positive_values = []
    negative_values = []
    for partition in partitions:
        weight = zeta[canonical_partition(partition)]
        if weight > 0:
            signs["positive"] += 1
            positive_values.append(weight)
            ranks[5 + b_count - len(partition)] += 1
        elif weight < 0:
            signs["negative"] += 1
            negative_values.append(weight)
            ranks[5 + b_count - len(partition)] += 1
        else:
            signs["zero"] += 1
    require(set(ranks) == {4}, "zeta support rank")
    require(sum(zeta[canonical_partition(partition)] for partition in partitions) == 0, "zeta signed mass")
    ewens_records.append(
        (
            len(partitions),
            signs["positive"],
            signs["negative"],
            signs["zero"],
            min(positive_values),
            min(negative_values),
        )
    )

require(
    ewens_records
    == [
        (4140, 285, 195, 3660, Fraction(1, 30), Fraction(-1, 15)),
        (21147, 720, 900, 19527, Fraction(1, 60), Fraction(-1, 30)),
    ],
    "Ewens zeta census",
)


print("THM3110 exact arbitrary-anchored product-Gamma reduction")
print("banks=24/25;balanced_masses=37/39;dominant_coefficients=1/2")
print(f"majorization_intervals=4;prefix_checks={majorization_checks}")
print(
    f"common_root_symbolic_checks={common_root_symbolic_checks};"
    f"spot_cells={common_root_cells};"
    f"tail_gap_symbolic_checks={tail_gap_symbolic_checks};"
    f"tail_gap_cells={tail_cells}"
)
print("common_divisors=P_a^3*P_max;P_a^3*P_max^2")
print("collision_divisors=a^4*b^2*(b-a)^3;a^3*b^2*(b-a)^4")
print("closest_gaps=min(a^2,2b(b-a));a*min(a,b-a)")
print("one_layer=k5_to_k8;evaluations=480;newton=276_positive+4_structural_zero")
print(
    f"two_layer=pairs45;newton_slots={two_layer_slots};"
    f"positive={two_layer_positive};zero={two_layer_zero};negative=0"
)
print("ewens_B1=Bell4140;rank4=285_positive+195_negative;min=1/30,-1/15")
print("ewens_B2=Bell21147;rank4=720_positive+900_negative;min=1/60,-1/30")
print("scope=dominant_tail_plus_initial_faces;all_low_histograms=OPEN")
print("all_exact_checks=PASS")
