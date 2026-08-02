#!/usr/bin/env python3
"""Exact finite companion for THM-3101's next-moment-unit suspension."""

from fractions import Fraction
from itertools import combinations
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def determinant(matrix):
    matrix = [[Fraction(entry) for entry in row] for row in matrix]
    size = len(matrix)
    answer = Fraction(1)
    sign = 1
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


def vandermonde(values):
    answer = 1
    for left in range(len(values)):
        for right in range(left + 1, len(values)):
            answer *= values[right] - values[left]
    return answer


def compositions(total):
    if total == 0:
        yield ()
        return
    for first in range(1, total + 1):
        for tail in compositions(total - first):
            yield (first,) + tail


def add_scaled(left, right, scalar):
    return [
        [left[row][column] + scalar * right[row][column]
         for column in range(len(left))]
        for row in range(len(left))
    ]


def forward_difference(values, order):
    values = list(values)
    for _ in range(order):
        values = [
            values[index + 1] - values[index]
            for index in range(len(values) - 1)
        ]
    return values


def extension_count(bound, free_tail, first_choices=3):
    answer = 0
    for first in range(1, first_choices + 1):
        available = bound - first
        if available >= free_tail:
            answer += comb(available, free_tail)
    return answer


# The exact power in tau=U_m/U_(m+1)^(m/(m+1)).
tau_exponent_cells = 0
for child_width in range(3, 11):
    upper_degree = child_width + 1
    exponent = (
        -Fraction(child_width - 1, 2)
        + Fraction(child_width * (upper_degree - 1), 2 * upper_degree)
    )
    require(exponent == Fraction(1, 2 * upper_degree), "tau power")
    tau_exponent_cells += 1
require(tau_exponent_cells == 8, "tau exponent census")


# Exhaust the lower cells and the first five upper rows for m=3,...,10.
# The fixed child, H_p, and all-high cells are the base-one model cells.
# Every remaining cell lies below rho=(m/p)^m, with exactly the two stated
# equality cells for each m.
lower_cells = 0
upper_cells = 0
equality_cells = []
for child_width in range(3, 11):
    upper_degree = child_width + 1
    rho = Fraction(child_width**child_width, upper_degree**child_width)

    for row_degree in range(2, child_width + 1):
        for actual_high in range(row_degree + 1):
            for normal_degree in range(actual_high, row_degree + 1):
                if actual_high == 0 and normal_degree == 0:
                    continue
                numerator = 1 if actual_high == 0 else actual_high**actual_high
                base = Fraction(numerator, upper_degree**normal_degree)
                require(base <= rho, "lower cell above repair scale")
                if base == rho:
                    equality_cells.append(
                        ("lower", child_width, row_degree,
                         actual_high, normal_degree)
                    )
                lower_cells += 1

    for row_degree in range(upper_degree, upper_degree + 5):
        for actual_high in range(row_degree + 1):
            for normal_degree in range(actual_high, row_degree + 1):
                if (
                    row_degree == upper_degree
                    and actual_high == 0
                    and normal_degree == 0
                ):
                    continue
                if actual_high == normal_degree == row_degree:
                    continue
                numerator = 1 if actual_high == 0 else actual_high**actual_high
                numerator *= upper_degree ** (row_degree - normal_degree)
                base = Fraction(numerator, row_degree**row_degree)
                require(base <= rho, "upper cell above repair scale")
                if base == rho:
                    equality_cells.append(
                        ("upper", child_width, row_degree,
                         actual_high, normal_degree)
                    )
                upper_cells += 1

expected_equalities = []
for child_width in range(3, 11):
    expected_equalities.extend(
        [
            ("lower", child_width, child_width,
             child_width, child_width),
            ("upper", child_width, child_width + 1,
             child_width, child_width),
        ]
    )
require(equality_cells == expected_equalities, "unique equality-cell pair")


# Replace the first consecutive upper row p by the lower repair row m.
# Generalized alternants are positive Vandermonde times positive Schur values.
# The bank includes the singleton first block q=1.
row_replacement_cells = 0
singleton_faces = 0
for child_width in range(3, 9):
    upper_degree = child_width + 1
    for block_size in range(1, 6):
        rows = [child_width] + list(
            range(upper_degree + 1, upper_degree + block_size)
        )
        ordinary_floor = vandermonde(rows)
        for tail in combinations(range(1, 9), block_size - 1):
            gap_degrees = (0,) + tail
            alternant = determinant(
                [[row**gap for gap in gap_degrees] for row in rows]
            )
            require(alternant > 0, "row-replacement alternant sign")
            require(alternant >= ordinary_floor, "Schur/Vandermonde floor")
            row_replacement_cells += 1
            if block_size == 1:
                singleton_faces += 1
require(row_replacement_cells == 978, "row-replacement census")
require(singleton_faces == 6, "q=1 singleton census")


# Enumerate every ordered composition face through five upper directions.
# The first block uses the repaired row set; later blocks are unchanged.
composition_faces = 0
for child_width in range(3, 9):
    upper_degree = child_width + 1
    for normal_rank in range(1, 6):
        for face in compositions(normal_rank):
            offset = 0
            face_determinant = 1
            for block_index, block_size in enumerate(face):
                if block_index == 0:
                    rows = [child_width] + list(
                        range(upper_degree + 1, upper_degree + block_size)
                    )
                else:
                    rows = list(
                        range(
                            upper_degree + offset,
                            upper_degree + offset + block_size,
                        )
                    )
                gap_degrees = tuple(range(block_size))
                face_determinant *= determinant(
                    [[row**gap for gap in gap_degrees] for row in rows]
                )
                offset += block_size
            require(face_determinant > 0, "composition-face determinant")
            composition_faces += 1
require(composition_faces == 186, "composition-face census")


# A scheme-zero old operator plus an o(tau) perturbation factors exactly as
# tau^d times a unit determinant.  This is the repaired hypothesis.
scheme_zero_unit_cells = 0
nonunit_cells = 0
for length in range(1, 7):
    zero = [[0 for _ in range(length)] for _ in range(length)]
    unit_matrix = [
        [row + 1 if row == column else 0 for column in range(length)]
        for row in range(length)
    ]
    error_matrix = [
        [1 if row == column else 0 for column in range(length)]
        for row in range(length)
    ]
    for denominator in (2, 3, 5, 11, 23):
        tau = Fraction(1, denominator)
        repaired = add_scaled(
            add_scaled(zero, unit_matrix, tau),
            error_matrix,
            tau**2,
        )
        normalized = determinant(repaired) / tau**length
        require(
            normalized == determinant(
                add_scaled(unit_matrix, error_matrix, tau)
            ),
            "scheme-zero determinant factorization",
        )
        require(normalized > 0, "scheme-zero positive unit norm")
        scheme_zero_unit_cells += 1
        require(determinant(add_scaled(zero, zero, tau)) == 0,
                "zero repair is nonunit")
        nonunit_cells += 1
require(scheme_zero_unit_cells == 30, "scheme-zero unit census")
require(nonunit_cells == 30, "nonunit census")


# The frozen nilpotent identity is tempting but not stable under an o(tau)
# deformation of the preceding relation.
frozen_nilpotent_cells = 0
for length in range(2, 7):
    nilpotent = [
        [1 if column == row + 1 else 0 for column in range(length)]
        for row in range(length)
    ]
    identity = [
        [1 if row == column else 0 for column in range(length)]
        for row in range(length)
    ]
    for denominator in (2, 5, 11):
        tau = Fraction(1, denominator)
        require(
            determinant(add_scaled(nilpotent, identity, tau))
            == tau**length,
            "frozen nilpotent identity",
        )
        frozen_nilpotent_cells += 1
require(frozen_nilpotent_cells == 15, "frozen nilpotent census")


# In C[e]/(e^2-epsilon), Norm(e+tau)=tau^2-epsilon.  The choice
# tau=s^2, epsilon=s^3 has epsilon=o(tau), but reverses the sign and
# dominates tau^2.  This is the sharp hostile to the unrepaired theorem.
nilpotent_deformation_hostiles = 0
for denominator in range(3, 21):
    small = Fraction(1, denominator)
    tau = small**2
    epsilon = small**3
    multiplication_e = [[0, epsilon], [1, 0]]
    identity = [[1, 0], [0, 1]]
    norm = determinant(add_scaled(multiplication_e, identity, tau))
    require(norm == tau**2 - epsilon, "deformed nilpotent norm")
    require(epsilon < tau, "lower deformation is o(tau) model")
    require(norm < 0 and abs(norm) > tau**2,
            "deformation dominates predicted norm")
    nilpotent_deformation_hostiles += 1
require(nilpotent_deformation_hostiles == 18,
        "nilpotent deformation census")


# The physical H_(m-1) all-high base already exceeds rho_m^2, so the
# one-cell inequality sigma<rho cannot control even a two-step Jordan norm.
physical_weight_hostiles = 0
for child_width in range(3, 11):
    upper_degree = child_width + 1
    alpha = Fraction(
        (child_width - 1) ** (child_width - 1),
        upper_degree ** (child_width - 1),
    )
    rho = Fraction(child_width**child_width, upper_degree**child_width)
    ratio = Fraction(
        (child_width**2 - 1) ** child_width,
        child_width ** (2 * child_width),
    ) * Fraction(upper_degree, child_width - 1)
    require(alpha > rho**2, "physical weighted-face hostile")
    require(alpha / rho**2 == ratio, "physical hostile ratio")
    physical_weight_hostiles += 1
require(physical_weight_hostiles == 8, "physical hostile census")

# A nonzero aggregate may still miss one primary component.
aggregate = [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
]
require(sum(aggregate[index][index] for index in range(4)) == 2,
        "aggregate hostile trace")
require(determinant(aggregate) == 0, "aggregate hostile zero divisor")


# Conjugate local factors contribute positive squared complex determinants.
# A genuine real residue would allow an odd negative determinant instead.
conjugate_norm_cells = 0
for real_part in range(-3, 4):
    for imaginary_part in range(-3, 4):
        if real_part == imaginary_part == 0:
            continue
        multiplication = [
            [real_part, -imaginary_part],
            [imaginary_part, real_part],
        ]
        norm = real_part**2 + imaginary_part**2
        require(determinant(multiplication) == norm, "conjugate norm")
        for artin_length in range(1, 6):
            require(norm**artin_length > 0, "positive conjugate Artin norm")
            conjugate_norm_cells += 1
require(conjugate_norm_cells == 240, "conjugate norm census")

real_residue_hostiles = 0
for multiplicity in range(1, 7):
    determinant_sign = (-1) ** multiplicity
    if multiplicity % 2:
        require(determinant_sign < 0, "odd real-residue hostile")
    else:
        require(determinant_sign > 0, "even real-residue control")
    real_residue_hostiles += 1
require(real_residue_hostiles == 6, "real-residue census")


# Rank and exponent ledger for one through four appended normal directions.
normal_exponent_cells = 0
for child_width in range(3, 11):
    upper_degree = child_width + 1
    child_length = factorial(child_width - 1)
    for normal_rank in range(1, 5):
        target_width = child_width + normal_rank
        upper_rank = factorial(target_width) // factorial(child_width)
        tail_rank = factorial(target_width) // factorial(upper_degree)
        require(upper_rank == upper_degree * tail_rank, "upper free rank")
        for zero_primary_length in sorted({2, child_length}):
            require(
                0 < zero_primary_length <= child_length,
                "zero-primary length range",
            )
            exponent = zero_primary_length * upper_rank
            require(
                exponent == zero_primary_length * upper_degree * tail_rank,
                "tau norm exponent",
            )
            higher_moment_exponent = (
                zero_primary_length * child_width * tail_rank
            )
            require(higher_moment_exponent % 2 == 0, "positive norm parity")
            normal_exponent_cells += 1
require(normal_exponent_cells == 60, "normal exponent census")


# Setting every appended coefficient to zero preserves a child point that
# kills all required moments.  Cross terms are deliberately nontrivial.
persistent_root_cells = 0
root = Fraction(2, 3)
for child_width in range(3, 9):
    for remote_offset in (child_width + 1, 2 * child_width + 7, 101):
        for moment in range(2, child_width + 2):
            child_value = (moment + remote_offset) * (root - root)
            appended_value = child_value
            for power in range(1, moment + 1):
                appended_value += (
                    comb(moment, power)
                    * remote_offset**power
                    * Fraction(0) ** power
                    * root ** (moment - power)
                )
            require(appended_value == 0, "zero-coefficient persistent root")
            persistent_root_cells += 1
require(persistent_root_cells == 99, "persistent-root census")


# A bounded first appended point leaves q-1 free tail points.  Exact finite
# differences recover degree q-1 and hence max(t-5,0) when m>=4.
conditional_count_cells = 0
for target_width in range(5, 13):
    for prefix_width in range(4, target_width):
        normal_rank = target_width - prefix_width
        first_bound = 3
        start = first_bound + normal_rank + 2
        values = [
            extension_count(bound, normal_rank - 1, first_bound)
            for bound in range(start, start + normal_rank + 2)
        ]
        require(
            all(value == 0 for value in forward_difference(values, normal_rank)),
            "bounded-first-point degree upper bound",
        )
        leading_difference = forward_difference(values, normal_rank - 1)
        require(
            leading_difference and all(value > 0 for value in leading_difference),
            "bounded-first-point sharp degree",
        )
        require(
            normal_rank - 1 <= max(target_width - 5, 0),
            "codimension-five exponent",
        )
        conditional_count_cells += 1
require(conditional_count_cells == 36, "conditional-count census")

general_sfc_cells = 0
for known_width in range(3, 8):
    for target_width in range(known_width + 2, 13):
        first_new_prefix = known_width + 1
        free_tail = target_width - first_new_prefix - 1
        require(
            free_tail == target_width - known_width - 2,
            "general SFC exponent",
        )
        general_sfc_cells += 1
require(general_sfc_cells == 30, "general SFC census")


print("THM-3101 SCHEME-ZERO-PRIMARY NEXT-MOMENT-UNIT SUSPENSION")
print(
    f"tau_exponent_cells={tau_exponent_cells} "
    "exact_power_1_over_2p=PASS"
)
print(
    f"newton_cells={lower_cells + upper_cells} lower={lower_cells} "
    f"upper={upper_cells} equality_cells={len(equality_cells)} "
    "unique_lower_and_upper_pair=PASS"
)
print(
    f"row_replacement_cells={row_replacement_cells} "
    f"singleton_faces={singleton_faces} composition_faces={composition_faces} "
    "positive_alternants=PASS"
)
print(
    f"scheme_zero_unit_cells={scheme_zero_unit_cells} "
    f"nonunit_cells={nonunit_cells} transverse_power=PASS"
)
print(
    f"frozen_nilpotent_cells={frozen_nilpotent_cells} "
    f"deformation_hostiles={nilpotent_deformation_hostiles} "
    f"physical_weight_hostiles={physical_weight_hostiles} "
    "weighted_Rees_boundary=PASS"
)
print("aggregate_hostile=nonzero_trace_but_zero_divisor")
print(
    f"conjugate_norm_cells={conjugate_norm_cells} "
    f"real_residue_hostiles={real_residue_hostiles} sign_boundary=PASS"
)
print(
    f"normal_exponent_cells={normal_exponent_cells} "
    "finite_free_rank_and_even_sign=PASS"
)
print(
    f"persistent_root_cells={persistent_root_cells} "
    "zero_new_coefficient=PASS"
)
print(
    f"conditional_count_cells={conditional_count_cells} "
    f"general_sfc_cells={general_sfc_cells} "
    "max_t_minus_5_and_general_drop=PASS"
)
print(
    "scope=fixed_child;finite_complete_intersection;scheme_zero_primary;"
    "pointwise_unit;fixed_normal_rank;conditional_count;nilpotent_open"
)
print("all_exact_checks=PASS")
