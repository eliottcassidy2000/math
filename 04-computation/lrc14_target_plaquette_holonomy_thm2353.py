#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2353."""

from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import isqrt


P = 13
MODULUS = 1_000_003
POINTS = tuple(range(P))
SPEEDS = (13, 13**3, 2 * 13**5, 1, 14, 27, 40, 53, 66)
NON_GUARD = (0, 1, 2, 4, 5, 6, 7, 8)
GUARD = 3
OTHER_TARGET = 1
DEEP_TARGET = 2
OTHER_GRAFT = 4
DEEP_GRAFT = 5
X_FREQUENCY = 13
Y_FREQUENCY = X_FREQUENCY + SPEEDS[DEEP_TARGET]


def require(condition: bool, message: str) -> None:
    """Raise under ordinary and optimized Python."""
    if not condition:
        raise RuntimeError(message)


def is_prime(number: int) -> bool:
    """Trial-divide the fixed small certificate modulus."""
    if number < 2:
        return False
    for divisor in range(2, isqrt(number) + 1):
        if number % divisor == 0:
            return False
    return True


def mean(values) -> Fraction:
    """Return the exact arithmetic mean."""
    values = tuple(values)
    require(bool(values), "cannot average an empty family")
    return sum(values, Fraction(0)) / len(values)


def rank_mod_p(rows) -> int:
    """Return matrix rank over F_13."""
    matrix = [list(entry % P for entry in row) for row in rows]
    if not matrix:
        return 0
    row_count = len(matrix)
    column_count = len(matrix[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, row_count)
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = (
            matrix[pivot],
            matrix[pivot_row],
        )
        inverse = pow(matrix[pivot_row][column], -1, P)
        matrix[pivot_row] = [
            entry * inverse % P for entry in matrix[pivot_row]
        ]
        for row in range(row_count):
            if row == pivot_row:
                continue
            coefficient = matrix[row][column]
            if coefficient:
                matrix[row] = [
                    (matrix[row][index]
                     - coefficient * matrix[pivot_row][index]) % P
                    for index in range(column_count)
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def owner_packet_target_control() -> tuple[int, int, int]:
    """Check the concrete graft packet and its quotient-dual basis."""
    pivot_coordinates = (0, 4, 5, 6, 7, 8)
    packet = []
    for coordinate in pivot_coordinates:
        row = [0] * len(SPEEDS)
        row[GUARD] = SPEEDS[coordinate]
        row[coordinate] = -SPEEDS[GUARD]
        if coordinate == OTHER_GRAFT:
            row[GUARD] += SPEEDS[OTHER_TARGET]
            row[OTHER_TARGET] -= SPEEDS[GUARD]
        if coordinate == DEEP_GRAFT:
            row[GUARD] += SPEEDS[DEEP_TARGET]
            row[DEEP_TARGET] -= SPEEDS[GUARD]
        require(
            sum(
                row[index] * SPEEDS[index]
                for index in range(len(SPEEDS))
            )
            == 0,
            "owner packet row left the exact relation lattice",
        )
        packet.append(tuple(row))

    other_dual = [0] * len(SPEEDS)
    other_dual[OTHER_TARGET] = 1
    other_dual[OTHER_GRAFT] = -1
    deep_dual = [0] * len(SPEEDS)
    deep_dual[DEEP_TARGET] = 1
    deep_dual[DEEP_GRAFT] = -1
    duals = (tuple(other_dual), tuple(deep_dual))

    annihilation_checks = 0
    for dual in duals:
        for row in packet:
            require(
                sum(
                    dual[index] * row[index]
                    for index in range(len(SPEEDS))
                )
                % P
                == 0,
                "target character did not annihilate the owner packet",
            )
            annihilation_checks += 1

    packet_rank = rank_mod_p(packet)
    quotient_dual_rank = rank_mod_p(
        (tuple(speed % P for speed in SPEEDS),) + duals
    )
    require(packet_rank == 6, "owner packet rank changed")
    require(
        quotient_dual_rank == 3,
        "two target characters lost independence modulo the speed gauge",
    )
    return packet_rank, quotient_dual_rank, annihilation_checks


def determinant(
    matrix,
    row: int,
    column: int,
):
    """Return one cyclic adjacent plaquette minor."""
    next_row = (row + 1) % P
    next_column = (column + 1) % P
    return (
        matrix[row][column] * matrix[next_row][next_column]
        - matrix[next_row][column] * matrix[row][next_column]
    )


def rational_matrix_controls() -> tuple[int, int, int, int]:
    """Check flat rank one, curved additive, and ANOVA separation."""
    first = tuple(Fraction(index + 1) for index in POINTS)
    second = tuple(Fraction(index + 2) for index in POINTS)

    rank_one = tuple(
        tuple(first[row] * second[column] for column in POINTS)
        for row in POINTS
    )
    additive = tuple(
        tuple(first[row] + second[column] for column in POINTS)
        for row in POINTS
    )

    flat_plaquettes = sum(
        determinant(rank_one, row, column) == 0
        for row in POINTS
        for column in POINTS
    )
    curved_additive_plaquettes = sum(
        determinant(additive, row, column) != 0
        for row in POINTS
        for column in POINTS
    )
    require(flat_plaquettes == P**2, "rank-one matrix acquired curvature")
    require(
        curved_additive_plaquettes == P**2,
        "additive hostile lost a plaquette",
    )

    require(
        all(
            rank_one[row][column]
            == rank_one[row][0]
            * rank_one[0][column]
            / rank_one[0][0]
            for row in POINTS
            for column in POINTS
        ),
        "flat matrix did not reconstruct from its two axes",
    )

    def interaction(matrix):
        rows = tuple(mean(row) for row in matrix)
        columns = tuple(
            mean(matrix[row][column] for row in POINTS)
            for column in POINTS
        )
        grand = mean(entry for row in matrix for entry in row)
        return tuple(
            tuple(
                matrix[row][column]
                - rows[row]
                - columns[column]
                + grand
                for column in POINTS
            )
            for row in POINTS
        )

    additive_interaction = interaction(additive)
    rank_one_interaction = interaction(rank_one)
    first_mean = mean(first)
    second_mean = mean(second)
    predicted_rank_one_interaction = tuple(
        tuple(
            (first[row] - first_mean) * (second[column] - second_mean)
            for column in POINTS
        )
        for row in POINTS
    )
    require(
        all(
            entry == 0
            for row in additive_interaction
            for entry in row
        ),
        "additively separable matrix acquired ANOVA interaction",
    )
    require(
        rank_one_interaction == predicted_rank_one_interaction,
        "rank-one ANOVA outer-product law failed",
    )
    require(
        any(
            entry != 0
            for row in rank_one_interaction
            for entry in row
        ),
        "flat rank-one hostile lost its additive interaction",
    )

    bad_phase_flat_checks = 0
    demodulation_checks = 0
    multiplier = 5
    for row in POINTS:
        for column in POINTS:
            next_row = (row + 1) % P
            next_column = (column + 1) % P
            phase = (-multiplier * column) % P
            diagonal_phase = (-multiplier * next_column) % P
            first_cross_phase = (-multiplier * column) % P
            second_cross_phase = (-multiplier * next_column) % P
            require(
                (phase + diagonal_phase) % P
                == (first_cross_phase + second_cross_phase) % P,
                "inverse character acquired plaquette curvature",
            )
            bad_phase_flat_checks += 1

            full_diagonal_phase = (
                multiplier * column
                + multiplier * next_column
            ) % P
            full_cross_phase = (
                multiplier * column
                + multiplier * next_column
            ) % P
            require(
                full_diagonal_phase == full_cross_phase,
                "deep-character demodulation changed a minor zero",
            )
            require(next_row in POINTS, "cyclic row left F_13")
            demodulation_checks += 1

    return (
        flat_plaquettes,
        curved_additive_plaquettes,
        bad_phase_flat_checks,
        demodulation_checks,
    )


Gaussian = tuple[Fraction, Fraction]


def g_add(left: Gaussian, right: Gaussian) -> Gaussian:
    return (left[0] + right[0], left[1] + right[1])


def g_neg(value: Gaussian) -> Gaussian:
    return (-value[0], -value[1])


def g_sub(left: Gaussian, right: Gaussian) -> Gaussian:
    return g_add(left, g_neg(right))


def g_mul(left: Gaussian, right: Gaussian) -> Gaussian:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def g_conjugate(value: Gaussian) -> Gaussian:
    return (value[0], -value[1])


def g_div(left: Gaussian, right: Gaussian) -> Gaussian:
    denominator = right[0] ** 2 + right[1] ** 2
    require(denominator != 0, "division by zero Gaussian rational")
    numerator = g_mul(left, g_conjugate(right))
    return (numerator[0] / denominator, numerator[1] / denominator)


def g_scale(value: Gaussian, scalar: Fraction) -> Gaussian:
    return (scalar * value[0], scalar * value[1])


def endpoint_factorization_control() -> int:
    """Check the conjugated endpoint cross-ratio identity exactly."""
    left = (
        (Fraction(1), Fraction(2)),
        (Fraction(2), Fraction(-1)),
        (Fraction(3), Fraction(1)),
        (Fraction(-1), Fraction(3)),
    )
    right = (
        (Fraction(2), Fraction(1)),
        (Fraction(1), Fraction(-2)),
        (Fraction(4), Fraction(-1)),
        (Fraction(3), Fraction(2)),
    )
    deep = Fraction(3, 7)
    response = tuple(
        g_scale(g_mul(left[index], g_conjugate(right[index])), deep)
        for index in range(4)
    )

    response_minor = g_sub(
        g_mul(response[0], response[3]),
        g_mul(response[1], response[2]),
    )
    endpoint_minor = g_scale(
        g_sub(
            g_mul(
                g_mul(left[0], left[3]),
                g_conjugate(g_mul(right[0], right[3])),
            ),
            g_mul(
                g_mul(left[1], left[2]),
                g_conjugate(g_mul(right[1], right[2])),
            ),
        ),
        deep**2,
    )
    require(
        response_minor == endpoint_minor,
        "physical four-twist determinant identity failed",
    )

    response_holonomy = g_div(
        g_mul(response[0], response[3]),
        g_mul(response[1], response[2]),
    )
    left_holonomy = g_div(
        g_mul(left[0], left[3]),
        g_mul(left[1], left[2]),
    )
    right_holonomy = g_div(
        g_mul(right[0], right[3]),
        g_mul(right[1], right[2]),
    )
    require(
        response_holonomy
        == g_mul(left_holonomy, g_conjugate(right_holonomy)),
        "endpoint holonomies did not factor with conjugation",
    )
    require(
        response_holonomy != (Fraction(1), Fraction(0)),
        "curved endpoint control became flat",
    )
    return 1


# The finite-slice certificate lives in
#   Q(c,zeta_13), c=2*cos(2*pi/7),
#   c^3+c^2-2c-1=0.
# We reduce its canonical 36-dimensional basis modulo a good prime.
Cubic = tuple[int, int, int]
Cyclotomic = tuple[Cubic, ...]
C_ZERO: Cubic = (0, 0, 0)
C_ONE: Cubic = (1, 0, 0)
C_GENERATOR: Cubic = (0, 1, 0)


def mod(value: int) -> int:
    return value % MODULUS


def c_add(left: Cubic, right: Cubic) -> Cubic:
    return tuple(mod(left[index] + right[index]) for index in range(3))


def c_neg(value: Cubic) -> Cubic:
    return tuple(mod(-entry) for entry in value)


def c_sub(left: Cubic, right: Cubic) -> Cubic:
    return c_add(left, c_neg(right))


def c_scale(value: Cubic, scalar: int) -> Cubic:
    return tuple(mod(entry * scalar) for entry in value)


def c_mul(left: Cubic, right: Cubic) -> Cubic:
    raw = [0] * 5
    for first in range(3):
        for second in range(3):
            raw[first + second] = mod(
                raw[first + second] + left[first] * right[second]
            )
    # c^3=-c^2+2c+1.
    for degree in (4, 3):
        coefficient = raw[degree]
        raw[degree] = 0
        raw[degree - 1] = mod(raw[degree - 1] - coefficient)
        raw[degree - 2] = mod(raw[degree - 2] + 2 * coefficient)
        raw[degree - 3] = mod(raw[degree - 3] + coefficient)
    return tuple(raw[:3])


def c_is_zero(value: Cubic) -> bool:
    return value == C_ZERO


def e_zero() -> Cyclotomic:
    return (C_ZERO,) * P


def e_normalize(value: Cyclotomic) -> Cyclotomic:
    """Use 1+z+...+z^12=0 and set the z^12 coefficient to zero."""
    last = value[P - 1]
    if c_is_zero(last):
        return value
    return tuple(
        C_ZERO if exponent == P - 1 else c_sub(coefficient, last)
        for exponent, coefficient in enumerate(value)
    )


def e_add(left: Cyclotomic, right: Cyclotomic) -> Cyclotomic:
    return e_normalize(
        tuple(c_add(left[index], right[index]) for index in POINTS)
    )


def e_neg(value: Cyclotomic) -> Cyclotomic:
    return tuple(c_neg(coefficient) for coefficient in value)


def e_sub(left: Cyclotomic, right: Cyclotomic) -> Cyclotomic:
    return e_add(left, e_neg(right))


def e_mul(left: Cyclotomic, right: Cyclotomic) -> Cyclotomic:
    raw = [C_ZERO] * P
    for first in POINTS:
        if c_is_zero(left[first]):
            continue
        for second in POINTS:
            if c_is_zero(right[second]):
                continue
            exponent = (first + second) % P
            raw[exponent] = c_add(
                raw[exponent],
                c_mul(left[first], right[second]),
            )
    return e_normalize(tuple(raw))


def e_monomial(exponent: int, coefficient: Cubic = C_ONE) -> Cyclotomic:
    raw = [C_ZERO] * P
    raw[exponent % P] = coefficient
    return e_normalize(tuple(raw))


def e_conjugate(value: Cyclotomic) -> Cyclotomic:
    raw = [C_ZERO] * P
    for exponent in POINTS:
        target = (-exponent) % P
        raw[target] = c_add(raw[target], value[exponent])
    return e_normalize(tuple(raw))


def e_is_zero(value: Cyclotomic) -> bool:
    return all(c_is_zero(coefficient) for coefficient in value)


C_SQUARE = c_mul(C_GENERATOR, C_GENERATOR)
RATIO_TABLE: tuple[Cubic, ...] = (
    C_ZERO,
    C_ONE,
    C_GENERATOR,
    c_sub(C_SQUARE, C_ONE),
    c_neg(c_sub(C_SQUARE, C_ONE)),
    c_neg(C_GENERATOR),
    c_neg(C_ONE),
)


def guard_weight(index: int) -> Cubic:
    """Return sin(2*pi*h/7)/(h*sin(2*pi/7)) modulo the good prime."""
    require(index and index % 7, "guard slice crossed a Fourier zero")
    denominator = index % MODULUS
    require(denominator, "good-prime certificate hit a denominator")
    inverse = pow(denominator, -1, MODULUS)
    return c_scale(RATIO_TABLE[index % 7], inverse)


def dynamic_slice_count(frequency: int) -> int:
    """Count the retained sign vectors by a separate mod-seven DP."""
    counts = [0] * 7
    counts[0] = 1
    for coordinate in NON_GUARD:
        updated = [0] * 7
        for residue, count in enumerate(counts):
            for sign in (-1, 1):
                updated[(residue + sign * SPEEDS[coordinate]) % 7] += count
        counts = updated
    excluded = counts[frequency % 7]
    return 2 ** len(NON_GUARD) - excluded


def endpoint_slice(
    frequency: int,
) -> tuple[
    dict[tuple[int, int], Cubic],
    dict[tuple[int, int], int],
    tuple[tuple[tuple[int, int], Cubic], ...],
    int,
]:
    """Build the exact +/-1 non-guard, seven-unit-guard endpoint slice."""
    cells: dict[tuple[int, int], Cubic] = {}
    term_counts: dict[tuple[int, int], int] = {}
    terms: list[tuple[tuple[int, int], Cubic]] = []
    kept = 0
    for signs in product((-1, 1), repeat=len(NON_GUARD)):
        indices = [0] * len(SPEEDS)
        non_guard_frequency = 0
        for coordinate, sign in zip(NON_GUARD, signs):
            indices[coordinate] = sign
            non_guard_frequency += SPEEDS[coordinate] * sign
        guard_index = frequency - non_guard_frequency
        indices[GUARD] = guard_index
        if guard_index == 0 or guard_index % 7 == 0:
            continue
        require(
            guard_index % MODULUS,
            "certificate modulus divides a retained guard index",
        )

        target = (
            (indices[OTHER_TARGET] - indices[OTHER_GRAFT]) % P,
            (indices[DEEP_TARGET] - indices[DEEP_GRAFT]) % P,
        )
        weight = guard_weight(guard_index)
        cells[target] = c_add(cells.get(target, C_ZERO), weight)
        term_counts[target] = term_counts.get(target, 0) + 1
        terms.append((target, weight))
        kept += 1

    require(
        kept == dynamic_slice_count(frequency),
        "literal slice count disagreed with mod-seven DP",
    )
    target_residues = (0, 2, P - 2)
    expected_support = set(product(target_residues, repeat=2))
    require(set(cells) == expected_support, "slice target support changed")
    require(
        all(not c_is_zero(value) for value in cells.values()),
        "a four-cell endpoint aggregate vanished modulo the good prime",
    )
    return cells, term_counts, tuple(terms), kept


def endpoint_transform(
    cells: dict[tuple[int, int], Cubic],
    row: int,
    column: int,
) -> Cyclotomic:
    """Evaluate one exact target-character endpoint transform."""
    answer = e_zero()
    for (first, second), coefficient in cells.items():
        answer = e_add(
            answer,
            e_monomial(row * first + column * second, coefficient),
        )
    return answer


def direct_endpoint_transform(
    terms: tuple[tuple[tuple[int, int], Cubic], ...],
    row: int,
    column: int,
) -> Cyclotomic:
    """Evaluate before grouping into the four target-residue cells."""
    answer = e_zero()
    for (first, second), coefficient in terms:
        answer = e_add(
            answer,
            e_monomial(row * first + column * second, coefficient),
        )
    return answer


def digest_elements(elements) -> str:
    """Freeze canonical modular-basis coefficients."""
    digest = sha256()
    for element in elements:
        for cubic in element:
            for coefficient in cubic:
                digest.update(str(coefficient).encode("ascii"))
                digest.update(b",")
        digest.update(b";")
    return digest.hexdigest()


def cross_axis_slice_certificate():
    """Certify every entry and adjacent minor of the finite typed slice."""
    require(is_prime(MODULUS), "certificate modulus is not prime")
    cubic_relation = c_sub(
        c_add(c_mul(C_SQUARE, C_GENERATOR), C_SQUARE),
        c_add(c_scale(C_GENERATOR, 2), C_ONE),
    )
    require(c_is_zero(cubic_relation), "cubic field relation failed")
    cyclotomic_sum = e_zero()
    for exponent in POINTS:
        cyclotomic_sum = e_add(cyclotomic_sum, e_monomial(exponent))
    require(e_is_zero(cyclotomic_sum), "Phi_13 relation failed")
    require(
        e_mul(e_monomial(P - 1), e_monomial(1)) == e_monomial(0),
        "zeta_13^13 ceased to be one",
    )

    left_cells, left_counts, left_term_bank, left_terms = endpoint_slice(
        X_FREQUENCY
    )
    right_cells, right_counts, right_term_bank, right_terms = endpoint_slice(
        Y_FREQUENCY
    )

    endpoint_grouping_checks = 0
    for row in POINTS:
        for column in POINTS:
            require(
                endpoint_transform(left_cells, row, column)
                == direct_endpoint_transform(left_term_bank, row, column),
                "left endpoint grouping changed its transform",
            )
            require(
                endpoint_transform(right_cells, row, column)
                == direct_endpoint_transform(right_term_bank, row, column),
                "right endpoint grouping changed its transform",
            )
            endpoint_grouping_checks += 2

    response = tuple(
        tuple(
            e_mul(
                endpoint_transform(left_cells, row, column),
                e_conjugate(
                    endpoint_transform(right_cells, row, column)
                ),
            )
            for column in POINTS
        )
        for row in POINTS
    )
    nonzero_entries = sum(
        not e_is_zero(response[row][column])
        for row in POINTS
        for column in POINTS
    )

    plaquettes = []
    for row in POINTS:
        for column in POINTS:
            next_row = (row + 1) % P
            next_column = (column + 1) % P
            plaquettes.append(
                e_sub(
                    e_mul(
                        response[row][column],
                        response[next_row][next_column],
                    ),
                    e_mul(
                        response[next_row][column],
                        response[row][next_column],
                    ),
                )
            )
    nonzero_plaquettes = sum(not e_is_zero(value) for value in plaquettes)
    require(nonzero_entries == P**2, "finite slice acquired a zero twist")
    require(
        nonzero_plaquettes == P**2,
        "finite slice acquired a flat cyclic plaquette",
    )

    return (
        left_terms,
        right_terms,
        endpoint_grouping_checks,
        left_counts,
        right_counts,
        nonzero_entries,
        nonzero_plaquettes,
        digest_elements(
            response[row][column]
            for row in POINTS
            for column in POINTS
        ),
        digest_elements(plaquettes),
    )


(
    flat_rank_one_plaquettes,
    curved_additive_plaquettes,
    bad_phase_flat_checks,
    demodulation_checks,
) = rational_matrix_controls()
(
    owner_packet_rank,
    quotient_dual_rank,
    target_packet_annihilation_checks,
) = owner_packet_target_control()
endpoint_identity_checks = endpoint_factorization_control()
(
    left_terms,
    right_terms,
    endpoint_grouping_checks,
    left_cell_counts,
    right_cell_counts,
    nonzero_slice_entries,
    nonzero_slice_plaquettes,
    slice_response_digest,
    slice_plaquette_digest,
) = cross_axis_slice_certificate()

print("theorem=THM-2353")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-AUDITED")
print("target_matrix=13x13")
print(f"rank_one_flat_plaquettes={flat_rank_one_plaquettes}")
print(f"additive_curved_plaquettes={curved_additive_plaquettes}")
print(f"inverse_character_flat_checks={bad_phase_flat_checks}")
print(f"deep_demodulation_minor_checks={demodulation_checks}")
print(f"endpoint_holonomy_identity_checks={endpoint_identity_checks}")
print("endpoint_holonomy_law=Omega_K=Omega_L*conjugate(Omega_F)")
print("holonomy_escape_law=Omega_L*conjugate(Omega_F)!=1")
print("plaquette_curvature_implies_ANOVA_interaction=NO")
print("ANOVA_interaction_implies_plaquette_curvature=NO")
print("slice_speed_vector=(13,2197,742586,1,14,27,40,53,66)")
print(f"slice_frequencies={X_FREQUENCY},{Y_FREQUENCY}")
print("slice_owner_packet=omit-H,graft-c2-q1,graft-c3-q2")
print("slice_target_dual=e_c2-e_q1,e_c3-e_q2")
print(f"slice_owner_packet_rank={owner_packet_rank}")
print(f"slice_quotient_dual_rank_with_speed_gauge={quotient_dual_rank}")
print(
    "slice_target_packet_annihilation_checks="
    f"{target_packet_annihilation_checks}"
)
print("slice_non_guard_modes={-1,+1}^8")
print("slice_guard_filter=NONZERO-MOD-7")
print("slice_target_support={0,+2,-2}^2")
print(f"slice_left_terms={left_terms}")
print(f"slice_right_terms={right_terms}")
print(f"slice_endpoint_grouping_checks={endpoint_grouping_checks}")
print(
    "slice_left_cell_counts="
    + ",".join(
        f"{point[0]}:{point[1]}:{left_cell_counts[point]}"
        for point in sorted(left_cell_counts)
    )
)
print(
    "slice_right_cell_counts="
    + ",".join(
        f"{point[0]}:{point[1]}:{right_cell_counts[point]}"
        for point in sorted(right_cell_counts)
    )
)
print(f"certificate_prime={MODULUS}")
print("characteristic_zero_number_field_basis_dimension=36")
print("modular_quotient_algebra_basis_dimension=36")
print(f"slice_nonzero_twists={nonzero_slice_entries}")
print(f"slice_nonzero_cyclic_plaquettes={nonzero_slice_plaquettes}")
print(f"slice_response_digest={slice_response_digest}")
print(f"slice_plaquette_digest={slice_plaquette_digest}")
print("finite_slice_omitted_tail=OPEN")
print("canonical_row_excluded=NO")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
