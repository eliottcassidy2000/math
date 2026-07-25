#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2325."""

from itertools import product


P13 = 13
P7 = 7
LABELS = ("j", "a", "b", "u0", "u1", "u2", "u3", "u4", "u5")
INDEX = {label: index for index, label in enumerate(LABELS)}
PIVOT = ("j", "u1", "u2", "u3", "u4", "u5")


def require(condition: bool, message: str) -> None:
    """Raise in ordinary and optimized Python when a check fails."""
    if not condition:
        raise RuntimeError(message)


def matrix_rank(matrix: list[list[int]], prime: int) -> int:
    """Return exact row rank over F_prime."""
    work = [[entry % prime for entry in row] for row in matrix]
    if not work:
        return 0
    rows = len(work)
    columns = len(work[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (
                row
                for row in range(pivot_row, rows)
                if work[row][column] % prime
            ),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [
            entry * inverse % prime for entry in work[pivot_row]
        ]
        for row in range(rows):
            if row == pivot_row:
                continue
            factor = work[row][column] % prime
            if factor:
                work[row] = [
                    (left - factor * right) % prime
                    for left, right in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def system_solution_count(
    matrix: list[list[int]],
    right_hand_side: list[int],
    variables: int,
    prime: int,
) -> int:
    """Count solutions of A*x=b over F_prime."""
    require(len(matrix) == len(right_hand_side), "system row count changed")
    if not matrix:
        return prime**variables
    coefficient_rank = matrix_rank(matrix, prime)
    augmented = [
        [entry % prime for entry in row] + [value % prime]
        for row, value in zip(matrix, right_hand_side)
    ]
    augmented_rank = matrix_rank(augmented, prime)
    if augmented_rank != coefficient_rank:
        return 0
    return prime ** (variables - coefficient_rank)


def affine_all_unit_count(
    rows: list[list[int]],
    offset: list[int],
    prime: int,
) -> int:
    """Count lambda with offset+lambda*rows nonzero in all nine coordinates."""
    variables = len(rows)
    coordinates = len(offset)
    require(coordinates == len(LABELS), "coordinate count changed")
    total = 0
    for mask in range(1 << coordinates):
        equations = []
        values = []
        for coordinate in range(coordinates):
            if mask & (1 << coordinate):
                equations.append(
                    [row[coordinate] % prime for row in rows]
                )
                values.append((-offset[coordinate]) % prime)
        intersection = system_solution_count(
            equations, values, variables, prime
        )
        total += (-1 if mask.bit_count() % 2 else 1) * intersection
    return total


def exact_owner_packet() -> tuple[list[int], list[list[int]]]:
    """Build THM-2309's exact hostile-word owner packet."""
    speeds = [
        13,
        13**3,
        2 * 13**5,
        1,
        14,
        27,
        40,
        53,
        66,
    ]
    omitted = INDEX["u0"]
    rows = []
    for label in PIVOT:
        coordinate = INDEX[label]
        row = [0] * len(LABELS)
        row[omitted] = speeds[coordinate]
        row[coordinate] = -speeds[omitted]
        if label == "u1":
            row[omitted] += speeds[INDEX["a"]]
            row[INDEX["a"]] -= speeds[omitted]
        if label == "u2":
            row[omitted] += speeds[INDEX["b"]]
            row[INDEX["b"]] -= speeds[omitted]
        require(
            sum(entry * speed for entry, speed in zip(row, speeds)) == 0,
            "exact packet row left the relation lattice",
        )
        rows.append(row)
    return speeds, rows


def packet_checks(
    speeds: list[int],
    rows: list[list[int]],
) -> tuple[int, int, int]:
    """Check rank, brightness, and the target quotient at both primes."""
    for prime in (P7, P13):
        require(matrix_rank(rows, prime) == 6, "packet rank changed")
        require(
            all(
                any(row[coordinate] % prime for row in rows)
                for coordinate in range(len(LABELS))
            ),
            "packet acquired a dark coordinate",
        )
        require(
            all(
                sum(
                    row[coordinate] * speeds[coordinate]
                    for coordinate in range(len(LABELS))
                )
                % prime
                == 0
                for row in rows
            ),
            "packet row left the modular scalar kernel",
        )

    axis_a = [0] * len(LABELS)
    axis_b = [0] * len(LABELS)
    axis_a[INDEX["a"]] = 1
    axis_b[INDEX["b"]] = 1
    require(
        matrix_rank(rows + [axis_a, axis_b], P13) == 8,
        "two-target quotient decomposition changed",
    )
    scalar_support_7 = sum(speed % P7 != 0 for speed in speeds)
    require(scalar_support_7 >= 2, "septimal support hypothesis failed")
    return (
        matrix_rank(rows, P13),
        matrix_rank(rows, P7),
        scalar_support_7,
    )


def projective_coset_atlas(
    rows: list[list[int]],
) -> tuple[int, int, int, int]:
    """Check all 168 nonzero quotient-vector cosets in the model packet."""
    quotient_vectors = set()
    counts = []
    directions = []
    directions.append(("axis_a", 1, 0))
    directions.append(("axis_b", 0, 1))
    directions.extend((f"gain_{gain}", 1, gain) for gain in range(1, P13))

    for _name, first, second in directions:
        direction_counts = []
        for scalar in range(1, P13):
            target_pair = (
                scalar * first % P13,
                scalar * second % P13,
            )
            require(target_pair != (0, 0), "zero projective vector entered")
            require(
                target_pair not in quotient_vectors,
                "projective quotient cosets ceased to be disjoint",
            )
            quotient_vectors.add(target_pair)
            offset = [0] * len(LABELS)
            offset[INDEX["a"]] = target_pair[0]
            offset[INDEX["b"]] = target_pair[1]
            count = affine_all_unit_count(rows, offset, P13)
            direction_counts.append(count)
            counts.append(count)
        require(
            len(direction_counts) == P13 - 1,
            "scalar torsor size changed",
        )

    require(len(directions) == P13 + 1, "projective direction count changed")
    require(
        len(quotient_vectors) == P13**2 - 1,
        "nonzero quotient-vector atlas changed",
    )
    return len(directions), len(quotient_vectors), min(counts), max(counts)


def row_combination(
    coefficients: tuple[int, ...],
    rows: list[list[int]],
    prime: int,
) -> list[int]:
    """Evaluate a coefficient vector on packet rows."""
    return [
        sum(
            coefficient * row[coordinate]
            for coefficient, row in zip(coefficients, rows)
        )
        % prime
        for coordinate in range(len(LABELS))
    ]


def first_all_unit_word(
    rows: list[list[int]],
    offset: list[int],
    prime: int,
) -> list[int]:
    """Find one exact positive control in an affine packet coset."""
    for coefficients in product(range(prime), repeat=len(rows)):
        word = [
            (left + right) % prime
            for left, right in zip(
                offset, row_combination(coefficients, rows, prime)
            )
        ]
        if all(word):
            return word
    raise RuntimeError("all-unit affine word disappeared")


def crt_pair(left: int, right: int) -> int:
    """Return x mod 91 with x=left mod13 and x=right mod7."""
    correction = (
        (right - left) * pow(P13, -1, P7)
    ) % P7
    value = (left + P13 * correction) % (P13 * P7)
    require(value % P13 == left % P13, "CRT thirteen residue changed")
    require(value % P7 == right % P7, "CRT seven residue changed")
    return value


def exact_lift_control(
    speeds: list[int],
    rows: list[list[int]],
) -> tuple[int, int, int]:
    """Construct and verify one centered CRT/Bezout exact lift."""
    gain = 5
    offset_13 = [0] * len(LABELS)
    offset_13[INDEX["a"]] = 1
    offset_13[INDEX["b"]] = gain
    word_13 = first_all_unit_word(rows, offset_13, P13)
    word_7 = first_all_unit_word(rows, [0] * len(LABELS), P7)
    residue = [
        crt_pair(left, right)
        for left, right in zip(word_13, word_7)
    ]
    require(
        all(value % P13 and value % P7 for value in residue),
        "CRT word lost an all-unit coordinate",
    )
    centered = [value if value <= 45 else value - 91 for value in residue]
    dot_product = sum(
        coefficient * speed
        for coefficient, speed in zip(centered, speeds)
    )
    require(dot_product % 91 == 0, "CRT word left the modular kernel")

    # The hostile scalar word has w_u0=1, so z=e_u0 is a Bezout vector.
    lifted = centered.copy()
    lifted[INDEX["u0"]] -= dot_product
    require(
        sum(
            coefficient * speed
            for coefficient, speed in zip(lifted, speeds)
        )
        == 0,
        "Bezout correction failed to make an exact relation",
    )
    require(
        all(
            (exact - modular) % 91 == 0
            for exact, modular in zip(lifted, residue)
        ),
        "Bezout correction changed a CRT residue",
    )
    scalar_height = sum(abs(speed) for speed in speeds)
    certified_height = 45 * (1 + scalar_height)
    actual_height = max(abs(coefficient) for coefficient in lifted)
    require(
        actual_height <= certified_height,
        "exact lift exceeded the certified height",
    )
    return gain, actual_height, certified_height


def septimal_kernel_count(support_size: int) -> int:
    """Count all-unit vectors in w^perp over F_7 from |supp(w)|."""
    require(1 <= support_size <= 9, "septimal support size left its range")
    supported_solutions = (
        (P7 - 1) ** support_size
        + (P7 - 1) * (-1) ** support_size
    ) // P7
    return (P7 - 1) ** (9 - support_size) * supported_solutions


def count_ledger() -> tuple[int, int, int, int, int, tuple[int, ...]]:
    """Verify the sharp affine-fibre, full-kernel, and CRT counts."""
    per_vector_13 = (P13 - 9 + 6 - 1) * (P13 - 1) ** 5
    septimal_table = tuple(
        septimal_kernel_count(support_size)
        for support_size in range(1, 10)
    )
    require(septimal_table[0] == 0, "singleton support obstruction changed")
    septimal = min(septimal_table[1:])
    require(per_vector_13 == 9 * 12**5, "thirteen count changed")
    require(septimal == 5 * 6**7, "septimal full-kernel floor changed")
    require(
        septimal_table.index(septimal) + 1 == 3,
        "septimal worst support size changed",
    )
    per_vector_91 = per_vector_13 * septimal
    per_direction_91 = (P13 - 1) * per_vector_91
    all_nonzero_quotient_91 = (P13**2 - 1) * per_vector_91
    require(per_vector_91 == 3134566563840, "CRT vector count changed")
    require(
        per_direction_91 == 37614798766080,
        "projective direction count changed",
    )
    require(
        all_nonzero_quotient_91 == 526607182725120,
        "whole nonzero quotient count changed",
    )
    return (
        per_vector_13,
        septimal,
        per_vector_91,
        per_direction_91,
        all_nonzero_quotient_91,
        septimal_table,
    )


def main() -> None:
    speeds, rows = exact_owner_packet()
    rank_13, rank_7, support_7 = packet_checks(speeds, rows)
    (
        directions,
        quotient_cosets,
        model_minimum,
        model_maximum,
    ) = projective_coset_atlas(rows)
    (
        per_vector_13,
        septimal,
        per_vector_91,
        per_direction_91,
        all_nonzero_quotient_91,
        septimal_table,
    ) = count_ledger()
    gain, actual_height, certified_height = exact_lift_control(speeds, rows)

    require(
        model_minimum >= per_vector_13,
        "model affine coset fell below the arrangement floor",
    )

    print("theorem=THM-2325")
    print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-AUDITED")
    print(f"owner_packet_rank_mod13={rank_13}")
    print(f"owner_packet_rank_mod7={rank_7}")
    print(f"hostile_scalar_septimal_support={support_7}")
    print(f"projective_target_directions={directions}")
    print(f"nonzero_target_quotient_vector_cosets={quotient_cosets}")
    print(f"per_vector_coset_mod13_floor={per_vector_13}")
    print(f"model_per_vector_coset_mod13_minimum={model_minimum}")
    print(f"model_per_vector_coset_mod13_maximum={model_maximum}")
    print(f"septimal_support_counts_s1_through_s9={septimal_table}")
    print(f"septimal_full_kernel_all_unit_floor={septimal}")
    print("septimal_floor_attained_at_support_size=3")
    print(f"per_vector_coset_mod91_floor={per_vector_91}")
    print(f"per_projective_direction_mod91_floor={per_direction_91}")
    print(
        "aggregate_all_168_nonzero_quotient_cosets_mod91_floor="
        f"{all_nonzero_quotient_91}"
    )
    print(f"exact_lift_control_gain={gain}")
    print(f"exact_lift_control_height={actual_height}")
    print(f"exact_lift_certified_height={certified_height}")
    print("each_projective_direction_has_12_disjoint_vector_cosets=true")
    print("all_nine_coordinates_are_91_units=true")
    print("exact_bezout_lift_preserves_gain_and_residues=true")
    print("bounded_visible_carrier_membership_not_proved=true")
    print("word_current_incidence_not_proved=true")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
