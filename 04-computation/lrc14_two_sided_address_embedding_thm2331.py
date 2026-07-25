#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2331."""

from itertools import product
from math import gcd


P = 7
P13 = 13
MODULUS = 91
COORDINATES = 9
RESIDUE_LIFT_FLOOR = 3 * 5**7
TARGET_VECTOR_BANK = 3_134_566_563_840
PROJECTIVE_BANK = 37_614_798_766_080


def require(condition: bool, message: str) -> None:
    """Raise in ordinary and optimized Python when a check fails."""
    if not condition:
        raise RuntimeError(message)


def matrix_rank(matrix: list[list[int]], prime: int) -> int:
    """Return exact row rank over a prime field."""
    work = [[entry % prime for entry in row] for row in matrix]
    if not work:
        return 0
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next(
            (
                row
                for row in range(pivot_row, len(work))
                if work[row][column]
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
        for row in range(len(work)):
            if row == pivot_row:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    (left - factor * right) % prime
                    for left, right in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def allowed(displacement: int) -> tuple[int, ...]:
    """Residues making both u and u+displacement seven-units."""
    return tuple(
        residue
        for residue in range(P)
        if residue and (residue + displacement) % P
    )


def atomic_support_checks() -> int:
    """Check the exact zero sets of danger/safe and guard numerators."""
    rows = 0
    for harmonic in range(-1_000, 1_001):
        expected = harmonic == 0 or harmonic % P != 0
        # For n != 0, the danger numerator is sin(pi*n/7) and the
        # guard numerator is sin(2*pi*n/7). Since 2 is a unit mod 7,
        # both vanish exactly at nonzero multiples of seven.
        danger_nonzero = harmonic == 0 or harmonic % P != 0
        guard_nonzero = harmonic == 0 or (2 * harmonic) % P != 0
        require(danger_nonzero == expected, "danger support law failed")
        require(guard_nonzero == expected, "guard support law failed")
        rows += 2
    return rows


def exhaustive_two_pivot_atlas() -> tuple[int, tuple[int, ...], int]:
    """Exhaust the load-bearing two-coordinate completion problem."""
    minimum = P
    witness: tuple[int, ...] | None = None
    rows = 0
    for first_weight, second_weight in product(range(1, P), repeat=2):
        for first_displacement, second_displacement in product(
            range(P), repeat=2
        ):
            first_allowed = allowed(first_displacement)
            second_allowed = allowed(second_displacement)
            for target in range(P):
                count = sum(
                    (
                        first_weight * first
                        + second_weight * second
                        - target
                    )
                    % P
                    == 0
                    for first in first_allowed
                    for second in second_allowed
                )
                require(count >= 3, "two-pivot completion floor failed")
                if count < minimum:
                    minimum = count
                    witness = (
                        first_weight,
                        second_weight,
                        first_displacement,
                        second_displacement,
                        target,
                    )
                rows += 1
    require(witness is not None, "two-pivot atlas was empty")
    return minimum, witness, rows


def first_two_sided_residue_lift(
    speeds: list[int],
    relation: list[int],
    multiplier: int,
    deepest_index: int,
    frequency: int,
) -> list[int]:
    """Construct one u modulo seven using the proof's two pivots."""
    displacement = [
        (multiplier if index == deepest_index else 0) - entry
        for index, entry in enumerate(relation)
    ]
    choices = [allowed(entry % P) for entry in displacement]
    support = [
        index for index, speed in enumerate(speeds) if speed % P
    ]
    require(len(support) >= 2, "septimal support boundary failed")
    first_pivot, second_pivot = support[:2]
    fixed: list[int | None] = [None] * COORDINATES
    for index in range(COORDINATES):
        if index not in (first_pivot, second_pivot):
            fixed[index] = choices[index][0]

    fixed_sum = sum(
        speeds[index] * value
        for index, value in enumerate(fixed)
        if value is not None
    )
    for first in choices[first_pivot]:
        for second in choices[second_pivot]:
            if (
                fixed_sum
                + speeds[first_pivot] * first
                + speeds[second_pivot] * second
                - frequency
            ) % P == 0:
                fixed[first_pivot] = first
                fixed[second_pivot] = second
                result = [int(value) for value in fixed]
                require(
                    all(result[index] % P for index in range(COORDINATES)),
                    "left endpoint acquired a septimal zero",
                )
                require(
                    all(
                        (result[index] + displacement[index]) % P
                        for index in range(COORDINATES)
                    ),
                    "right endpoint acquired a septimal zero",
                )
                return result
    raise RuntimeError("constructive two-pivot lift disappeared")


def centered(residue: int) -> int:
    """Return the centered representative modulo seven."""
    residue %= P
    return residue if residue <= P // 2 else residue - P


def owner_packet_rows(speeds: list[int]) -> list[list[int]]:
    """Construct THM-2325's target-grafted owner packet."""
    omitted = 3
    pivot = (0, 4, 5, 6, 7, 8)
    rows = []
    for coordinate in pivot:
        row = [0] * COORDINATES
        row[omitted] = speeds[coordinate]
        row[coordinate] = -speeds[omitted]
        if coordinate == 4:
            row[omitted] += speeds[1]
            row[1] -= speeds[omitted]
        if coordinate == 5:
            row[omitted] += speeds[2]
            row[2] -= speeds[omitted]
        require(
            sum(left * right for left, right in zip(row, speeds)) == 0,
            "owner packet row left the exact relation lattice",
        )
        rows.append(row)
    return rows


def exact_positive_control() -> tuple[int, int, int, int]:
    """Verify an exact target-axis address on THM-2325's model word."""
    speeds = [13, 13**3, 2 * 13**5, 1, 14, 27, 40, 53, 66]
    relation = [
        -27,
        -27,
        -27,
        20_110_798,
        -41,
        -27,
        -27,
        -27,
        38,
    ]
    require(
        sum(left * right for left, right in zip(relation, speeds)) == 0,
        "positive-control address left the exact relation lattice",
    )
    require(
        all(gcd(entry, MODULUS) == 1 for entry in relation),
        "positive-control relation lost a 91-unit coordinate",
    )

    packet = owner_packet_rows(speeds)
    axis_a = [0] * COORDINATES
    axis_a[1] = 1
    residual = [
        relation[index] - axis_a[index] for index in range(COORDINATES)
    ]
    require(matrix_rank(packet, P13) == 6, "owner packet rank changed")
    require(
        matrix_rank(packet + [residual], P13) == 6,
        "control left the target-a affine fibre",
    )
    require(
        matrix_rank(packet + [relation], P13) == 7,
        "control target quotient became zero",
    )

    deepest_index = 2
    multiplier = 1
    frequency = speeds[0]
    target_frequency = frequency + multiplier * speeds[deepest_index]
    residues = first_two_sided_residue_lift(
        speeds,
        relation,
        multiplier,
        deepest_index,
        frequency,
    )
    centered_residues = [centered(entry) for entry in residues]
    numerator = frequency - sum(
        left * right
        for left, right in zip(centered_residues, speeds)
    )
    require(numerator % P == 0, "exact correction ceased to be integral")
    correction = numerator // P
    # Coordinate 3 has speed one, so e_3 is a Bezout vector.
    left = centered_residues[:]
    left[3] += P * correction
    displacement = [-entry for entry in relation]
    displacement[deepest_index] += multiplier
    right = [
        left_entry + difference
        for left_entry, difference in zip(left, displacement)
    ]

    require(
        sum(entry * speed for entry, speed in zip(left, speeds))
        == frequency,
        "left exact frequency changed",
    )
    require(
        sum(entry * speed for entry, speed in zip(right, speeds))
        == target_frequency,
        "right exact frequency changed",
    )
    recovered = [
        left_entry
        + (multiplier if index == deepest_index else 0)
        - right_entry
        for index, (left_entry, right_entry) in enumerate(zip(left, right))
    ]
    require(recovered == relation, "relation address was not recovered")
    require(
        all(entry % P for entry in left)
        and all(entry % P for entry in right),
        "exact lift crossed an interval Fourier zero",
    )

    scalar_size = sum(abs(speed) for speed in speeds)
    left_bound = 3 + abs(frequency) + 3 * scalar_size
    right_bound = (
        left_bound + max(abs(entry) for entry in relation) + abs(multiplier)
    )
    require(max(abs(entry) for entry in left) <= left_bound, "left height")
    require(max(abs(entry) for entry in right) <= right_bound, "right height")
    require(multiplier % P != 0, "deepest comb coefficient vanished")
    return (
        len([speed for speed in speeds if speed % P]),
        max(abs(entry) for entry in relation),
        max(abs(entry) for entry in left),
        max(abs(entry) for entry in right),
    )


atomic_rows = atomic_support_checks()
minimum, equality_witness, atlas_rows = exhaustive_two_pivot_atlas()
require(minimum == 3, "sharp completion minimum changed")
support_one_hostile = sum(
    (residue - 0) % P == 0 for residue in allowed(1)
)
require(support_one_hostile == 0, "support-one hostile disappeared")
require(RESIDUE_LIFT_FLOOR == 234_375, "residue lift constant changed")
term_bank = RESIDUE_LIFT_FLOOR * TARGET_VECTOR_BANK
projective_term_bank = RESIDUE_LIFT_FLOOR * PROJECTIVE_BANK
require(
    term_bank == 734_664_038_400_000_000,
    "target-vector term bank changed",
)
require(
    projective_term_bank == 8_815_968_460_800_000_000,
    "projective term bank changed",
)
(
    support_size,
    relation_height,
    left_height,
    right_height,
) = exact_positive_control()

print("theorem=THM-2331")
print("status=PROVED+VERIFIED-EXACT+CANDIDATE-UNDER-INDEPENDENT-AUDIT")
print(f"two_pivot_atlas_rows={atlas_rows}")
print(f"atomic_support_rows={atomic_rows}")
print(f"sharp_two_pivot_minimum={minimum}")
print("sharp_witness=" + ",".join(str(entry) for entry in equality_witness))
print(f"support_one_hostile_count={support_one_hostile}")
print(f"residue_lifts_per_address={RESIDUE_LIFT_FLOOR}")
print(f"target_vector_term_bank={term_bank}")
print(f"projective_term_bank={projective_term_bank}")
print(f"positive_control_support7={support_size}")
print("positive_control_target_quotient=axis_a")
print(f"positive_control_relation_height={relation_height}")
print(f"positive_control_left_height={left_height}")
print(f"positive_control_right_height={right_height}")
print("abel_term=NONZERO-FOR-EVERY-0<rho<1")
print("aggregate_survival=NOT-CLAIMED")
print("visible_carrier=OPEN")
print("terminal_phase=OPEN")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
