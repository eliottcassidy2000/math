#!/usr/bin/env python3
"""Exact companion for THM-3468's Rule 30 carrier boundaries.

This file deliberately separates four finite objects attached to the
distinguished single-seed center column:

* two independent packed-row generators;
* Boolean query functions of the binary time index;
* ordinary (unary-time) Hankel matrices; and
* low/high binary-address flattenings.

The checks are finite certificates only.  They are not extrapolated to any
of the three Rule 30 prize questions.
"""

from __future__ import annotations

import hashlib
from functools import lru_cache


CENTER_HORIZON = 1 << 18
CENTER_SHA256 = "75be8407c265798fa046baa3ba0f9336e4cbe27bff0be21aeba3e87a7681bea4"
HANKEL_SIZES = tuple(1 << exponent for exponent in range(4, 14))
HANKEL_RANKS = (15, 31, 64, 128, 255, 511, 1023, 2047, 4096, 8192)
FLATTENING_RANKS_GF2 = (
    2,
    4,
    8,
    16,
    32,
    64,
    128,
    256,
    510,
    256,
    128,
    64,
    32,
    16,
    8,
    4,
    2,
)
GREEN_H1_CAP = 4095
GREEN_RECURRENCE_M_CAP = 127
RADIAL_RESTART_CAP = 80
INNOVATION_DEPTH_CAP = 30
CONE_CYLINDER_CASES = ((4, 8), (8, 32), (12, 64), (16, 256))


def check(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


def center_from_right_edge(horizon: int) -> bytearray:
    """Pack each row inward from its right edge, as in THM-3458."""

    center = bytearray(horizon)
    row = 1
    for time in range(horizon):
        center[time] = (row >> time) & 1
        row ^= (row << 1) | (row << 2)
    return center


def center_from_left_edge(horizon: int) -> bytearray:
    """Independent opposite-edge polynomial recurrence."""

    center = bytearray(horizon)
    row = 1
    for time in range(horizon):
        center[time] = (row >> time) & 1
        row = (row << 2) ^ (row << 1) ^ row ^ ((row << 1) & row)
    return center


def rule30(left: int, center: int, right: int) -> int:
    return left ^ center ^ right ^ (center & right)


def centered_rows(horizon: int) -> list[frozenset[int]]:
    """Direct centered-coordinate rows, independent of both packed charts."""

    rows: list[frozenset[int]] = []
    row = frozenset({0})
    for time in range(horizon + 1):
        rows.append(row)
        following = {
            site
            for site in range(-time - 1, time + 2)
            if rule30(
                int(site - 1 in row),
                int(site in row),
                int(site + 1 in row),
            )
        }
        row = frozenset(following)
    return rows


def centered_green_powers(cap: int) -> list[int]:
    """Translate Q^n, Q=y^-1+1+y, to (1+x+x^2)^n."""

    powers = [1]
    for _ in range(cap):
        current = powers[-1]
        powers.append(current ^ (current << 1) ^ (current << 2))
    return powers


def centered_green(powers: list[int], distance: int, time: int) -> int:
    """H_distance(time)=[y^distance](y^-1+1+y)^time."""

    if time < 0 or distance < 0 or distance > time:
        return 0
    return (powers[time] >> (time + distance)) & 1


def two_adic_valuation(value: int) -> int:
    check(value > 0, "2-adic valuation input is positive")
    return (value & -value).bit_length() - 1


def audit_radial_green_restart() -> tuple[int, int, int, tuple[int, ...]]:
    powers = centered_green_powers(GREEN_H1_CAP)

    for middle in range(GREEN_RECURRENCE_M_CAP + 1):
        for half_distance in range(GREEN_RECURRENCE_M_CAP + 2):
            check(
                centered_green(powers, 2 * half_distance, 2 * middle)
                == centered_green(powers, half_distance, middle),
                f"H even/even m={middle} e={half_distance}",
            )
            check(
                centered_green(powers, 2 * half_distance + 1, 2 * middle) == 0,
                f"H odd/even m={middle} e={half_distance}",
            )
            check(
                centered_green(powers, 2 * half_distance, 2 * middle + 1)
                == centered_green(powers, half_distance, middle),
                f"H even/odd m={middle} e={half_distance}",
            )
            check(
                centered_green(powers, 2 * half_distance + 1, 2 * middle + 1)
                == (
                    centered_green(powers, half_distance, middle)
                    ^ centered_green(powers, half_distance + 1, middle)
                ),
                f"H odd/odd m={middle} e={half_distance}",
            )

    h1_word = []
    for time in range(GREEN_H1_CAP + 1):
        observed = centered_green(powers, 1, time)
        expected = two_adic_valuation(time + 1) & 1
        check(observed == expected, f"H1 valuation parity n={time}")
        if time < 32:
            h1_word.append(observed)

    for distance in range(RADIAL_RESTART_CAP + 1):
        check(
            centered_green(powers, distance, distance) == 1,
            f"radial shell first arrival d={distance}",
        )
        for time in range(distance):
            check(
                centered_green(powers, distance, time) == 0,
                f"radial shell premature arrival d={distance} n={time}",
            )

    rows = centered_rows(RADIAL_RESTART_CAP)
    check(rows[2] == frozenset({-2, -1, 2}), "centered row-two restart support")
    for time in range(2, RADIAL_RESTART_CAP + 1):
        reconstructed = centered_green(powers, 1, time - 2)
        for source_time in range(2, time):
            propagation_time = time - 1 - source_time
            source_row = rows[source_time]
            for site in range(-source_time, source_time):
                if site in source_row and site + 1 in source_row:
                    reconstructed ^= centered_green(
                        powers,
                        abs(site),
                        propagation_time,
                    )
        check(
            reconstructed == int(0 in rows[time]),
            f"radial restart center t={time}",
        )

    return (
        GREEN_RECURRENCE_M_CAP,
        GREEN_H1_CAP,
        RADIAL_RESTART_CAP,
        tuple(h1_word),
    )


def packed_edge_rows(horizon: int) -> list[int]:
    rows = []
    row = 1
    for _ in range(horizon + 1):
        rows.append(row)
        row ^= (row << 1) | (row << 2)
    return rows


def packed_phi_width(row: int, width: int) -> int:
    return (row ^ ((row << 1) | (row << 2))) & ((1 << width) - 1)


def seed_period(width: int) -> int:
    row = 1
    period = 0
    while True:
        row = packed_phi_width(row, width)
        period += 1
        if row == 1:
            return period


def period_lift_bit(width: int, period: int) -> int:
    row = 1
    for _ in range(period):
        row = packed_phi_width(row, width + 1)
    return (row >> width) & 1


def periodic_edge_bit(
    rows: list[int],
    periods: list[int],
    depth: int,
    time: int,
) -> int:
    period = periods[depth]
    return (rows[time % period] >> depth) & 1


def arrival_trace_bit(
    rows: list[int],
    periods: list[int],
    phase: int,
    depth: int,
) -> int:
    return periodic_edge_bit(rows, periods, depth, phase + depth)


def audit_innovation_cocycle_and_cylinders(
    center: bytearray,
) -> tuple[str, list[tuple[int, int, int, int]], list[tuple[int, int]]]:
    periods = [seed_period(width) for width in range(1, INNOVATION_DEPTH_CAP + 2)]
    lift_bits = [
        period_lift_bit(width, periods[width - 1])
        for width in range(1, INNOVATION_DEPTH_CAP + 1)
    ]
    lift_word = "".join(str(bit) for bit in lift_bits)
    check(
        lift_word == "101101101000001100000001101010",
        "frozen innovation word eps1..30",
    )

    maximum_period = max(periods)
    maximum_cylinder_q = max(q for _, q in CONE_CYLINDER_CASES)
    rows = packed_edge_rows(max(2 * maximum_period + INNOVATION_DEPTH_CAP, maximum_cylinder_q))

    for depth in range(INNOVATION_DEPTH_CAP + 1):
        period = periods[depth]
        for time in range(period):
            check(
                ((rows[time] >> depth) & 1)
                == ((rows[time + period] >> depth) & 1),
                f"fixed edge-column period k={depth} t={time}",
            )

    for depth in range(2, INNOVATION_DEPTH_CAP + 1):
        lower_period = periods[depth - 1]
        full_period = periods[depth]
        epsilon = lift_bits[depth - 1]
        current_parity = 0
        for phase in range(full_period):
            left = arrival_trace_bit(rows, periods, phase + 1, depth)
            right = arrival_trace_bit(rows, periods, phase, depth)
            parent_current = (
                arrival_trace_bit(rows, periods, phase + 1, depth - 1)
                | arrival_trace_bit(rows, periods, phase + 2, depth - 2)
            )
            check(
                (left ^ right) == parent_current,
                f"arrival transition current k={depth} h={phase}",
            )

        for phase in range(lower_period):
            current = (
                arrival_trace_bit(rows, periods, phase + 1, depth - 1)
                | arrival_trace_bit(rows, periods, phase + 2, depth - 2)
            )
            current_parity ^= current
            check(
                arrival_trace_bit(rows, periods, phase + lower_period, depth)
                == (arrival_trace_bit(rows, periods, phase, depth) ^ epsilon),
                f"innovation cocycle k={depth} h={phase}",
            )
            check(
                current
                == (
                    arrival_trace_bit(rows, periods, phase + lower_period + 1, depth - 1)
                    | arrival_trace_bit(rows, periods, phase + lower_period + 2, depth - 2)
                ),
                f"current period P_k k={depth} h={phase}",
            )
        check(current_parity == epsilon, f"current parity epsilon k={depth}")

    expected_motifs = {
        5: (8, 6, 3),
        8: (32, 24, 16),
        17: (256, 192, 128),
        26: (1024, 768, 516),
    }
    observed_motif_depths = {
        depth
        for depth in range(3, INNOVATION_DEPTH_CAP + 1)
        if lift_bits[depth - 3 : depth] == [1, 1, 0]
    }
    check(
        observed_motif_depths == set(expected_motifs),
        "complete innovation motif-110 depths through cap",
    )
    motif_rows = []
    for depth, expected in expected_motifs.items():
        check(
            lift_bits[depth - 3 : depth] == [1, 1, 0],
            f"innovation motif 110 ending k={depth}",
        )
        period = periods[depth - 1]
        current_weight = sum(
            arrival_trace_bit(rows, periods, phase + 1, depth - 1)
            | arrival_trace_bit(rows, periods, phase + 2, depth - 2)
            for phase in range(period)
        )
        trace_weight = sum(
            arrival_trace_bit(rows, periods, phase, depth)
            for phase in range(period)
        )
        observed = (period, current_weight, trace_weight)
        check(observed == expected, f"110 motif census k={depth}")
        check(
            3 * period <= 8 * trace_weight <= 5 * period,
            f"110 motif discrepancy window k={depth}",
        )
        motif_rows.append((depth, period, current_weight, trace_weight))

    cylinder_rows = []
    for cutoff, q in CONE_CYLINDER_CASES:
        check(q > cutoff, f"cylinder q exceeds cutoff K={cutoff}")
        check(q % periods[cutoff] == 0, f"cylinder q is multiple of P_(K+1) K={cutoff}")
        phase = -q
        for depth in range(cutoff + 1):
            check(
                arrival_trace_bit(rows, periods, phase, depth) == center[depth],
                f"cylinder center prefix K={cutoff} q={q} k={depth}",
            )
        cone_block = [
            (rows[phase + depth] >> depth) & 1
            for depth in range(q, 2 * q + 1)
        ]
        check(
            cone_block == [0] * q + [1],
            f"cylinder light-cone block K={cutoff} q={q}",
        )
        cylinder_rows.append((cutoff, q))

    return lift_word, motif_rows, cylinder_rows


def anf_degree(values: bytes) -> int:
    """Degree of a Boolean truth table in increasing mask order."""

    size = len(values)
    variables = size.bit_length() - 1
    check(size == 1 << variables, "ANF truth-table size is a power of two")
    coefficients = bytearray(values)
    for bit in range(variables):
        stride = 1 << bit
        for base in range(0, size, 2 * stride):
            for offset in range(stride):
                coefficients[base + stride + offset] ^= coefficients[base + offset]
    return max(
        (mask.bit_count() for mask, value in enumerate(coefficients) if value),
        default=-1,
    )


def restrict_truth_table(values: bytes, bit: int, value: int) -> bytes:
    variables = len(values).bit_length() - 1
    check(0 <= bit < variables, "restriction bit range")
    check(value in (0, 1), "restriction value is Boolean")
    return bytes(
        values[mask]
        for mask in range(1 << variables)
        if ((mask >> bit) & 1) == value
    )


@lru_cache(maxsize=None)
def has_full_query_depth(values: bytes) -> bool:
    """Certify deterministic decision-tree depth equal to all variables.

    For an n-variable nonconstant function f, D(f)=n exactly when, for every
    possible first query i, at least one of f|x_i=0 and f|x_i=1 has depth
    n-1.  Full ANF degree supplies a cheap terminal certificate.
    """

    variables = len(values).bit_length() - 1
    if all(value == values[0] for value in values):
        return False
    if anf_degree(values) == variables:
        return True
    for bit in range(variables):
        zero = restrict_truth_table(values, bit, 0)
        one = restrict_truth_table(values, bit, 1)
        if not (has_full_query_depth(zero) or has_full_query_depth(one)):
            return False
    return True


def audit_fixed_seed_queries(center: bytearray) -> tuple[list[int], list[int], list[int]]:
    padded_degrees = []
    padded_depths = []
    for variables in range(2, 19):
        values = bytes(center[: 1 << variables])
        degree = anf_degree(values)
        check(has_full_query_depth(values), f"padded query depth m={variables}")
        padded_degrees.append(degree)
        padded_depths.append(variables)

    standard_degrees = []
    standard_depths = []
    for length in range(5, 17):
        free_bits = length - 1
        values = bytes(center[1 << free_bits : 1 << length])
        degree = anf_degree(values)
        check(has_full_query_depth(values), f"standard positive query depth m={length}")
        standard_degrees.append(degree)
        standard_depths.append(free_bits)

    check(
        padded_degrees == [2, 3, 4, 4, 6, 7, 8, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18],
        "padded ANF degrees",
    )
    check(
        standard_degrees == [4, 5, 5, 6, 8, 8, 9, 10, 11, 13, 13, 14],
        "standard positive ANF degrees",
    )
    return padded_degrees, padded_depths, standard_depths


def pack_center_bits(center: bytearray, count: int) -> int:
    packed = 0
    for index in range(count):
        packed |= int(center[index]) << index
    return packed


def gf2_rank_high_pivot(rows: list[int]) -> int:
    basis: dict[int, int] = {}
    for original in rows:
        row = original
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def gf2_rank_low_pivot(rows: list[int], width: int) -> int:
    work = rows[:]
    rank = 0
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(work)) if (work[index] >> column) & 1),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_row = work[rank]
        for index in range(rank + 1, len(work)):
            if (work[index] >> column) & 1:
                work[index] ^= pivot_row
        rank += 1
        if rank == len(work):
            break
    return rank


def hankel_rows(packed_center: int, start: int, size: int) -> list[int]:
    mask = (1 << size) - 1
    return [((packed_center >> (start + row)) & mask) for row in range(size)]


def audit_hankel_ranks(center: bytearray) -> list[int]:
    packed = pack_center_bits(center, 2 * HANKEL_SIZES[-1])
    ranks = []
    for size in HANKEL_SIZES:
        rows = hankel_rows(packed, 0, size)
        high_rank = gf2_rank_high_pivot(rows)
        low_rank = gf2_rank_low_pivot(rows, size)
        check(high_rank == low_rank, f"independent Hankel ranks d={size}")
        ranks.append(high_rank)
    check(tuple(ranks) == HANKEL_RANKS, "frozen ordinary Hankel rank list")
    return ranks


def flattening_rows(center: bytearray, low_bits: int, total_bits: int = 18) -> list[int]:
    low_count = 1 << low_bits
    high_count = 1 << (total_bits - low_bits)
    rows = []
    for residue in range(low_count):
        row = 0
        for quotient in range(high_count):
            row |= int(center[residue + low_count * quotient]) << quotient
        rows.append(row)
    return rows


def ternary_add(
    left: tuple[int, int],
    right: tuple[int, int],
    mask: int,
) -> tuple[int, int]:
    """Add two vectors over F_3 encoded by disjoint one/two bitsets."""

    left_one, left_two = left
    right_one, right_two = right
    left_zero = mask ^ left_one ^ left_two
    right_zero = mask ^ right_one ^ right_two
    result_one = (
        (left_zero & right_one)
        | (left_one & right_zero)
        | (left_two & right_two)
    )
    result_two = (
        (left_zero & right_two)
        | (left_two & right_zero)
        | (left_one & right_one)
    )
    check((result_one & result_two) == 0, "ternary bitsets remain disjoint")
    return result_one, result_two


def gf3_rank_and_square_determinant(binary_rows: list[int], width: int) -> tuple[int, int]:
    """Exact bit-parallel Gaussian elimination over F_3."""

    mask = (1 << width) - 1
    work = [(row & mask, 0) for row in binary_rows]
    rank = 0
    determinant = 1
    for column in range(width):
        pivot = next(
            (
                index
                for index in range(rank, len(work))
                if ((work[index][0] | work[index][1]) >> column) & 1
            ),
            None,
        )
        if pivot is None:
            continue
        if pivot != rank:
            work[rank], work[pivot] = work[pivot], work[rank]
            determinant = (-determinant) % 3

        pivot_one, pivot_two = work[rank]
        if (pivot_two >> column) & 1:
            determinant = (2 * determinant) % 3
            pivot_one, pivot_two = pivot_two, pivot_one
            work[rank] = (pivot_one, pivot_two)

        for index in range(rank + 1, len(work)):
            one, two = work[index]
            if (one >> column) & 1:
                work[index] = ternary_add((one, two), (pivot_two, pivot_one), mask)
            elif (two >> column) & 1:
                work[index] = ternary_add((one, two), (pivot_one, pivot_two), mask)
        rank += 1
        if rank == len(work):
            break

    if len(work) != width or rank != width:
        determinant = 0
    return rank, determinant


def audit_flattening_ranks(
    center: bytearray,
) -> tuple[list[int], tuple[int, int], tuple[int, int, int]]:
    ranks = []
    central_rows: list[int] | None = None
    for low_bits in range(1, 18):
        rows = flattening_rows(center, low_bits)
        rank = gf2_rank_high_pivot(rows)
        ranks.append(rank)
        if low_bits == 9:
            central_rows = rows
    check(tuple(ranks) == FLATTENING_RANKS_GF2, "18-bit GF2 flattening ranks")
    check(central_rows is not None, "central flattening retained")
    central_gf3 = gf3_rank_and_square_determinant(central_rows, 1 << 9)
    check(central_gf3 == (512, 1), "central flattening rank and determinant mod 3")

    block = 1 << 9
    standard_rows = []
    for quotient in range(block // 2, block):
        row = 0
        for residue in range(block):
            row |= int(center[residue + block * quotient]) << residue
        standard_rows.append(row)
    standard_rank = gf2_rank_high_pivot(standard_rows)
    standard_row_profiles = len(set(standard_rows))

    standard_columns = []
    for residue in range(block):
        column = 0
        for quotient in range(block // 2, block):
            column |= int(center[residue + block * quotient]) << (quotient - block // 2)
        standard_columns.append(column)
    standard_column_profiles = len(set(standard_columns))
    standard = (standard_rank, standard_row_profiles, standard_column_profiles)
    check(standard == (256, 256, 512), "standard 18-bit restricted carrier ranks")
    return ranks, central_gf3, standard


def main() -> None:
    center_right = center_from_right_edge(CENTER_HORIZON)
    center_left = center_from_left_edge(CENTER_HORIZON)
    check(center_right == center_left, "independent center generators through 2^18")
    center_hash = hashlib.sha256(center_right).hexdigest()
    check(center_hash == CENTER_SHA256, "frozen 2^18 center hash")

    green_audit = audit_radial_green_restart()
    lift_word, motif_rows, cylinder_rows = audit_innovation_cocycle_and_cylinders(center_right)
    padded_degrees, padded_depths, standard_depths = audit_fixed_seed_queries(center_right)
    hankel_ranks = audit_hankel_ranks(center_right)
    flattening_ranks, central_gf3, standard_carrier = audit_flattening_ranks(center_right)

    print("THM-3468 RULE 30 RADIAL/INNOVATION/CARRIER EXACT AUDIT")
    print(f"center_generators=right_edge,left_edge agreement_t=0..{CENTER_HORIZON - 1}")
    print(f"center_sha256_0..{CENTER_HORIZON - 1}={center_hash}")
    print(
        "green_recurrence_universe="
        f"m=0..{green_audit[0]},e=0..{green_audit[0] + 1} "
        f"H1_n=0..{green_audit[1]} shell_d=0..{green_audit[2]}"
    )
    print(f"green_H1_prefix_n0..31={''.join(map(str, green_audit[3]))}")
    print(f"radial_restart_center_t=2..{RADIAL_RESTART_CAP}")
    print(f"innovation_lift_bits_eps1..{INNOVATION_DEPTH_CAP}={lift_word}")
    print(f"innovation_110_rows_k_p_qweight_tweight={motif_rows}")
    print(f"cone_cylinder_cases_K_q={cylinder_rows}")
    print(f"padded_query_anf_degrees_m2..18={padded_degrees}")
    print(f"padded_query_depths_m2..18={padded_depths}")
    print(f"standard_positive_query_depths_m5..16={standard_depths}")
    print(f"ordinary_hankel_sizes={list(HANKEL_SIZES)}")
    print(f"ordinary_hankel_ranks_gf2={hankel_ranks}")
    print(f"lookup18_low_high_ranks_gf2={flattening_ranks}")
    print(f"lookup18_central_rank_det_mod3={central_gf3}")
    print(f"lookup18_standard_rank_rows_columns_gf2={standard_carrier}")
    print("scope=FINITE-EXACT; no Rule30 prize conclusion")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
