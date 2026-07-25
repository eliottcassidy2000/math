#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2333."""

from fractions import Fraction
from itertools import product


P = 13
GROUP = tuple(product(range(P), repeat=2))
ZERO = (0, 0)
GROUP_SIZE = P**2


def require(condition: bool, message: str) -> None:
    """Raise in ordinary and optimized Python when a check fails."""
    if not condition:
        raise RuntimeError(message)


def add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    """Add in F_13^2."""
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def sub(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    """Subtract in F_13^2."""
    return ((left[0] - right[0]) % P, (left[1] - right[1]) % P)


def dot(left: tuple[int, int], right: tuple[int, int]) -> int:
    """Dot product in F_13."""
    return (left[0] * right[0] + left[1] * right[1]) % P


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


def owner_packet_quotient_check() -> tuple[int, int]:
    """Reconstruct THM-2309's rank-six packet and two-target quotient."""
    speeds = [13, 13**3, 2 * 13**5, 1, 14, 27, 40, 53, 66]
    omitted = 3
    pivot = (0, 4, 5, 6, 7, 8)
    rows = []
    for coordinate in pivot:
        row = [0] * 9
        row[omitted] = speeds[coordinate]
        row[coordinate] = -speeds[omitted]
        if coordinate == 4:
            row[omitted] += speeds[1]
            row[1] -= speeds[omitted]
        if coordinate == 5:
            row[omitted] += speeds[2]
            row[2] -= speeds[omitted]
        rows.append(row)
    axis_a = [0] * 9
    axis_b = [0] * 9
    axis_a[1] = 1
    axis_b[2] = 1
    packet_rank = matrix_rank(rows, P)
    kernel_rank = matrix_rank(rows + [axis_a, axis_b], P)
    require(packet_rank == 6, "owner packet rank changed")
    require(kernel_rank == 8, "two-target quotient changed")
    return packet_rank, kernel_rank - packet_rank


def orthogonality_atlas() -> tuple[int, int]:
    """Check the exact exponent distribution behind character projection."""
    zero_rows = 0
    balanced_rows = 0
    for difference in GROUP:
        counts = [0] * P
        for character in GROUP:
            counts[dot(character, difference)] += 1
        if difference == ZERO:
            require(
                counts == [GROUP_SIZE] + [0] * (P - 1),
                "trivial character row changed",
            )
            zero_rows += 1
        else:
            require(
                counts == [P] * P,
                "nontrivial character row stopped balancing",
            )
            balanced_rows += 1
    return zero_rows, balanced_rows


def hostile_weights(
    shift: tuple[int, int],
) -> tuple[dict[tuple[int, int], Fraction], dict[tuple[int, int], Fraction]]:
    """Return exact full-support endpoint weights U,V."""
    one_over = Fraction(1, GROUP_SIZE + 1)
    left = {
        point: Fraction(2 if point == ZERO else 1)
        for point in GROUP
    }
    right = {}
    for point in GROUP:
        unshifted = sub(point, shift)
        right[point] = (
            Fraction(GROUP_SIZE, GROUP_SIZE + 1)
            if unshifted == ZERO
            else -one_over
        )
    return left, right


def hostile_fibre_atlas() -> tuple[int, int, int]:
    """Verify full termwise support with only the zero aggregate surviving."""
    shift = (1, 4)
    left, right = hostile_weights(shift)
    require(all(left.values()), "left hostile lost support")
    require(all(right.values()), "right hostile lost support")
    fibre_sums = {}
    term_pairs = {}
    for target in GROUP:
        total = Fraction(0)
        count = 0
        for point in GROUP:
            partner = add(sub(point, target), shift)
            term = left[point] * right[partner]
            require(term != 0, "hostile term vanished")
            total += term
            count += 1
        fibre_sums[target] = total
        term_pairs[target] = count

    require(
        fibre_sums[ZERO] == 1,
        "zero hostile fibre changed",
    )
    require(
        all(
            value == 0
            for target, value in fibre_sums.items()
            if target != ZERO
        ),
        "a nonzero hostile fibre survived",
    )
    require(
        all(count == GROUP_SIZE for count in term_pairs.values()),
        "hostile termwise fibre support changed",
    )
    require(
        sum(fibre_sums.values(), Fraction(0)) == 1,
        "hostile total changed",
    )
    return (
        sum(value != 0 for value in fibre_sums.values()),
        min(term_pairs.values()),
        sum(term_pairs.values()),
    )


packet_rank, quotient_dimension = owner_packet_quotient_check()
zero_character_rows, balanced_character_rows = orthogonality_atlas()
surviving_fibres, terms_per_fibre, total_terms = hostile_fibre_atlas()

require(GROUP_SIZE == 169, "target quotient size changed")
require(quotient_dimension == 2, "target quotient dimension changed")
require(zero_character_rows == 1, "zero orthogonality row changed")
require(balanced_character_rows == 168, "balanced row count changed")
require(surviving_fibres == 1, "hostile survivor count changed")
require(terms_per_fibre == 169, "hostile fibre term count changed")
require(total_terms == 169**2, "hostile total term count changed")

print("theorem=THM-2333")
print("status=PROVED+VERIFIED-EXACT+CANDIDATE-UNDER-INDEPENDENT-AUDIT")
print(f"owner_packet_rank={packet_rank}")
print(f"target_quotient_dimension={quotient_dimension}")
print(f"target_quotient_size={GROUP_SIZE}")
print(f"zero_character_rows={zero_character_rows}")
print(f"balanced_character_rows={balanced_character_rows}")
print(f"hostile_terms_per_fibre={terms_per_fibre}")
print(f"hostile_total_terms={total_terms}")
print(f"hostile_surviving_fibres={surviving_fibres}")
print("hostile_surviving_label=zero")
print("abel_fibre_limits=EXIST")
print("some_fibre_aggregate=NONZERO")
print("nonzero_target_fibre=OPEN")
print("semantic_word_mode=OPEN")
print("all91_unit_aggregate=OPEN")
print("visible_carrier=OPEN")
print("terminal_phase=OPEN")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
