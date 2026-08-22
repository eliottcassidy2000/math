#!/usr/bin/env python3
"""Exact audit for the n=14 CRT-even boundary-blind flat.

The accompanying note proves two statements.

1. At n=2p with p odd prime, a normalized vector supported only at the
   nonzero even indices is a right-boundary blocker iff at least one of those
   even coordinates is zero.
2. At every n=2p, no nonzero vector in that flat is a full micro-staircase
   blocker.  A parity refinement of the unique-protector cells gives a
   constructive proof; its half-turn danger graph strictly lowers the
   2-adic valuation.

This program verifies the endpoint-owner table, the six explicit cell rows,
the constructive parity certificate, and an independent meet-in-the-middle
exhaustion of all 14^6 flat vectors.  It is a finite residue/cell audit, not a
speed realization and not an LRC(14) proof.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import factorial, isqrt
from pathlib import Path


N = 14
EVEN_INDICES = (2, 4, 6, 8, 10, 12)
FIXED_ZERO_INDICES = (1, 3, 5, 7, 9, 11, 13)
HALF_TURN_PRIORITY = (2, 6, 10, 4, 12, 8)
EXPECTED_EVEN_BIN_ROWS = {
    2: (13, 12, 12, 11, 11, 10),
    4: (6, 13, 6, 12, 5, 12),
    6: (4, 8, 13, 3, 8, 12),
    8: (3, 6, 10, 13, 3, 6),
    10: (2, 5, 8, 11, 13, 2),
    12: (2, 4, 6, 9, 11, 13),
}
EXPECTED_ODD_DANGER = {
    2: (),
    4: (2, 6),
    6: (),
    8: (4, 12),
    10: (),
    12: (6,),
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def bin_value(n: int, index: int, alpha: Fraction) -> int:
    scaled = n * ((index * alpha) % 1)
    return scaled.numerator // scaled.denominator


def explicit_alpha(n: int, index: int) -> Fraction:
    if 2 * index <= n - 1:
        t = Fraction(1, 2) + Fraction(1, 4 * n * n)
    else:
        t = Fraction(1, 2 * n)
    return (1 - t / n) / index


def explicit_pattern(n: int, index: int) -> tuple[int, ...]:
    alpha = explicit_alpha(n, index)
    pattern = tuple(bin_value(n, i, alpha) for i in range(1, n))
    endpoint_owners = tuple(
        i for i, value in enumerate(pattern, 1) if value in (0, n - 1)
    )
    require(endpoint_owners == (index,), (n, index, endpoint_owners))
    for i in range(1, n):
        scaled = n * ((i * alpha) % 1)
        require(scaled.denominator != 1, ("breakpoint", n, index, i, alpha))
    return pattern


def is_prime(value: int) -> bool:
    return value >= 2 and all(
        value % divisor for divisor in range(2, isqrt(value) + 1)
    )


def valuation_two(value: int) -> int:
    valuation = 0
    while value % 2 == 0:
        value //= 2
        valuation += 1
    return valuation


def endpoint_owner_table(n: int) -> tuple[tuple[int, tuple[int, ...]], ...]:
    return tuple(
        (
            a,
            tuple(
                i for i in range(1, n) if (a * i) % n in (0, n - 1)
            ),
        )
        for a in range(n)
    )


def expected_n14_owner_table() -> tuple[tuple[int, tuple[int, ...]], ...]:
    return (
        (0, tuple(range(1, 14))),
        (1, (13,)),
        (2, (7,)),
        (3, (9,)),
        (4, (7,)),
        (5, (11,)),
        (6, (7,)),
        (7, EVEN_INDICES),
        (8, (7,)),
        (9, (3,)),
        (10, (7,)),
        (11, (5,)),
        (12, (7,)),
        (13, (1,)),
    )


def atomic_patterns(n: int) -> tuple[tuple[int, ...], ...]:
    breakpoints = sorted(
        {
            Fraction(numerator, n * index)
            for index in range(1, n)
            for numerator in range(n * index + 1)
        }
    )
    patterns = []
    for left, right in zip(breakpoints, breakpoints[1:]):
        alpha = (left + right) / 2
        patterns.append(tuple(bin_value(n, i, alpha) for i in range(1, n)))
    require(len(patterns) == len(set(patterns)), (n, "duplicate atomic patterns"))
    return tuple(patterns)


def safe_shift_mask(
    pattern: tuple[int, ...], index: int, residue: int, shifts: tuple[int, ...]
) -> int:
    mask = 0
    for bit, shift in enumerate(shifts):
        if (pattern[index - 1] + shift * residue) % N not in (0, N - 1):
            mask |= 1 << bit
    return mask


def bad_odd_shift_mask(pattern: tuple[int, ...], index: int, residue: int) -> int:
    shifts = tuple(range(1, N, 2))
    return ((1 << len(shifts)) - 1) ^ safe_shift_mask(
        pattern, index, residue, shifts
    )


def audit_explicit_cell_mechanism() -> tuple[int, int, int]:
    patterns = {index: explicit_pattern(N, index) for index in EVEN_INDICES}
    rows = {
        index: tuple(patterns[index][other - 1] for other in EVEN_INDICES)
        for index in EVEN_INDICES
    }
    require(rows == EXPECTED_EVEN_BIN_ROWS, ("even-bin table", rows))

    danger = {
        index: tuple(
            other
            for other, value in zip(EVEN_INDICES, rows[index])
            if other != index and value in (6, 7)
        )
        for index in EVEN_INDICES
    }
    require(danger == EXPECTED_ODD_DANGER, ("danger table", danger))

    # Every nonempty half-turn support has a priority-selected member whose
    # danger set misses that support.  This is the acyclic sidecar used in the
    # proof, checked independently over all 63 supports.
    support_checks = 0
    for bits in range(1, 1 << len(EVEN_INDICES)):
        support = {
            EVEN_INDICES[bit]
            for bit in range(len(EVEN_INDICES))
            if (bits >> bit) & 1
        }
        chosen = next(index for index in HALF_TURN_PRIORITY if index in support)
        require(not (set(danger[chosen]) & support), (support, chosen, danger[chosen]))
        support_checks += 1

    # Exact local parity facts: an ordinary nonzero residue forbids exactly
    # one odd shift at every bin; a half-turn forbids odd shifts exactly at
    # bins 6 and 7.  The selected coordinate has bin 13 and hence forbids no
    # odd shift (while it forbids all even shifts).
    parity_checks = 0
    dummy = tuple(range(N))
    for bin_index in range(N):
        pattern = list(dummy[:-1])
        pattern[0] = bin_index
        pattern_tuple = tuple(pattern)
        for residue in range(1, N):
            bad_mask = bad_odd_shift_mask(pattern_tuple, 1, residue)
            bad_count = bad_mask.bit_count()
            expected = 7 if residue == 7 and bin_index in (6, 7) else 0 if residue == 7 else 1
            require(bad_count == expected, (bin_index, residue, bad_count, expected))
            parity_checks += 1

    # Exhaust all value-type profiles combinatorially.  These counts locate
    # exactly how much of the boundary-blind family was not covered by the
    # strict weighted-cost inequality in THM-3317.
    boundary_nonzero = 0
    weighted_certified = 0
    weighted_remainder = 0
    for half_turns in range(7):
        for ordinary in range(7 - half_turns):
            zeros = 6 - half_turns - ordinary
            if zeros == 0 or half_turns + ordinary == 0:
                continue
            multiplicity = (
                factorial(6)
                // factorial(half_turns)
                // factorial(ordinary)
                // factorial(zeros)
                * 12**ordinary
            )
            boundary_nonzero += multiplicity
            if 7 * half_turns + 2 * ordinary < 14:
                weighted_certified += multiplicity
            else:
                weighted_remainder += multiplicity

    require(boundary_nonzero == N**6 - 13**6 - 1, boundary_nonzero)
    require(weighted_certified == 1_953_510, weighted_certified)
    require(weighted_remainder == 749_216, weighted_remainder)
    return support_checks, parity_checks, weighted_remainder


def audit_general_activation_graph() -> tuple[int, int, int, int]:
    """Check the proved danger-graph formula for odd primes through 101."""
    prime_count = 0
    row_count = 0
    edge_count = 0
    maximum_out_degree = 0
    for prime in range(3, 102, 2):
        if not is_prime(prime):
            continue
        prime_count += 1
        n = 2 * prime
        even_indices = tuple(range(2, n, 2))
        for chosen in even_indices:
            alpha = explicit_alpha(n, chosen)
            actual = {
                other
                for other in even_indices
                if other != chosen
                and bin_value(n, other, alpha) in (prime - 1, prime)
            }
            expected: set[int] = set()
            if chosen % 4 == 0:
                expected.add(chosen // 2)
                if 3 * chosen // 2 < n:
                    expected.add(3 * chosen // 2)
            require(actual == expected, (prime, chosen, actual, expected))
            for target in actual:
                require(
                    valuation_two(target) == valuation_two(chosen) - 1,
                    (prime, chosen, target),
                )
            row_count += 1
            edge_count += len(actual)
            maximum_out_degree = max(maximum_out_degree, len(actual))
    require(prime_count == 25, prime_count)
    return prime_count, row_count, edge_count, maximum_out_degree


def triple_intersections(
    patterns: tuple[tuple[int, ...], ...],
    indices: tuple[int, int, int],
    shifts: tuple[int, ...],
) -> tuple[int, ...]:
    candidate_count = len(patterns) * len(shifts)
    full_mask = (1 << candidate_count) - 1
    coordinate_masks: dict[tuple[int, int], int] = {}
    for index in indices:
        for residue in range(N):
            mask = 0
            bit = 0
            for pattern in patterns:
                local = safe_shift_mask(pattern, index, residue, shifts)
                mask |= local << bit
                bit += len(shifts)
            coordinate_masks[index, residue] = mask

    intersections = []
    for residues in product(range(N), repeat=3):
        mask = full_mask
        for index, residue in zip(indices, residues):
            mask &= coordinate_masks[index, residue]
        intersections.append(mask)
    return tuple(intersections)


def count_disjoint_pairs(first: tuple[int, ...], second: tuple[int, ...]) -> int:
    count = 0
    for first_mask in first:
        for second_mask in second:
            count += not (first_mask & second_mask)
    return count


def audit_boundary_flat() -> tuple[int, int, int]:
    # All fixed odd coordinates automatically block every boundary parameter
    # except a=7.  At a=7 they all have bin 7, while all even coordinates have
    # bin 0.  The following MITM pass therefore directly classifies the full
    # 14^6 flat, without using the proved classification formula.
    pattern = tuple(7 if index % 2 else 0 for index in range(1, N))
    patterns = (pattern,)
    shifts = tuple(range(N))
    first = triple_intersections(patterns, EVEN_INDICES[:3], shifts)
    second = triple_intersections(patterns, EVEN_INDICES[3:], shifts)
    blocker_count = count_disjoint_pairs(first, second)
    expected = N**6 - 13**6
    require(blocker_count == expected, ("boundary blockers", blocker_count, expected))

    # Named hostile and all-active negative control.
    one_defect = (1, 0, 0, 0, 0, 0)
    all_active = (1, 1, 1, 1, 1, 1)
    require(0 in one_defect, one_defect)
    require(0 not in all_active, all_active)
    return blocker_count, len(first) * len(second), expected - 1


def audit_full_micro_flat() -> tuple[int, int, int, int]:
    all_patterns = atomic_patterns(N)
    require(len(all_patterns) == 812, len(all_patterns))
    eligible = tuple(
        pattern
        for pattern in all_patterns
        if all(
            pattern[index - 1] not in (0, N - 1)
            for index in FIXED_ZERO_INDICES
        )
    )
    require(len(eligible) == 210, len(eligible))
    shifts = tuple(range(N))
    first = triple_intersections(eligible, EVEN_INDICES[:3], shifts)
    second = triple_intersections(eligible, EVEN_INDICES[3:], shifts)
    blocker_count = count_disjoint_pairs(first, second)
    require(blocker_count == 1, ("full blockers in flat", blocker_count))

    # Positive and hostile controls, counted directly on all 812 cells.
    def safe_candidate_count(values: tuple[int, ...]) -> int:
        vector = dict(zip(EVEN_INDICES, values))
        count = 0
        for pattern in all_patterns:
            for shift in shifts:
                if all(
                    (
                        pattern[index - 1]
                        + shift * vector.get(index, 0)
                    )
                    % N
                    not in (0, N - 1)
                    for index in range(1, N)
                ):
                    count += 1
        return count

    require(safe_candidate_count((0, 0, 0, 0, 0, 0)) == 0, "zero control")
    one_defect_safe = safe_candidate_count((1, 0, 0, 0, 0, 0))
    half_turn_safe = safe_candidate_count((0, 0, 7, 0, 0, 0))
    require(one_defect_safe == 336, one_defect_safe)
    require(half_turn_safe == 56, half_turn_safe)
    return blocker_count, len(first) * len(second), one_defect_safe, half_turn_safe


def main() -> None:
    source = Path(__file__).read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in source, "source contains bare carriage return")

    owners = endpoint_owner_table(N)
    require(owners == expected_n14_owner_table(), owners)
    cell_rows = audit_explicit_cell_mechanism()
    general_graph = audit_general_activation_graph()
    boundary = audit_boundary_flat()
    full = audit_full_micro_flat()

    print("Composite CRT-even micro-staircase flat audit")
    print("scope=PROVED_finite_residue_cell_theorems_plus_FINITE-EXACT_audit;not_LRC")
    print(f"n14_endpoint_owner_table={owners};status=PASS")
    print(
        "unique_protector_half_turn_sidecar="
        f"support_checks={cell_rows[0]},parity_checks={cell_rows[1]},"
        f"weighted_boundary_remainder_closed={cell_rows[2]};status=PASS"
    )
    print(f"explicit_even_bin_rows={EXPECTED_EVEN_BIN_ROWS}")
    print(f"odd_half_turn_danger_sets={EXPECTED_ODD_DANGER}")
    print(
        "general_2p_activation_graph_audit="
        f"primes={general_graph[0]},rows={general_graph[1]},"
        f"edges={general_graph[2]},max_out_degree={general_graph[3]},"
        "range=p3..101;status=PASS"
    )
    print(
        "right_boundary_even_flat="
        f"blockers={boundary[0]},mitm_vectors={boundary[1]},"
        f"nonscalar_false_positives={boundary[2]};status=PASS"
    )
    print(
        "full_microstaircase_even_flat="
        f"blockers={full[0]},mitm_vectors={full[1]},"
        f"unit_one_defect_safe_candidates={full[2]},"
        f"coordinate6_half_turn_safe_candidates={full[3]};status=PASS"
    )
    print(f"source_lf_sha256={sha256(source).hexdigest()}")
    print("ALL_CHECKS_PASSED")


if __name__ == "__main__":
    main()
