#!/usr/bin/env python3
"""Dependency-free exact referee for THM-1256.

The audit has four independent layers.

1.  The THM-1254 residual invoice is reduced to the endpoint-address order
    ``P >= nr``.
2.  A binary relative digit records the order of the two centered phases.
3.  Ordered minimal interval chains are checked for the alignment-or-adjacent-
    overlap dichotomy.
4.  Exact danger teeth are checked for the backtrack detuning law; abstract
    owner words are checked for the resulting ABAB exclusion and turn floor.

Only integer arithmetic and ``fractions.Fraction`` are used.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import ceil


def determinant_audit() -> int:
    rows = 0
    for s0 in range(1, 21):
        for sr in range(1, 21):
            for n0 in range(1, 16):
                for nr in range(1, 16):
                    endpoint = nr * s0 - n0 * sr
                    if endpoint <= 0:
                        continue
                    delta_recip = Fraction(s0, n0) - Fraction(sr, nr)
                    assert Fraction(endpoint, n0 * nr) == delta_recip
                    for k in range(11):
                        for digit in (0, 1):
                            P = nr + k + digit
                            residual = P * s0 - n0 * sr
                            assert residual == endpoint + s0 * (k + digit)
                            assert residual >= endpoint
                            assert (residual >= endpoint) == (P >= nr)
                            rows += 1
    return rows


def binary_phase_audit() -> int:
    rows = 0
    for Q in range(2, 151):
        for ell in range(1, Q):
            if not 5 * Q < 28 * ell < 23 * Q:
                continue
            theta = Fraction(ell, Q)
            for digit in (0, 1):
                phase_difference = Fraction(theta - digit, Q)
                if digit == 0:
                    assert phase_difference > 0
                else:
                    assert phase_difference < 0
                assert Fraction(5, 28 * Q) < abs(phase_difference)
                assert abs(phase_difference) < Fraction(23, 28 * Q)
                rows += 1
    return rows


def interval_chain_audit() -> tuple[int, int, int, int]:
    chains = 0
    marked_pairs = 0
    mismatch_witnesses = 0
    marked_phase_words = 0
    grid = range(10)
    half_grid = [Fraction(z, 2) for z in range(19)]

    for length in range(2, 6):
        for lefts in combinations(grid, length):
            for rights in combinations(range(1, 11), length):
                if any(lefts[a] >= rights[a] for a in range(length)):
                    continue
                if any(lefts[a + 1] >= rights[a] for a in range(length - 1)):
                    continue
                if any(rights[a] > lefts[a + 2] for a in range(length - 2)):
                    continue
                chains += 1

                for p in range(length):
                    for q in range(p + 1, length):
                        marked_pairs += 1
                        if lefts[q] < rights[p]:
                            assert q == p + 1
                        if q >= p + 2:
                            assert rights[p] <= lefts[q]

                        xs = [x for x in half_grid if lefts[p] < x < rights[p]]
                        ys = [y for y in half_grid if lefts[q] < y < rights[q]]
                        for x in xs:
                            for y in ys:
                                if x <= y:
                                    # The consecutive-overlap chain has union
                                    # (left_p,right_q), hence covers [x,y].
                                    for z in half_grid:
                                        if x <= z <= y:
                                            assert any(
                                                lefts[a] < z < rights[a]
                                                for a in range(p, q + 1)
                                            ) or z in (x, y)
                                else:
                                    mismatch_witnesses += 1
                                    assert q == p + 1
                                    assert lefts[q] < rights[p]

                # Two canonical interior choices per interval audit the
                # global consequence: every inversion is an adjacent swap,
                # and the inversion edges form a matching.
                mark_choices = [
                    (
                        Fraction(3 * lefts[a] + rights[a], 4),
                        Fraction(lefts[a] + 3 * rights[a], 4),
                    )
                    for a in range(length)
                ]
                for marks in product(*mark_choices):
                    inversions = []
                    for p in range(length):
                        for q in range(p + 1, length):
                            if marks[p] > marks[q]:
                                assert q == p + 1
                                inversions.append((p, q))
                    used = [v for edge in inversions for v in edge]
                    assert len(used) == len(set(used))
                    marked_phase_words += 1

    return chains, marked_pairs, mismatch_witnesses, marked_phase_words


def tooth(speed: int, address: int) -> tuple[Fraction, Fraction]:
    return (
        Fraction(14 * address - 1, 14 * speed),
        Fraction(14 * address + 1, 14 * speed),
    )


def backtrack_tooth_audit() -> int:
    rows = 0
    for middle_speed in range(1, 31):
        for outer_speed in range(2, 251):
            for rung in range(1, 6):
                for outer_address in range(1, 31):
                    midpoint = Fraction(
                        middle_speed * (2 * outer_address + rung),
                        2 * outer_speed,
                    )
                    base = midpoint.numerator // midpoint.denominator
                    for middle_address in range(max(1, base - 2), base + 3):
                        left = tooth(outer_speed, outer_address)
                        middle = tooth(middle_speed, middle_address)
                        right = tooth(outer_speed, outer_address + rung)
                        if not (left[0] < middle[0] < right[0]):
                            continue
                        if not (left[1] < middle[1] < right[1]):
                            continue
                        if not (middle[0] < left[1] and right[0] < middle[1]):
                            continue
                        assert left[1] <= right[0]

                        omega_left = left[1] - middle[0]
                        omega_right = middle[1] - right[0]
                        detuning = outer_speed - (7 * rung - 1) * middle_speed
                        assert omega_left > 0 and omega_right > 0
                        assert detuning > 0
                        assert outer_speed > 6 * middle_speed
                        assert omega_left + omega_right == Fraction(
                            detuning, 7 * middle_speed * outer_speed
                        )
                        rows += 1
    return rows


def word_audit() -> tuple[int, int, int]:
    rows = 0
    max_edge_run = 0
    minimum_turn_slack = 10**9
    alphabet = range(4)

    for length in range(3, 11):
        for word in product(alphabet, repeat=length):
            if any(word[a] == word[a + 1] for a in range(length - 1)):
                continue
            abab = any(
                word[a] == word[a + 2] and word[a + 1] == word[a + 3]
                for a in range(length - 3)
            )
            if abab:
                continue

            edges = [frozenset((word[a], word[a + 1])) for a in range(length - 1)]
            backtracks = sum(edges[a - 1] == edges[a] for a in range(1, len(edges)))
            turns = len(edges) - 1 - backtracks
            assert backtracks + turns == length - 2
            assert 2 * backtracks <= length - 1
            floor = ceil(Fraction(length - 1, 2)) - 1
            assert turns >= floor
            minimum_turn_slack = min(minimum_turn_slack, turns - floor)

            run = 1
            local_max = 1
            for a in range(1, len(edges)):
                if edges[a] == edges[a - 1]:
                    run += 1
                    local_max = max(local_max, run)
                else:
                    run = 1
            assert local_max <= 2
            max_edge_run = max(max_edge_run, local_max)
            rows += 1

    assert max_edge_run == 2
    assert minimum_turn_slack == 0
    return rows, max_edge_run, minimum_turn_slack


def main() -> None:
    determinant_rows = determinant_audit()
    phase_rows = binary_phase_audit()
    chains, marked_pairs, mismatches, phase_words = interval_chain_audit()
    backtrack_rows = backtrack_tooth_audit()
    word_rows, max_run, min_slack = word_audit()

    print("THM-1256 BINARY PHASE / TOOTH-WORD LANDING EXACT AUDIT")
    print(f"endpoint residual/address rows = {determinant_rows}")
    print(f"binary centered-phase rows = {phase_rows}")
    print(f"minimal interval chains = {chains}")
    print(f"marked interval pairs = {marked_pairs}")
    print(f"explicit phase/order mismatch witnesses = {mismatches}")
    print(f"marked phase words / adjacent-swap matching checks = {phase_words}")
    print(f"exact same-provider backtrack teeth = {backtrack_rows}")
    print(f"ABAB-free owner words = {word_rows}")
    print(f"maximum unordered-edge run = {max_run}")
    print(f"minimum turn-floor slack = {min_slack}")
    assert 6**4 < 2345 < 6**5
    print("compact-box backtrack-DAG height ceiling = 4")
    print("invoice interpretation: R-endpoint = s0*(P-nr) = s0*(k+delta)")
    print("word law: backtracks + nonbacktracking turns = N-2")
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
