#!/usr/bin/env python3
"""Exact controls for THM-2624's two-clock signed root tomography.

The canonical matrix W_(s,ell) has one row for each active THM-2614 target
section q and twelve columns for the nonzero deep-probe roots r.  Entries are
globally primitive integer weights, summed over the two retained rails and all
seven future clocks.  Exact Fraction RREF certifies every stated Q-rank and
kernel; no floating-point or optimized-mode assert is used.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from math import gcd

import lrc14_punctured_target_root_cosupport_thm2614 as old


P = 13
Q7 = old.Q7
ROOTS = tuple(range(1, P))
GENERIC_S = (1, 3, 4, 5, 8, 9, 10, 12)
EXCEPTIONAL_A = (2, 7)
EXCEPTIONAL_B = (6, 11)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def primitive_integer_vector(vector):
    denominator = 1
    for value in vector:
        denominator = denominator * value.denominator // gcd(
            denominator, value.denominator
        )
    values = [value.numerator * (denominator // value.denominator)
              for value in vector]
    content = 0
    for value in values:
        content = gcd(content, abs(value))
    require(content > 0, "zero vector cannot be projectivized")
    values = [value // content for value in values]
    if next(value for value in values if value) < 0:
        values = [-value for value in values]
    return tuple(values)


def exact_rref(matrix):
    """Return exact rank, pivot columns, and reduced rows over Q."""
    rows = [[Fraction(value) for value in row] for row in matrix]
    rank = 0
    pivots = []
    for col in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows))
                      if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][col]
        rows[rank] = [value / scale for value in rows[rank]]
        for i in range(len(rows)):
            if i == rank or not rows[i][col]:
                continue
            factor = rows[i][col]
            rows[i] = [left - factor * right
                       for left, right in zip(rows[i], rows[rank])]
        pivots.append(col)
        rank += 1
        if rank == min(len(rows), len(rows[0])):
            break
    return rank, tuple(pivots), tuple(tuple(row) for row in rows)


def exact_rank_nullspace(matrix):
    rank, pivots, rows = exact_rref(matrix)
    free = tuple(col for col in range(len(matrix[0])) if col not in pivots)
    null = []
    for col in free:
        vector = [Fraction(0)] * len(matrix[0])
        vector[col] = Fraction(1)
        for i, pivot in enumerate(pivots):
            vector[pivot] = -rows[i][col]
        integer = primitive_integer_vector(vector)
        require(all(sum(value * coordinate
                        for value, coordinate in zip(row, integer)) == 0
                    for row in matrix), "exact null vector failed")
        null.append(integer)
    require(len(null) == len(matrix[0]) - rank,
            "rank-nullity failed")
    return rank, tuple(null)


def canonical_dihedral_projective(vector):
    orbit = []
    for direction in (vector, tuple(reversed(vector))):
        for shift in range(len(vector)):
            rotated = direction[shift:] + direction[:shift]
            if next(value for value in rotated if value) < 0:
                rotated = tuple(-value for value in rotated)
            orbit.append(rotated)
    return min(orbit)


def cyclotomic_nonzero(matrix, qlabels, alpha, beta):
    """Exact Q(zeta_13) nonvanishing by reduction modulo Phi_13."""
    coefficients = [0] * P
    for row, q in zip(matrix, qlabels):
        for r, value in zip(ROOTS, row):
            coefficients[(alpha * q + beta * r) % P] += value
    # A degree-at-most-12 polynomial vanishes at zeta_13 iff all thirteen
    # coefficients are equal, since Phi_13=1+...+x^12.
    return len(set(coefficients)) > 1


def stack(matrices):
    return [row for matrix in matrices for row in matrix]


def rebuild_matrices():
    rails, joint, unit, by_cell, _, _ = old.rebuild_bank()
    require(len(rails) == 162 and len(by_cell) == 84,
            "THM-2614 carrier changed")
    matrices = {}
    qlabels = {}
    for cell, edges in sorted(by_cell.items()):
        qs = tuple(q for q in range(P)
                   if any(unit[j][q] for j, _ in edges))
        matrix = tuple(
            tuple(
                sum(joint[j][q][ell5][r]
                    for j, _ in edges if unit[j][q]
                    for ell5 in range(Q7)) // old.GLOBAL_CONTENT
                for r in ROOTS
            )
            for q in qs
        )
        require(all(all(value > 0 for value in row) for row in matrix),
                "an active q row lost full nonzero-root support")
        matrices[cell] = matrix
        qlabels[cell] = qs
    return matrices, qlabels


def main():
    matrices, qlabels = rebuild_matrices()

    ranks = {}
    nullspaces = {}
    rank_hist = Counter()
    mixed_hist = Counter()
    rank11_orbits = defaultdict(list)
    for cell, matrix in matrices.items():
        rank, null = exact_rank_nullspace(matrix)
        ranks[cell] = rank
        nullspaces[cell] = null
        rank_hist[rank] += 1
        mixed = sum(
            cyclotomic_nonzero(matrix, qlabels[cell], alpha, beta)
            for alpha in range(1, P) for beta in range(1, P)
        )
        mixed_hist[mixed] += 1
        if rank == 11:
            require(len(null) == 1, "rank-eleven nullity changed")
            rank11_orbits[canonical_dihedral_projective(null[0])].append(cell)

    expected_rank_hist = Counter({5: 6, 6: 12, 7: 4, 8: 22,
                                  9: 6, 10: 2, 11: 32})
    require(rank_hist == expected_rank_hist, "one-clock rank atlas changed")
    require(mixed_hist == Counter({144: 84}),
            "a mixed cyclotomic character vanished")
    require(len(rank11_orbits) == 2 and
            sorted(map(len, rank11_orbits.values())) == [16, 16],
            "rank-eleven null orbits changed")

    signatures = defaultdict(list)
    for s in range(1, P):
        signatures[tuple(ranks[s, ell] for ell in range(Q7))].append(s)

    pair_rank_hist = Counter()
    full_pairs = {}
    pair_patterns = defaultdict(list)
    for s in range(1, P):
        winners = []
        for left, right in combinations(range(Q7), 2):
            rank, _ = exact_rank_nullspace(
                stack((matrices[s, left], matrices[s, right]))
            )
            pair_rank_hist[rank] += 1
            if rank == len(ROOTS):
                winners.append((left, right))
        require(winners, f"source shift {s} has no two-clock tomography")
        full_pairs[s] = tuple(winners)
        pair_patterns[tuple(winners)].append(s)

    require(len(pair_patterns) == 3, "two-clock graph types changed")
    require({tuple(values) for values in pair_patterns.values()} ==
            {EXCEPTIONAL_A, EXCEPTIONAL_B, GENERIC_S},
            "two-clock graph symmetry classes changed")

    adjacent_rank_hist = Counter()
    for s in range(1, P):
        for ell in range(Q7):
            rank, _ = exact_rank_nullspace(
                stack((matrices[s, ell], matrices[s, (ell + 1) % Q7]))
            )
            adjacent_rank_hist[rank] += 1
    require(adjacent_rank_hist == Counter({9: 8, 10: 6, 11: 10, 12: 60}),
            "adjacent-clock rank atlas changed")

    thresholds = {}
    for s in range(1, P):
        for ell in range(Q7):
            threshold = None
            for width in range(1, Q7 + 1):
                cells = tuple(matrices[s, (ell + offset) % Q7]
                              for offset in range(width))
                if exact_rank_nullspace(stack(cells))[0] == len(ROOTS):
                    threshold = width
                    break
            require(threshold is not None, "seven clocks lost full rank")
            thresholds[s, ell] = threshold
    threshold_hist = Counter(thresholds.values())
    require(threshold_hist == Counter({2: 60, 3: 22, 4: 2}),
            "consecutive-clock threshold changed")
    max_cells = tuple(cell for cell, value in thresholds.items() if value == 4)
    require(max_cells == ((2, 5), (11, 0)),
            "maximal consecutive-window cells changed")

    # Positivity boundary: every available row is strictly positive in all
    # twelve root columns.  Hence any nonzero nonnegative row combination is
    # strictly positive in all columns and cannot be a coordinate row e_r.
    # Every rational left inverse supplied by a full stack is therefore signed.
    positive_rows = sum(len(matrix) for matrix in matrices.values())
    require(positive_rows == 1_080,
            "active-row census changed")

    orbit_summary = tuple(sorted(
        (tuple(cells), vector, max(map(abs, vector)))
        for vector, cells in rank11_orbits.items()
    ))
    signature_summary = tuple(sorted(
        (tuple(values), signature)
        for signature, values in signatures.items()
    ))
    pair_pattern_summary = tuple(sorted(
        (tuple(values), pattern)
        for pattern, values in pair_patterns.items()
    ))

    print("THM-2624 exact two-clock root tomography controls")
    print(f"one_clock_exact_rank_hist={sorted(rank_hist.items())} "
          f"mixed_character_nonzero_hist={sorted(mixed_hist.items())}")
    print(f"rank_signatures_by_s={signature_summary}")
    print(f"rank11_projective_null_orbits={orbit_summary}")
    print(f"fixed_s_pair_rank_hist={sorted(pair_rank_hist.items())}")
    print(f"full_rank_pair_graph_types={pair_pattern_summary}")
    print(f"adjacent_pair_rank_hist={sorted(adjacent_rank_hist.items())}")
    print(f"consecutive_window_threshold_hist={sorted(threshold_hist.items())} "
          f"width4_cells={max_cells}")
    print(f"strictly_positive_active_rows={positive_rows} "
          "nonnegative_left_inverse_exists=0")
    print("verdict=PASS: two clocks give signed root tomography; disjoint clock strata do not give physical descent")


if __name__ == "__main__":
    main()
