#!/usr/bin/env python3
"""Exact clause-state Bellman referee for THM-2224.

The computation uses the uniform root-sheet cap 26/169 for every
13^2-stripped danger comb.  It enumerates exact primal and dual LP vertices
for all 25 compressed 169-adic offset schedules.
"""

from fractions import Fraction
from itertools import combinations


A = 169
Q = Fraction(26, A)
CLAUSES = 4
ALL = (1 << CLAUSES) - 1


def solve_square(matrix, rhs):
    """Solve a square rational system, returning None when singular."""
    n = len(rhs)
    aug = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if aug[r][col]), None)
        if pivot is None:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [x / scale for x in aug[col]]
        for r in range(n):
            if r == col:
                continue
            scale = aug[r][col]
            if scale:
                aug[r] = [
                    aug[r][j] - scale * aug[col][j] for j in range(n + 1)
                ]
    return tuple(aug[i][-1] for i in range(n))


def maximize_marginal_lp(values, active_caps):
    """Maximize E values[Z] with upper bounds on active marginals.

    The LP is put in standard form by adjoining one slack variable for each
    marginal inequality, then every basic feasible solution is enumerated.
    `values` is indexed by local active-subset masks.
    """
    m = len(active_caps)
    if m == 0:
        return values[0]
    patterns = 1 << m
    columns = []
    objective = []
    # Probability variables p_Z.
    for z in range(patterns):
        columns.append(
            (Fraction(1),)
            + tuple(Fraction((z >> i) & 1) for i in range(m))
        )
        objective.append(values[z])
    # Marginal slack variables.
    for i in range(m):
        columns.append(
            (Fraction(0),)
            + tuple(Fraction(int(i == j)) for j in range(m))
        )
        objective.append(Fraction(0))
    rhs = (Fraction(1),) + tuple(active_caps)
    rank = m + 1
    best = None
    for basis in combinations(range(len(columns)), rank):
        matrix = tuple(
            tuple(columns[basis[col]][row] for col in range(rank))
            for row in range(rank)
        )
        solution = solve_square(matrix, rhs)
        if solution is None or any(x < 0 for x in solution):
            continue
        candidate = sum(
            objective[basis[i]] * solution[i] for i in range(rank)
        )
        if best is None or candidate > best:
            best = candidate
    if best is None:
        raise RuntimeError("marginal LP has no basic feasible solution")
    dual = minimize_marginal_dual(values, active_caps)
    if best != dual:
        raise RuntimeError(f"primal/dual drift: {best} != {dual}")
    return best


def minimize_marginal_dual(values, active_caps):
    """Independently enumerate vertices of the exact marginal-LP dual.

    Dual variables are one free alpha and nonnegative beta_i:

        minimize alpha + sum_i cap_i beta_i
        subject to alpha + sum_(i in Z) beta_i >= values[Z].
    """
    m = len(active_caps)
    if m == 0:
        return values[0]
    patterns = 1 << m
    # Equality rows for active pattern constraints, followed by beta_i=0.
    rows = []
    rhs = []
    for z in range(patterns):
        rows.append(
            (Fraction(1),)
            + tuple(Fraction((z >> i) & 1) for i in range(m))
        )
        rhs.append(values[z])
    for i in range(m):
        rows.append(
            (Fraction(0),)
            + tuple(Fraction(int(i == j)) for j in range(m))
        )
        rhs.append(Fraction(0))
    best = None
    rank = m + 1
    for active in combinations(range(len(rows)), rank):
        matrix = tuple(rows[i] for i in active)
        solution = solve_square(matrix, tuple(rhs[i] for i in active))
        if solution is None:
            continue
        alpha, *beta = solution
        if any(b < 0 for b in beta):
            continue
        if any(
            alpha
            + sum(beta[i] for i in range(m) if (z >> i) & 1)
            < values[z]
            for z in range(patterns)
        ):
            continue
        candidate = alpha + sum(
            active_caps[i] * beta[i] for i in range(m)
        )
        if best is None or candidate < best:
            best = candidate
    if best is None:
        raise RuntimeError("marginal dual has no enumerated vertex")
    return best


def eliminate_time(value, clause_caps):
    """Eliminate one root digit, allowing arbitrary coupling of active clauses."""
    active = tuple(sorted(clause_caps))
    local_clause_bits = tuple(1 << r for r in active)
    out = []
    for later in range(1 << CLAUSES):
        local_values = []
        for z in range(1 << len(active)):
            newly_satisfied = 0
            for i, bit in enumerate(local_clause_bits):
                if (z >> i) & 1:
                    newly_satisfied |= bit
            local_values.append(value[later | newly_satisfied])
        out.append(
            maximize_marginal_lp(
                tuple(local_values),
                tuple(clause_caps[r] for r in active),
            )
        )
    return tuple(out)


def checkpoint_bound(offsets):
    """Clause-LP bound for stripped 169-adic offsets h_1,h_2,h_3."""
    value = tuple(Fraction(int(s == ALL)) for s in range(1 << CLAUSES))
    for t in range(min(offsets), max(offsets) + CLAUSES):
        clause_caps = {}
        for h in offsets:
            r = t - h
            if 0 <= r < CLAUSES:
                clause_caps[r] = clause_caps.get(r, Fraction(0)) + Q
        value = eliminate_time(value, clause_caps)
    return value[0]


def capped_offsets(gap1, gap2):
    return (0, gap1, gap1 + gap2)


def main():
    target = Fraction(961, 6930)
    rows = []
    for gap1 in range(5):
        for gap2 in range(5):
            offsets = capped_offsets(gap1, gap2)
            rows.append((checkpoint_bound(offsets), offsets))
    rows.sort(reverse=True)
    for bound, offsets in rows:
        print(f"offsets={offsets} bound={bound} decimal={float(bound):.15f}")
    best, offsets = rows[0]
    if offsets != (0, 1, 2):
        raise RuntimeError(f"unexpected worst offset schedule: {offsets}")
    if best != Fraction(3272, 28561):
        raise RuntimeError(f"worst exact bound drift: {best}")
    if sum(bound == best for bound, _ in rows) != 1:
        raise RuntimeError("worst schedule stopped being unique")
    expected_gap = Fraction(4_772_161, 197_927_730)
    if target - best != expected_gap:
        raise RuntimeError("target gap drift")
    # The canonical geometric-chain event is a strict positive control inside
    # the relaxed worst schedule.
    geometric_actual = Fraction(3_385_513, 33_787_663)
    if not geometric_actual < best < target:
        raise RuntimeError("geometric/relaxed/target ordering drift")
    # Empty time slices are identity: hostile spot checks of both independent
    # gap compressions used in the 25-schedule reduction.
    compression_controls = (
        ((0, 17, 18), (0, 4, 5)),
        ((0, 3, 25), (0, 3, 7)),
        ((0, 12, 29), (0, 4, 8)),
    )
    for original, compressed in compression_controls:
        if checkpoint_bound(original) != checkpoint_bound(compressed):
            raise RuntimeError(
                f"gap compression drift: {original} -> {compressed}"
            )
    print(f"best_offsets={offsets}")
    print(f"best_bound={best}")
    print(f"target={target}")
    print(f"target_minus_best={target-best}")
    print(f"geometric_actual={geometric_actual}")
    print("gap_compression_controls=PASS")
    print("exact_primal_dual_vertex_parity=PASS")
    if not best < target:
        raise RuntimeError("uniform clause LP does not close the target")


if __name__ == "__main__":
    main()
