#!/usr/bin/env python3
"""Exact three-checkpoint clause-state Bellman referee for THM-2226.

The computation uses the uniform root-sheet cap 26/169 for every
13^2-stripped danger comb.  It enumerates exact primal and dual LP vertices
for all 16 compressed 169-adic offset schedules and audits the valuation
profile consequence at first depth four and five.
"""

from fractions import Fraction
from itertools import combinations


A = 169
Q = Fraction(26, A)
CLAUSES = 3
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


def minimize_marginal_dual(values, active_caps):
    """Independently enumerate vertices of the exact marginal-LP dual."""
    m = len(active_caps)
    if m == 0:
        return values[0]
    patterns = 1 << m
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


def maximize_marginal_lp(values, active_caps):
    """Maximize E values[Z] subject only to active marginal upper bounds."""
    m = len(active_caps)
    if m == 0:
        return values[0]
    patterns = 1 << m
    columns = []
    objective = []
    for z in range(patterns):
        columns.append(
            (Fraction(1),)
            + tuple(Fraction((z >> i) & 1) for i in range(m))
        )
        objective.append(values[z])
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


def eliminate_time(value, clause_caps):
    """Eliminate one root digit, allowing arbitrary active-clause coupling."""
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


def compressed_offsets(offsets):
    """Delete empty slices by capping both sorted start gaps at three."""
    h0, h1, h2 = sorted(offsets)
    gap1 = min(h1 - h0, CLAUSES)
    gap2 = min(h2 - h1, CLAUSES)
    return (0, gap1, gap1 + gap2)


def valuation_profiles(rows_by_offset):
    """Audit all profiles with first depth four or five."""
    closed = []
    open_profiles = []
    closed_by_first = {4: 0, 5: 0}
    total_by_first = {4: 0, 5: 0}
    for lambda1 in (4, 5):
        for lambda2 in range(lambda1, 19):
            for lambda3 in range(lambda2 + 1, 20):
                total_by_first[lambda1] += 1
                offsets = (
                    0,
                    (lambda2 - lambda1) // 2,
                    (lambda3 - lambda1) // 2,
                )
                bound = rows_by_offset[compressed_offsets(offsets)]
                profile = (lambda1, lambda2, lambda3)
                if bound < Fraction(961, 6930):
                    closed.append(profile)
                    closed_by_first[lambda1] += 1
                else:
                    open_profiles.append(profile)
    return (
        tuple(closed),
        tuple(open_profiles),
        closed_by_first,
        total_by_first,
    )


def main():
    target = Fraction(961, 6930)
    rows = []
    for gap1 in range(CLAUSES + 1):
        for gap2 in range(CLAUSES + 1):
            offsets = (0, gap1, gap1 + gap2)
            rows.append((checkpoint_bound(offsets), offsets))
    rows.sort(reverse=True)
    for bound, offsets in rows:
        verdict = "PASS" if bound < target else "FAIL"
        print(
            f"offsets={offsets} bound={bound} "
            f"decimal={float(bound):.15f} verdict={verdict}"
        )

    failing = tuple((bound, offsets) for bound, offsets in rows if bound >= target)
    expected_fail = ((Fraction(5934, 28561), (0, 1, 2)),)
    if failing != expected_fail:
        raise RuntimeError(f"three-checkpoint failure schedule drift: {failing}")
    passing = tuple((bound, offsets) for bound, offsets in rows if bound < target)
    worst_pass = max(passing)
    expected_worst_pass = (Fraction(7_444_136, 62_748_517), (0, 2, 4))
    if worst_pass != expected_worst_pass:
        raise RuntimeError(f"worst passing schedule drift: {worst_pass}")
    if target - worst_pass[0] != Fraction(8_713_462_357, 434_847_222_810):
        raise RuntimeError("worst passing margin drift")
    if failing[0][0] - target != Fraction(13_675_499, 197_927_730):
        raise RuntimeError("failing schedule excess drift")

    geometric_actual = Fraction(916_159, 4_826_809)
    if not target < geometric_actual < failing[0][0]:
        raise RuntimeError("geometric hostile ordering drift")

    rows_by_offset = {offsets: bound for bound, offsets in rows}
    closed, open_profiles, closed_by_first, total_by_first = valuation_profiles(
        rows_by_offset
    )
    expected_open = (
        (4, 6, 8),
        (4, 6, 9),
        (4, 7, 8),
        (4, 7, 9),
        (5, 7, 9),
        (5, 7, 10),
        (5, 8, 9),
        (5, 8, 10),
    )
    if open_profiles != expected_open:
        raise RuntimeError(f"remaining profile drift: {open_profiles}")
    if len(closed) != 217:
        raise RuntimeError("closed profile count drift")
    if closed_by_first != {4: 116, 5: 101}:
        raise RuntimeError("per-first-depth closure drift")
    if total_by_first != {4: 120, 5: 105}:
        raise RuntimeError("per-first-depth total drift")
    if 675 - len(closed) != 458:
        raise RuntimeError("current scalar ledger drift")

    compression_controls = (
        ((0, 17, 18), (0, 3, 4)),
        ((0, 2, 25), (0, 2, 5)),
        ((0, 12, 29), (0, 3, 6)),
    )
    for original, compressed in compression_controls:
        if checkpoint_bound(original) != checkpoint_bound(compressed):
            raise RuntimeError(
                f"gap compression drift: {original} -> {compressed}"
            )

    print(f"unique_failing_offsets={failing[0][1]}")
    print(f"failing_bound={failing[0][0]}")
    print(f"target={target}")
    print(f"failing_bound_minus_target={failing[0][0]-target}")
    print(f"worst_passing_offsets={worst_pass[1]}")
    print(f"worst_passing_bound={worst_pass[0]}")
    print(f"target_minus_worst_passing={target-worst_pass[0]}")
    print(f"geometric_three_checkpoint_actual={geometric_actual}")
    print(f"closed_profiles={len(closed)}")
    print("closed_lambda1_4=116/120")
    print("closed_lambda1_5=101/105")
    print(f"remaining_profiles={open_profiles}")
    print("combined_current_scalar_ledger=458")
    print("gap_compression_controls=PASS")
    print("exact_primal_dual_vertex_parity=PASS")


if __name__ == "__main__":
    main()
