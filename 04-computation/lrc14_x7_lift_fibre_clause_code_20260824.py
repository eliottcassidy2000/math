#!/usr/bin/env python3
"""Exact first-x7 lift-fibre clause code for the LRC(14) paper sieve.

For a parent ``v`` modulo ``q`` with gcd(q, 7) = 1, its first x7 children are

    x_i = v_i + q a_i,          a_i in F_7.

At a genuinely new clock ``h/(7q)`` (h nonzero modulo 7), varying ``a_i``
moves the i-th position through the seven q-spaced lifts of one residue modulo
q.  The open danger interval at radius 1/14 contains exactly one of those
lifts, except at an exact endpoint equality, where it contains none.  Hence
every new clock is a positive clause: an improper child must agree with the
clock's forbidden word in at least one coordinate.

This file verifies that reduction directly, builds its 91-literal grouped CNF,
and compares two independent exact solvers on the canonical AP parent.  The
UNSAT rows are finite-exact solver results, not a theorem about every parent or
every prime.
"""

from __future__ import annotations

from collections import Counter
from functools import reduce
from hashlib import sha256
from itertools import combinations, product
from math import gcd
from typing import Iterable, Optional, Sequence

from ortools.sat.python import cp_model
from pysat.solvers import Solver


C = 7
K = 13
AP = tuple(range(1, K + 1))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def pattern_key(pattern: Sequence[Optional[int]]) -> tuple[int, ...]:
    return tuple(-1 if value is None else value for value in pattern)


def forbidden_pattern(
    parent: Sequence[int], q: int, h: int
) -> tuple[Optional[int], ...]:
    """Return the unique forbidden digit per coordinate at h/(7q).

    ``None`` records the endpoint equality 2*(h*v_i mod q) = q, where all
    seven lift digits are safe.  The function is used only for h not divisible
    by 7.
    """

    require(gcd(q, C) == 1, "q must be coprime to 7")
    residue = h % C
    require(residue != 0, "old clocks do not define an x7 digit clause")
    inverse = pow(residue, -1, C)
    modulus = C * q
    word: list[Optional[int]] = []
    for value in parent:
        lifted_residue = (h * value) % modulus
        layer, base_residue = divmod(lifted_residue, q)
        if 2 * base_residue < q:
            danger_layer: Optional[int] = 0
        elif 2 * base_residue > q:
            danger_layer = C - 1
        else:
            danger_layer = None
        word.append(
            None
            if danger_layer is None
            else ((danger_layer - layer) * inverse) % C
        )
    return tuple(word)


def raw_clock_words(
    parent: Sequence[int], q: int
) -> tuple[tuple[Optional[int], ...], ...]:
    return tuple(
        forbidden_pattern(parent, q, h)
        for h in range(C * q)
        if h % C
    )


def distinct_clock_words(
    parent: Sequence[int], q: int
) -> tuple[tuple[Optional[int], ...], ...]:
    return tuple(sorted(set(raw_clock_words(parent, q)), key=pattern_key))


def child(parent: Sequence[int], q: int, digits: Sequence[int]) -> tuple[int, ...]:
    require(len(parent) == len(digits), "digit length drift")
    return tuple(value + q * digit for value, digit in zip(parent, digits))


def clock_clause_survivor(
    parent: Sequence[int], q: int, digits: Sequence[int]
) -> bool:
    """True iff no genuinely new 1/(7q)-grid clock witnesses the child."""

    return all(
        any(
            forbidden is not None and digit == forbidden
            for digit, forbidden in zip(digits, word)
        )
        for word in raw_clock_words(parent, q)
    )


def direct_new_clock_survivor(
    parent: Sequence[int], q: int, digits: Sequence[int]
) -> bool:
    """Direct modular scan, independent of the forbidden-word construction."""

    values = child(parent, q, digits)
    modulus = C * q
    for h in range(modulus):
        if h % C == 0:
            continue
        safe = True
        for value in values:
            residue = (h * value) % modulus
            if 2 * min(residue, modulus - residue) < q:
                safe = False
                break
        if safe:
            return False
    return True


def parent_grid_survivor(parent: Sequence[int], q: int) -> bool:
    """True iff the parent has no 1/q-grid lonely clock at radius 1/14."""

    for h in range(q):
        if all(
            14 * min((h * value) % q, q - ((h * value) % q)) >= q
            for value in parent
        ):
            return False
    return True


def direct_gcd_proper(values: Sequence[int], lift_modulus: int) -> bool:
    """The exact S--T gcd alternative at the displayed lift modulus."""

    for omitted in range(len(values)):
        common = reduce(
            gcd,
            (lift_modulus,) + tuple(
                value for index, value in enumerate(values) if index != omitted
            ),
        )
        if common > 1:
            return True
    return False


def seven_divisibility_digits(parent: Sequence[int], q: int) -> tuple[int, ...]:
    inverse = pow(q, -1, C)
    return tuple((-value * inverse) % C for value in parent)


def new_seven_gcd_proper(
    parent: Sequence[int], q: int, digits: Sequence[int]
) -> bool:
    special = seven_divisibility_digits(parent, q)
    return sum(a == b for a, b in zip(digits, special)) >= len(parent) - 1


def variable(index: int, digit: int) -> int:
    return 1 + C * index + digit


def augmented_cnf(parent: Sequence[int], q: int) -> list[list[int]]:
    """One-hot digits + all clock clauses + the new 7-gcd exclusion."""

    clauses: list[list[int]] = []
    for index in range(len(parent)):
        clauses.append([variable(index, digit) for digit in range(C)])
        for left, right in combinations(range(C), 2):
            clauses.append([-variable(index, left), -variable(index, right)])

    for word in distinct_clock_words(parent, q):
        clauses.append(
            [
                variable(index, forbidden)
                for index, forbidden in enumerate(word)
                if forbidden is not None
            ]
        )

    special = seven_divisibility_digits(parent, q)
    for omitted in range(len(parent)):
        clauses.append(
            [
                -variable(index, special[index])
                for index in range(len(parent))
                if index != omitted
            ]
        )
    return clauses


def solve_pysat(parent: Sequence[int], q: int) -> bool:
    with Solver(name="cadical195", bootstrap_with=augmented_cnf(parent, q)) as solver:
        return bool(solver.solve())


def solve_ortools(parent: Sequence[int], q: int) -> bool:
    model = cp_model.CpModel()
    digits = [
        [model.new_bool_var(f"x_{index}_{digit}") for digit in range(C)]
        for index in range(len(parent))
    ]
    for row in digits:
        model.add_exactly_one(row)
    for word in distinct_clock_words(parent, q):
        model.add_bool_or(
            digits[index][forbidden]
            for index, forbidden in enumerate(word)
            if forbidden is not None
        )
    special = seven_divisibility_digits(parent, q)
    for omitted in range(len(parent)):
        model.add_bool_or(
            digits[index][special[index]].Not()
            for index in range(len(parent))
            if index != omitted
        )

    solver = cp_model.CpSolver()
    solver.parameters.num_workers = 1
    solver.parameters.random_seed = 0
    status = solver.solve(model)
    require(
        status in (cp_model.OPTIMAL, cp_model.INFEASIBLE),
        f"unexpected OR-Tools status {solver.status_name(status)}",
    )
    return status == cp_model.OPTIMAL


def semantic_digest(words: Iterable[Sequence[Optional[int]]]) -> str:
    payload = "\n".join(
        ",".join("*" if value is None else str(value) for value in word)
        for word in sorted(words, key=pattern_key)
    )
    return sha256(payload.encode("ascii")).hexdigest()


def length_histogram(words: Iterable[Sequence[Optional[int]]]) -> tuple[tuple[int, int], ...]:
    counts = Counter(sum(value is not None for value in word) for word in words)
    return tuple(sorted(counts.items()))


def scalar_capacity(words: Sequence[Sequence[Optional[int]]]) -> tuple[int, int]:
    """The level-one-style union capacity, included as a hostile control."""

    total = 0
    for index in range(K):
        counts = Counter(
            word[index] for word in words if word[index] is not None
        )
        total += max(counts.values())
    return total, total - len(words)


def exhaustive_reduction_control() -> int:
    """Exhaust four free lift digits, comparing clauses with raw clocks."""

    q = 22
    require(parent_grid_survivor(AP, q), "q=22 AP parent unexpectedly proper")
    tested = 0
    for prefix in product(range(C), repeat=4):
        digits = prefix + (0,) * (K - 4)
        clause = clock_clause_survivor(AP, q, digits)
        direct = direct_new_clock_survivor(AP, q, digits)
        require(clause == direct, f"clause/direct drift at {digits}")
        values = child(AP, q, digits)
        require(
            direct_gcd_proper(values, 14) == new_seven_gcd_proper(AP, q, digits),
            f"gcd sidecar drift at {digits}",
        )
        tested += 1
    return tested


def main() -> None:
    require(C**K == 96_889_010_407, "7^13 child count drift")
    exhaustive = exhaustive_reduction_control()

    q_hostile = 86  # 2*43
    hostile_digits = (3, 6, 2, 0, 1, 4, 0, 3, 6, 3, 5, 4, 4)
    hostile_child = child(AP, q_hostile, hostile_digits)
    require(parent_grid_survivor(AP, q_hostile), "q=86 parent drift")
    require(clock_clause_survivor(AP, q_hostile, hostile_digits), "hostile clause drift")
    require(direct_new_clock_survivor(AP, q_hostile, hostile_digits), "hostile direct drift")
    require(not direct_gcd_proper(hostile_child, 14), "hostile killed by gcd")
    require(not new_seven_gcd_proper(AP, q_hostile, hostile_digits), "hostile 7-gcd drift")

    q_target = 382  # 2*191, the smallest k=13 prime in the finite-check ledger
    require(parent_grid_survivor(AP, q_target), "q=382 parent drift")
    ghost_digits = seven_divisibility_digits(AP, q_target)
    ghost_child = child(AP, q_target, ghost_digits)
    require(clock_clause_survivor(AP, q_target, ghost_digits), "divisibility ghost drift")
    require(direct_new_clock_survivor(AP, q_target, ghost_digits), "ghost direct drift")
    require(all(value % 7 == 0 for value in ghost_child), "ghost not all divisible by 7")
    require(direct_gcd_proper(ghost_child, 14), "ghost gcd condition missing")

    solver_rows = {}
    for q in (q_hostile, q_target):
        pysat_status = solve_pysat(AP, q)
        ortools_status = solve_ortools(AP, q)
        require(pysat_status == ortools_status, f"two-engine drift at q={q}")
        solver_rows[q] = pysat_status
    require(solver_rows[q_hostile], "q=86 hostile should survive")
    require(not solver_rows[q_target], "q=382 AP fibre should be empty")

    print("LRC14_FIRST_X7_LIFT_FIBRE_CLAUSE_CODE")
    print("status=FINITE-EXACT_REDUCTION_PLUS_TWO_ENGINE_PARENT_CONTROLS")
    print(f"raw_child_states=7^13={C**K}")
    print("digit_literals=13*7=91")
    print(f"exhaustive_clause_vs_direct_slice=q22_free4={exhaustive}")

    for q in (q_hostile, q_target):
        raw = raw_clock_words(AP, q)
        distinct = distinct_clock_words(AP, q)
        capacity, excess = scalar_capacity(distinct)
        print(
            f"q={q} raw_new_clocks={len(raw)} distinct_clock_words={len(distinct)} "
            f"raw_length_hist={length_histogram(raw)} "
            f"distinct_length_hist={length_histogram(distinct)}"
        )
        print(
            f"q={q} scalar_capacity={capacity} scalar_excess={excess} "
            f"augmented_cnf_clauses={len(augmented_cnf(AP, q))} "
            f"semantic_sha256={semantic_digest(distinct)}"
        )
        print(
            f"q={q} pysat_cadical195={'SAT' if solver_rows[q] else 'UNSAT'} "
            f"ortools_cpsat={'SAT' if solver_rows[q] else 'UNSAT'}"
        )

    print(f"q=86_exact_improper_child_digits={hostile_digits}")
    print(f"q=86_exact_improper_child={hostile_child}")
    print(
        "q=382_clock_only_divisibility_ghost="
        f"{ghost_digits} all_13_children_divisible_by_7=1 gcd_proper=1"
    )
    print("boundary_hostile=q86_and_q382_each_have_one_distinct_length6_clause")
    print("source_scope=ST_section7_names_I13p1;_x7_cost_is_a_repo_inference")
    print("lrc14_status=OPEN")
    print("PASS")


if __name__ == "__main__":
    main()
