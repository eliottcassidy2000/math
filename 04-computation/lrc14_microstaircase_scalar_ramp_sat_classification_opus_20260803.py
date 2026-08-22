#!/usr/bin/env python3
"""Exact SAT classification of the n=14 full micro-staircase blockers.

The finite object is the one introduced in the historical S364/S371 probes.
For every open-cell pattern ``b`` and shift ``s mod n``, a residue vector
``v`` blocks the candidate when some coordinate satisfies

    s*v_i + b_i in {0, n-1} (mod n).

THM-363 proves that adding a scalar ramp ``m*i`` merely reindexes cells.  We
therefore fix ``v_1=0`` and ask whether a nonzero normalized full blocker
exists.  Two independent SAT engines answer UNSAT.  A separate all-model
enumeration, without fixing the gauge, finds exactly the n scalar ramps.

The script rebuilds the rational breakpoint arrangement and the clauses
directly, then cross-checks every pattern and every coverage bit against the
older S364 implementation.  All arithmetic before the SAT call is exact.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from pathlib import Path
from platform import python_version

import pysat
from pysat.solvers import Solver


ROOT = Path(__file__).resolve().parents[1]
S364_PATH = ROOT / "04-computation" / "lonely_runner_feedback_loop_s364.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def cell_pattern(n: int, alpha: Fraction) -> tuple[int, ...]:
    values = []
    for i in range(1, n):
        scaled = n * ((i * alpha) % 1)
        values.append(scaled.numerator // scaled.denominator)
    return tuple(values)


def cell_patterns(n: int) -> tuple[tuple[int, ...], ...]:
    breaks = {Fraction(0), Fraction(1)}
    for i in range(1, n):
        for a in range(n * i + 1):
            breaks.add(Fraction(a, n * i))

    patterns: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    ordered = sorted(breaks)
    for lo, hi in zip(ordered, ordered[1:]):
        if lo == hi:
            continue
        pattern = cell_pattern(n, (lo + hi) / 2)
        if pattern not in seen:
            seen.add(pattern)
            patterns.append(pattern)
    return tuple(patterns)


def variable(n: int, coordinate: int, residue: int) -> int:
    return coordinate * n + residue + 1


def scalar_multiplier(vector: tuple[int, ...], n: int) -> int | None:
    for multiplier in range(n):
        if all(
            value == multiplier * index % n
            for index, value in enumerate(vector, start=1)
        ):
            return multiplier
    return None


def blocks_candidate(
    vector: tuple[int, ...], n: int, shift: int, pattern: tuple[int, ...]
) -> bool:
    return any(
        (shift * residue + pattern[index]) % n in (0, n - 1)
        for index, residue in enumerate(vector)
    )


def cover_clause(
    n: int, shift: int, pattern: tuple[int, ...]
) -> list[int]:
    return [
        variable(n, index, residue)
        for index in range(n - 1)
        for residue in range(n)
        if (shift * residue + pattern[index]) % n in (0, n - 1)
    ]


def build_cnf(
    n: int,
    patterns: tuple[tuple[int, ...], ...],
    *,
    fix_gauge: bool,
    forbid_zero: bool,
) -> list[list[int]]:
    clauses: list[list[int]] = []
    for index in range(n - 1):
        clauses.append([variable(n, index, residue) for residue in range(n)])
        for left in range(n):
            for right in range(left + 1, n):
                clauses.append(
                    [
                        -variable(n, index, left),
                        -variable(n, index, right),
                    ]
                )

    if fix_gauge:
        clauses.append([variable(n, 0, 0)])
    if forbid_zero:
        require(fix_gauge, "the nonzero normalized clause requires a fixed gauge")
        clauses.append(
            [
                variable(n, index, residue)
                for index in range(1, n - 1)
                for residue in range(1, n)
            ]
        )

    for shift in range(n):
        for pattern in patterns:
            clause = cover_clause(n, shift, pattern)
            require(bool(clause), "a candidate has no possible blocker")
            clauses.append(clause)
    return clauses


def decode_model(model: list[int], n: int) -> tuple[int, ...]:
    positive = {literal for literal in model if literal > 0}
    vector = []
    for index in range(n - 1):
        chosen = [
            residue
            for residue in range(n)
            if variable(n, index, residue) in positive
        ]
        require(len(chosen) == 1, "SAT model violates the exact-one encoding")
        vector.append(chosen[0])
    return tuple(vector)


def verify_full_blocker(
    vector: tuple[int, ...], n: int, patterns: tuple[tuple[int, ...], ...]
) -> None:
    require(len(vector) == n - 1, "decoded vector has the wrong length")
    require(all(0 <= residue < n for residue in vector), "residue outside Z/n")
    for shift in range(n):
        for pattern in patterns:
            require(
                blocks_candidate(vector, n, shift, pattern),
                f"model misses shift={shift}, pattern={pattern}",
            )


def dimacs_bytes(clauses: list[list[int]], variable_count: int) -> bytes:
    rows = [f"p cnf {variable_count} {len(clauses)}\n"]
    rows.extend(" ".join(map(str, clause)) + " 0\n" for clause in clauses)
    return "".join(rows).encode("ascii")


def solve_verdict(name: str, clauses: list[list[int]]) -> bool:
    with Solver(name=name, bootstrap_with=clauses) as solver:
        return solver.solve()


def enumerate_models(
    clauses: list[list[int]], n: int, patterns: tuple[tuple[int, ...], ...]
) -> tuple[tuple[int, ...], ...]:
    models: list[tuple[int, ...]] = []
    with Solver(name="cadical195", bootstrap_with=clauses) as solver:
        while solver.solve():
            vector = decode_model(solver.get_model(), n)
            verify_full_blocker(vector, n, patterns)
            models.append(vector)
            solver.add_clause(
                [
                    -variable(n, index, residue)
                    for index, residue in enumerate(vector)
                ]
            )
    return tuple(sorted(models))


def cross_check_s364(
    n: int, patterns: tuple[tuple[int, ...], ...]
) -> None:
    old = SourceFileLoader("s364_cross_check", str(S364_PATH)).load_module()
    system = old.build_pattern_system(n)
    require(system.patterns == patterns, "direct and S364 pattern orders differ")
    require(system.candidate_count == n * len(patterns), "candidate count mismatch")

    candidate = 0
    direct_masks = [[0 for _ in range(n)] for _ in range(n - 1)]
    for shift in range(n):
        for pattern in patterns:
            bit = 1 << candidate
            direct_clause = set(cover_clause(n, shift, pattern))
            old_clause: set[int] = set()
            for index in range(n - 1):
                for residue in range(n):
                    if (shift * residue + pattern[index]) % n in (0, n - 1):
                        direct_masks[index][residue] |= bit
                    if (system.masks[index][residue] >> candidate) & 1:
                        old_clause.add(variable(n, index, residue))
            require(direct_clause == old_clause, f"coverage clause mismatch at {candidate}")
            candidate += 1
    require(
        tuple(tuple(row) for row in direct_masks) == system.masks,
        "direct and S364 coverage masks differ",
    )


def run_dimension(n: int, *, enumerate_all: bool) -> None:
    patterns = cell_patterns(n)
    cross_check_s364(n, patterns)
    variables = n * (n - 1)

    normalized = build_cnf(
        n, patterns, fix_gauge=True, forbid_zero=True
    )
    verdicts = {
        name: solve_verdict(name, normalized)
        for name in ("cadical195", "glucose4")
    }
    require(
        verdicts == {"cadical195": False, "glucose4": False},
        f"unexpected normalized verdicts: {verdicts}",
    )

    widths = Counter(map(len, normalized))
    digest = sha256(dimacs_bytes(normalized, variables)).hexdigest()
    print(f"n={n}")
    print(f"  patterns={len(patterns)} candidates={n * len(patterns)}")
    print(f"  variables={variables} normalized_clauses={len(normalized)}")
    print(f"  normalized_dimacs_sha256={digest}")
    print(f"  normalized_clause_widths={tuple(sorted(widths.items()))}")
    print("  cadical195_normalized_nonzero=UNSAT")
    print("  glucose4_normalized_nonzero=UNSAT")

    if enumerate_all:
        ungauged = build_cnf(
            n, patterns, fix_gauge=False, forbid_zero=False
        )
        models = enumerate_models(ungauged, n, patterns)
        multipliers = tuple(scalar_multiplier(vector, n) for vector in models)
        require(len(models) == n, f"expected {n} ungauged models, found {len(models)}")
        require(
            set(multipliers) == set(range(n)),
            f"ungauged models are not exactly the scalar ramps: {multipliers}",
        )
        require(
            all(multiplier is not None for multiplier in multipliers),
            "a nonscalar model survived enumeration",
        )
        print(f"  ungauged_models={len(models)}")
        print(f"  ungauged_scalar_multipliers={tuple(sorted(multipliers))}")
        print("  ungauged_nonscalar_models=0")


def main() -> None:
    print("Full micro-staircase scalar-ramp SAT classification")
    print("scope=FINITE-EXACT residue-cell model; not LRC(14)")
    print(f"runtime=python:{python_version()};pysat:{pysat.__version__}")
    run_dimension(14, enumerate_all=True)
    run_dimension(15, enumerate_all=False)
    print("ALL_CHECKS_PASSED")


if __name__ == "__main__":
    main()
