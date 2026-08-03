#!/usr/bin/env python3
"""Independent binary-CNF audit of the S364 micro-staircase classification.

This script deliberately does not import either historical S364/S371 script or
the one-hot SAT classifier under audit.  It rebuilds the exact open-cell
arrangement, represents every residue by four Boolean bits, derives equality
atoms from those bits, and checks normalized nonzero infeasibility with two SAT
backends not used by the classifier (MiniSat22 and Lingeling).

For n=14 it also enumerates every ungauged model and verifies directly, against
all 11,368 (shift, open-cell) candidates, that the models are exactly the 14
scalar ramps.  The optional --skip-enumeration flag retains both normalized
UNSAT checks and is useful for a quick -O replay.

Logical reductions audited here:

* The discontinuities of floor(n*{i alpha}) in [0,1] are exactly a/(n*i).
  Consecutive points in their union are therefore all atomic open cells.
* Adding the scalar ramp m*i with m=-v_1 makes the first coordinate zero;
  the normalized vector is zero iff the original vector is v_i=v_1*i.
* Shift zero is universally blocked on every atomic cell.
* If b is a cell pattern, then (n-1-b_i)_i is the pattern at 1-alpha.
  The blocker predicate is invariant under (s,b)->(-s,n-1-b), so only
  shifts 1,...,floor(n/2) need be encoded.

Scope: a FINITE-EXACT residue-cell classification.  This is not LRC(n).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import platform

import pysat
from pysat.solvers import Solver


BITS = 4


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


@dataclass(frozen=True)
class CellSystem:
    n: int
    breakpoints: tuple[Fraction, ...]
    intervals: tuple[tuple[Fraction, Fraction], ...]
    patterns: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class BinaryEncoding:
    n: int
    variable_count: int
    base_clauses: tuple[tuple[int, ...], ...]
    representative_cover_clauses: tuple[tuple[int, ...], ...]
    representative_shifts: tuple[int, ...]
    raw_representative_candidate_count: int


def floor_bin(n: int, index: int, alpha: Fraction) -> int:
    fractional = (index * alpha) % 1
    scaled = n * fractional
    return scaled.numerator // scaled.denominator


def pattern_at(n: int, alpha: Fraction) -> tuple[int, ...]:
    return tuple(floor_bin(n, index, alpha) for index in range(1, n))


def build_cell_system(n: int) -> CellSystem:
    discontinuities = {Fraction(0), Fraction(1)}
    for index in range(1, n):
        for numerator in range(n * index + 1):
            discontinuities.add(Fraction(numerator, n * index))

    breakpoints = tuple(sorted(discontinuities))
    intervals = tuple(zip(breakpoints, breakpoints[1:]))
    patterns: list[tuple[int, ...]] = []
    for lo, hi in intervals:
        require(lo < hi, "nonpositive atomic interval")
        left_sample = (2 * lo + hi) / 3
        right_sample = (lo + 2 * hi) / 3
        left_pattern = pattern_at(n, left_sample)
        right_pattern = pattern_at(n, right_sample)
        require(
            left_pattern == right_pattern,
            f"pattern changed inside claimed atomic interval {(lo, hi)}",
        )
        patterns.append(left_pattern)

    require(
        len(set(patterns)) == len(patterns),
        "two distinct atomic intervals have the same pattern",
    )
    return CellSystem(n, breakpoints, intervals, tuple(patterns))


def blocks_candidate(
    vector: tuple[int, ...], n: int, shift: int, pattern: tuple[int, ...]
) -> bool:
    return any(
        (shift * residue + pattern[index]) % n in (0, n - 1)
        for index, residue in enumerate(vector)
    )


def support(
    n: int, shift: int, pattern: tuple[int, ...]
) -> tuple[tuple[int, int], ...]:
    return tuple(
        (index, residue)
        for index in range(n - 1)
        for residue in range(n)
        if (shift * residue + pattern[index]) % n in (0, n - 1)
    )


def complement_pattern(pattern: tuple[int, ...], n: int) -> tuple[int, ...]:
    return tuple(n - 1 - value for value in pattern)


def audit_cell_universe_and_shift_reduction(system: CellSystem) -> None:
    n = system.n
    pattern_set = set(system.patterns)
    require(
        all(any(value in (0, n - 1) for value in pattern) for pattern in system.patterns),
        "shift zero is not universally blocked",
    )
    for pattern in system.patterns:
        reflected = complement_pattern(pattern, n)
        require(reflected in pattern_set, "cell patterns are not reflection-closed")
        for shift in range(1, n):
            require(
                support(n, shift, pattern)
                == support(n, (-shift) % n, reflected),
                f"shift-sign support mismatch at shift={shift}, pattern={pattern}",
            )


def bit_variable(index: int, bit: int) -> int:
    return 1 + BITS * index + bit


def equality_variable(n: int, index: int, residue: int) -> int:
    bit_count = BITS * (n - 1)
    return bit_count + 1 + index * n + residue


def mismatch_clause(index: int, value: int) -> tuple[int, ...]:
    return tuple(
        -bit_variable(index, bit)
        if (value >> bit) & 1
        else bit_variable(index, bit)
        for bit in range(BITS)
    )


def equality_equivalence_clauses(
    n: int, index: int, residue: int
) -> tuple[tuple[int, ...], ...]:
    atom = equality_variable(n, index, residue)
    forward = tuple(
        (
            -atom,
            bit_variable(index, bit)
            if (residue >> bit) & 1
            else -bit_variable(index, bit),
        )
        for bit in range(BITS)
    )
    reverse = ((atom,) + mismatch_clause(index, residue),)
    return forward + reverse


def clause_value(clause: tuple[int, ...], assignment: dict[int, bool]) -> bool:
    return any(
        assignment[abs(literal)] if literal > 0 else not assignment[abs(literal)]
        for literal in clause
    )


def audit_local_binary_semantics(n: int) -> None:
    # Exhaustively audit both directions of E(i,r) <-> (bits(i)=r).
    index = 0
    for residue in range(n):
        clauses = equality_equivalence_clauses(n, index, residue)
        atom = equality_variable(n, index, residue)
        for value in range(1 << BITS):
            for atom_value in (False, True):
                assignment = {
                    bit_variable(index, bit): bool((value >> bit) & 1)
                    for bit in range(BITS)
                }
                assignment[atom] = atom_value
                cnf_value = all(clause_value(clause, assignment) for clause in clauses)
                require(
                    cnf_value == (atom_value == (value == residue)),
                    f"bad equality encoding for n={n}, r={residue}, value={value}",
                )

    # The domain clauses accept precisely 0,...,n-1.
    domain = tuple(mismatch_clause(index, value) for value in range(n, 1 << BITS))
    for value in range(1 << BITS):
        assignment = {
            bit_variable(index, bit): bool((value >> bit) & 1)
            for bit in range(BITS)
        }
        accepted = all(clause_value(clause, assignment) for clause in domain)
        require(accepted == (value < n), f"bad binary domain encoding at value={value}")


def build_encoding(system: CellSystem) -> BinaryEncoding:
    n = system.n
    clauses: list[tuple[int, ...]] = []
    for index in range(n - 1):
        for residue in range(n):
            clauses.extend(equality_equivalence_clauses(n, index, residue))
        clauses.extend(
            mismatch_clause(index, invalid)
            for invalid in range(n, 1 << BITS)
        )

    representative_shifts = tuple(range(1, n // 2 + 1))
    cover_clauses: list[tuple[int, ...]] = []
    seen_cover_clauses: set[tuple[int, ...]] = set()
    raw_representative_candidate_count = 0
    for shift in representative_shifts:
        for pattern in system.patterns:
            raw_representative_candidate_count += 1
            clause = tuple(
                equality_variable(n, index, residue)
                for index, residue in support(n, shift, pattern)
            )
            require(bool(clause), "empty representative cover clause")
            if clause not in seen_cover_clauses:
                seen_cover_clauses.add(clause)
                cover_clauses.append(clause)

    variable_count = BITS * (n - 1) + n * (n - 1)
    return BinaryEncoding(
        n=n,
        variable_count=variable_count,
        base_clauses=tuple(clauses),
        representative_cover_clauses=tuple(cover_clauses),
        representative_shifts=representative_shifts,
        raw_representative_candidate_count=raw_representative_candidate_count,
    )


def normalized_clauses(encoding: BinaryEncoding) -> tuple[tuple[int, ...], ...]:
    n = encoding.n
    gauge = tuple((-bit_variable(0, bit),) for bit in range(BITS))
    nonzero = (
        tuple(
            bit_variable(index, bit)
            for index in range(1, n - 1)
            for bit in range(BITS)
        ),
    )
    return (
        encoding.base_clauses
        + gauge
        + nonzero
        + encoding.representative_cover_clauses
    )


def dimacs_digest(clauses: tuple[tuple[int, ...], ...], variable_count: int) -> str:
    rows = [f"p cnf {variable_count} {len(clauses)}\n"]
    rows.extend(" ".join(map(str, clause)) + " 0\n" for clause in clauses)
    return sha256("".join(rows).encode("ascii")).hexdigest()


def solve_unsat(
    solver_name: str, clauses: tuple[tuple[int, ...], ...]
) -> None:
    with Solver(name=solver_name, bootstrap_with=clauses) as solver:
        verdict = solver.solve()
    require(not verdict, f"{solver_name} found a normalized nonzero model")


def decode_bits(model: list[int], n: int) -> tuple[int, ...]:
    positive = {literal for literal in model if literal > 0}
    vector = tuple(
        sum(
            1 << bit
            for bit in range(BITS)
            if bit_variable(index, bit) in positive
        )
        for index in range(n - 1)
    )
    require(all(value < n for value in vector), "solver model escaped binary domain")
    for index, value in enumerate(vector):
        for residue in range(n):
            atom_true = equality_variable(n, index, residue) in positive
            require(
                atom_true == (value == residue),
                "decoded model violates equality-atom equivalence",
            )
    return vector


def verify_full_unreduced(system: CellSystem, vector: tuple[int, ...]) -> None:
    require(len(vector) == system.n - 1, "wrong vector length")
    for shift in range(system.n):
        for pattern in system.patterns:
            require(
                blocks_candidate(vector, system.n, shift, pattern),
                f"vector misses shift={shift}, pattern={pattern}",
            )


def scalar_ramp(n: int, multiplier: int) -> tuple[int, ...]:
    return tuple(multiplier * index % n for index in range(1, n))


def normalized_vector(vector: tuple[int, ...], n: int) -> tuple[int, ...]:
    multiplier = (-vector[0]) % n
    return tuple(
        (value + multiplier * index) % n
        for index, value in enumerate(vector, start=1)
    )


def audit_gauge_and_scalar_direction(system: CellSystem) -> None:
    n = system.n
    zero = (0,) * (n - 1)
    for multiplier in range(n):
        ramp = scalar_ramp(n, multiplier)
        require(normalized_vector(ramp, n) == zero, "scalar ramp did not normalize to zero")
        verify_full_unreduced(system, ramp)

    # Algebraic converse used by the classification:
    # normalize(v)=0 gives v_i-v_1*i=0 coordinatewise, hence v_i=v_1*i.
    # The following finite controls check both a scalar and a nonscalar vector.
    control = list(zero)
    if n > 2:
        control[1] = 1
        require(
            normalized_vector(tuple(control), n) != zero,
            "nonscalar hostile control normalized to zero",
        )


def enumerate_ungauged_models(
    system: CellSystem, encoding: BinaryEncoding
) -> tuple[tuple[int, ...], ...]:
    clauses = encoding.base_clauses + encoding.representative_cover_clauses
    models: list[tuple[int, ...]] = []
    with Solver(name="minisat22", bootstrap_with=clauses) as solver:
        while solver.solve():
            model = solver.get_model()
            vector = decode_bits(model, system.n)
            verify_full_unreduced(system, vector)
            models.append(vector)
            positive = {literal for literal in model if literal > 0}
            solver.add_clause(
                [
                    -bit_variable(index, bit)
                    if bit_variable(index, bit) in positive
                    else bit_variable(index, bit)
                    for index in range(system.n - 1)
                    for bit in range(BITS)
                ]
            )
    return tuple(sorted(models))


def run_dimension(n: int, *, enumerate_models: bool) -> None:
    system = build_cell_system(n)
    audit_cell_universe_and_shift_reduction(system)
    audit_local_binary_semantics(n)
    audit_gauge_and_scalar_direction(system)
    encoding = build_encoding(system)
    normalized = normalized_clauses(encoding)

    print(f"n={n}")
    print(
        "  "
        f"breakpoints={len(system.breakpoints)} "
        f"atomic_open_intervals={len(system.intervals)} "
        f"distinct_patterns={len(set(system.patterns))} "
        f"full_candidates={n * len(system.patterns)}"
    )
    print("  open_cell_samples=one_third,two_thirds agree_on_every_interval")
    print("  shift_zero=universally_blocked")
    print(
        "  "
        f"reflection_shift_representatives={encoding.representative_shifts} "
        f"raw_representative_candidates={encoding.raw_representative_candidate_count} "
        f"unique_cover_clauses={len(encoding.representative_cover_clauses)}"
    )
    print("  reflection_rule=(s,b)->(-s,n-1-b); blocker_set_preserved")
    print("  gauge_rule=w_i=v_i-v_1*i; w=0 iff v_i=v_1*i")
    print("  residue_encoding=4_bits+equality_atoms_iff; no_one_hot_clauses")
    print("  binary_truth_table_audit=PASSED")
    print(
        "  "
        f"variables={encoding.variable_count} base_clauses={len(encoding.base_clauses)} "
        f"normalized_clauses={len(normalized)}"
    )
    print(
        "  "
        f"normalized_reduced_dimacs_sha256={dimacs_digest(normalized, encoding.variable_count)}"
    )

    for solver_name in ("minisat22", "lingeling"):
        solve_unsat(solver_name, normalized)
        print(f"  {solver_name}_normalized_nonzero=UNSAT")

    print(f"  directly_verified_scalar_ramps={n}/{n}")
    if enumerate_models:
        models = enumerate_ungauged_models(system, encoding)
        ramps = {scalar_ramp(n, multiplier) for multiplier in range(n)}
        require(set(models) == ramps, "ungauged models are not exactly scalar ramps")
        print(f"  ungauged_models={len(models)} nonscalar_models=0")
        print(f"  ungauged_scalar_multipliers={tuple(range(n))}")
    else:
        print("  ungauged_enumeration=SKIPPED")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip-enumeration",
        action="store_true",
        help="skip the slower direct enumeration of all n=14 ungauged models",
    )
    args = parser.parse_args()

    source_bytes = Path(__file__).read_bytes()
    require(
        b"\r" not in source_bytes.replace(b"\r\n", b""),
        "script contains a bare carriage return",
    )
    source_lf = source_bytes.replace(b"\r\n", b"\n")
    source_digest = sha256(source_lf).hexdigest()
    print("Independent binary-CNF micro-staircase audit")
    print("scope=FINITE-EXACT residue-cell model; not LRC(14) or LRC(15)")
    print(f"python={platform.python_version()} implementation={platform.python_implementation()}")
    print(f"python_sat={pysat.__version__} platform={platform.platform()}")
    print(f"script_sha256={source_digest}")
    run_dimension(14, enumerate_models=not args.skip_enumeration)
    run_dimension(15, enumerate_models=False)
    print("ALL_CHECKS_PASSED")


if __name__ == "__main__":
    main()
