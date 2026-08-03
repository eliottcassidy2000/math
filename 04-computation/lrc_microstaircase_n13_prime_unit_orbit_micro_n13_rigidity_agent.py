#!/usr/bin/env python3
"""Exact n=13 micro-staircase rigidity via prime-unit orbit splitting.

For n=13, this program classifies full blockers in the residue/open-cell
micro-staircase model.  It is deliberately standalone: it does not import the
older one-hot classifier, the n=14 binary audit, S364, or S371.

The key extra reduction is specific to a prime modulus.  After the scalar
gauge fixes v_1=0, multiplication v -> c*v by any nonzero c in F_13 preserves
the full-blocker predicate because it merely reindexes the complete shift
set.  Every nonzero normalized vector can therefore be put uniquely into one
of eleven branches: its first nonzero coordinate is j and has value 1.

Each branch is encoded with four residue bits plus audited equality atoms and
checked by CaDiCaL 1.9.5 and MapleSAT.  The full unreduced candidate system is
also used to verify all thirteen scalar ramps and a nonscalar hostile control.

Scope: FINITE-EXACT in the n=13 residue-cell model.  This is not LRC(13), is
not LRC(14), and is not an all-n theorem.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import platform

import pysat
from pysat.solvers import Solver


N = 13
BITS = 4
BRANCH_CONFLICT_BUDGET = 250_000
EXPECTED_BREAKPOINTS = 599
EXPECTED_PATTERNS = 598
EXPECTED_FULL_CANDIDATES = 7_774
EXPECTED_BASE_CLAUSES = 816
EXPECTED_REPRESENTATIVE_COVERS = 3_588
EXPECTED_VARIABLES = 204
EXPECTED_BOUNDARY_SUPPORT_ONE = 132
EXPECTED_BOUNDARY_SUPPORT_TWO = 7_920


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def floor_bin(index: int, alpha: Fraction) -> int:
    fractional = (index * alpha) % 1
    scaled = N * fractional
    return scaled.numerator // scaled.denominator


def pattern_at(alpha: Fraction) -> tuple[int, ...]:
    return tuple(floor_bin(index, alpha) for index in range(1, N))


def build_cells() -> tuple[
    tuple[Fraction, ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[int, ...], ...],
]:
    points = {Fraction(0), Fraction(1)}
    for index in range(1, N):
        for numerator in range(N * index + 1):
            points.add(Fraction(numerator, N * index))

    breakpoints = tuple(sorted(points))
    intervals = tuple(zip(breakpoints, breakpoints[1:]))
    patterns: list[tuple[int, ...]] = []
    for lo, hi in intervals:
        require(lo < hi, ("nonpositive interval", lo, hi))
        one_third = (2 * lo + hi) / 3
        two_thirds = (lo + 2 * hi) / 3
        left = pattern_at(one_third)
        right = pattern_at(two_thirds)
        require(left == right, ("pattern changes inside interval", lo, hi))
        patterns.append(left)

    require(len(breakpoints) == EXPECTED_BREAKPOINTS, len(breakpoints))
    require(len(intervals) == EXPECTED_PATTERNS, len(intervals))
    require(len(set(patterns)) == EXPECTED_PATTERNS, "patterns are not distinct")
    return breakpoints, intervals, tuple(patterns)


def blocks(
    vector: tuple[int, ...], shift: int, pattern: tuple[int, ...]
) -> bool:
    return any(
        (shift * vector[index] + pattern[index]) % N in (0, N - 1)
        for index in range(N - 1)
    )


def support(shift: int, pattern: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(
        (index, residue)
        for index in range(N - 1)
        for residue in range(N)
        if (shift * residue + pattern[index]) % N in (0, N - 1)
    )


def bit_variable(index: int, bit: int) -> int:
    return 1 + BITS * index + bit


def equality_variable(index: int, residue: int) -> int:
    return BITS * (N - 1) + 1 + index * N + residue


def mismatch_clause(index: int, value: int) -> tuple[int, ...]:
    return tuple(
        -bit_variable(index, bit)
        if (value >> bit) & 1
        else bit_variable(index, bit)
        for bit in range(BITS)
    )


def equality_clauses(index: int, residue: int) -> tuple[tuple[int, ...], ...]:
    atom = equality_variable(index, residue)
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


def audit_local_encoding() -> None:
    index = 0
    for residue in range(N):
        clauses = equality_clauses(index, residue)
        atom = equality_variable(index, residue)
        for value in range(1 << BITS):
            for atom_value in (False, True):
                assignment = {
                    bit_variable(index, bit): bool((value >> bit) & 1)
                    for bit in range(BITS)
                }
                assignment[atom] = atom_value
                accepted = all(clause_value(clause, assignment) for clause in clauses)
                require(
                    accepted == (atom_value == (value == residue)),
                    ("bad equality iff", residue, value, atom_value),
                )

    domain = tuple(
        mismatch_clause(index, invalid) for invalid in range(N, 1 << BITS)
    )
    for value in range(1 << BITS):
        assignment = {
            bit_variable(index, bit): bool((value >> bit) & 1)
            for bit in range(BITS)
        }
        accepted = all(clause_value(clause, assignment) for clause in domain)
        require(accepted == (value < N), ("bad residue domain", value))


def build_base_clauses() -> tuple[tuple[int, ...], ...]:
    clauses: list[tuple[int, ...]] = []
    for index in range(N - 1):
        for residue in range(N):
            clauses.extend(equality_clauses(index, residue))
        clauses.extend(
            mismatch_clause(index, invalid) for invalid in range(N, 1 << BITS)
        )
    require(len(clauses) == EXPECTED_BASE_CLAUSES, len(clauses))
    return tuple(clauses)


def audit_reflection(patterns: tuple[tuple[int, ...], ...]) -> None:
    pattern_set = set(patterns)
    for pattern in patterns:
        require(
            any(value in (0, N - 1) for value in pattern),
            ("shift zero is not blocked", pattern),
        )
        reflected = tuple(N - 1 - value for value in pattern)
        require(reflected in pattern_set, ("missing reflected pattern", pattern))
        for shift in range(1, N):
            require(
                support(shift, pattern) == support((-shift) % N, reflected),
                ("reflection support mismatch", shift, pattern),
            )


def build_representative_cover_clauses(
    patterns: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], ...]:
    clauses: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    for shift in range(1, (N - 1) // 2 + 1):
        for pattern in patterns:
            pairs = support(shift, pattern)
            require(len(pairs) == 2 * (N - 1), ("bad support width", shift, pattern))
            clause = tuple(
                equality_variable(index, residue) for index, residue in pairs
            )
            require(clause not in seen, ("duplicate representative cover", shift))
            seen.add(clause)
            clauses.append(clause)
    require(len(clauses) == EXPECTED_REPRESENTATIVE_COVERS, len(clauses))
    return tuple(clauses)


def scalar_ramp(multiplier: int) -> tuple[int, ...]:
    return tuple(multiplier * index % N for index in range(1, N))


def normalize_scalar_gauge(vector: tuple[int, ...]) -> tuple[int, ...]:
    multiplier = (-vector[0]) % N
    return tuple(
        (value + multiplier * index) % N
        for index, value in enumerate(vector, start=1)
    )


def verify_positive_and_hostile_controls(
    intervals: tuple[tuple[Fraction, Fraction], ...],
    patterns: tuple[tuple[int, ...], ...],
) -> None:
    zero = (0,) * (N - 1)
    for multiplier in range(N):
        ramp = scalar_ramp(multiplier)
        require(normalize_scalar_gauge(ramp) == zero, ("bad scalar gauge", multiplier))
        for shift in range(N):
            for pattern in patterns:
                require(blocks(ramp, shift, pattern), ("scalar ramp miss", multiplier))

    hostile = (0, 1) + (0,) * 10
    misses = tuple(
        (shift, cell_index)
        for shift in range(N)
        for cell_index, pattern in enumerate(patterns)
        if not blocks(hostile, shift, pattern)
    )
    require(len(misses) == 242, ("hostile miss count", len(misses)))
    require(misses[0] == (1, 311), ("hostile first miss", misses[0]))
    require(
        intervals[311] == (Fraction(27, 52), Fraction(61, 117)),
        ("hostile interval", intervals[311]),
    )
    require(
        patterns[311] == (6, 0, 7, 1, 7, 1, 8, 2, 8, 2, 9, 3),
        ("hostile pattern", patterns[311]),
    )

    # Scaling by a nonzero field element only reindexes the complete shifts.
    hostile_counts = []
    for scalar in range(1, N):
        scaled = tuple(scalar * value % N for value in hostile)
        hostile_counts.append(
            sum(
                blocks(scaled, shift, pattern)
                for shift in range(N)
                for pattern in patterns
            )
        )
    require(len(set(hostile_counts)) == 1, ("scaling changed cover count", hostile_counts))


def boundary_blocks(vector: tuple[int, ...]) -> bool:
    """Test only cells immediately right of alpha=a/13.

    Those cells have b_i=a*i mod 13, so this is the rank-two affine-cover
    shadow <(a,s),(i,v_i)> in {0,-1}.  It is necessary, not assumed
    sufficient for the full-cell system.
    """
    return all(
        any(
            (a * index + shift * vector[index - 1]) % N in (0, N - 1)
            for index in range(1, N)
        )
        for a in range(N)
        for shift in range(N)
    )


def audit_boundary_low_support() -> tuple[int, int]:
    eligible = tuple(range(1, N - 1))  # v_1 is gauge-fixed; these are v_2,...,v_12.
    support_one = 0
    support_one_survivors = 0
    for index in eligible:
        for residue in range(1, N):
            vector = [0] * (N - 1)
            vector[index] = residue
            support_one += 1
            support_one_survivors += int(boundary_blocks(tuple(vector)))

    support_two = 0
    support_two_survivors = 0
    for left, right in combinations(eligible, 2):
        for left_residue, right_residue in product(range(1, N), repeat=2):
            vector = [0] * (N - 1)
            vector[left] = left_residue
            vector[right] = right_residue
            support_two += 1
            support_two_survivors += int(boundary_blocks(tuple(vector)))

    require(support_one == EXPECTED_BOUNDARY_SUPPORT_ONE, support_one)
    require(support_two == EXPECTED_BOUNDARY_SUPPORT_TWO, support_two)
    require(support_one_survivors == 0, ("boundary support-one survivor",))
    require(support_two_survivors == 0, ("boundary support-two survivor",))
    return support_one, support_two


def value_literals(index: int, value: int) -> tuple[int, ...]:
    return tuple(
        bit_variable(index, bit)
        if (value >> bit) & 1
        else -bit_variable(index, bit)
        for bit in range(BITS)
    )


def branch_assumptions(first_nonzero_coordinate: int) -> tuple[int, ...]:
    require(2 <= first_nonzero_coordinate <= N - 1, first_nonzero_coordinate)
    literals = list(value_literals(0, 0))
    for coordinate in range(2, first_nonzero_coordinate):
        literals.extend(value_literals(coordinate - 1, 0))
    literals.extend(value_literals(first_nonzero_coordinate - 1, 1))
    return tuple(literals)


def audit_prime_unit_orbits() -> None:
    # Every possible first nonzero residue has a unique multiplicative inverse
    # sending it to 1.  Zero-prefix positions do not change under scaling.
    for residue in range(1, N):
        inverses = tuple(scalar for scalar in range(1, N) if scalar * residue % N == 1)
        require(len(inverses) == 1, ("nonunique inverse", residue, inverses))
    for coordinate in range(2, N):
        assumptions = branch_assumptions(coordinate)
        require(len(assumptions) == BITS * coordinate, (coordinate, len(assumptions)))


def dimacs_digest(
    clauses: tuple[tuple[int, ...], ...], variable_count: int
) -> str:
    rows = [f"p cnf {variable_count} {len(clauses)}\n"]
    rows.extend(" ".join(map(str, clause)) + " 0\n" for clause in clauses)
    return sha256("".join(rows).encode("ascii")).hexdigest()


def branch_digest(
    common_clauses: tuple[tuple[int, ...], ...], coordinate: int
) -> str:
    units = tuple((literal,) for literal in branch_assumptions(coordinate))
    return dimacs_digest(common_clauses + units, EXPECTED_VARIABLES)


def solve_all_branches(
    solver_name: str,
    common_clauses: tuple[tuple[int, ...], ...],
) -> tuple[int, ...]:
    closed: list[int] = []
    with Solver(name=solver_name, bootstrap_with=common_clauses) as solver:
        for coordinate in range(2, N):
            solver.conf_budget(BRANCH_CONFLICT_BUDGET)
            verdict = solver.solve_limited(
                assumptions=branch_assumptions(coordinate),
                expect_interrupt=True,
            )
            require(verdict is not None, (solver_name, coordinate, "budget exhausted"))
            require(verdict is False, (solver_name, coordinate, "SAT model"))
            closed.append(coordinate)
    return tuple(closed)


def main() -> None:
    source_bytes = Path(__file__).read_bytes()
    source_lf = source_bytes.replace(b"\r\n", b"\n")
    require(b"\r" not in source_lf, "source contains a bare carriage return")

    breakpoints, intervals, patterns = build_cells()
    require(N * len(patterns) == EXPECTED_FULL_CANDIDATES, "bad candidate count")
    audit_local_encoding()
    audit_reflection(patterns)
    audit_prime_unit_orbits()
    verify_positive_and_hostile_controls(intervals, patterns)
    boundary_support_one, boundary_support_two = audit_boundary_low_support()
    base_clauses = build_base_clauses()
    cover_clauses = build_representative_cover_clauses(patterns)
    common_clauses = base_clauses + cover_clauses

    require(
        equality_variable(N - 2, N - 1) == EXPECTED_VARIABLES,
        "bad variable count",
    )

    print("n=13 prime-unit orbit micro-staircase rigidity audit")
    print("scope=FINITE-EXACT residue/open-cell model; not LRC")
    print(f"python={platform.python_version()} python_sat={pysat.__version__}")
    print(f"script_lf_sha256={sha256(source_lf).hexdigest()}")
    print(
        f"breakpoints={len(breakpoints)} atomic_open_intervals={len(intervals)} "
        f"distinct_patterns={len(set(patterns))} full_candidates={N * len(patterns)}"
    )
    print("open_cell_samples=one_third,two_thirds agree_on_every_interval")
    print("shift_zero=universally_blocked")
    print("reflection=(s,b)->(-s,12-b); representative_shifts=(1,2,3,4,5,6)")
    print(
        f"residue_encoding=4_bits+equality_atoms_iff variables={EXPECTED_VARIABLES} "
        f"base_clauses={len(base_clauses)} representative_covers={len(cover_clauses)}"
    )
    print("binary_domain_and_equality_truth_tables=PASS")
    print("gauge=v_i-v_1*i; prime_unit_scaling=v->c*v reindexes_nonzero_shifts")
    print("orbit_branches=first_nonzero_coordinate_2_through_12,value_1")
    digest_rows = tuple(
        (coordinate, branch_digest(common_clauses, coordinate))
        for coordinate in range(2, N)
    )
    for coordinate, digest in digest_rows:
        print(f"branch_first_nonzero={coordinate};dimacs_sha256={digest}")
    print(
        "branch_digest_bank_sha256="
        + sha256(repr(digest_rows).encode("ascii")).hexdigest()
    )

    for solver_name in ("cadical195", "maplesat"):
        closed = solve_all_branches(solver_name, common_clauses)
        print(
            f"{solver_name}_budget_per_branch={BRANCH_CONFLICT_BUDGET};"
            f"unsat_first_nonzero_coordinates={closed}"
        )

    print("positive_control=13_of_13_scalar_ramps_block_all_7774_candidates")
    print(
        "hostile_control=(0,1,0,0,0,0,0,0,0,0,0,0);misses=242;"
        "first_miss=shift_1_cell_311_[27/52,61/117)"
    )
    print(
        f"boundary_affine_cover_support_1_scanned={boundary_support_one};survivors=0;"
        f"support_2_scanned={boundary_support_two};survivors=0"
    )
    print("classification=exactly_13_scalar_ramps;normalized_nonzero=UNSAT")
    print("ALL_CHECKS_PASSED")


if __name__ == "__main__":
    main()
