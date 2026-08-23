#!/usr/bin/env python3
"""Independent SymPy audit and exact 2x2-cell optimality certificate.

Unlike primary.py, this works directly in the original (a,b) coordinates,
does not import the maintained AK engine, and proves the two-seed no-go by a
generic-rank case split over rational slope parameters.
"""

from itertools import combinations, combinations_with_replacement
import hashlib
import json
import sys

import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


V = 4
COORDS = ((1, 1), (1, 2), (2, 2), (2, 1))  # A,B,C,D
RAW_ROWS = (
    ("hAD", 0, 3, (5, -3)),
    ("hBC", 1, 2, (5, -3)),
    ("vAB", 0, 1, (0, 2)),
    ("vDC", 3, 2, (4, -2)),
    ("seedA", 0, None, (-2, 4)),
    ("seedB", 1, None, (2, 0)),
    ("seedC", 2, None, (47, -27)),
)


def original_matrix(spec=RAW_ROWS):
    rows = []
    for _, u, v, label in spec:
        row = [sp.Rational(0)] * (2 * V)
        row[2*u], row[2*u+1] = label
        if v is not None:
            row[2*v], row[2*v+1] = -label[0], -label[1]
        rows.append(row)
    return sp.Matrix(rows)


def fire_set(matrix, live):
    live = tuple(sorted(live))
    cols = tuple(2*v+j for v in live for j in (0, 1))
    projected = matrix[:, cols]
    base_rank = projected.rank()
    fired = []
    for i, v in enumerate(live):
        target = sp.zeros(1, len(cols))
        target[0, 2*i], target[0, 2*i+1] = 1, -1
        if projected.col_join(target).rank() == base_rank:
            fired.append(v)
    return tuple(fired)


def closure(matrix):
    live = set(range(V))
    trace = []
    while live:
        fired = fire_set(matrix, live)
        if not fired:
            break
        trace.append(fired)
        live.difference_update(fired)
    return not live, tuple(trace)


def loose_matrix():
    # The literal suffix-unconstrained first-axis rule adds B-A and C-D
    # transporters but needs only one extra basis row beyond the strict pair.
    spec = list(RAW_ROWS[4:])
    h = (5, -3)
    spec += [
        ("loose_AD", 0, 3, h),
        ("loose_BA", 1, 0, h),
        ("loose_CD", 2, 3, h),
        RAW_ROWS[2],
        RAW_ROWS[3],
    ]
    return original_matrix(spec)


# Symbolic normalized-slope two-seed no-go.
a, b, c, s0, s1 = sp.symbols("a b c s0 s1")
PARAMS = (a, b, c, s0, s1)


def q_vectors():
    return ((0, 0, 0, 0), (b, -b, 0, 0),
            (b, a-b, -a, 0), (b, a-b, c-a, -c))


def seed_form(pos, rho):
    q = q_vectors()[pos]
    return tuple(q[j] + (rho if j == pos else 0) for j in range(4))


def cycle_constraints(positions):
    f0, f1 = seed_form(positions[0], s0), seed_form(positions[1], s1)
    return sp.Matrix(((a-b, b-a, a-c, c-a),
                      tuple(f1[j]-f0[j] for j in range(4))))


def normalized_generator_matrix(positions):
    rows = []
    for (u, v), rho in zip(((0, 3), (1, 2), (0, 1), (3, 2)), (a, a, b, c)):
        row = [0] * 8
        row[2*u], row[2*u+1], row[2*v], row[2*v+1] = 1, rho, -1, -rho
        rows.append(row)
    for v, rho in zip(positions, (s0, s1)):
        row = [0] * 8
        row[2*v], row[2*v+1] = 1, rho
        rows.append(row)
    return sp.Matrix(rows)


def two_seed_case_split():
    cases = support_families = rank_possible = 0
    rank_pairs = {}
    for positions in combinations_with_replacement(range(4), 2):
        W = cycle_constraints(positions)
        A = normalized_generator_matrix(positions)
        for fire in combinations(range(4), 2):
            rest = tuple(v for v in range(4) if v not in fire)
            support_equations = [W[i, v] for i in range(2) for v in rest]
            coefficient_matrix = sp.Matrix([
                [sp.expand(eq).coeff(x) for x in PARAMS]
                for eq in support_equations
            ])
            basis = coefficient_matrix.nullspace()
            cases += 1
            if not basis:
                continue
            support_families += 1
            z = sp.symbols(f"z0:{len(basis)}")
            substitutions = {
                PARAMS[i]: sum(z[j]*basis[j][i] for j in range(len(basis)))
                for i in range(5)
            }
            specialized = A.subs(substitutions)
            terminal_cols = tuple(2*v+j for v in rest for j in (0, 1))
            pair = (specialized.rank(), specialized[:, terminal_cols].rank())
            rank_pairs[pair] = rank_pairs.get(pair, 0) + 1
            if pair == (6, 4):
                rank_possible += 1
    return cases, support_families, rank_possible, tuple(sorted(rank_pairs.items()))


def main():
    gates = 0
    strict = original_matrix()
    require(strict.shape == (7, 8) and strict.rank() == 7, "strict raw matrix rank mismatch")
    require(closure(strict) == (True, ((3,), (0, 1, 2))), "raw-coordinate strict closure mismatch")
    gates += 2

    coefficients = sp.Matrix([[-6, -14, 21, -6, -15, 35, -2]])
    first = coefficients * strict
    target = sp.zeros(1, 8)
    target[0, 6], target[0, 7] = 6, -6
    require(first == target, "independent first witness mismatch")
    gates += 1

    live_after_D = (0, 1, 2)
    terminal_cols = tuple(2*v+j for v in live_after_D for j in (0, 1))
    require(strict[:, terminal_cols].rank() == 6, "terminal raw matrix not full")
    gates += 1

    loose = loose_matrix()
    require(loose.shape == (8, 8), "loose hostile row accounting mismatch")
    require(closure(loose) == (True, ((0, 1, 2, 3),)), "loose hostile trace mismatch")
    gates += 2

    deletion_traces = []
    for i, row in enumerate(RAW_ROWS):
        reduced = strict.copy()
        reduced.row_del(i)
        ok, trace = closure(reduced)
        require(not ok, f"row deletion {row[0]} unexpectedly succeeds")
        deletion_traces.append((row[0], trace))
        gates += 1

    cases, families, possible, rank_pairs = two_seed_case_split()
    require((cases, families, possible) == (60, 60, 0), "two-seed symbolic no-go changed")
    require(sum(count for _, count in rank_pairs) == 60, "case census mismatch")
    gates += 2

    semantic = {
        "status": "STRICT_2x2_OPTIMUM_7_OVER_4_PROVED_EXACT",
        "scope": "Q_NORMALIZED_SLOPES_PUBLIC_STRICT_PAID_ROWS_ONLY",
        "raw_rank": 7,
        "trace": ((3,), (0, 1, 2)),
        "first_target": (6, -6),
        "two_seed_cases": cases,
        "generic_rank_pairs": rank_pairs,
        "deletion_traces": deletion_traces,
        "loose_trace": ((0, 1, 2, 3),),
    }
    semantic_hash = hashlib.sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()
    print("AK_STRICT_C4_7_OVER_4_INDEPENDENT_AUDIT_20260823")
    print("status=INDEPENDENT_PASS;STRICT_2x2_OPTIMUM_7_OVER_4")
    print("raw_coordinate_trace=((D,),(A,B,C));raw_rank=7;terminal_rank=6")
    print("first_raw_combination=(-6,-14,21,-6,-15,35,-2)->6*(1,-1)_D")
    print(f"two_seed_cases={cases};support_families={families};rank_possible={possible}")
    print(f"generic_rank_pairs={rank_pairs}")
    print("hostile=loose_unpaid_extra_row_changes_trace_to_one_round")
    print(f"semantic_sha256={semantic_hash}")
    print(f"gates={gates}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
