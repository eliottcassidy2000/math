#!/usr/bin/env python3
"""Independent raw-coordinate SymPy audit of the strict [2,3] 11/6 control."""

import hashlib
import json
import sys
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def row(u, v, label, vertices=6):
    result = [sp.Integer(0)] * (2 * vertices)
    result[2 * u], result[2 * u + 1] = label
    if v is not None:
        result[2 * v], result[2 * v + 1] = -label[0], -label[1]
    return result


def fire_set(matrix, live):
    live = tuple(sorted(live))
    columns = tuple(2 * v + j for v in live for j in (0, 1))
    projected = matrix[:, columns]
    rank = projected.rank()
    fired = []
    for i, v in enumerate(live):
        target = sp.zeros(1, len(columns))
        target[0, 2 * i], target[0, 2 * i + 1] = 1, -1
        if projected.col_join(target).rank() == rank:
            fired.append(v)
    return tuple(fired)


def closure(matrix, vertices):
    live = set(range(vertices))
    trace = []
    while live:
        fired = fire_set(matrix, live)
        if not fired:
            break
        trace.append(fired)
        live.difference_update(fired)
    return not live, tuple(trace)


def strict_23():
    h = (5, -3)
    return sp.Matrix([
        row(0, 3, h), row(1, 4, h), row(2, 5, h),
        row(0, 1, (0, 2)), row(1, 2, (-7, 9)),
        row(3, 4, (4, -2)), row(4, 5, (-7, 9)),
        row(0, None, (-2, 4)), row(1, None, (2, 0)),
        row(4, None, (47, -27)), row(2, None, (-6, 8)),
    ])


def loose_23():
    # A-D plus within-layer transporters B-A,E-D,C-A,F-D span the literal
    # complete-bipartite first-axis bank.  The remaining eight rows are the
    # four verticals and four seeds from the strict matrix.
    h = (5, -3)
    loose_horizontal = [
        row(0, 3, h), row(1, 0, h), row(4, 3, h),
        row(2, 0, h), row(5, 3, h),
    ]
    strict = strict_23()
    return sp.Matrix(loose_horizontal + [list(strict.row(i)) for i in range(3, 11)])


def strict_22():
    h = (5, -3)
    return sp.Matrix([
        row(0, 2, h, 4), row(1, 3, h, 4),
        row(0, 1, (0, 2), 4), row(2, 3, (4, -2), 4),
        row(0, None, (-2, 4), 4), row(1, None, (2, 0), 4),
        row(3, None, (47, -27), 4),
    ])


def main():
    gates = 0
    c4 = strict_22()
    require(c4.shape == (7, 8) and c4.rank() == 7, "C4 raw matrix changed")
    require(closure(c4, 4) == (True, ((2,), (0, 1, 3))), "C4 positive trace changed")
    gates += 2

    strict = strict_23()
    require(strict.shape == (11, 12) and strict.rank() == 11, "strict matrix changed")
    trace = closure(strict, 6)
    require(trace == (True, ((3,), (0, 1, 2, 4, 5))), "strict closure changed")
    after_D = tuple(c for c in range(12) if c not in (6, 7))
    require(strict[:, after_D].rank() == 10, "terminal restriction not full")
    gates += 3

    coeff = sp.Matrix([[-6, -14, 0, 21, 0, -6, 0, -15, 35, -2, 0]])
    target = sp.zeros(1, 12)
    target[0, 6], target[0, 7] = 6, -6
    require(coeff * strict == target, "explicit first witness changed")
    gates += 1

    loose = loose_23()
    require(loose.shape == (13, 12) and loose.rank() == 12, "loose matrix changed")
    require(closure(loose, 6) == (True, ((0, 1, 2, 3, 4, 5),)),
            "loose hostile closure changed")
    gates += 2

    deletion_traces = []
    # Physical-row deletion is stronger than legal-slot deletion for the
    # three tied horizontal rows.
    for i in range(11):
        reduced = strict.copy()
        reduced.row_del(i)
        ok, waves = closure(reduced, 6)
        require(not ok, f"physical row {i} dispensable")
        deletion_traces.append((i, waves))
        gates += 1

    semantic = {
        "status": "INDEPENDENT_STRICT_2x3_11_OVER_6_CONTROL_PASS",
        "scope": "RAW_Q_MATRIX_NO_ENGINE_NO_OPTIMALITY_NO_AK",
        "strict_trace": trace[1],
        "ranks": (strict.rank(), strict[:, after_D].rank()),
        "loose_shape_rank": (loose.shape, loose.rank()),
        "deletion_traces": deletion_traces,
    }
    digest = hashlib.sha256(json.dumps(
        semantic, sort_keys=True, separators=(",", ":"), default=str
    ).encode()).hexdigest()
    print("AK_STRICT_2x3_11_OVER_6_INDEPENDENT_AUDIT_20260823")
    print("status=INDEPENDENT_FINITE_EXACT_PASS;POSITIVE_CONTROL_ONLY;NO_AK_CLAIM")
    print("raw_trace=((D,),(A,B,C,E,F));raw_rank=11;terminal_rank=10")
    print("first_raw_combination=(-6,-14,0,21,0,-6,0,-15,35,-2,0)->6*(1,-1)_D")
    print("positive_control_C4=raw_rank_7_trace_1+3")
    print("hostile=loose_raw_basis_13x12_rank12_fires_all_once")
    print(f"physical_deletion_traces={tuple(deletion_traces)}")
    print(f"semantic_sha256={digest}")
    print(f"gates={gates}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
