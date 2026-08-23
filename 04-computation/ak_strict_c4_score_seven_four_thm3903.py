#!/usr/bin/env python3
"""Exact primary audit of the strict 2x2 Arithmetic-Kakeya 7/4 cell.

This audits only the repository's public paid-row forcing semantics.  It does
not audit the private verifier or prove an Arithmetic-Kakeya bound.
"""

from fractions import Fraction
import hashlib
import json
import sys

sys.path.insert(0, "04-computation")
from ak_forcing_engine import AKInstance, force

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


VERTICES = ((1, 1), (1, 2), (2, 2), (2, 1))  # A,B,C,D cyclically
N = 4
RAW = {
    "hAD": (5, -3),       # rho=4; the same f1 slot also creates hBC
    "vAB": (0, 2),        # rho=-1
    "vDC": (4, -2),       # rho=3
    "seedA": (-2, 4),     # rho=-3
    "seedB": (2, 0),      # rho=1
    "seedC": (47, -27),   # rho=37/10
}


def sd(raw):
    x, y = raw
    return Fraction(x + y), Fraction(x - y)


def edge_row(u, v, raw):
    s, d = sd(raw)
    row = [Fraction(0)] * (2 * N)
    row[2*u], row[2*u+1] = s, d
    row[2*v], row[2*v+1] = -s, -d
    return tuple(row)


def seed_row(v, raw):
    s, d = sd(raw)
    row = [Fraction(0)] * (2 * N)
    row[2*v], row[2*v+1] = s, d
    return tuple(row)


NAMES = ("hAD", "hBC", "vAB", "vDC", "seedA", "seedB", "seedC")
ROWS = (
    edge_row(0, 3, RAW["hAD"]),
    edge_row(1, 2, RAW["hAD"]),
    edge_row(0, 1, RAW["vAB"]),
    edge_row(3, 2, RAW["vDC"]),
    seed_row(0, RAW["seedA"]),
    seed_row(1, RAW["seedB"]),
    seed_row(2, RAW["seedC"]),
)


def rank(rows, columns):
    matrix = [[Fraction(row[c]) for c in columns] for row in rows]
    rr = 0
    for col in range(len(columns)):
        pivot = next((i for i in range(rr, len(matrix)) if matrix[i][col]), None)
        if pivot is None:
            continue
        matrix[rr], matrix[pivot] = matrix[pivot], matrix[rr]
        q = matrix[rr][col]
        matrix[rr] = [x/q for x in matrix[rr]]
        for i in range(len(matrix)):
            if i != rr and matrix[i][col]:
                q = matrix[i][col]
                matrix[i] = [x-q*y for x, y in zip(matrix[i], matrix[rr])]
        rr += 1
        if rr == len(matrix):
            break
    return rr


def fire_set(rows, live):
    columns = [2*v+j for v in sorted(live) for j in (0, 1)]
    full_rank = rank(rows, columns)
    return tuple(v for v in sorted(live)
                 if rank(rows, [c for c in columns if c != 2*v+1]) == full_rank-1)


def closure(rows):
    live = set(range(N))
    trace = []
    while live:
        fired = fire_set(rows, live)
        if not fired:
            break
        trace.append(fired)
        live.difference_update(fired)
    return not live, tuple(trace)


def combination(coefficients):
    return tuple(sum(Fraction(c)*row[j] for c, row in zip(coefficients, ROWS)) for j in range(2*N))


def live_part(vector, live):
    return tuple(vector[2*v+j] for v in live for j in (0, 1))


def raw_for_rho(numerator, denominator=1):
    return (denominator + numerator, denominator - numerator)


def main():
    gates = 0
    for label, raw in RAW.items():
        require(raw != (0, 0) and sum(raw) != 0, f"invalid label {label}")
        gates += 1

    X = [(0, 0)] + list(dict.fromkeys(RAW.values()))
    fs = [
        {(1,): RAW["hAD"]},
        {(1, 1): RAW["vAB"], (2, 1): RAW["vDC"]},
    ]
    seeds = [
        (VERTICES[0], RAW["seedA"]),
        (VERTICES[1], RAW["seedB"]),
        (VERTICES[2], RAW["seedC"]),
    ]
    inst = AKInstance(X, [2, 2], fs, [], seeds)
    require(inst.validate() == [], "benchmark input validation failed")
    require((inst.m(), inst.r(), inst.n, inst.t(), inst.score()) == (4, 3, 4, 0, Fraction(7, 4)), "cost mismatch")
    require(len(inst.generators("strict")) == 7, "strict generator accounting mismatch")
    require(len(inst.generators("loose")) == 8, "loose hostile bank did not add its uncharged row")
    gates += 4

    success, trace = closure(ROWS)
    expected_trace = ((3,), (0, 1, 2))
    require(success and trace == expected_trace, "independent Fraction closure mismatch")
    ok_engine, _, engine_trace = force(inst, "strict")
    engine_indices = tuple(tuple(VERTICES.index(v) for v in wave) for wave in engine_trace)
    require(ok_engine and engine_indices == expected_trace, "maintained strict engine mismatch")
    loose_ok, _, loose_trace = force(inst, "loose")
    require(loose_ok and len(loose_trace) == 1, "loose hostile control changed")
    gates += 3

    # Exact first-round strict witness.  In raw rows it yields 12*d_D.
    first_coeff = (-6, -14, 21, -6, -15, 35, -2)
    first_vector = combination(first_coeff)
    target = [Fraction(0)] * 8
    target[7] = 12
    require(first_vector == tuple(target), "first-round combination mismatch")
    gates += 1

    # Second round: D is already a wildcard, so only A,B,C live coordinates matter.
    second = {
        0: ((1, 0, 0, 0, -1, 0, 0), 14),
        1: ((-2, 0, 7, 0, -5, 7, 0), 28),
        2: ((-6, -14, 21, 14, -15, 35, 0), 28),
    }
    for vertex, (coeff, scale) in second.items():
        vector = combination(coeff)
        desired = [Fraction(0)] * 6
        desired[2*vertex+1] = scale
        require(live_part(vector, (0, 1, 2)) == tuple(desired), f"second witness {vertex} mismatch")
        gates += 1

    full_cols = tuple(range(8))
    without_dD = tuple(c for c in full_cols if c != 7)
    without_D = tuple(range(6))
    require((rank(ROWS, full_cols), rank(ROWS, without_dD), rank(ROWS, without_D)) == (7, 6, 6), "round rank/p profile mismatch")
    require(rank(ROWS, tuple(range(6))) == 6, "terminal restriction not full")
    gates += 2

    # The three normalized cycle equations explain the unique first firing.
    W = (
        (Fraction(5), Fraction(-5), Fraction(1), Fraction(-1)),
        (Fraction(2), Fraction(2), Fraction(0), Fraction(0)),
        (Fraction(2), Fraction(5), Fraction(-3, 10), Fraction(0)),
    )
    kernel = (Fraction(-1, 10), Fraction(1, 10), Fraction(1), Fraction(0))
    require(all(sum(row[j]*kernel[j] for j in range(4)) == 0 for row in W), "cycle-kernel mismatch")
    require(rank(W, tuple(range(4))) == 3, "cycle constraints lost rank")
    require(tuple(i for i, x in enumerate(kernel) if x == 0) == (3,), "cycle kernel does not isolate D")
    gates += 3

    # Every one-row deletion destroys success; every legal paid-slot deletion
    # does too (the f1 slot deletes both shared horizontal rows).
    physical_deletions = []
    for i, name in enumerate(NAMES):
        ok, tr = closure(ROWS[:i] + ROWS[i+1:])
        require(not ok, f"physical row {name} was dispensable")
        physical_deletions.append((name, tr))
        gates += 1
    slots = ((0, 1), (2,), (3,), (4,), (5,), (6,))
    for slot in slots:
        reduced = tuple(row for i, row in enumerate(ROWS) if i not in slot)
        require(not closure(reduced)[0], f"paid slot {slot} was dispensable")
        gates += 1

    # Codimension-one hostile perturbation: changing only rho(seed C) from
    # 37/10 to 47/10 destroys the first firing.
    perturbed = ROWS[:-1] + (seed_row(2, raw_for_rho(47, 10)),)
    require(closure(perturbed) == (False, ()), "perturbation control unexpectedly forces")
    gates += 1

    semantic = {
        "status": "STRICT_2x2_SCORE_7_OVER_4_PROVED_FINITE_EXACT",
        "scope": "PUBLIC_PAID_ROW_STRICT_SEMANTICS_ONLY_NO_AK_CLAIM",
        "score": "7/4",
        "cost": (4, 3, 4, 0),
        "raw_labels": RAW,
        "trace": expected_trace,
        "rank_profile": (7, 6, 6),
        "first_coeff": first_coeff,
        "strict_generators": 7,
        "loose_generators": 8,
        "deletion_traces": physical_deletions,
    }
    semantic_hash = hashlib.sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()
    print("AK_STRICT_C4_7_OVER_4_PRIMARY_20260823")
    print("status=PROVED_FINITE_EXACT_STRICT_CELL;NO_ARITHMETIC_KAKEYA_CLAIM")
    print("dims=(2,2);m=4;r=3;n=4;t=0;score=7/4")
    print("raw_slopes=(h:4,left:-1,right:3,seeds:-3,1,37/10)")
    print("trace=((D,),(A,B,C));round_ranks=(7,6,0);p=(0,3)")
    print(f"first_raw_coefficients={first_coeff};first_target=12*d_D")
    print(f"physical_deletion_traces={tuple(physical_deletions)}")
    print("hostile=loose_bank_has_8_rows_for_paid_cost_7;strict_proof_uses_exactly_7")
    print(f"semantic_sha256={semantic_hash}")
    print(f"gates={gates}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
