#!/usr/bin/env python3
"""Exact engine/Fraction audit of the saturated strict [2,3] 11/6 control.

The control is not an improvement on the strict [2,2] optimum 7/4.  It only
shows that the next saturated ladder is forceable with one additional seed.
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


def raw(rho):
    rho = Fraction(rho)
    return rho.denominator + rho.numerator, rho.denominator - rho.numerator


def matrix_rank(rows, columns):
    a = [[Fraction(row[c]) for c in columns] for row in rows]
    rr = 0
    for col in range(len(columns)):
        pivot = next((i for i in range(rr, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rr], a[pivot] = a[pivot], a[rr]
        q = a[rr][col]
        a[rr] = [x / q for x in a[rr]]
        for i in range(len(a)):
            if i != rr and a[i][col]:
                q = a[i][col]
                a[i] = [x - q * y for x, y in zip(a[i], a[rr])]
        rr += 1
    return rr


def dense_rows(instance, mode):
    rows = []
    for sparse in instance.generators(mode):
        row = [Fraction(0)] * (2 * instance.n)
        for vertex, label in sparse.items():
            s, d = label[0] + label[1], label[0] - label[1]
            row[2 * vertex], row[2 * vertex + 1] = s, d
        rows.append(tuple(row))
    return tuple(rows)


def instance_23():
    vertices = tuple((i, j) for i in (1, 2) for j in (1, 2, 3))
    fs = [
        {(1,): raw(4)},
        {(1, 1): raw(-1), (1, 2): raw(-8),
         (2, 1): raw(3), (2, 2): raw(-8)},
    ]
    seeds = [
        (vertices[0], raw(-3)),
        (vertices[1], raw(1)),
        (vertices[4], raw(Fraction(37, 10))),
        (vertices[2], raw(-7)),
    ]
    X = [(0, 0)] + list(dict.fromkeys(
        [x for f in fs for x in f.values()] + [x for _, x in seeds]
    ))
    return AKInstance(X, [2, 3], fs, [], seeds)


def instance_22():
    vertices = ((1, 1), (1, 2), (2, 1), (2, 2))
    fs = [
        {(1,): raw(4)},
        {(1, 1): raw(-1), (2, 1): raw(3)},
    ]
    seeds = [
        (vertices[0], raw(-3)),
        (vertices[1], raw(1)),
        (vertices[3], raw(Fraction(37, 10))),
    ]
    X = [(0, 0)] + list(dict.fromkeys(
        [x for f in fs for x in f.values()] + [x for _, x in seeds]
    ))
    return AKInstance(X, [2, 2], fs, [], seeds)


def main():
    gates = 0
    c4 = instance_22()
    require(c4.validate() == [] and c4.score() == Fraction(7, 4), "C4 control invalid")
    c4_ok, _, c4_trace = force(c4, "strict")
    require(c4_ok and c4_trace == [[(2, 1)], [(1, 1), (1, 2), (2, 2)]],
            "C4 strict positive trace changed")
    require(len(c4.generators("strict")) == 7 and len(c4.generators("loose")) == 8,
            "C4 strict/loose accounting changed")
    gates += 3

    inst = instance_23()
    require(inst.validate() == [], "[2,3] input invalid")
    require((inst.m(), inst.r(), inst.n, inst.t(), inst.score()) ==
            (7, 4, 6, 0, Fraction(11, 6)), "[2,3] cost mismatch")
    strict_rows = dense_rows(inst, "strict")
    loose_rows = dense_rows(inst, "loose")
    require((len(strict_rows), len(loose_rows)) == (11, 13), "row accounting mismatch")
    gates += 3

    ok, _, trace = force(inst, "strict")
    expected = [[(2, 1)], [(1, 1), (1, 2), (1, 3), (2, 2), (2, 3)]]
    require(ok and trace == expected, "strict [2,3] trace changed")
    loose_ok, _, loose_trace = force(inst, "loose")
    require(loose_ok and loose_trace == [[(1, 1), (1, 2), (1, 3),
                                          (2, 1), (2, 2), (2, 3)]],
            "loose hostile trace changed")
    gates += 2

    # Engine row order: four seeds, three shared horizontals, four verticals.
    coeff = (-15, 35, -2, 0, -6, -14, 0, 21, 0, -6, 0)
    vector = tuple(sum(Fraction(q) * row[j] for q, row in zip(coeff, strict_rows))
                   for j in range(12))
    target = [Fraction(0)] * 12
    target[7] = 12
    require(vector == tuple(target), "explicit D witness changed")
    full_columns = tuple(range(12))
    after_D = tuple(c for c in full_columns if c not in (6, 7))
    require((matrix_rank(strict_rows, full_columns),
             matrix_rank(strict_rows, tuple(c for c in full_columns if c != 7)),
             matrix_rank(strict_rows, after_D)) == (11, 10, 10),
            "rank profile changed")
    gates += 2

    # Every legal slot deletion destroys eventual success.
    slot_failures = []
    for axis, key in ((0, (1,)), (1, (1, 1)), (1, (1, 2)),
                      (1, (2, 1)), (1, (2, 2))):
        fs = [dict(f) for f in inst.fs]
        del fs[axis][key]
        reduced = AKInstance(inst.X, inst.dims, fs, [], inst.R)
        rok, _, rtrace = force(reduced, "strict")
        require(not rok, f"slot {(axis, key)} dispensable")
        slot_failures.append((axis, key, tuple(tuple(x) for x in rtrace)))
        gates += 1
    for i in range(4):
        reduced = AKInstance(inst.X, inst.dims, inst.fs, [], inst.R[:i] + inst.R[i + 1:])
        rok, _, rtrace = force(reduced, "strict")
        require(not rok, f"seed {i} dispensable")
        slot_failures.append(("seed", i, tuple(tuple(x) for x in rtrace)))
        gates += 1

    semantic = {
        "status": "STRICT_2x3_11_OVER_6_POSITIVE_CONTROL",
        "scope": "PUBLIC_STRICT_PAID_ROWS_ONLY_NOT_OPTIMALITY_NOT_AK",
        "cost": (7, 4, 6, 0),
        "trace": tuple(tuple(x) for x in trace),
        "ranks": (11, 10, 10),
        "coeff": coeff,
        "strict_loose_rows": (11, 13),
        "slot_failures": slot_failures,
    }
    digest = hashlib.sha256(json.dumps(
        semantic, sort_keys=True, separators=(",", ":"), default=str
    ).encode()).hexdigest()
    print("AK_STRICT_2x3_11_OVER_6_POSITIVE_CONTROL_20260823")
    print("status=FINITE_EXACT_POSITIVE_CONTROL;NOT_OPTIMALITY;NO_AK_CLAIM")
    print("dims=(2,3);m=7;r=4;n=6;t=0;score=11/6")
    print("slopes=(cross:4;vertical:-1,-8,3,-8;seeds:A:-3,B:1,E:37/10,C:-7)")
    print("strict_trace=((D,),(A,B,C,E,F));ranks=(11,10,0);p=(0,5)")
    print(f"first_raw_coefficients={coeff};first_target=12*d_D")
    print("positive_control_C4=score_7/4_trace_1+3")
    print("hostile=loose_bank_has_13_rows_for_paid_cost_11_and_fires_all_once")
    print(f"slot_deletion_failures={tuple(slot_failures)}")
    print(f"semantic_sha256={digest}")
    print(f"gates={gates}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
