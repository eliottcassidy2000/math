#!/usr/bin/env python3
import csv
import hashlib
import sys

from ortools.sat.python import cp_model

EXPECTED_SHA256 = "cfe73a7d08795dde301af82629b3a6003b5afd6c889ca8da017826fd6afd4e21"


def load(path):
    raw = open(path, "rb").read()
    assert hashlib.sha256(raw).hexdigest() == EXPECTED_SHA256
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            p = int(r["w0"], 16) | (int(r["w1"], 16) << 64) | (int(r["w2"], 16) << 128)
            assert p.bit_count() == int(r["cover"])
            rows.append((p, int(r["least_mask_hex"], 16)))
    assert len(rows) == 2620
    assert all(p and p >> 176 == 0 for p, _ in rows)
    return rows


def maximal(rows):
    out = []
    for i, x in enumerate(rows):
        p = x[0]
        if not any(i != j and p != y[0] and p & ~y[0] == 0
                   for j, y in enumerate(rows)):
            out.append(x)
    return out


def configure(solver):
    solver.parameters.num_search_workers = 8
    solver.parameters.random_seed = 0
    solver.parameters.max_time_in_seconds = 300


def full_cover(rows):
    m = cp_model.CpModel()
    x = [m.new_bool_var(f"f{i}") for i in range(len(rows))]
    for bit in range(176):
        m.add(sum(x[i] for i, (p, _) in enumerate(rows) if p >> bit & 1) >= 1)
    m.minimize(sum(x))
    s = cp_model.CpSolver()
    configure(s)
    st = s.solve(m)
    chosen = [i for i, v in enumerate(x) if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) and s.value(v)]
    return s, st, chosen


def exact14_residual(rows):
    m = cp_model.CpModel()
    x = [m.new_bool_var(f"r{i}") for i in range(len(rows))]
    for bit in range(101):
        m.add(sum(x[i] for i, (p, _) in enumerate(rows) if p >> bit & 1) >= 1)
    m.add(sum(x) <= 14)
    missed = []
    for bit in range(101, 176):
        y = m.new_bool_var(f"miss{bit}")
        m.add(sum(x[i] for i, (p, _) in enumerate(rows) if p >> bit & 1) + y >= 1)
        missed.append(y)
    m.minimize(sum(missed))
    s = cp_model.CpSolver()
    configure(s)
    st = s.solve(m)
    chosen = [i for i, v in enumerate(x) if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) and s.value(v)]
    missed_bits = [101 + i for i, y in enumerate(missed) if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) and s.value(y)]
    return s, st, chosen, missed_bits


def print_witness(tag, rows, chosen):
    union = 0
    for i in chosen:
        union |= rows[i][0]
    print(tag, "COUNT", len(chosen), "COVER636", (union & ((1 << 101) - 1)).bit_count(),
          "COVER632", (union >> 101).bit_count())
    print(tag, "CLASSES", ",".join(map(str, chosen)))
    print(tag, "MASKS", ",".join(f"{rows[i][1]:08x}" for i in chosen))


def main():
    assert len(sys.argv) == 2
    all_rows = load(sys.argv[1])
    union = 0
    for p, _ in all_rows:
        union |= p
    absent = [bit for bit in range(176) if not (union >> bit) & 1]
    assert absent == [158]
    rows = maximal(all_rows)
    print("LRC14_176_EXACT_COVER_V1", flush=True)
    print("ATLAS_SHA256", EXPECTED_SHA256, "ALL", len(all_rows), "MAXIMAL", len(rows), flush=True)
    print("RESPONSE_UNION_COVER", union.bit_count(), "ABSENT_BITS", ",".join(map(str, absent)), flush=True)
    s, st, chosen, misses = exact14_residual(rows)
    print("EXACT14_STATUS", s.status_name(st), "MISSES", len(misses), "BOUND", s.best_objective_bound, flush=True)
    print_witness("EXACT14", rows, chosen)
    print("EXACT14_MISSED_BITS", ",".join(map(str, misses)))
    print("EXACT14_STATS", s.response_stats().strip().replace("\n", " | "), flush=True)
    print("FULL_STATUS IMPOSSIBLE ZERO_RESPONSE_BIT_158", flush=True)
    assert st == cp_model.OPTIMAL and len(chosen) == 14 and len(misses) == 37
    print("VERDICT PASS FINITE_EXACT_ROWAWARE_INTEGER_SEARCH_AND_ZERO_RESPONSE_OBSTRUCTION", flush=True)


if __name__ == "__main__":
    main()
