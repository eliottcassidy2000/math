#!/usr/bin/env python3
"""Exact CP-SAT search for arbitrary-rank literal-wall responders to THM-4303.

Geometry and every positive witness are re-evaluated directly with integer
wall coordinates; no rank truncation is used.
"""

from __future__ import annotations

import argparse
import collections
import math
import pathlib
import sys
from ortools.sat.python import cp_model


POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
PAIRS = ((96, 595), (100, 595), (210, 595))
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def fnv_add(state: int, word: int) -> int:
    for shift in range(0, 64, 8):
        state ^= (word >> shift) & 0xFF
        state = (state * FNV_PRIME) & MASK64
    return state


def safe_midpoint(speed: int, grid: int, left: int, right: int) -> bool:
    residue = speed * (left + right) % (2 * grid)
    return 7 * residue >= grid and 7 * residue <= 13 * grid


def geometry(q: int, r: int) -> tuple[int, dict[int, int], int, int]:
    grid = 1
    for speed in POOL + (q, r):
        grid = math.lcm(grid, 14 * speed)
    walls = {0, grid}
    for speed in POOL + (q, r):
        unit = grid // (14 * speed)
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    walls = sorted(walls)
    classes: dict[int, int] = collections.defaultdict(int)
    pair_ticks = 0
    for left, right in zip(walls, walls[1:]):
        if not safe_midpoint(q, grid, left, right):
            continue
        if not safe_midpoint(r, grid, left, right):
            continue
        failure = 0
        for bit, speed in enumerate(POOL):
            if not safe_midpoint(speed, grid, left, right):
                failure |= 1 << bit
        width = right - left
        pair_ticks += width
        classes[failure] += width
    return grid, dict(classes), len(walls) - 1, pair_ticks


def load_failures(path: pathlib.Path) -> list[tuple[int, int]]:
    lines = path.read_text().splitlines()
    assert lines[0] == "q,r,body_hex"
    failures = []
    for line in lines[1:]:
        q, r, body = line.split(",")
        assert int(r) == 595
        body_int = int(body, 16)
        assert body_int.bit_count() == 9
        failures.append((int(q), body_int))
    assert collections.Counter(q for q, _ in failures) == {96: 116, 100: 13, 210: 16}
    return failures


def exact_margin(mask: int, grid: int, classes: dict[int, int]) -> tuple[int, int]:
    mass = sum(width for failure, width in classes.items() if failure & ~mask == 0)
    return mass, 63 * mass - 4 * grid


def bits(mask: int):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def make_and(model: cp_model.CpModel, xs, name: str):
    """Return boolean exactly equal to conjunction xs."""
    xs = list(xs)
    if not xs:
        return model.new_constant(1)
    y = model.new_bool_var(name)
    for x in xs:
        model.add(y <= x)
    model.add(y >= sum(xs) - len(xs) + 1)
    return y


def solve(path: pathlib.Path, mode: str, seconds: float, workers: int,
          hint: list[int], num_masks: int, selected_q: set[int],
          force_q210_all: bool, first_solution: bool, presolve: bool,
          symmetry: int, pure_feasibility: bool, unsat_core: bool,
          quiet: bool):
    all_failures = load_failures(path)
    failures = [(q, body) for q, body in all_failures if q in selected_q]
    geometries = [geometry(*pair) for pair in PAIRS]
    expected = (
        (36482318832960, 2453, 8515, 26803336285440),
        (91205797082400, 2519, 8523, 67008340713600),
        (18241159416480, 2543, 8743, 13412617218000),
    )
    for p, ((grid, classes, cells, ticks), exp) in enumerate(zip(geometries, expected)):
        assert (grid, len(classes), cells, ticks) == exp
        print(f"GEOMETRY {PAIRS[p][0]},595 GRID {grid} CELLS {cells} CLASSES {len(classes)} PAIR_TICKS {ticks}", flush=True)

    all_classes = sorted(set().union(*(classes for _, classes, _, _ in geometries)))
    assert len(all_classes) == 2783
    class_index = {f: i for i, f in enumerate(all_classes)}

    model = cp_model.CpModel()
    x = [[model.new_bool_var(f"x_{i}_{j}") for j in range(30)] for i in range(num_masks)]
    values = [sum((1 << j) * x[i][j] for j in range(30)) for i in range(num_masks)]
    for i in range(num_masks - 1):
        model.add(values[i] <= values[i + 1])

    # Exact subset indicators for every literal-wall failure class.
    contained = []
    for i in range(num_masks):
        row = []
        for ci, failure in enumerate(all_classes):
            row.append(make_and(model, (x[i][j] for j in bits(failure)), f"y_{i}_{ci}"))
        contained.append(row)

    active = [[model.new_bool_var(f"a_{i}_{p}") for p in range(3)] for i in range(num_masks)]
    for i in range(num_masks):
        for p, (grid, classes, _, _) in enumerate(geometries):
            threshold = (4 * grid + 62) // 63
            mass = sum(width * contained[i][class_index[failure]] for failure, width in classes.items())
            model.add(mass >= threshold).only_enforce_if(active[i][p])
            # Exact activity, useful for auditing signatures and propagation.
            model.add(mass <= threshold - 1).only_enforce_if(active[i][p].Not())

    if mode == "all":
        for i in range(num_masks):
            for p in range(3):
                model.add(active[i][p] == 1)
    elif mode == "general":
        # A scoped all-rank two-responder obstruction for q=210 implies that
        # every member of a three-cover must be active and useful there.
        if force_q210_all:
            for i in range(num_masks):
                model.add(active[i][2] == 1)
    else:
        raise ValueError(mode)

    hits = []
    assumptions = []
    assumption_rows = {}
    for ordinal, (q, body) in enumerate(failures):
        p = next(i for i, pair in enumerate(PAIRS) if pair[0] == q)
        body_bits = list(bits(body))
        body_hits = []
        for i in range(num_masks):
            disjoint = make_and(model, (x[i][j].Not() for j in body_bits), f"d_{ordinal}_{i}")
            hit = make_and(model, (active[i][p], disjoint), f"h_{ordinal}_{i}")
            body_hits.append(hit)
        cover_constraint = model.add(sum(body_hits) >= 1)
        if unsat_core:
            assumption = model.new_bool_var(f"cover_assumption_{ordinal}")
            cover_constraint.only_enforce_if(assumption)
            model.add_assumption(assumption)
            assumptions.append(assumption)
            assumption_rows[assumption.index] = (ordinal, q, body)
        hits.append(body_hits)

    # Search compact responders first; this is not a rank restriction.
    if not pure_feasibility:
        model.minimize(sum(x[i][j] for i in range(num_masks) for j in range(30)))
    if hint:
        assert len(hint) == num_masks
        for i, mask in enumerate(sorted(hint)):
            for j in range(30):
                model.add_hint(x[i][j], (mask >> j) & 1)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.log_search_progress = not quiet
    solver.parameters.cp_model_presolve = presolve
    solver.parameters.symmetry_level = symmetry
    solver.parameters.random_seed = 4303
    solver.parameters.stop_after_first_solution = first_solution
    print(
        f"SEARCH MASKS {num_masks} FAMILIES {','.join(map(str, sorted(selected_q)))} "
        f"PURE_FEASIBILITY {int(pure_feasibility)} WORKERS {workers} "
        f"PRESOLVE {int(presolve)} SYMMETRY {symmetry}", flush=True)
    print(f"MODEL VARIABLES {len(model.proto.variables)} CONSTRAINTS {len(model.proto.constraints)}", flush=True)
    status = solver.solve(model)
    print(f"STATUS {solver.status_name(status)} WALL {solver.wall_time:.6f} BRANCHES {solver.num_branches} CONFLICTS {solver.num_conflicts}", flush=True)
    if status == cp_model.INFEASIBLE and unsat_core:
        core = solver.sufficient_assumptions_for_infeasibility()
        print(f"UNSAT_CORE SIZE {len(core)}", flush=True)
        for literal in core:
            ordinal, q, body = assumption_rows[literal]
            print(f"CORE {ordinal} {q},595 {body:08x}", flush=True)
    if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        return 2 if status == cp_model.UNKNOWN else 1

    masks = [sum((1 << j) for j in range(30) if solver.value(x[i][j])) for i in range(num_masks)]
    print("MASKS " + " ".join(f"{m:08x}" for m in masks), flush=True)
    joined = collections.Counter()
    witness_fnv = FNV_OFFSET
    response_fnv = FNV_OFFSET
    for i, mask in enumerate(masks):
        witness_fnv = fnv_add(witness_fnv, mask)
        sig = []
        margins = []
        cover = []
        for p, (grid, classes, _, _) in enumerate(geometries):
            mass, ticks = exact_margin(mask, grid, classes)
            sig.append(ticks >= 0)
            margins.append((mass, ticks, grid))
            q = PAIRS[p][0]
            cover.append(sum(q0 == q and not (mask & body) for q0, body in failures) if sig[-1] else 0)
        print(
            f"MASK {mask:08x} RANK {mask.bit_count()} LABELS {','.join(map(str,bits(mask)))} "
            f"SIG {''.join('1' if z else '0' for z in sig)} COVER {cover[0]},{cover[1]},{cover[2]} "
            + " ".join(f"MARGIN_{PAIRS[p][0]} {margins[p][1]}/{63*margins[p][2]} MASS {margins[p][0]}" for p in range(3)),
            flush=True,
        )
    uncovered = []
    incidences = 0
    for ordinal, (q, body) in enumerate(all_failures):
        p = next(i for i, pair in enumerate(PAIRS) if pair[0] == q)
        responders = []
        for i, mask in enumerate(masks):
            _, ticks = exact_margin(mask, geometries[p][0], geometries[p][1])
            if ticks >= 0 and not (mask & body):
                responders.append(i)
                incidences += 1
                response_fnv = fnv_add(response_fnv, mask)
                response_fnv = fnv_add(response_fnv, q)
                response_fnv = fnv_add(response_fnv, 595)
                response_fnv = fnv_add(response_fnv, ordinal)
                response_fnv = fnv_add(response_fnv, body)
        if responders:
            joined[q] += 1
        else:
            uncovered.append((ordinal, q, body))
    print(f"JOINED {joined[96]},{joined[100]},{joined[210]} TOTAL {sum(joined.values())} UNCOVERED {len(uncovered)} INCIDENCES {incidences}", flush=True)
    print(f"WITNESS_FNV {witness_fnv:016x} RESPONSE_FNV {response_fnv:016x}", flush=True)
    for ordinal, q, body in uncovered:
        print(f"UNCOVERED {ordinal} {q},595 {body:08x}", flush=True)
    assert not [row for row in uncovered if row[1] in selected_q]
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("failures", type=pathlib.Path)
    parser.add_argument("--mode", choices=("all", "general"), default="all")
    parser.add_argument("--seconds", type=float, default=600)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--masks", type=int, default=3)
    parser.add_argument("--families", nargs="+", type=int, choices=(96, 100, 210), default=[96, 100, 210])
    parser.add_argument("--force-q210-all", action="store_true")
    parser.add_argument("--first-solution", action="store_true")
    parser.add_argument("--no-presolve", action="store_true")
    parser.add_argument("--symmetry", type=int, choices=(0, 1, 2, 3), default=2)
    parser.add_argument("--pure-feasibility", action="store_true")
    parser.add_argument("--unsat-core", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--hint", nargs="+", type=lambda s: int(s, 16), default=[])
    args = parser.parse_args()
    return solve(args.failures, args.mode, args.seconds, args.workers,
                 args.hint, args.masks, set(args.families),
                 args.force_q210_all, args.first_solution,
                 not args.no_presolve, args.symmetry,
                 args.pure_feasibility, args.unsat_core, args.quiet)


if __name__ == "__main__":
    sys.exit(main())
