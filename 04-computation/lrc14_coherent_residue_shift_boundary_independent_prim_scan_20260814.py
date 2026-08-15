#!/usr/bin/env python3
"""Independent all-packet audit of the coherent-shift boundary census.

This path uses Prim rather than Kruskal and computes pair intersections by
|A intersect B|=|A|+|B|-|A union B| rather than a two-pointer intersection.
"""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from multiprocessing import get_context
from pathlib import Path
import argparse


ROOT = Path(__file__).resolve().parents[1]
MEDIAN = ROOT / "04-computation/lrc14_j7_reflected_median_star_chord_scout_20260801.py"
EXPECTED_MEDIAN = "b0eba7785cd73ecf76d2887f3d36316eae6980a76652058daa2587ecbe031276"
SHIFTS = tuple(range(1, 85))
ASSIGNMENTS = tuple(combinations(range(6), 3))
EXPECTED_FAILURE_SHIFTS = (67, 69, 80, 82, 84)
EXPECTED_FAILURE_MARGINS = (
    Q(-3_689_469_617_499, 1_523_807_086_440_275),
    Q(-109_470_469, 68_745_941_264),
    Q(-390_660_252_142, 18_345_869_001_369),
    Q(-47_427_700_563, 4_633_736_015_480),
    Q(-52_458_412_717, 61_674_995_021_640),
)
EXPECTED_SEMANTIC = "ee616f1066203446c9672ede2490e17d46b2626faf9ddad01ca5287b2c9c67ef"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def file_sha(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")).hexdigest()


def load(name, path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha(MEDIAN) == EXPECTED_MEDIAN, (MEDIAN, file_sha(MEDIAN)))
M = load("coherent_shift_boundary_prim_median", MEDIAN)
R = M.R


def mass(rows):
    return sum((right - left for left, right in R.merge_intervals(tuple(rows))), Q())


def prim_credit(vertices, weights):
    selected = {vertices[0]}
    credit = Q()
    while len(selected) < len(vertices):
        candidates = (
            (weights[a, b], a, b)
            for a in selected for b in vertices
            if b not in selected and (a, b) in weights
        )
        weight, a, b = max(candidates)
        credit += weight
        selected.add(b)
    return credit


def scan_shift(shift):
    tested = 0
    minimum = None
    failures = []
    digest = sha256()
    bodies = M.body_universe()
    require(len(bodies) == 649, len(bodies))
    for body, L in bodies:
        ruler, ranges = R.safe_cell_ranges(body)
        require(ruler == L, (body, ruler, L))
        cells = tuple(j for left, right in ranges for j in range(left, right))
        cell = cells[len(cells) // 2]
        residues = tuple(e + shift for e in body)
        events = {
            (q, a): R.direct_multiplier_arcs(L, q * L - a, cell)
            for q in (3, 5) for a in residues
        }
        singles = {vertex: mass(event) for vertex, event in events.items()}
        weights = {}
        for a in residues:
            for b in residues:
                left_vertex, right_vertex = (3, a), (5, b)
                weight = (
                    singles[left_vertex] + singles[right_vertex]
                    - mass(events[left_vertex] + events[right_vertex])
                )
                require(weight >= 0, (shift, body, a, b, weight))
                weights[left_vertex, right_vertex] = weight
                weights[right_vertex, left_vertex] = weight
        for indices in ASSIGNMENTS:
            left = tuple(residues[index] for index in indices)
            right = tuple(residues[index] for index in range(6) if index not in indices)
            vertices = tuple((3, a) for a in left) + tuple((5, b) for b in right)
            credit = prim_credit(vertices, weights)
            debt = sum((singles[vertex] for vertex in vertices), Q()) - Q(6, 7)
            margin = credit - debt
            row = (margin, body, L, cell, shift, left, right, credit, debt)
            digest.update((repr(row) + "\n").encode())
            tested += 1
            if minimum is None or row < minimum:
                minimum = row
            if margin <= 0:
                failures.append(row)
    require(tested == 12_980, (shift, tested))
    return shift, minimum, tuple(failures), digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=6)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    if args.processes == 1:
        results = tuple(map(scan_shift, SHIFTS))
    else:
        with ProcessPoolExecutor(
            max_workers=args.processes, mp_context=get_context("fork")
        ) as pool:
            results = tuple(pool.map(scan_shift, SHIFTS))
    require(tuple(row[0] for row in results) == SHIFTS, "shift order")
    all_failures = tuple(
        failure for shift, minimum, failures, semantic in results
        for failure in failures
    )
    failure_shifts = tuple(row[4] for row in all_failures)
    require(failure_shifts == EXPECTED_FAILURE_SHIFTS, failure_shifts)
    require(tuple(row[0] for row in all_failures) == EXPECTED_FAILURE_MARGINS,
            all_failures)
    safe_minimum = min(
        (minimum for shift, minimum, failures, semantic in results if shift < 67),
        key=lambda row: row[0],
    )
    require(safe_minimum[0] == Q(6_721_694_301_305, 502_164_203_501_961),
            safe_minimum)
    digest = sha256()
    for result in results:
        digest.update((repr(result) + "\n").encode())
    semantic = digest.hexdigest()
    require(semantic == EXPECTED_SEMANTIC, semantic)
    print("LRC14 COHERENT-SHIFT BOUNDARY INDEPENDENT PRIM/UNION AUDIT")
    print("bodies", 649, "shifts", (SHIFTS[0], SHIFTS[-1]),
          "assignments", len(ASSIGNMENTS), "packets", len(SHIFTS) * 12_980)
    print("safe_prefix_packets", 66 * 12_980, "safe_prefix_minimum", safe_minimum)
    print("failure_shifts", failure_shifts)
    for failure in all_failures:
        print("failure", failure)
    print("semantic_sha256", semantic)
    print("status=INDEPENDENTLY VERIFIED-EXACT full census; Prim and union-identity overlap path")


if __name__ == "__main__":
    main()
