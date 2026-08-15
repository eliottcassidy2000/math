#!/usr/bin/env python3
"""Exact census for the first coherent-shift cross-K3,3 obstruction.

For shifts c=1,...,84, all 649 upper-median
bodies, and all C(6,3)=20 assignments of levels (3,3,3,5,5,5), compute
the exact Hunter margin of the maximum spanning tree in the restricted
cross-K3,3 overlap graph.  At every restricted failure, independently
compute the unrestricted full-K6 maximum-tree margin and literal union.
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
ASSIGNMENTS = tuple(combinations(range(6), 3))
LAST_SHIFT = 84
EXPECTED_CENSUS = "100f7ca9f7a07d0e36685483958092768f6b3422ff25a01de2fc6ef8b7d6ce52"
EXPECTED_FAILURE_SHIFTS = (67, 69, 80, 82, 84)
EXPECTED_SAFE_PREFIX_MINIMUM = Q(6_721_694_301_305, 502_164_203_501_961)


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
M = load("coherent_shift_boundary_median", MEDIAN)
R = M.R


def mass(rows):
    return sum((right - left for left, right in R.merge_intervals(tuple(rows))), Q())


def overlap(first, second):
    i = j = 0
    answer = Q()
    while i < len(first) and j < len(second):
        answer += max(Q(), min(first[i][1], second[j][1])
                      - max(first[i][0], second[j][0]))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return answer


def maximum_tree(vertices, admissible_edges):
    edges = sorted(admissible_edges, reverse=True)
    parent = {vertex: vertex for vertex in vertices}

    def root(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    credit = Q()
    tree = []
    for weight, a, b in edges:
        x, y = root(a), root(b)
        if x == y:
            continue
        parent[x] = y
        credit += weight
        tree.append((a, b, weight))
        if len(tree) == len(vertices) - 1:
            return credit, tuple(tree)
    raise RuntimeError("admissible graph disconnected")


def scan_shift(shift):
    bodies = M.body_universe()
    require(len(bodies) == 649, len(bodies))
    tested = 0
    minimum = None
    failures = []
    digest = sha256()
    for body, L in bodies:
        ruler, ranges = R.safe_cell_ranges(body)
        require(ruler == L, (body, ruler, L))
        cells = tuple(j for left, right in ranges for j in range(left, right))
        cell = cells[len(cells) // 2]
        residues = tuple(e + shift for e in body)
        require(min(residues) > 0, (shift, body, residues))
        require(max(residues) < 3 * L, (shift, body, L, residues))
        events = {
            (q, a): R.direct_multiplier_arcs(L, q * L - a, cell)
            for q in (3, 5)
            for a in residues
        }
        singles = {key: mass(value) for key, value in events.items()}
        cross_weights = {
            (a, b): overlap(events[3, a], events[5, b])
            for a in residues for b in residues
        }
        for indices in ASSIGNMENTS:
            left = tuple(residues[index] for index in indices)
            right = tuple(residues[index] for index in range(6) if index not in indices)
            vertices = tuple((3, a) for a in left) + tuple((5, b) for b in right)
            cross_edges = tuple(
                (cross_weights[a, b], (3, a), (5, b))
                for a in left for b in right
            )
            cross_credit, cross_tree = maximum_tree(vertices, cross_edges)
            singleton_sum = sum((singles[vertex] for vertex in vertices), Q())
            debt = singleton_sum - Q(6, 7)
            cross_margin = cross_credit - debt
            row = (
                cross_margin, body, L, cell, shift, left, right,
                singleton_sum, debt, cross_credit, cross_tree,
            )
            tested += 1
            digest.update((repr(row) + "\n").encode())
            if minimum is None or row < minimum:
                minimum = row
            if cross_margin <= 0:
                all_edges = tuple(
                    (overlap(events[a], events[b]), a, b)
                    for a, b in combinations(vertices, 2)
                )
                full_credit, full_tree = maximum_tree(vertices, all_edges)
                full_margin = full_credit - debt
                union = mass(tuple(
                    interval for vertex in vertices for interval in events[vertex]
                ))
                failures.append((
                    row, full_credit, full_margin, full_tree,
                    union, union - Q(6, 7), all_edges,
                ))
    require(tested == 12_980, (shift, tested))
    require(minimum is not None, shift)
    return shift, tested, minimum, tuple(failures), digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=6)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    shifts = tuple(range(1, LAST_SHIFT + 1))
    if args.processes == 1:
        results = tuple(map(scan_shift, shifts))
    else:
        with ProcessPoolExecutor(
            max_workers=args.processes, mp_context=get_context("fork")
        ) as pool:
            results = tuple(pool.map(scan_shift, shifts))
    require(tuple(row[0] for row in results) == shifts, "shift order")
    census_digest = sha256()
    total_packets = 0
    all_failures = []
    for result in results:
        shift, tested, minimum, failures, semantic = result
        census_digest.update((repr(result) + "\n").encode())
        total_packets += tested
        all_failures.extend((shift, failure) for failure in failures)
    failure_shifts = tuple(shift for shift, _ in all_failures)
    first_failure = min(failure_shifts, default=None)
    semantic = census_digest.hexdigest()
    require(total_packets == 1_090_320, total_packets)
    require(failure_shifts == EXPECTED_FAILURE_SHIFTS, failure_shifts)
    require(first_failure == 67, first_failure)
    require(semantic == EXPECTED_CENSUS, semantic)
    safe_prefix = tuple(
        (minimum[0], shift, minimum)
        for shift, tested, minimum, failures, shift_semantic in results
        if shift < first_failure
    )
    safe_prefix_minimum, safe_prefix_shift, safe_prefix_row = min(safe_prefix)
    require(safe_prefix_minimum == EXPECTED_SAFE_PREFIX_MINIMUM,
            (safe_prefix_shift, safe_prefix_row))
    require(safe_prefix_shift == 56, (safe_prefix_shift, safe_prefix_row))
    require(safe_prefix_row[1:7] == (
        (1, 2, 3, 4, 6, 12), 168, 90, 56,
        (58, 60, 62), (57, 59, 68),
    ), safe_prefix_row)
    failure_summaries = []
    for shift, failure in all_failures:
        row, full_credit, full_margin, full_tree, union, union_margin, all_edges = failure
        cross_margin, body, L, cell = row[:4]
        debt, cross_credit, cross_tree = row[8], row[9], row[10]
        within_edges = tuple(edge for edge in full_tree if edge[0][0] == edge[1][0])
        bridge_edges = tuple(edge for edge in full_tree if edge[0][0] != edge[1][0])
        require(len(within_edges) == 4 and len(bridge_edges) == 1,
                (shift, full_tree))
        within_credit = sum((edge[2] for edge in within_edges), Q())
        bridge_credit = bridge_edges[0][2]
        require(within_credit + bridge_credit == full_credit,
                (shift, within_credit, bridge_credit, full_credit))
        require(within_credit > debt, (shift, within_credit, debt))
        require(full_margin > 0, (shift, full_margin))
        require(union_margin < 0, (shift, union_margin))
        phase = (-shift * cell) % L
        if 2 * phase > L:
            phase -= L
        failure_summaries.append((
            shift, phase, row[1], row[2], row[3], row[5], row[6],
            cross_margin, debt, cross_credit, cross_tree,
            within_credit, bridge_credit, full_margin, full_tree,
            union, union_margin,
        ))
    first = all_failures[0][1]
    require(first[0][0:7] == (
        Q(-3_689_469_617_499, 1_523_807_086_440_275),
        (1, 2, 3, 4, 6, 12), 168, 90, 67,
        (68, 70, 79), (69, 71, 73),
    ), first[0][:7])
    require(first[2] == Q(876_026_208, 3_635_859_955), first[2])
    require(first[5] == Q(-4_849_929_096, 18_179_299_775), first[5])
    print("LRC14 COHERENT RESIDUE-SHIFT FIRST-BOUNDARY EXACT CENSUS")
    print("dependency", MEDIAN.relative_to(ROOT), EXPECTED_MEDIAN)
    print("bodies", 649, "assignments_per_body", len(ASSIGNMENTS),
          "shifts", (shifts[0], shifts[-1]), "packets", total_packets)
    print("safe_prefix", (1, first_failure - 1), "packets",
          (first_failure - 1) * 12_980, "cross_failures", 0)
    print("safe_prefix_minimum", (safe_prefix_shift, safe_prefix_row))
    print("failure_shifts", failure_shifts)
    print("postboundary_failure_counts", tuple(
        (shift, len(failures))
        for shift, tested, minimum, failures, shift_semantic in results
        if shift >= first_failure
    ))
    for summary in failure_summaries:
        print("failure_summary", summary)
    print("census_semantic_sha256", semantic)
    print("conclusion=c=67 is the unique first restricted failure; every c<=84 packet closes by cross-K3,3 or, at exactly five failures, full-K6 Hunter")
    print("status=FINITE-EXACT census; separate literal audits diagnose every restricted failure")


if __name__ == "__main__":
    main()
