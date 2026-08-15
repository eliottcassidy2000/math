#!/usr/bin/env python3
"""Exact multi-cell rescue audit for selected-cell coherent-shift failures.

Discover every nonpositive full-K6 Hunter margin at the deterministic
upper-median cell for c=156..491 through the pinned primary census.  For each
such packet, independently rebuild all pair overlaps through the union
identity on every body-safe cell and optimize the full K6 tree by Prim.

This is a finite packet/cell audit.  It does not change the body bank, choose
different cells for different edges, or claim an arbitrary-residue theorem.
"""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from multiprocessing import get_context
from pathlib import Path
import argparse
import sys


ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "04-computation/lrc14_coherent_residue_shift_boundary_census_20260814.py"
EXPECTED_PRIMARY = "5c69fbd209c332d76cc51814bc1ef126d8cc3dea0e9c9d16f00ccf095d6b9a7e"
MEDIAN = ROOT / "04-computation/lrc14_j7_reflected_median_star_chord_scout_20260801.py"
EXPECTED_MEDIAN = "b0eba7785cd73ecf76d2887f3d36316eae6980a76652058daa2587ecbe031276"
SHIFTS = tuple(range(156, 492))
EXPECTED_SELECTED_FULL_FAILURES = 269
EXPECTED_SELECTED_LITERAL_FAILURES = 187
EXPECTED_FAILURE_BODY_COUNTS = (((1, 2, 3, 4, 6, 12), 269),)
EXPECTED_WEAKEST_BEST_FULL = Q(31_091_887, 85_183_035)
EXPECTED_WEAKEST_BEST_LITERAL = Q(-8_062_511, 21_715_855)
EXPECTED_MAXIMUM_BAD_CELLS = 46
EXPECTED_SEMANTIC = "85622d480d328d4f3aa84efbe99ae84ec43d7acfe033f8512e50c7f655ac77c6"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def file_sha(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")).hexdigest()


def load(name, path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


require(file_sha(PRIMARY) == EXPECTED_PRIMARY, file_sha(PRIMARY))
require(file_sha(MEDIAN) == EXPECTED_MEDIAN, file_sha(MEDIAN))
B = load("coherent_shift_multicell_primary", PRIMARY)
M = load("coherent_shift_multicell_median", MEDIAN)
R = M.R


def mass(rows):
    return sum((right - left for left, right in R.merge_intervals(tuple(rows))), Q())


def prim_credit(vertices, weights):
    selected = {vertices[0]}
    credit = Q()
    tree = []
    while len(selected) < len(vertices):
        weight, a, b = max(
            (weights[a, b], a, b)
            for a in selected for b in vertices if b not in selected
        )
        selected.add(b)
        credit += weight
        tree.append((a, b, weight))
    return credit, tuple(tree)


def discover_shift(shift):
    result = B.scan_shift(shift)
    packets = []
    for failure in result[3]:
        row, full_credit, full_margin, full_tree, union, union_margin, all_edges = failure
        if full_margin <= 0:
            packets.append((
                shift, row[1], row[2], row[3], row[5], row[6],
                row[0], full_margin, union_margin,
            ))
    return tuple(packets)


def audit_packet(packet):
    (
        shift, body, L, selected_cell, left, right,
        selected_cross_margin, selected_full_margin, selected_union_margin,
    ) = packet
    ruler, ranges = R.safe_cell_ranges(body)
    require(ruler == L, (body, ruler, L))
    cells = tuple(cell for first, last in ranges for cell in range(first, last))
    require(selected_cell == cells[len(cells) // 2], (body, selected_cell, cells))
    vertices = tuple((3, residue) for residue in left) + tuple(
        (5, residue) for residue in right
    )
    rows = []
    for cell in cells:
        events = {
            vertex: tuple(
                R.direct_multiplier_arcs(L, vertex[0] * L - vertex[1], cell)
            )
            for vertex in vertices
        }
        singles = {vertex: mass(event) for vertex, event in events.items()}
        weights = {}
        for a, b in combinations(vertices, 2):
            overlap = singles[a] + singles[b] - mass(events[a] + events[b])
            require(overlap >= 0, (packet, cell, a, b, overlap))
            weights[a, b] = weights[b, a] = overlap
        debt = sum(singles.values(), Q()) - Q(6, 7)
        full_credit, full_tree = prim_credit(vertices, weights)
        full_margin = full_credit - debt
        union_margin = mass(
            tuple(interval for event in events.values() for interval in event)
        ) - Q(6, 7)
        rows.append((cell, full_margin, union_margin, full_tree))

    selected = next(row for row in rows if row[0] == selected_cell)
    require(selected[1] == selected_full_margin, (packet, selected))
    require(selected[2] == selected_union_margin, (packet, selected))
    alternatives = tuple(row for row in rows if row[0] != selected_cell)
    best_full = max(alternatives, key=lambda row: (row[1], -row[0]))
    best_literal = min(alternatives, key=lambda row: (row[2], row[0]))
    full_bad_cells = tuple(row[0] for row in rows if row[1] <= 0)
    union_bad_cells = tuple(row[0] for row in rows if row[2] >= 0)
    return (
        packet,
        len(cells),
        full_bad_cells,
        union_bad_cells,
        best_full,
        best_literal,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    context = get_context("fork")
    if args.processes == 1:
        by_shift = tuple(map(discover_shift, SHIFTS))
    else:
        with ProcessPoolExecutor(max_workers=args.processes, mp_context=context) as pool:
            by_shift = tuple(pool.map(discover_shift, SHIFTS))
    packets = tuple(packet for group in by_shift for packet in group)
    require(len(packets) == EXPECTED_SELECTED_FULL_FAILURES, len(packets))
    if args.processes == 1:
        audits = tuple(map(audit_packet, packets))
    else:
        with ProcessPoolExecutor(max_workers=args.processes, mp_context=context) as pool:
            audits = tuple(pool.map(audit_packet, packets))

    full_unrescued = tuple(row for row in audits if row[4][1] <= 0)
    literal_unrescued = tuple(row for row in audits if row[5][2] >= 0)
    require(not full_unrescued, full_unrescued)
    require(not literal_unrescued, literal_unrescued)
    selected_literal_failures = tuple(
        row for row in audits if row[0][8] >= 0
    )
    weakest_best_full = min(
        (row[4][1], row[0], row[4][0], row[1]) for row in audits
    )
    weakest_best_literal = max(
        (row[5][2], row[0], row[5][0], row[1]) for row in audits
    )
    maximum_full_bad = max(
        (len(row[2]), row[0], row[2], row[1]) for row in audits
    )
    maximum_union_bad = max(
        (len(row[3]), row[0], row[3], row[1]) for row in audits
    )
    ruler, ranges = R.safe_cell_ranges((1, 2, 3, 4, 6, 12))
    require(ruler == 168, ruler)
    all_cells = tuple(cell for first, last in ranges for cell in range(first, last))
    require(len(all_cells) == 88, len(all_cells))
    full_bad_by_cell = Counter(cell for row in audits for cell in row[2])
    union_bad_by_cell = Counter(cell for row in audits for cell in row[3])
    universal_full_rescue_cells = tuple(
        cell for cell in all_cells if full_bad_by_cell[cell] == 0
    )
    universal_literal_rescue_cells = tuple(
        cell for cell in all_cells if union_bad_by_cell[cell] == 0
    )
    digest = sha256()
    for row in audits:
        digest.update((repr(row) + "\n").encode())
    semantic = digest.hexdigest()
    failure_body_counts = tuple(sorted(Counter(row[0][1] for row in audits).items()))
    require(
        len(selected_literal_failures) == EXPECTED_SELECTED_LITERAL_FAILURES,
        len(selected_literal_failures),
    )
    require(failure_body_counts == EXPECTED_FAILURE_BODY_COUNTS, failure_body_counts)
    require(weakest_best_full[0] == EXPECTED_WEAKEST_BEST_FULL, weakest_best_full)
    require(
        weakest_best_literal[0] == EXPECTED_WEAKEST_BEST_LITERAL,
        weakest_best_literal,
    )
    require(maximum_full_bad[0] == EXPECTED_MAXIMUM_BAD_CELLS, maximum_full_bad)
    require(maximum_union_bad[0] == EXPECTED_MAXIMUM_BAD_CELLS, maximum_union_bad)
    require(semantic == EXPECTED_SEMANTIC, semantic)

    print("LRC14 COHERENT SHIFT COMMON-POSITIVITY MULTI-CELL RESCUE AUDIT")
    print("primary_dependency", PRIMARY.relative_to(ROOT), EXPECTED_PRIMARY)
    print("event_dependency", MEDIAN.relative_to(ROOT), EXPECTED_MEDIAN)
    print("shifts", (SHIFTS[0], SHIFTS[-1]), "selected_full_failures", len(packets))
    print("failure_shift_counts", tuple(sorted(Counter(row[0][0] for row in audits).items())))
    print("failure_body_counts", failure_body_counts)
    print("selected_literal_failures", len(selected_literal_failures))
    print("alternative_full_unrescued", len(full_unrescued))
    print("alternative_literal_unrescued", len(literal_unrescued))
    print("weakest_best_alternative_full", weakest_best_full)
    print("weakest_best_alternative_literal", weakest_best_literal)
    print("maximum_full_bad_cells", maximum_full_bad)
    print("maximum_union_bad_cells", maximum_union_bad)
    print("universal_full_rescue_cells", universal_full_rescue_cells)
    print("universal_literal_rescue_cells", universal_literal_rescue_cells)
    print("maximum_packet_failures_at_one_cell", max(full_bad_by_cell.values()),
          max(union_bad_by_cell.values()))
    print("c392_selected_literal_audits", tuple(
        row for row in selected_literal_failures if row[0][0] == 392
    ))
    print("semantic_sha256", semantic)
    print("status=FINITE-EXACT scratch audit; one common alternative cell per packet")


if __name__ == "__main__":
    main()
