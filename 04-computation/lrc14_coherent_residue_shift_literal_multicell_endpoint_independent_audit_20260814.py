#!/usr/bin/env python3
"""Standalone endpoint-sweep audit of all selected-cell literal failures.

This path independently rebuilds the minimum body's 88 body-safe cells and
every danger interval from the defining inequality.  It scans all twenty
assignments at c=390..491, selects the 187 packets whose upper-median-cell
literal union is at least 6/7, and tests every safe cell by an endpoint-slab
union sweep.  It imports no interval, overlap, or spanning-tree engine.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "04-computation/lrc14_coherent_residue_shift_boundary_census_20260814.py"
EXPECTED_PRIMARY = "5c69fbd209c332d76cc51814bc1ef126d8cc3dea0e9c9d16f00ccf095d6b9a7e"
BODY = (1, 2, 3, 4, 6, 12)
L = 168
SELECTED_CELL = 90
SHIFTS = tuple(range(390, 492))
ASSIGNMENTS = tuple(combinations(range(6), 3))
EXPECTED_FAILURES = 187
EXPECTED_MAXIMUM_BAD_CELLS = 46
EXPECTED_WEAKEST_BEST_LITERAL = Q(-8_062_511, 21_715_855)
EXPECTED_SEMANTIC = "db0d29f421c2533d4303d37418f5a9b853353b2cd351d3dd0b67e97459d53881"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def file_sha(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")).hexdigest()


require(file_sha(PRIMARY) == EXPECTED_PRIMARY, file_sha(PRIMARY))


def clip(left, right):
    left = max(Q(), left)
    right = min(Q(1), right)
    return (left, right) if left < right else None


def direct_arcs(multiplier, cell):
    """Solve ||multiplier*(cell+u)/L||<1/14 directly for 0<=u<=1."""
    require(multiplier > 0, (multiplier, cell))
    first = (multiplier * cell) // L - 1
    last = (multiplier * (cell + 1)) // L + 1
    rows = []
    for integer in range(first, last + 1):
        left = Q(L * (14 * integer - 1), 14 * multiplier) - cell
        right = Q(L * (14 * integer + 1), 14 * multiplier) - cell
        row = clip(left, right)
        if row is not None:
            rows.append(row)
    return tuple(sorted(set(rows)))


def union_mass_slab(arcs):
    """Literal union mass from independently classified endpoint slabs."""
    endpoints = sorted({Q(), Q(1), *(point for arc in arcs for point in arc)})
    answer = Q()
    for left, right in zip(endpoints, endpoints[1:]):
        midpoint = (left + right) / 2
        if any(a < midpoint < b for a, b in arcs):
            answer += right - left
    return answer


def union_margin(left, right, cell):
    arcs = []
    for residue in left:
        arcs.extend(direct_arcs(3 * L - residue, cell))
    for residue in right:
        arcs.extend(direct_arcs(5 * L - residue, cell))
    return union_mass_slab(tuple(arcs)) - Q(6, 7)


def safe_cells():
    cells = tuple(
        cell for cell in range(L)
        if all(not direct_arcs(label, cell) for label in BODY)
    )
    require(len(cells) == 88, len(cells))
    require(cells[len(cells) // 2] == SELECTED_CELL, cells)
    return cells


def main():
    cells = safe_cells()
    failures = []
    for shift in SHIFTS:
        residues = tuple(shift + label for label in BODY)
        require(max(residues) < 3 * L, (shift, residues))
        for indices in ASSIGNMENTS:
            left = tuple(residues[index] for index in indices)
            right = tuple(
                residues[index] for index in range(6) if index not in indices
            )
            selected_margin = union_margin(left, right, SELECTED_CELL)
            if selected_margin >= 0:
                failures.append((shift, left, right, selected_margin))
    require(len(failures) == EXPECTED_FAILURES, len(failures))

    audits = []
    for packet in failures:
        shift, left, right, selected_margin = packet
        rows = tuple(
            (cell, union_margin(left, right, cell)) for cell in cells
        )
        require(
            next(margin for cell, margin in rows if cell == SELECTED_CELL)
            == selected_margin,
            (packet, rows),
        )
        bad_cells = tuple(cell for cell, margin in rows if margin >= 0)
        alternatives = tuple(row for row in rows if row[0] != SELECTED_CELL)
        best = min(alternatives, key=lambda row: (row[1], row[0]))
        require(best[1] < 0, (packet, best, bad_cells))
        audits.append((packet, bad_cells, best))
    audits = tuple(audits)

    maximum_bad = max((len(row[1]), row) for row in audits)
    weakest_best = max((row[2][1], row[0], row[2][0]) for row in audits)
    require(maximum_bad[0] == EXPECTED_MAXIMUM_BAD_CELLS, maximum_bad)
    require(weakest_best[0] == EXPECTED_WEAKEST_BEST_LITERAL, weakest_best)
    digest = sha256()
    for row in audits:
        digest.update((repr(row) + "\n").encode())
    semantic = digest.hexdigest()
    require(semantic == EXPECTED_SEMANTIC, semantic)

    print("LRC14 COHERENT SHIFT LITERAL MULTI-CELL ENDPOINT AUDIT")
    print("pinned_universe", PRIMARY.relative_to(ROOT), EXPECTED_PRIMARY)
    print("body", BODY, "L", L, "safe_cells", len(cells), "selected_cell", SELECTED_CELL)
    print("shifts", (SHIFTS[0], SHIFTS[-1]), "assignments", len(ASSIGNMENTS))
    print("selected_literal_failures", len(failures))
    print("failure_shift_counts", tuple(sorted(Counter(row[0][0] for row in audits).items())))
    print("alternative_unrescued", 0)
    print("weakest_best_alternative_literal", weakest_best)
    print("maximum_bad_cells", maximum_bad)
    print("c392", tuple(row for row in audits if row[0][0] == 392))
    print("semantic_sha256", semantic)
    print("status=INDEPENDENT FINITE-EXACT; defining arcs plus endpoint-slab unions")


if __name__ == "__main__":
    main()
