#!/usr/bin/env python3
"""Exact multi-cell repair for the two c=392 one-cell union obstructions."""
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MEDIAN = ROOT / "04-computation/lrc14_j7_reflected_median_star_chord_scout_20260801.py"
EXPECTED_MEDIAN = "b0eba7785cd73ecf76d2887f3d36316eae6980a76652058daa2587ecbe031276"
BODY = (1, 2, 3, 4, 6, 12)
L = 168
CASES = (
    ((393, 394, 398), (395, 396, 404)),
    ((395, 396, 404), (393, 394, 398)),
)
EXPECTED_RANGES = (
    (12, 13), (15, 26), (30, 39), (45, 52), (60, 69), (71, 78),
    (90, 97), (99, 108), (116, 123), (129, 138), (142, 153), (155, 156),
)
EXPECTED_FULL_FAILURE_CELLS = ((77, 90), (15, 60, 77, 90, 107, 152))
EXPECTED_UNION_FAILURE_CELLS = ((77, 90), (77, 90))
EXPECTED_SEMANTIC = "9cae9193ff51a4a1bf4581e2e5b4168f1ae577291eab1cabf854d0bd9669bc8c"


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


require(file_sha(MEDIAN) == EXPECTED_MEDIAN, file_sha(MEDIAN))
M = load("c392_multicell_median", MEDIAN)
R = M.R


def mass(rows):
    return sum((right - left for left, right in R.merge_intervals(tuple(rows))), Q())


def maximum_tree(vertices, weights):
    edges = sorted(
        ((weight, a, b) for (a, b), weight in weights.items() if a < b),
        reverse=True,
    )
    parent = {vertex: vertex for vertex in vertices}

    def root(vertex):
        while parent[vertex] != vertex:
            vertex = parent[vertex]
        return vertex

    credit = Q()
    for weight, a, b in edges:
        x, y = root(a), root(b)
        if x == y:
            continue
        parent[x] = y
        credit += weight
    return credit


def main():
    ruler, ranges = R.safe_cell_ranges(BODY)
    require(ruler == L, ruler)
    require(ranges == EXPECTED_RANGES, ranges)
    cells = tuple(cell for left, right in ranges for cell in range(left, right))
    require(cells[len(cells) // 2] == 90, cells[len(cells) // 2])
    digest = sha256()
    print("LRC14 c392 MULTI-CELL EXACT RESELECTION")
    print("body", BODY, "L", L, "safe_cells", len(cells), "ranges", ranges)
    for case_index, (left, right) in enumerate(CASES):
        vertices = tuple((3, residue) for residue in left) + tuple(
            (5, residue) for residue in right
        )
        rows = []
        for cell in cells:
            events = {
                vertex: tuple(R.direct_multiplier_arcs(L, vertex[0] * L - vertex[1], cell))
                for vertex in vertices
            }
            singles = {vertex: mass(event) for vertex, event in events.items()}
            weights = {}
            for a, b in combinations(vertices, 2):
                overlap = singles[a] + singles[b] - mass(events[a] + events[b])
                weights[a, b] = weights[b, a] = overlap
            debt = sum(singles.values(), Q()) - Q(6, 7)
            full_margin = maximum_tree(vertices, weights) - debt
            union_margin = mass(tuple(row for event in events.values() for row in event)) - Q(6, 7)
            row = (cell, full_margin, union_margin)
            rows.append(row)
            digest.update((repr((left, right, row)) + "\n").encode())
        union_failures = tuple(row for row in rows if row[2] >= 0)
        full_failures = tuple(row for row in rows if row[1] <= 0)
        require(
            tuple(row[0] for row in full_failures)
            == EXPECTED_FULL_FAILURE_CELLS[case_index],
            full_failures,
        )
        require(
            tuple(row[0] for row in union_failures)
            == EXPECTED_UNION_FAILURE_CELLS[case_index],
            union_failures,
        )
        print("case", left, right)
        print("full_failures", len(full_failures), full_failures)
        print("union_failures", len(union_failures), union_failures)
        print("best", max(rows, key=lambda row: row[1]))
        print("worst", min(rows, key=lambda row: row[1]))
    semantic = digest.hexdigest()
    require(semantic == EXPECTED_SEMANTIC, semantic)
    print("semantic_sha256", semantic)
    print("status=FINITE-EXACT cell-reselection repair of both c392 upper-median union failures")


if __name__ == "__main__":
    main()
