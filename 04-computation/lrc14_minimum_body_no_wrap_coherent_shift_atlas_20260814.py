#!/usr/bin/env python3
"""Exact coherent-shift atlas on the unique minimum-ruler body.

Scope: E=(1,2,3,4,6,12), (L,j)=(168,90), every positive coherent
shift c for which all shifted residues stay below L (1 <= c <= 155),
and all 20 assignments of levels (3,3,3,5,5,5).
"""
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
CELL = 90
SHIFTS = tuple(range(1, L - max(BODY)))
ASSIGNMENTS = tuple(combinations(range(6), 3))
TRIPLE_A = (1, 3, 12)
TRIPLE_B = (2, 4, 6)
EXPECTED_FAILURE_SHIFTS = (
    67, 69, 80, 82, 84, 91, 93, 97, 104, 106, 110, 112, 115,
    117, 119, 123, 125, 128, 130, 136, 138, 141, 143, 149, 151, 154,
)
EXPECTED_SEMANTIC = "52a7bab56606a50ebbdc71209a4d27e54e4ea325969f9cc9cceda77566d36c57"


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
M = load("minimal_body_no_wrap_median", MEDIAN)
R = M.R


def mass(rows):
    return sum((right - left for left, right in R.merge_intervals(tuple(rows))), Q())


def overlap(first, second):
    i = j = 0
    total = Q()
    while i < len(first) and j < len(second):
        total += max(Q(), min(first[i][1], second[j][1])
                     - max(first[i][0], second[j][0]))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def maximum_tree(vertices, edges):
    parent = {vertex: vertex for vertex in vertices}

    def root(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    credit = Q()
    tree = []
    for weight, a, b in sorted(edges, reverse=True):
        x, y = root(a), root(b)
        if x == y:
            continue
        parent[x] = y
        credit += weight
        tree.append((a, b, weight))
        if len(tree) == len(vertices) - 1:
            return credit, tuple(tree)
    raise RuntimeError("disconnected edge graph")


def main():
    digest = sha256()
    packets = 0
    failures = []
    global_minimum = None
    for shift in SHIFTS:
        residues = tuple(e + shift for e in BODY)
        require(0 < min(residues) <= max(residues) < L, (shift, residues))
        events = {
            (q, a): R.direct_multiplier_arcs(L, q * L - a, CELL)
            for q in (3, 5) for a in residues
        }
        singles = {vertex: mass(event) for vertex, event in events.items()}
        for indices in ASSIGNMENTS:
            left = tuple(residues[index] for index in indices)
            right = tuple(residues[index] for index in range(6) if index not in indices)
            vertices = tuple((3, a) for a in left) + tuple((5, b) for b in right)
            debt = sum((singles[vertex] for vertex in vertices), Q()) - Q(6, 7)
            cross_edges = tuple(
                (overlap(events[3, a], events[5, b]), (3, a), (5, b))
                for a in left for b in right
            )
            cross_credit, cross_tree = maximum_tree(vertices, cross_edges)
            cross_margin = cross_credit - debt
            row = (
                cross_margin, shift, left, right, debt,
                cross_credit, cross_tree,
            )
            packets += 1
            digest.update((repr(row) + "\n").encode())
            if global_minimum is None or row < global_minimum:
                global_minimum = row
            if cross_margin <= 0:
                base_left = tuple(a - shift for a in left)
                base_right = tuple(a - shift for a in right)
                require((base_left, base_right) in (
                    (TRIPLE_A, TRIPLE_B), (TRIPLE_B, TRIPLE_A),
                ), (shift, base_left, base_right))
                full_edges = tuple(
                    (overlap(events[a], events[b]), a, b)
                    for a, b in combinations(vertices, 2)
                )
                full_credit, full_tree = maximum_tree(vertices, full_edges)
                full_margin = full_credit - debt
                union = mass(tuple(
                    interval for vertex in vertices for interval in events[vertex]
                ))
                phase = (-shift * CELL) % L
                if 2 * phase > L:
                    phase -= L
                within = tuple(edge for edge in full_tree if edge[0][0] == edge[1][0])
                bridge = tuple(edge for edge in full_tree if edge[0][0] != edge[1][0])
                require(len(within) == 4 and len(bridge) == 1,
                        (shift, full_tree))
                require(full_margin > 0, (shift, full_margin))
                require(union < Q(6, 7), (shift, union))
                failures.append((
                    shift, phase, left, right, cross_margin,
                    full_margin, union - Q(6, 7), full_tree,
                ))
    require(packets == 3_100, packets)
    failure_shifts = tuple(row[0] for row in failures)
    require(len(failure_shifts) == len(set(failure_shifts)), failure_shifts)
    require(failure_shifts == EXPECTED_FAILURE_SHIFTS, failure_shifts)
    semantic = digest.hexdigest()
    require(semantic == EXPECTED_SEMANTIC, semantic)
    print("LRC14 MINIMUM-BODY POSITIVE NO-WRAP COHERENT-SHIFT ATLAS")
    print("body", BODY, "L", L, "cell", CELL,
          "shifts", (SHIFTS[0], SHIFTS[-1]), "packets", packets)
    print("restricted_failure_count", len(failures),
          "failure_shifts", failure_shifts)
    print("global_cross_minimum", global_minimum[:4])
    for row in failures:
        shift, phase, left, right, cross_margin, full_margin, union_margin, full_tree = row
        print("failure", (
            shift, phase, left, right, cross_margin, full_margin, union_margin,
        ))
    print("semantic_sha256", semantic)
    print("conclusion=every no-wrap packet closes by restricted cross-K3,3 or full-K6 Hunter; every restricted failure has literal union below 6/7")
    print("status=FINITE-EXACT atlas on one body; not an all-649 extension beyond c=84")


if __name__ == "__main__":
    main()
