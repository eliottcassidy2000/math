#!/usr/bin/env python3
"""Isolate every Hunter obstruction at the coherent residue shift c=84.

Universe: all 649 upper-median bodies; every 3+3 assignment of base levels
3 and 5; literal events with residues e+84.  Report cross-K3,3 and full-K6
maximum-tree margins, together with the exact six-event union mass.
"""
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path


HERE_ROOT = Path(__file__).resolve().parents[1]
ROOT = HERE_ROOT if (HERE_ROOT / "04-computation").is_dir() else Path.cwd()
MEDIAN = ROOT / "04-computation/lrc14_j7_reflected_median_star_chord_scout_20260801.py"
EXPECTED_MEDIAN = "b0eba7785cd73ecf76d2887f3d36316eae6980a76652058daa2587ecbe031276"
SHIFT = 84
EXPECTED_SEMANTIC = "4c2629fe164b49309639f2267540299749990e6df1e61a107396591f3bbdd173"


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
M = load("reflected_shift84_median", MEDIAN)
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
    raise RuntimeError("edge graph disconnected")


def main():
    bodies = M.body_universe()
    require(len(bodies) == 649, len(bodies))
    tested = 0
    cross_failures = []
    full_failures = []
    digest = sha256()
    for body, L in bodies:
        ruler, ranges = R.safe_cell_ranges(body)
        require(ruler == L, (body, ruler, L))
        cells = tuple(j for left, right in ranges for j in range(left, right))
        cell = cells[len(cells) // 2]
        residues = tuple(e + SHIFT for e in body)
        events = {
            (q, a): R.direct_multiplier_arcs(L, q * L - a, cell)
            for q in (3, 5)
            for a in residues
        }
        singles = {key: mass(value) for key, value in events.items()}
        for indices in combinations(range(6), 3):
            left = tuple(residues[index] for index in indices)
            right = tuple(residues[index] for index in range(6) if index not in indices)
            vertices = tuple((3, a) for a in left) + tuple((5, b) for b in right)
            cross_edges = tuple(
                (overlap(events[3, a], events[5, b]), (3, a), (5, b))
                for a in left for b in right
            )
            all_edges = tuple(
                (overlap(events[a], events[b]), a, b)
                for a, b in combinations(vertices, 2)
            )
            cross_credit, cross_tree = maximum_tree(vertices, cross_edges)
            full_credit, full_tree = maximum_tree(vertices, all_edges)
            singleton_sum = sum((singles[vertex] for vertex in vertices), Q())
            debt = singleton_sum - Q(6, 7)
            cross_margin = cross_credit - debt
            full_margin = full_credit - debt
            union = mass(tuple(interval for vertex in vertices for interval in events[vertex]))
            row = (
                body, L, cell, left, right, singleton_sum, debt,
                cross_credit, cross_margin, cross_tree,
                full_credit, full_margin, full_tree, union,
            )
            tested += 1
            digest.update((repr(row) + "\n").encode())
            if cross_margin <= 0:
                cross_failures.append(row)
            if full_margin <= 0:
                full_failures.append(row)
    require(tested == 12_980, tested)
    require(len(cross_failures) == 1, len(cross_failures))
    require(not full_failures, full_failures)
    semantic = digest.hexdigest()
    require(semantic == EXPECTED_SEMANTIC, semantic)
    witness = cross_failures[0]
    require(witness[:5] == ((1, 2, 3, 4, 6, 12), 168, 90,
                            (86, 88, 90), (85, 87, 96)), witness[:5])
    require(witness[8] == -Q(52_458_412_717, 61_674_995_021_640), witness[8])
    require(witness[11] == Q(129_336_503, 523_645_738), witness[11])
    require(witness[-1] - Q(6, 7) == -Q(209_657_717, 785_468_607), witness[-1])
    print("LRC14 COHERENT RESIDUE-SHIFT c=84 TREE OBSTRUCTION")
    print("bodies", len(bodies), "packets", tested)
    print("cross_failures", len(cross_failures), "full_tree_failures", len(full_failures))
    for row in cross_failures:
        print("cross_failure", row)
        print("union_minus_6/7", row[-1] - Q(6, 7))
    print("semantic_sha256", semantic)
    print("status=FINITE-EXACT; failure is restricted to the cross-K3,3 Hunter certificate")


if __name__ == "__main__":
    main()
