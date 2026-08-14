#!/usr/bin/env python3
"""Finite-exact all-649 probe of small uniform noncanonical residue shifts.

For every upper-median body, split its six residues into three level-3 and
three level-5 clauses. Then replace every canonical residue e by e+c for
c in {-2,-1,1,2}, retaining only positive residues. The audit restricts to
cross edges between the two equal-level classes, hence to K_3,3. This script
checks the literal cross-tree Hunter margin using actual singleton masses.
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
SHIFTS = (-2, -1, 1, 2)
EXPECTED_COUNTS = {-2: 1_700, -1: 5_180, 1: 12_980, 2: 12_980}
EXPECTED_SEMANTIC = "e165e13d5c2dd190ae3873e22b5775c9eead61d8939f5c1734057851f48d5c1f"
EXPECTED_MINIMA = {
    -2: Q(704_886_090_693_279_025_178, 10_402_065_960_840_068_379_865),
    -1: Q(179_516_485_647_379_506_696_761, 2_706_434_410_235_514_048_893_873),
    1: Q(2_423_380_281_636_460_862_201, 43_432_524_347_730_591_299_828),
    2: Q(64_566_844_993_819_835, 1_119_348_893_518_526_181),
}


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
M = load("reflected_shift_median", MEDIAN)
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


def cross_tree(left, right, weights):
    edges = sorted(
        ((weights[a, b], a, b) for a in left for b in right), reverse=True
    )
    parent = {vertex: vertex for vertex in left + right}

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
        if len(tree) == 5:
            return credit, tuple(tree)
    raise RuntimeError("K3,3 disconnected")


def main():
    bodies = M.body_universe()
    require(len(bodies) == 649, len(bodies))
    counts = {shift: 0 for shift in SHIFTS}
    minima = {}
    digest = sha256()
    failures = 0
    for body, L in bodies:
        ruler, ranges = R.safe_cell_ranges(body)
        require(ruler == L, (body, ruler, L))
        cells = tuple(j for left, right in ranges for j in range(left, right))
        cell = cells[len(cells) // 2]
        for shift in SHIFTS:
            residues = tuple(e + shift for e in body)
            if min(residues) <= 0:
                continue
            clauses = {
                (q, a): R.direct_multiplier_arcs(L, q * L - a, cell)
                for q in (3, 5)
                for a in residues
            }
            singles = {key: mass(value) for key, value in clauses.items()}
            weights = {
                (a, b): overlap(clauses[3, a], clauses[5, b])
                for a in residues
                for b in residues
            }
            for indices in combinations(range(6), 3):
                left = tuple(residues[index] for index in indices)
                right = tuple(
                    residues[index] for index in range(6) if index not in indices
                )
                credit, tree = cross_tree(left, right, weights)
                debt = sum((singles[3, a] for a in left), Q())
                debt += sum((singles[5, b] for b in right), Q()) - Q(6, 7)
                margin = credit - debt
                row = (margin, body, L, cell, shift, left, right, credit, debt, tree)
                counts[shift] += 1
                failures += margin <= 0
                if shift not in minima or row < minima[shift]:
                    minima[shift] = row
                digest.update((repr(row) + "\n").encode())
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, counts)
    require(failures == 0, failures)
    require({shift: row[0] for shift, row in minima.items()} == EXPECTED_MINIMA,
            minima)
    semantic = digest.hexdigest()
    require(semantic == EXPECTED_SEMANTIC, semantic)
    print("LRC14 UNIFORM NONCANONICAL RESIDUE-SHIFT PROBE")
    print("bodies", len(bodies), "counts", counts, "total", sum(counts.values()))
    print("failures", failures)
    for shift in SHIFTS:
        print("shift", shift, "minimum", minima[shift])
    print("semantic_sha256", semantic)
    print("conclusion=all tested shifted 3+3 repeated-level packets close by literal K3,3 Hunter tree")
    print("status=FINITE-EXACT positive probe; not an all-shift or arbitrary-residue theorem")


if __name__ == "__main__":
    main()
