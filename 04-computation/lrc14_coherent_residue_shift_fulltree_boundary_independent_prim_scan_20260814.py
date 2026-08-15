#!/usr/bin/env python3
"""Independent exact census for the first full-tree/coherent-shift boundary.

This path scans c=85..392.  It uses Prim and computes every overlap through
|A intersect B|=|A|+|B|-|A union B|.  The earlier independent census covers
c=1..84, so together the two paths test the entire prefix through the first
literal-union obstruction.
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
SHIFTS = tuple(range(85, 393))
ASSIGNMENTS = tuple(combinations(range(6), 3))
EXPECTED_RESTRICTED_FAILURES = 719
EXPECTED_SAFE_PREFIX_USED_CERTIFICATE = (
    Q(227605554049, 3492993217078365),
    315,
    (1, 2, 3, 4, 6, 8),
    336,
    180,
    (317, 318, 319),
    (316, 321, 323),
)
EXPECTED_FULL_FAILURE_SHIFTS = (390, 390, 392, 392, 392, 392)
EXPECTED_UNION_FAILURE_SHIFTS = (392, 392)
EXPECTED_SEMANTIC = "b14ec22a230fc34f78c9abf639cd29932eea12f7ca2972ee00cc7c73c6b4eb29"


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
M = load("coherent_shift_fulltree_independent_median", MEDIAN)
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
            if (a, b) in weights
        )
        credit += weight
        tree.append((a, b, weight))
        selected.add(b)
    return credit, tuple(tree)


def scan_shift(shift):
    tested = 0
    restricted_failures = 0
    minimum_safe = None
    full_failures = []
    union_failures = []
    digest = sha256()
    bodies = M.body_universe()
    require(len(bodies) == 649, len(bodies))
    for body, L in bodies:
        ruler, ranges = R.safe_cell_ranges(body)
        require(ruler == L, (body, ruler, L))
        cells = tuple(j for left, right in ranges for j in range(left, right))
        cell = cells[len(cells) // 2]
        residues = tuple(e + shift for e in body)
        require(max(residues) < 3 * L, (shift, body, L, residues))
        events = {
            (q, a): tuple(R.direct_multiplier_arcs(L, q * L - a, cell))
            for q in (3, 5) for a in residues
        }
        singles = {vertex: mass(event) for vertex, event in events.items()}

        def weight(a, b):
            return singles[a] + singles[b] - mass(events[a] + events[b])

        cross_weights = {}
        for a in residues:
            for b in residues:
                x, y = (3, a), (5, b)
                w = weight(x, y)
                cross_weights[x, y] = cross_weights[y, x] = w
        for indices in ASSIGNMENTS:
            left = tuple(residues[index] for index in indices)
            right = tuple(residues[index] for index in range(6) if index not in indices)
            vertices = tuple((3, a) for a in left) + tuple((5, b) for b in right)
            cross_credit, cross_tree = prim_credit(vertices, cross_weights)
            debt = sum((singles[vertex] for vertex in vertices), Q()) - Q(6, 7)
            cross_margin = cross_credit - debt
            full_credit = cross_credit
            full_tree = cross_tree
            union = None
            union_margin = None
            if cross_margin <= 0:
                restricted_failures += 1
                full_weights = dict(cross_weights)
                for a, b in combinations(vertices, 2):
                    if (a, b) not in full_weights:
                        w = weight(a, b)
                        full_weights[a, b] = full_weights[b, a] = w
                full_credit, full_tree = prim_credit(vertices, full_weights)
                union = mass(tuple(interval for vertex in vertices for interval in events[vertex]))
                union_margin = union - Q(6, 7)
            full_margin = full_credit - debt
            row = (
                shift, body, L, cell, left, right, cross_margin, full_margin,
                union, union_margin, cross_tree, full_tree,
            )
            digest.update((repr(row) + "\n").encode())
            tested += 1
            if full_margin <= 0:
                full_failures.append(row)
            elif minimum_safe is None or full_margin < minimum_safe[0]:
                minimum_safe = (full_margin, row)
            if union_margin is not None and union_margin >= 0:
                union_failures.append(row)
    require(tested == 12_980, (shift, tested))
    return (
        shift, tested, restricted_failures, minimum_safe,
        tuple(full_failures), tuple(union_failures), digest.hexdigest(),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    args = parser.parse_args()
    if args.processes == 1:
        results = tuple(map(scan_shift, SHIFTS))
    else:
        with ProcessPoolExecutor(
            max_workers=args.processes, mp_context=get_context("fork")
        ) as pool:
            results = tuple(pool.map(scan_shift, SHIFTS))
    require(tuple(row[0] for row in results) == SHIFTS, "shift order")
    full_failures = tuple(row for result in results for row in result[4])
    union_failures = tuple(row for result in results for row in result[5])
    safe_prefix_used_certificate_minimum = min(
        result[3] for result in results if result[0] < 390
    )
    digest = sha256()
    for result in results:
        digest.update((repr(result) + "\n").encode())
    restricted_failures = sum(result[2] for result in results)
    safe_control = (
        safe_prefix_used_certificate_minimum[0],
        safe_prefix_used_certificate_minimum[1][0],
        safe_prefix_used_certificate_minimum[1][1],
        safe_prefix_used_certificate_minimum[1][2],
        safe_prefix_used_certificate_minimum[1][3],
        safe_prefix_used_certificate_minimum[1][4],
        safe_prefix_used_certificate_minimum[1][5],
    )
    require(restricted_failures == EXPECTED_RESTRICTED_FAILURES, restricted_failures)
    require(
        safe_control == EXPECTED_SAFE_PREFIX_USED_CERTIFICATE,
        safe_control,
    )
    require(
        tuple(row[0] for row in full_failures) == EXPECTED_FULL_FAILURE_SHIFTS,
        tuple(row[0] for row in full_failures),
    )
    require(
        tuple(row[0] for row in union_failures) == EXPECTED_UNION_FAILURE_SHIFTS,
        tuple(row[0] for row in union_failures),
    )
    require(digest.hexdigest() == EXPECTED_SEMANTIC, digest.hexdigest())
    print("LRC14 COHERENT SHIFT FULL-TREE BOUNDARY INDEPENDENT CENSUS")
    print("dependency", MEDIAN.relative_to(ROOT), EXPECTED_MEDIAN)
    print("shifts", (SHIFTS[0], SHIFTS[-1]), "packets", len(SHIFTS) * 12_980)
    print("restricted_failures", restricted_failures)
    # Before the first failure, the restricted cross tree already certifies
    # every row.  This is the minimum margin of the certificate actually used,
    # not a census of the (possibly larger) unrestricted K6 optimum.
    print(
        "safe_prefix_used_certificate_minimum",
        safe_prefix_used_certificate_minimum,
    )
    print("full_failure_shifts", tuple(row[0] for row in full_failures))
    for row in full_failures:
        print("full_failure", row)
    print("union_failure_shifts", tuple(row[0] for row in union_failures))
    for row in union_failures:
        print("union_failure", row)
    print("semantic_sha256", digest.hexdigest())
    print("status=independent exact census beyond c=84; Prim plus union-identity overlaps")


if __name__ == "__main__":
    main()
