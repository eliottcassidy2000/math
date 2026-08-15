#!/usr/bin/env python3
"""Independent exact audit of all c<=84 restricted-tree failures.

Unlike the census, this verifier enumerates every spanning tree (81 in
K3,3 and 1296 in K6) and computes literal union mass by an endpoint sweep.
"""
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MEDIAN = ROOT / "04-computation/lrc14_j7_reflected_median_star_chord_scout_20260801.py"
EXPECTED_MEDIAN = "b0eba7785cd73ecf76d2887f3d36316eae6980a76652058daa2587ecbe031276"
EXPECTED_SEMANTIC = "8f08fa06c0b1cebea6335cb51b49a2f05028579f8c4d58012907418af25ca8d9"
CASES = (
    (67, (68, 70, 79), (69, 71, 73),
     Q(-3_689_469_617_499, 1_523_807_086_440_275),
     Q(876_026_208, 3_635_859_955),
     Q(10_732_327_854, 18_179_299_775)),
    (69, (71, 73, 75), (70, 72, 81),
     Q(-109_470_469, 68_745_941_264),
     Q(79_844_187, 319_006_688),
     Q(52_553_773, 91_144_768)),
    (80, (81, 83, 92), (82, 84, 86),
     Q(-390_660_252_142, 18_345_869_001_369),
     Q(10_505_347_270, 43_370_848_703),
     Q(3_685_017_575, 6_195_835_529)),
    (82, (84, 86, 88), (83, 85, 94),
     Q(-47_427_700_563, 4_633_736_015_480),
     Q(52_954_831, 205_015_720),
     Q(14_777_818, 25_626_965)),
    (84, (86, 88, 90), (85, 87, 96),
     Q(-52_458_412_717, 61_674_995_021_640),
     Q(129_336_503, 523_645_738),
     Q(66_228_727, 112_209_801)),
)


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
M = load("coherent_shift_failure_independent_median", MEDIAN)
R = M.R


def sweep_mass(families):
    """Union mass by sorted endpoints and midpoint membership."""
    rows = tuple(row for family in families for row in family)
    endpoints = sorted({endpoint for row in rows for endpoint in row})
    total = Q()
    for left, right in zip(endpoints, endpoints[1:]):
        midpoint = (left + right) / 2
        if any(a < midpoint < b for a, b in rows):
            total += right - left
    return total


def pair_mass(first, second):
    """Intersection mass by a two-pointer path independent of sweep_mass."""
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


def brute_maximum_tree(vertices, edges):
    maximum = None
    tree_count = 0
    for chosen in combinations(edges, len(vertices) - 1):
        parent = {vertex: vertex for vertex in vertices}

        def root(vertex):
            while parent[vertex] != vertex:
                vertex = parent[vertex]
            return vertex

        acyclic = True
        for weight, a, b in chosen:
            x, y = root(a), root(b)
            if x == y:
                acyclic = False
                break
            parent[x] = y
        if not acyclic:
            continue
        tree_count += 1
        credit = sum((edge[0] for edge in chosen), Q())
        if maximum is None or credit > maximum:
            maximum = credit
    require(maximum is not None, vertices)
    return maximum, tree_count


def main():
    body = (1, 2, 3, 4, 6, 12)
    L = 168
    cell = 90
    digest = sha256()
    rows = []
    for shift, left, right, expected_cross_margin, expected_full_margin, expected_union in CASES:
        require(tuple(sorted(left + right)) == tuple(e + shift for e in body),
                (shift, left, right))
        vertices = tuple((3, a) for a in left) + tuple((5, b) for b in right)
        events = {
            vertex: tuple(R.direct_multiplier_arcs(L, vertex[0] * L - vertex[1], cell))
            for vertex in vertices
        }
        singles = {vertex: sweep_mass((event,)) for vertex, event in events.items()}
        debt = sum(singles.values(), Q()) - Q(6, 7)
        cross_edges = tuple(
            (pair_mass(events[3, a], events[5, b]), (3, a), (5, b))
            for a in left for b in right
        )
        full_edges = tuple(
            (pair_mass(events[a], events[b]), a, b)
            for a, b in combinations(vertices, 2)
        )
        cross_credit, cross_tree_count = brute_maximum_tree(vertices, cross_edges)
        full_credit, full_tree_count = brute_maximum_tree(vertices, full_edges)
        cross_margin = cross_credit - debt
        full_margin = full_credit - debt
        union = sweep_mass(tuple(events.values()))
        require(cross_tree_count == 81, (shift, cross_tree_count))
        require(full_tree_count == 1_296, (shift, full_tree_count))
        require(cross_margin == expected_cross_margin, (shift, cross_margin))
        require(full_margin == expected_full_margin, (shift, full_margin))
        require(union == expected_union, (shift, union))
        require(union < Q(6, 7), (shift, union))
        row = (
            shift, cross_tree_count, full_tree_count, debt,
            cross_credit, cross_margin, full_credit, full_margin,
            union, union - Q(6, 7),
        )
        digest.update((repr(row) + "\n").encode())
        rows.append(row)
    print("LRC14 COHERENT-SHIFT FAILURE INDEPENDENT TREE/UNION AUDIT")
    print("body", body, "L", L, "cell", cell, "cases", len(rows))
    for row in rows:
        print("case", row)
    semantic = digest.hexdigest()
    require(semantic == EXPECTED_SEMANTIC, semantic)
    print("semantic_sha256", semantic)
    print("status=INDEPENDENT VERIFIED-EXACT failure audit; exhaustive tree enumeration and endpoint-sweep union")


if __name__ == "__main__":
    main()
