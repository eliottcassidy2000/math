#!/usr/bin/env python3
"""Seven-private-point Moser extension: exact metric and full-boundary controls."""
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import hashlib
import importlib.util
import json
import subprocess
import sys


def need(ok, detail):
    if not ok:
        raise RuntimeError(detail)


def partitions(n):
    P = [(0,)]
    for _ in range(1, n):
        P = [a + (j,) for a in P for j in range(max(a) + 2)]
    return P


def falling(q, k):
    return factorial(q) // factorial(q - k) if k <= q else 0


def proper(words, edges):
    return all(words[u] != words[v] for u, v in edges)


def edge_mask(edges):
    pairs = list(combinations(range(7), 2))
    return sum(1 << pairs.index(tuple(sorted(e))) for e in edges)


def iso_map(source_edges, target_mask):
    for permutation in permutations(range(7)):
        if edge_mask([(permutation[u], permutation[v]) for u, v in source_edges]) == target_mask:
            return permutation
    return None


def main():
    here = Path(__file__).resolve().parent
    root = here.parents[3]
    source = Path(__file__).read_bytes()
    pins = json.loads((here / "SOURCE_PINS.json").read_text(encoding="utf-8"))
    frozen_sources = {}
    for pin in pins:
        b = (root / pin["path"]).read_bytes()
        need(hashlib.sha256(b).hexdigest() == pin["sha256"], ("source pin", pin["path"]))
        frozen_sources[pin["path"]] = b
    expected_census = json.loads((here / "census_receipt.json").read_text(encoding="utf-8"))
    census_command = [sys.executable] + (["-O"] if sys.flags.optimize else []) + [str(here / "census.py")]
    actual_census = json.loads(subprocess.check_output(census_command, cwd=here))
    need(actual_census == expected_census, "complete labelled seven-core census replay")
    need([r["mask"] for r in actual_census["representatives"]] == [120383, 120444, 120568, 122616, 378106, 380278],
         "six canonical residual graphs")

    old_path = root / "05-knowledge/results/two_port_beyond_nine_20261103/geometry/verify.py"
    spec = importlib.util.spec_from_file_location("seven_core_old_geometry", old_path)
    old = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = old
    spec.loader.exec_module(old)
    inherited = json.loads(old_path.with_name("coordinates.json").read_text(encoding="utf-8"))
    P = [tuple(old.K(tuple(old.F(s) for s in a)) for a in point) for point in inherited["spindle"]]
    E = old.graph(P)
    need(len(E) == 11 and len(P) == 7, "actual Moser carrier")
    mapping = iso_map(E, 378106)
    need(mapping is not None, "Moser canonical mask")
    need(all(any(tuple(sorted(e)) in E for e in combinations(T, 2)) for T in combinations(range(7), 3)),
         "Moser independence number at most two")
    need(any(tuple(sorted(e)) not in E for e in combinations(range(7), 2)), "Moser independence number exactly two")

    # Both exact formulations of the exterior triple-contact certificate.
    unit_triples = 0
    for T in combinations(range(7), 3):
        center = old.circle_center(*(P[t] for t in T))
        d, e = old.minus(P[T[1]], P[T[0]]), old.minus(P[T[2]], P[T[0]])
        determinant = d[0] * e[1] - d[1] * e[0]
        radius1 = determinant != old.ZERO and old.norm(d) * old.norm(e) * old.distance(P[T[1]], P[T[2]]) == 4 * determinant * determinant
        need(radius1 == (center is not None), ("circumradius routes", T))
        need(center is None or center in P, ("no external triple contact", T))
        unit_triples += int(radius1)
    need(unit_triples == 10, "ten old unit-triple centres")

    # The two non-three-colourable false metric survivors.
    pairs = list(combinations(range(7), 2))
    E380 = [e for j, e in enumerate(pairs) if 380278 >> j & 1 and e != (2, 3)]
    need(len(E380) == 11 and iso_map(E380, 378106) is not None, "380278 has spanning Moser after one deletion")
    alpha = (old.scalar(old.F(1, 2)), old.R3 / 2)
    central = [(old.ZERO, old.ZERO), (old.ONE, old.ZERO), alpha]
    caps = [tuple(a + b - c for a, b, c in zip(central[i], central[j], central[k]))
            for i, j, k in [(0, 1, 2), (0, 2, 1), (1, 2, 0)]]
    need(all(old.distance(caps[i], caps[j]) == old.scalar(4) for i, j in combinations(range(3), 2)),
         "122616 forces side-two cap triangle")
    need(old.circle_center(*caps) is None, "side-two triangle has no common unit centre")

    # Rigidity ingredients: diamond far endpoints and closing-angle data.
    need(old.distance(P[0], P[1]) == old.distance(P[0], P[4]) == old.scalar(3), "diamond far lengths")
    d, e = old.minus(P[1], P[0]), old.minus(P[4], P[0])
    dot = d[0] * e[0] + d[1] * e[1]
    need(dot == old.scalar(old.F(5, 2)) and old.distance(P[1], P[4]) == old.ONE, "rigid hinge cosine5/6")

    # Complete boundary patterns, including all states rather than a selected pin.
    patterns = partitions(7)
    valid = [beta for beta in patterns if max(beta) < 5 and proper(beta, E)]
    need(len(valid) == 42 and sum(len(set(b)) == 4 for b in valid) == 16
         and sum(len(set(b)) == 5 for b in valid) == 26, "full five-colour boundary partition spectrum")
    capacity_table = []
    for beta in valid:
        for a, b in combinations(range(5), 2):
            forbidden_vertices = [v for v, c in enumerate(beta) if c in (a, b)]
            need(len(forbidden_vertices) <= 4 and 2 * len(forbidden_vertices) < 14,
                 ("reciprocal contact capacity", beta, a, b))
            capacity_table.append(dict(partition=list(beta), forbidden_colours=[a, b],
                                       vertices=forbidden_vertices, capacity=2 * len(forbidden_vertices), required=14))

    # Connected exact14-point control with a genuine private Moser3-core.
    # Unit translation also creates every corresponding contact plus one extra.
    edge = old.minus(P[1], P[4])
    shift = ((edge[0] + old.R3 * edge[1]) / 2,
             (edge[1] - old.R3 * edge[0]) / 2)
    need(old.norm(shift) == old.ONE, "unit translation")
    Q = [tuple(a + b for a, b in zip(point, shift)) for point in P]
    all_points = P + Q
    full_edges = old.graph(all_points)
    cross = {(u, v - 7) for u, v in full_edges if u < 7 <= v}
    need(cross == {(j, j) for j in range(7)} | {(1, 4)}, "all eight cross contacts")
    need(len(full_edges) == 30, "complete fourteen-point incidence")
    need({(u - 7, v - 7) for u, v in full_edges if u >= 7} == E, "actual private Moser core")
    partial_tables = []
    whole_counts = {}
    witnesses = {}
    for q in (4, 5):
        private_colourings = [c for c in product(range(q), repeat=7) if proper(c, E)]
        need(len(private_colourings) == (384 if q == 4 else 5040), ("private literal counts", q))
        total = 0
        nonzero = 0
        for beta in valid:
            if max(beta) >= q:
                continue
            extensions = [c for c in private_colourings if all(beta[u] != c[v] for u, v in cross)]
            w = len(extensions)
            need(w > 0 if q == 5 else True, ("five-colour connected control", beta))
            nonzero += int(w > 0)
            total += falling(q, len(set(beta))) * w
            if extensions:
                need(proper(beta + extensions[0], full_edges), ("literal complete colouring", beta))
                witnesses.setdefault(q, list(beta + extensions[0]))
            partial_tables.append(dict(palette=q, partition=list(beta), weight=w,
                                        witness=None if not extensions else list(extensions[0])))
        whole_counts[q] = dict(colourings=total, supported_boundary_partitions=nonzero)
    need(whole_counts[4]["colourings"] > 0 and whole_counts[5]["supported_boundary_partitions"] == 42,
         "connected control remains4-chromatic and retains all42five-colour states")
    coordinates = dict(field="Q(sqrt3,sqrt11), basis1,sqrt3,sqrt11,sqrt33",
                       points=[[z.serial() for z in p] for p in all_points],
                       unit_edges=sorted(map(list, full_edges)),
                       cross_contacts=sorted(map(list, cross)),
                       translation=[z.serial() for z in shift],
                       colouring4=witnesses[4], colouring5=witnesses[5])
    outputs = {"coordinates.json": coordinates, "boundary_tables.json": partial_tables,
               "capacity_table.json": capacity_table}
    for name, data in outputs.items():
        if "--write" in sys.argv:
            (here / name).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
        else:
            need(json.loads((here / name).read_text(encoding="utf-8")) == data, ("frozen artifact", name))
    receipt = dict(seven_graph_universe=actual_census["universe"],
                   excluded_metric_masks=[122616, 380278],
                   three_choosable_masks=[120383, 120444, 120568],
                   residual_uniform_badlist_mask=378106,
                   standard_to_canonical_moser=list(mapping),
                   spindle_independence_number=2, checked_spindle_triples=35,
                   unit_triples_with_old_centres=unit_triples,
                   proper_fivecolour_boundary_patterns=42, fourblock_patterns=16, fiveblock_patterns=26,
                   reciprocal_capacity_cases=len(capacity_table), maximum_capacity=8, required_contacts=14,
                   connected_control_points=14, connected_control_unit_edges=30,
                   connected_control_cross_edges=8, connected_control_counts=whole_counts,
                   source_pins=len(pins), source_sha256=hashlib.sha256(source).hexdigest())
    need(Path(__file__).read_bytes() == source, "source changed")
    need(all((root / p).read_bytes() == b for p, b in frozen_sources.items()), "inherited source changed")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
