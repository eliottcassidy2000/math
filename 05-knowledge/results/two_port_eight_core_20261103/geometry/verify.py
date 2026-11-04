#!/usr/bin/env python3
"""Exact certificates for the eight-private-point Moser extension theorem.

The complete 2^28 graph enumeration is a pinned dependency in ../cores/.
This checker independently evaluates its 30 polynomial coefficients and
checks all eight geometric exclusions, plus a full actual cube control.
Stdlib only. Normal execution checks frozen artifacts; --write creates them.
"""
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import hashlib
import importlib.util
import json
import sys


def need(ok, detail):
    if not ok:
        raise RuntimeError(detail)


def edges(n, mask):
    return {pair for k, pair in enumerate(combinations(range(n), 2)) if mask >> k & 1}


def encode(n, E):
    return sum(1 << k for k, pair in enumerate(combinations(range(n), 2)) if pair in E)


def proper(colours, E):
    return all(colours[u] != colours[v] for u, v in E)


def partitions(n):
    out = [(0,)]
    for _ in range(1, n):
        out = [p + (j,) for p in out for j in range(max(p) + 2)]
    return out


def coefficient_by_endpoints(n, E, target):
    """Sum all 2^m signed endpoint choices, independently of truncated DP."""
    E = sorted(E)
    degrees = [0] * n
    for u, v in E:
        degrees[u] += 1
    answer = int(degrees == target)
    state = 0
    sign = 1
    for step in range(1, 1 << len(E)):
        bit = (step & -step).bit_length() - 1
        u, v = E[bit]
        if state >> bit & 1:
            degrees[u] += 1
            degrees[v] -= 1
        else:
            degrees[u] -= 1
            degrees[v] += 1
        state ^= 1 << bit
        sign = -sign
        if degrees == target:
            answer += sign
    return answer


def filters(n, E):
    A = [{v for v in range(n) if tuple(sorted((u, v))) in E} for u in range(n)]
    if any(len(A[u] & A[v]) > 2 for u, v in combinations(range(n), 2)):
        return False
    for v in range(n):
        unseen = set(A[v])
        while unseen:
            stack = [min(unseen)]
            component = set()
            while stack:
                u = stack.pop()
                if u in component:
                    continue
                component.add(u)
                stack.extend((A[u] & A[v]) - component)
            degrees = [len(A[u] & A[v]) for u in component]
            if max(degrees) > 2 or (min(degrees) == 2 and len(component) != 6):
                return False
            unseen -= component
    return True


def main():
    here = Path(__file__).resolve().parent
    root = here.parents[3]
    source = Path(__file__).read_bytes()
    pins = json.loads((here / "SOURCE_PINS.json").read_text(encoding="utf-8"))
    sources = {}
    for pin in pins:
        data = (root / pin["path"]).read_bytes()
        need(hashlib.sha256(data).hexdigest() == pin["sha256"], ("source pin", pin["path"]))
        sources[pin["path"]] = data
    census = json.loads((here.parent / "cores/order8.json").read_text(encoding="utf-8"))
    need((census["n"], census["all_graphs"], census["min_degree3_graphs"], census["metric_filter_graphs"], census["classes"])
         == (8, 1 << 28, 53470067, 724080, 38), "complete graph universe")
    for name, digest in census["source_sha256"].items():
        need(hashlib.sha256((here.parent / "cores" / name).read_bytes()).hexdigest() == digest,
             ("census producer digest", name))
    rows = census["representatives"]
    need(len(rows) == 38 and len({r["mask"] for r in rows}) == 38, "38 distinct representative masks")
    need(sum(r["orbit"] for r in rows) == 724080, "complete orbit cardinality")
    checked_coefficients = []
    unresolved = []
    for row in rows:
        E = edges(8, row["mask"])
        adjacency = [sum(1 << v for v in range(8) if tuple(sorted((u, v))) in E) for u in range(8)]
        degrees = sorted(a.bit_count() for a in adjacency)
        need(row["rows"] == adjacency and row["degrees"] == degrees and row["edges"] == len(E),
             ("graph decoding", row["mask"]))
        need(degrees[0] >= 3 and filters(8, E), ("necessary local conditions", row["mask"]))
        need(row["automorphisms"] * row["orbit"] == factorial(8), "orbit-stabilizer")
        cert = row["list_certificate"]
        if cert is None:
            unresolved.append(row["mask"])
            continue
        target = cert["exponents"]
        need(len(target) == 8 and all(0 <= d <= 2 for d in target) and sum(target) == len(E),
             ("degree-bounded homogeneous monomial", row["mask"]))
        coefficient = coefficient_by_endpoints(8, E, target)
        need(coefficient != 0 and coefficient == cert["coefficient"],
             ("independent signed endpoint coefficient", row["mask"], coefficient))
        checked_coefficients.append({"mask": row["mask"], "exponents": target,
                                     "coefficient": coefficient, "endpoint_terms": 1 << len(E)})
    need(len(checked_coefficients) == 30, "thirty independent three-list certificates")
    need(unresolved == [5675891, 5683827, 5691763, 5692274, 6215909, 6216436, 14587370, 48675698],
         "eight exact residual classes")

    old_path = root / "05-knowledge/results/two_port_beyond_nine_20261103/geometry/verify.py"
    spec = importlib.util.spec_from_file_location("eight_core_old_geometry", old_path)
    old = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = old
    spec.loader.exec_module(old)
    inherited = json.loads(old_path.with_name("coordinates.json").read_text(encoding="utf-8"))
    P = [tuple(old.K(tuple(old.F(s) for s in a)) for a in point) for point in inherited["spindle"]]
    old_edges = old.graph(P)
    need(len(old_edges) == 11 and len(P) == 7, "exact old Moser coordinates")
    unit_triples = []
    for T in combinations(range(7), 3):
        center = old.circle_center(*(P[t] for t in T))
        d, e = old.minus(P[T[1]], P[T[0]]), old.minus(P[T[2]], P[T[0]])
        det = d[0] * e[1] - d[1] * e[0]
        radius_one = det != old.ZERO and old.norm(d) * old.norm(e) * old.distance(P[T[1]], P[T[2]]) == 4 * det * det
        need(radius_one == (center is not None), ("two exact circumradius routes", T))
        need(center is None or center in P, ("no external triple contact", T))
        if center is not None:
            unit_triples.append({"triple": T, "old_center": P.index(center)})
    need(len(unit_triples) == 10, "all 35 triples tested, exactly ten old centres")

    # Every relabelled Moser embedding on seven labels, retaining a witness.
    moser_orbit = {}
    for p in permutations(range(7)):
        image = {tuple(sorted((p[u], p[v]))) for u, v in old_edges}
        moser_orbit.setdefault(encode(7, image), p)
    need(len(moser_orbit) == 630 and 378106 in moser_orbit, "full labelled Moser orbit")
    exclusions = []
    for mask in unresolved[2:]:
        E = edges(8, mask)
        found = None
        for omitted in range(8):
            keep = [u for u in range(8) if u != omitted]
            induced = encode(7, {(i, j) for i, j in combinations(range(7), 2) if (keep[i], keep[j]) in E})
            for moser, p in sorted(moser_orbit.items()):
                if moser & induced == moser:
                    embedding = [keep[p[i]] for i in range(7)]
                    need(len(set(embedding)) == 7 and omitted not in embedding, "injective Moser embedding")
                    need(all(tuple(sorted((embedding[u], embedding[v]))) in E for u, v in old_edges),
                         "all eleven required Moser edges")
                    neighbours = [v for v in range(8) if tuple(sorted((omitted, v))) in E]
                    need(len(neighbours) >= 3, "exterior vertex exceeds rigid Moser contact bound")
                    found = {"mask": mask, "reason": "rigid_Moser_exterior_capacity",
                             "moser_embedding_from_standard_coordinates": embedding,
                             "exterior_vertex": omitted, "exterior_neighbours": neighbours,
                             "maximum_possible_exterior_contacts": 2}
                    break
            if found is not None:
                break
        need(found is not None, ("Moser subgraph exclusion", mask))
        exclusions.append(found)

    # The remaining two masks contain the same six-vertex forbidden circle fan.
    hub, path, external = 0, [1, 7, 2, 5], 4
    fan_edges = {(min(hub, v), max(hub, v)) for v in path}
    fan_edges |= {tuple(sorted(e)) for e in zip(path, path[1:])}
    fan_edges |= {(1, 4), (4, 5)}
    alpha = (old.scalar(old.F(1, 2)), old.R3 / 2)
    circle = [(old.ONE, old.ZERO), alpha, (alpha[0] - 1, alpha[1]), (-old.ONE, old.ZERO)]
    need(len(set(circle)) == 4 and all(old.norm(p) == old.ONE for p in circle), "four distinct unit neighbours")
    need(all(old.distance(circle[i], circle[i + 1]) == old.ONE for i in range(3)), "consecutive hexagon path")
    need(old.distance(circle[0], circle[3]) == old.scalar(4), "antipodal fan endpoints at distance two")
    need(tuple(a + b for a, b in zip(circle[0], circle[3])) == (old.ZERO, old.ZERO),
         "the sole common unit centre is the original hub")
    for mask in unresolved[:2]:
        need(fan_edges <= edges(8, mask), ("required forbidden fan contacts", mask))
        exclusions.append({"mask": mask, "reason": "four_neighbour_circle_fan_forces_collision",
                           "hub": hub, "ordered_neighbour_path": path,
                           "antipodal_vertices": [1, 5], "forced_duplicate_vertex": external,
                           "required_edges": sorted(fan_edges)})
    exclusions.sort(key=lambda r: r["mask"])
    need([r["mask"] for r in exclusions] == unresolved, "all eight residuals excluded")

    # Actual connected positive control: old Moser plus an eight-point cube.
    generators = [(old.ONE, old.ZERO), (old.ZERO, old.ONE),
                  (old.scalar(old.F(3, 5)), old.scalar(old.F(4, 5)))]
    Q = [(sum((generators[i][0] for i in range(3) if mask >> i & 1), old.ZERO),
          old.ONE + sum((generators[i][1] for i in range(3) if mask >> i & 1), old.ZERO))
         for mask in range(8)]
    points = P + Q
    full_edges = old.graph(points)
    private_edges = {(u - 7, v - 7) for u, v in full_edges if u >= 7}
    cross_edges = sorted((u, v - 7) for u, v in full_edges if u < 7 <= v)
    cube_edges = {(u, v) for u, v in combinations(range(8), 2) if (u ^ v).bit_count() == 1}
    need(private_edges == cube_edges and len(private_edges) == 12, "exact planar unit cube with no extra private contacts")
    need(cross_edges == [(0, 0), (2, 0)] and len(full_edges) == 25, "all actual cross-contacts")
    need((0, 2) in old_edges, "every proper boundary removes two different colours at private vertex zero")
    all_patterns = [p for p in partitions(7) if proper(p, old_edges)]
    need(len(all_patterns) == 53, "complete old Moser equality partitions")
    old_counts = {q: sum(proper(c, old_edges) for c in product(range(q), repeat=7)) for q in (3, 4, 5)}
    need(old_counts == {3: 0, 4: 384, 5: 5040}, "literal old chromatic controls")
    boundary_tables = {}
    summaries = {}
    for q in (4, 5, 6):
        # Literal complete q^8 private assignment universe; no list theorem is used here.
        fixed_counts = [0] * q
        witnesses = [None] * q
        for word in product(range(q), repeat=8):
            if proper(word, private_edges):
                fixed_counts[word[0]] += 1
                if witnesses[word[0]] is None:
                    witnesses[word[0]] = word
        need(len(set(fixed_counts)) == 1, ("private palette permutation symmetry", q))
        table = []
        for pi in all_patterns:
            blocks = max(pi) + 1
            if blocks > q:
                continue
            forbidden = {pi[u] for u, v in cross_edges}
            allowed = [c for c in range(q) if c not in forbidden]
            count = sum(fixed_counts[c] for c in allowed)
            witness = witnesses[allowed[0]]
            need(count > 0 and proper(pi + witness, full_edges), ("full actual boundary extension", q, pi))
            named = factorial(q) // factorial(q - blocks)
            table.append({"partition": pi, "blocks": blocks, "extension_count": count,
                          "named_boundary_assignments": named, "private_witness": witness})
        need(len(table) == {4: 16, 5: 42, 6: 52}[q], ("full boundary spectrum", q))
        boundary_tables[str(q)] = table
        summaries[str(q)] = {"private_assignment_universe": q ** 8,
                              "proper_cube_colourings": sum(fixed_counts),
                              "boundary_partition_types": len(table),
                              "extensions_per_boundary": sorted({r["extension_count"] for r in table}),
                              "full_graph_colourings": sum(r["extension_count"] * r["named_boundary_assignments"] for r in table)}
    need(summaries["5"]["extensions_per_boundary"] == [17796], "all 42 proper five-colour states retained")

    artifacts = {
        "polynomial_certificates.json": checked_coefficients,
        "metric_exclusions.json": exclusions,
        "boundary_tables.json": boundary_tables,
        "coordinates.json": {"field_basis": ["1", "sqrt3", "sqrt11", "sqrt33"],
                             "old_spindle": [[a.serial() for a in p] for p in P],
                             "private_cube": [[a.serial() for a in p] for p in Q],
                             "full_unit_edges": sorted(full_edges), "cross_edges_old_private": cross_edges,
                             "all_unordered_pairs_checked": 105, "unit_pairs": 25, "nonunit_pairs": 80},
        "spindle_contact_audit.json": {"triples_checked": 35, "unit_triples": unit_triples,
                                       "external_unit_triple_centres": 0},
    }
    serialized = {name: json.dumps(obj, indent=2, sort_keys=True) + "\n" for name, obj in artifacts.items()}
    receipt = {"status": "PASS", "scope": "eight-private-point fixed-Moser extension certificates",
               "verifier_sha256": hashlib.sha256(source).hexdigest(),
               "complete_core_universe": 1 << 28, "metric_filter_graphs": 724080,
               "classes": 38, "independent_polynomial_certificates": 30,
               "moser_exterior_exclusions": 6, "circle_fan_exclusions": 2,
               "all_actual_pairs_checked": 105, "actual_unit_pairs": 25,
               "actual_nonunit_pairs": 80, "control_chromatic_number": 4,
               "boundary_controls": summaries,
               "artifacts_sha256": {name: hashlib.sha256(data.encode()).hexdigest() for name, data in serialized.items()},
               "source_pins_sha256": hashlib.sha256((here / "SOURCE_PINS.json").read_bytes()).hexdigest()}
    for path, data in sources.items():
        need((root / path).read_bytes() == data, ("dependency changed during run", path))
    need(Path(__file__).read_bytes() == source, "verifier changed during run")
    if "--write" in sys.argv:
        for name, data in serialized.items():
            (here / name).write_text(data, encoding="utf-8", newline="\n")
        (here / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    else:
        for name, data in serialized.items():
            need((here / name).read_bytes() == data.encode(), ("frozen artifact", name))
        need(json.loads((here / "receipt.json").read_text(encoding="utf-8")) == receipt, "frozen receipt")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
