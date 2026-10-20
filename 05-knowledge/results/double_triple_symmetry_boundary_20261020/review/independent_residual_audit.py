#!/usr/bin/env python3
"""Independent expanded-adjacency audit, without producer imports.

Universe: n=11,12,13 graphs with two distinguished independent-twin blocks
of sizes a,b>=5; every other vertex is a separate unrestricted vertex.
Every edge setting on this weighted presentation is enumerated.  This is
strictly larger than the four quotient patterns in the producer's audit.
Automorphism signs and both subgroup predicates are tested on expanded
vertex permutations, not from weighted-capacity/sign formulas.
"""
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
import hashlib
import json

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def expand(weights, bits):
    owners = [i for i, w in enumerate(weights) for _ in range(w)]
    k = len(weights)
    quotient = [[0] * k for _ in range(k)]
    for bit, (i, j) in enumerate(combinations(range(k), 2)):
        quotient[i][j] = quotient[j][i] = (bits >> bit) & 1
    return tuple(tuple(quotient[i][j] for j in owners) for i in owners)


def from_edges(weights, edges):
    pairs = list(combinations(range(len(weights)), 2))
    bits = sum(1 << pairs.index(tuple(sorted(edge))) for edge in edges)
    return expand(weights, bits)


def twin_classes(adjacency):
    classes = {}
    for i, row in enumerate(adjacency):
        classes.setdefault(row, []).append(i)
    return list(classes.values())


def canonical(adjacency):
    classes = twin_classes(adjacency)
    k = len(classes)
    records = []
    for q in permutations(range(k)):
        records.append((
            tuple(len(classes[i]) for i in q),
            tuple(adjacency[classes[q[i]][0]][classes[q[j]][0]]
                  for i, j in combinations(range(k), 2)),
        ))
    return min(records)


def is_automorphism(adjacency, p):
    return all(adjacency[i][j] == adjacency[p[i]][p[j]]
               for i, j in combinations(range(len(p)), 2))


def edge_sign(adjacency, p):
    # Direct definition on individual graph edges; no class-weight formula.
    return (-1) ** sum(adjacency[i][j] and p[i] > p[j]
                       for i, j in combinations(range(len(p)), 2))


def cycle(n, triple, reverse=False):
    p = list(range(n))
    a, b, c = triple
    if reverse:
        b, c = c, b
    p[a], p[b], p[c] = b, c, a
    return tuple(p)


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def predicates(adjacency, support):
    """Test actual 3+3 permutations and actual independent triple generators."""
    n = len(adjacency)
    remaining = tuple(i for i in range(n) if i not in support)
    triples = list(combinations(remaining, 3))
    cycles = {t: cycle(n, t) for t in triples}
    single_automorphism = {t: is_automorphism(adjacency, cycles[t])
                           for t in triples}
    synchronous = independent = False
    for index, s in enumerate(triples):
        for t in triples[index + 1:]:
            if set(s).intersection(t):
                continue
            if single_automorphism[s] and single_automorphism[t]:
                independent = synchronous = True
                return synchronous, independent
            if not synchronous:
                for reverse in (False, True):
                    p = compose(cycles[s], cycle(n, t, reverse))
                    if is_automorphism(adjacency, p):
                        synchronous = True
    return synchronous, independent


def analyze(adjacency):
    n = len(adjacency)
    classes = twin_classes(adjacency)
    k = len(classes)
    weights = [len(c) for c in classes]
    automorphisms = []
    for q in permutations(range(k)):
        if any(weights[i] != weights[q[i]] for i in range(k)):
            continue
        p = list(range(n))
        for i in range(k):
            for a, b in zip(classes[i], classes[q[i]]):
                p[a] = b
        if is_automorphism(adjacency, p):
            automorphisms.append((q, edge_sign(adjacency, p)))
    # All omitted internal symmetric-group generators are checked directly.
    for c in classes:
        for x, y in zip(c, c[1:]):
            p = list(range(n))
            p[x], p[y] = p[y], p[x]
            require(is_automorphism(adjacency, p), "twin generator not an automorphism")
            require(edge_sign(adjacency, p) == 1, "internal twin sign is not even")
    eligible_d, eligible_i = [], []
    for i, c in enumerate(classes):
        if len(c) < 5:
            continue
        d, independent = predicates(adjacency, set(c[:5]))
        if d:
            eligible_d.append(i)
        if independent:
            eligible_i.append(i)
    even = all(sign == 1 for _, sign in automorphisms)
    result = {}
    for name, eligible in (("synchronized", eligible_d), ("independent", eligible_i)):
        relative = [i for i in eligible
                    if all(sign == 1 for q, sign in automorphisms if q[i] == i)]
        orbits = []
        seen = set()
        for i in relative:
            if i not in seen:
                orbit = sorted({q[i] for q, _ in automorphisms})
                seen.update(orbit)
                orbits.append(orbit)
        result[name] = {
            "eligible_classes": eligible,
            "relative_even_classes": relative,
            "relative_copy_orbits": orbits,
            "fusion": max(0, len(orbits) - 1) if even else 0,
            "escape": len(orbits) if not even else 0,
        }
    return {
        "weights": weights,
        "quotient_edges": [[i, j] for i, j in combinations(range(k), 2)
                           if adjacency[classes[i][0]][classes[j][0]]],
        "degrees": [sum(adjacency[c[0]]) for c in classes],
        "globally_even": even,
        "quotient_automorphism_order": len(automorphisms),
        **result,
    }


def main():
    rows = []
    for n in (11, 12, 13):
        examples = {}
        raw = 0
        presentations = []
        for a in range(5, n - 4):
            for b in range(a, n - a + 1):
                weights = (a, b) + (1,) * (n - a - b)
                count = 1 << (len(weights) * (len(weights) - 1) // 2)
                presentations.append({"weights": weights, "edge_settings": count})
                for bits in range(count):
                    raw += 1
                    adjacency = expand(weights, bits)
                    examples.setdefault(canonical(adjacency), adjacency)
        anomalies = []
        totals = {name: Counter() for name in ("synchronized", "independent")}
        predicate_differences = 0
        for key, adjacency in sorted(examples.items()):
            record = analyze(adjacency)
            if record["synchronized"]["eligible_classes"] != record["independent"]["eligible_classes"]:
                predicate_differences += 1
            for name in totals:
                totals[name].update({key: record[name][key] for key in ("fusion", "escape")})
            if any(record[name]["fusion"] or record[name]["escape"] for name in totals):
                anomalies.append(record)
        expected = {"fusion": 4, "escape": 2} if n == 13 else {"fusion": 0, "escape": 0}
        for name in totals:
            require(dict(totals[name]) == expected, f"unexpected anomaly totals n={n}, {name}")
        require(predicate_differences == 0, "predicates differ in broad candidate universe")
        require(len(anomalies) == (6 if n == 13 else 0), "unexpected anomaly type count")
        rows.append({
            "n": n, "presentations": presentations,
            "raw_presentations": raw, "distinct_global_graph_types": len(examples),
            "anomaly_totals": totals, "predicate_differences": predicate_differences,
            "anomalies": anomalies,
        })
    hostile_controls = []
    for name, weights, edges, eligible_count, fusion, escape in (
        ("K55 plus I2", (5, 5, 2), ((0, 1),), 0, 0, 0),
        ("K55 plus I3", (5, 5, 3), ((0, 1),), 2, 0, 1),
        ("K55 plus K3", (5, 5, 1, 1, 1), ((0, 1), (2, 3), (2, 4), (3, 4)), 2, 0, 0),
        ("K56", (5, 6), ((0, 1),), 1, 0, 0),
        ("K57", (5, 7), ((0, 1),), 1, 0, 0),
        ("K66", (6, 6), ((0, 1),), 2, 0, 0),
    ):
        adjacency = from_edges(weights, edges)
        record = analyze(adjacency)
        for name_predicate in ("synchronized", "independent"):
            result = record[name_predicate]
            require(len(result["eligible_classes"]) == eligible_count,
                    f"wrong eligible-copy count for {name}")
            require((result["fusion"], result["escape"]) == (fusion, escape),
                    f"wrong anomaly count for {name}")
        if name == "K55 plus I2":
            n = len(adjacency)
            require(any(is_automorphism(adjacency, cycle(n, triple))
                        for triple in combinations(range(5, n), 3)),
                    "single-triple hostile is not admitted")
        hostile_controls.append({"name": name, **record})
    circulant_controls = []
    for kind, distances, multiplier in (
        ("five-cycle", {1, 4}, -1),
        ("complementary-five-cycle", {2, 3}, -1),
        ("complete-five-graph", {1, 2, 3, 4}, 2),
    ):
        adjacency = tuple(tuple(int((j - i) % 5 in distances)
                                for j in range(5)) for i in range(5))
        p = tuple(multiplier * i % 5 for i in range(5))
        require(is_automorphism(adjacency, p), "nonempty circulant witness failed")
        require(edge_sign(adjacency, p) == -1, "nonempty circulant witness not odd")
        circulant_controls.append({"kind": kind, "permutation": p, "edge_sign": -1})
    # Sharpness control: three regular cyclic modules, including m=5 used at n=15.
    sharp_controls = []
    for m in (3, 5, 7, 9):
        n = 3 * m
        adjacency = tuple(tuple(
            int((j // m - i // m) % 3 == 1) if i // m != j // m
            else int(0 < (j - i) % m <= (m - 1) // 2)
            for j in range(n)) for i in range(n))
        require(all(sum(row) == (n - 1) // 2 for row in adjacency), "sharp control not regular")
        p = tuple((i + 1) % m if i < m else i for i in range(n))
        require(is_automorphism(adjacency, p), "single-cycle sharp control failed")
        sharp_controls.append({"module_size": m, "degree": n, "outdegree": (n - 1) // 2})
    print(json.dumps({
        "status": "PASS",
        "scope": "All expanded graphs from two independent twin blocks >=5 plus unrestricted residual vertices, n=11,12,13; not all graphs or a subgroup census.",
        "method": "Direct expanded-vertex adjacency and edge-inversion sign tests; explicit double-three-cycle and independent-three-cycle tests; no producer imports.",
        "checks": CHECKS,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "rows": rows, "hostile_controls": hostile_controls,
        "nonempty_circulant_controls": circulant_controls,
        "sharp_module_controls": sharp_controls,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
