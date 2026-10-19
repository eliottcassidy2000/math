#!/usr/bin/env python3
"""Ten labelled coordinates, joint matching predicates, and marked Mathieu action.

Exact finite universe, no Collatz acceptance or geometric realization assumed.
Independent constructions: symplectic span enumeration versus H8 embeddings;
explicit matchings versus signature multiplicities; group orbit versus designs.
Checks remain active under python -O. Standard library only.
"""
import argparse
from collections import Counter
from functools import reduce
from itertools import combinations, permutations
import json
from math import prod
from operator import xor
from pathlib import Path

from cardinality_carriers_20261018m import carrier, matchings, relabel

CHECKS = 0


def check(ok, reason):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(reason)


def span(rows):
    result = {0}
    for row in rows:
        result |= {x ^ row for x in result}
    return frozenset(result)


def closure(gens):
    identity = tuple(range(len(gens[0])))
    seen, stack = {identity}, [identity]
    while stack:
        p = stack.pop()
        for g in gens:
            q = tuple(g[i] for i in p)
            if q not in seen:
                seen.add(q)
                stack.append(q)
    return seen


def mobius(p, a, b, c, d):
    def evaluate(x):
        if x == p:
            return a * pow(c, -1, p) % p if c % p else p
        num, den = (a*x+b) % p, (c*x+d) % p
        return num * pow(den, -1, p) % p if den else p
    return tuple(evaluate(x) for x in range(p+1))


def orbit_sets(base, gens):
    seen, stack = {frozenset(base)}, [frozenset(base)]
    while stack:
        old = stack.pop()
        for g in gens:
            new = frozenset(g[x] for x in old)
            if new not in seen:
                seen.add(new)
                stack.append(new)
    return seen


def signature_count(rows, n=10):
    counts = Counter(sum(((r >> i) & 1) << j for j, r in enumerate(rows))
                     for i in range(n))
    if any(v % 2 for v in counts.values()):
        return 0
    return prod(prod(range(1, v, 2)) for v in counts.values())


def matching_accepts(matching, rows):
    return all(((r >> i) & 1) == ((r >> j) & 1)
               for i, j in matching for r in rows)


def group_orbits(objects, gens):
    pending = set(objects)
    sizes = []
    while pending:
        base = min(pending, key=lambda s: tuple(sorted(s)))
        orbit, stack = {base}, [base]
        while stack:
            old = stack.pop()
            for g in gens:
                new = frozenset(relabel(v, g) for v in old)
                check(new in objects, "permutation preserves code universe")
                if new not in orbit:
                    orbit.add(new)
                    stack.append(new)
        pending -= orbit
        sizes.append(len(orbit))
    return sorted(sizes)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    canonical, spaces, images, census = carrier(4)
    check(census["isotropic_subspaces_by_dimension"] == [1, 255, 5355, 11475, 2295],
          "all isotropic subspace levels")
    codes = {frozenset(v for x in s for v in (x, x ^ 1023)) for s in spaces}
    matched = {frozenset(v for x in s for v in (x, x ^ 1023)) for s in images}
    _, spaces8, images8, _ = carrier(3)
    hamming8 = {frozenset(v for x in s for v in (x, x ^ 255))
                for s in spaces8 - images8.keys()}
    embedded = {}
    for pair in combinations(range(10), 2):
        pairmask = sum(1 << i for i in pair)
        rest = [i for i in range(10) if i not in pair]
        for h in hamming8:
            lift = {sum(1 << rest[i] for i in range(8) if v >> i & 1) for v in h}
            code = frozenset(lift | {v ^ pairmask for v in lift})
            check(code not in embedded, "unique distinguished pair")
            embedded[code] = pair
    check(set(embedded) == codes - matched, "independent classification")
    weights = Counter()
    for c in codes:
        check(len(c) == 32, "dimension five")
        check(all((a ^ b) in c and (a & b).bit_count() % 2 == 0
                  for a in c for b in c), "linear self-dual")
        weights[sum(v.bit_count() == 2 for v in c)] += 1
    check(weights == {5: 945, 1: 1350}, "weight two decoder")

    # Rank-one changes of controller: common 16-word code is held fixed.
    surgeries = {}
    incoming = Counter()
    for matching in matchings(tuple(range(10))):
        pairmasks = [sum(1 << i for i in pair) for pair in matching]
        old = span(pairmasks)
        for p in pairmasks:
            rest = [r for r in pairmasks if r != p]
            common = span([p] + [rest[0] ^ r for r in rest[1:]])
            extensions = [c for c in codes if common <= c]
            check(len(common) == 16 and len(extensions) == 3 and old in extensions,
                  "three completions of fixed common observations")
            new = frozenset(c for c in extensions if c != old)
            check(len(new) == 2 and all(c in embedded and p in c for c in new),
                  "two Hamming replacements retain the distinguished pair")
            surgeries[(old, p)] = new
            incoming.update(new)
    check(set(incoming) == set(embedded) and set(incoming.values()) == {7},
          "ten-to-seven incidence degrees")
    # Fixed pair 8,9: 30 Hamming codes form a 7-regular bipartite graph.
    fibre = {c for c, p in embedded.items() if p == (8, 9)}
    edges = {tuple(sorted(new, key=lambda c: tuple(sorted(c))))
             for (old, p), new in surgeries.items() if p == 768}
    adjacency = {c: set() for c in fibre}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    colours = {}
    for c in sorted(fibre, key=lambda c: tuple(sorted(c))):
        if c in colours:
            continue
        colours[c], stack = 0, [c]
        while stack:
            v = stack.pop()
            for w in adjacency[v]:
                if w not in colours:
                    colours[w] = 1-colours[v]
                    stack.append(w)
                check(colours[w] != colours[v], "Hamming extension graph is bipartite")
    check(Counter(colours.values()) == {0: 15, 1: 15} and len(edges) == 105,
          "15+15 extension incidence")
    check(all(len(v) == 7 for v in adjacency.values()), "seven directions")
    check(all(len(adjacency[a] & adjacency[b]) == 3
              for a, b in combinations(fibre, 2) if colours[a] == colours[b]),
          "symmetric 2-(15,7,3) incidence")
    left = sorted((c for c in fibre if colours[c] == 0), key=lambda c: tuple(sorted(c)))
    right = sorted((c for c in fibre if colours[c] == 1), key=lambda c: tuple(sorted(c)))
    nonincidence = [sum(1 << j for j,c in enumerate(right) if c not in adjacency[a]) for a in left]
    check(span(nonincidence) == frozenset([0]+nonincidence), "projective nonincidence row space")
    basis, current = [], frozenset([0])
    for v in nonincidence:
        if v not in current:
            basis.append(v)
            current = span(basis)
    check(len(basis) == 4, "actual PG(3,2) coordinate space")
    coordinates = {k: reduce(xor, (basis[i] for i in range(4) if k >> i & 1), 0)
                   for k in range(16)}
    inverse = {v:k for k,v in coordinates.items()}
    points = [inverse[v] for v in nonincidence]
    planes = [sum(((b >> j)&1) << i for i,b in enumerate(basis)) for j in range(15)]
    check(set(points) == set(planes) == set(range(1,16)), "all projective points and dual planes")
    check(all(((points[i]&planes[j]).bit_count() % 2 == 0) == (right[j] in adjacency[left[i]])
              for i in range(15) for j in range(15)), "explicit point-plane isomorphism")

    # Fixed visible labels 0..7; hidden labels 8,9. Retain actual fibre counts.
    projection_hist = {"matching": Counter(), "missing": Counter()}
    multiplicity = {"matching": [0] * 256, "missing": [0] * 256}
    for c in codes:
        kind = "matching" if c in matched else "missing"
        fibres = Counter(v & 255 for v in c)
        projection_hist[kind][len(fibres)] += 1
        check(set(fibres.values()) == {32 // len(fibres)}, "all marked fibres")
        for v in fibres:
            multiplicity[kind][v] += 1
    profile = []
    for w in range(9):
        values = {tuple(multiplicity[k][v] for k in ("matching", "missing"))
                  for v in range(256) if v.bit_count() == w}
        check(len(values) == 1, "visible permutation invariance")
        a, b = values.pop()
        expected_all = 2295 if w in (0, 8) else 270 if w % 2 else 255
        df = lambda t: prod(range(1, t + 1, 2))
        expected_match = 2 * df(w) * df(8 - w) if w % 2 else 9 * df(w - 1) * df(7 - w)
        check(a == expected_match and a + b == expected_all, "derived support profile")
        profile.append({"port_weight": w, "matching": a, "missing": b, "all": a + b})
    check(all(multiplicity["matching"]), "matching union already covers every port word")

    explicit = [170, 204, 240, 255]
    all_matchings = list(matchings(tuple(range(10))))
    joint = []
    for k in range(1, 5):
        rows = explicit[:k]
        actual = sum(matching_accepts(m, rows) for m in all_matchings)
        check(actual == signature_count(rows), "signature compiler")
        joint.append({"rows": rows, "matching_extensions": actual,
                      "all_code_extensions": sum(all(r in c for r in rows) for c in codes)})
    check([x["matching_extensions"] for x in joint] == [45, 3, 0, 0], "joint hostile")
    check([x["all_code_extensions"] for x in joint] == [135, 15, 3, 1], "isotropic extension")
    chart_linear = span([15, 240, 768])
    affine_charts = [{v ^ r for v in chart_linear} for r in (0,170,204,102)]
    check(all(len(c) == 8 for c in affine_charts)
          and len(set.union(*affine_charts)) == 32
          and set.union(*affine_charts) == set(span([170,204,240,255,768])),
          "optimal four affine partition charts for ten-label Hamming relation")
    for offset, chart in zip((0,170,204,102),affine_charts):
        direct = {x for x in range(1024) if all(
            (((x^offset)>>i)&1) == (((x^offset)>>j)&1)
            for block in ((0,1,2,3),(4,5,6,7),(8,9)) for i in block for j in block)}
        check(direct == chart, "binary equality/disequality chart with original ports")
    # Exhaust all row pairs among all 512 even words (ordered, repetitions allowed).
    even = [r for r in range(1024) if r.bit_count() % 2 == 0]
    pair_tests = 0
    for a in even:
        for b in even:
            check(bool(signature_count([a, b])) == ((a & b).bit_count() % 2 == 0),
                  "pairwise orthogonality iff matching compatibility")
            pair_tests += 1
    # All ordered triples within a representative H8+I2 code, including dependence.
    representative = span([170, 204, 240, 255, 768])
    triple_tests = 0
    for a in representative:
        for b in representative:
            for c in representative:
                check(bool(signature_count([a, b, c])) == ((a & b & c).bit_count() % 2 == 0),
                      "cubic intersection detects three-row obstruction")
                triple_tests += 1

    # Same audited M12 generators as the repo's September moonshine construction.
    gens12 = [mobius(11, 1, 1, 0, 1), mobius(11, 3, 0, 0, 1),
              mobius(11, 0, 10, 1, 0)]
    d = list(range(12))
    for a, b in [(2, 10), (3, 4), (5, 9), (6, 7)]:
        d[a], d[b] = b, a
    gens12.append(tuple(d))
    m12 = closure(gens12)
    check(len(m12) == 95040 and len({g[:5] for g in m12}) == 95040,
          "M12 sharply five-transitive")
    m10 = {g[:10] for g in m12 if g[10:] == (10, 11)}
    check(len(m10) == 720 and len({g[:3] for g in m10}) == 720,
          "marked two-point stabilizer sharply three-transitive")
    gens10, generated = [], {tuple(range(10))}
    for g in sorted(m10):
        if g not in generated:
            gens10.append(g)
            generated = closure(gens10)
        if len(generated) == 720:
            break
    check(generated == m10, "small generator set for marked stabilizer")
    orbit_census = {"matching": group_orbits(matched, gens10),
                    "missing": group_orbits(set(embedded), gens10)}
    hexads = None
    for base in combinations(range(12), 6):
        orb = orbit_sets(base, gens12)
        if len(orb) == 132:
            hexads = orb
            break
    check(hexads is not None, "Witt hexads")
    w10 = {frozenset(h - {10, 11}) for h in hexads if {10, 11} <= h}
    check(len(w10) == 30, "twice derived Witt design")
    for triple in combinations(range(10), 3):
        check(sum(set(triple) <= block for block in w10) == 1, "S(3,4,10)")
    w10masks = {sum(1 << i for i in b) for b in w10}
    derivations = Counter()
    selected_hamming = set()
    for pair in combinations(range(10), 2):
        pairs = [pair] + [tuple(sorted(b - set(pair))) for b in w10 if set(pair) <= b]
        check(len(pairs) == 5 and sorted(i for p in pairs for i in p) == list(range(10)),
              "derived Witt matching preserves marked pair")
        c = span([sum(1 << i for i in p) for p in pairs])
        check(c in matched, "Witt derived matching is in actual carrier")
        derivations[c] += 1
        selected_hamming.update(surgeries[(c, sum(1 << i for i in pair))])
    intersections = Counter((a & b).bit_count() for a in w10masks for b in w10masks if a != b)
    check(any(k % 2 for k in intersections), "hostile: W10 blocks do not form self-dual code")
    check(len(selected_hamming) == 90 and group_orbits(selected_hamming, gens10) == [90],
          "Witt anchored surgery selects precisely the small Hamming orbit")

    # The next primitive object: even unions of six pairs plus a transversal.
    pair_transversal = []
    for k in (4,5,6,8):
        ps = [3 << (2*i) for i in range(k)]
        t = sum(1 << (2*i) for i in range(k))
        dcode = span([ps[0]^v for v in ps[1:]]+[t])
        self_orthogonal = all((a&b).bit_count() % 2 == 0 for a in dcode for b in dcode)
        check(self_orthogonal == (k % 2 == 0), "transversal parity gate")
        distribution = Counter(v.bit_count() for v in dcode)
        if k == 6:
            check(distribution == {0:1,4:15,6:32,8:15,12:1}, "primitive length-twelve code")
            tetrads = [v for v in dcode if v.bit_count() == 4]
            signatures = [tuple((v>>i)&1 for v in tetrads) for i in range(12)]
            check(all((signatures[i] == signatures[j]) == (i//2 == j//2)
                      for i in range(12) for j in range(12)), "unique matching from tetrad incidence")
        pair_transversal.append({"pairs":k,"self_dual":self_orthogonal,
                                 "weights":dict(sorted(distribution.items()))})
    twelve = {"matching":prod(range(1,12,2)),"two_pairs_plus_H8":1485*30,
              "primitive_transversal":2*prod(range(1,12,2)),"all":prod(2**i+1 for i in range(1,6))}
    check(sum(twelve[k] for k in ("matching","two_pairs_plus_H8","primitive_transversal")) == twelve["all"],
          "twelve-label mass closes after the new primitive family")

    out = {"status": "FINITE-EXACT with independent constructions; no universal Collatz or geometric claim",
           "universe": "all labelled binary self-dual length-ten codes, equivalently V4 Lagrangians",
           "census": census, "weight_two_histogram": dict(sorted(weights.items())),
           "port_projection_sizes": {k: dict(sorted(v.items())) for k, v in projection_hist.items()},
           "port_support_multiplicities": profile, "joint_rows": joint,
           "pair_tests": pair_tests, "triple_tests": triple_tests,
           "optimal_partition_chart_cover": {"charts":4,"words_per_chart":8,
                                             "coset_representatives":[0,170,204,102]},
           "next_boundary":{"pair_transversal_representatives":pair_transversal,
                            "twelve_label_counts_from_proved_injections":twelve},
           "surgeries": {"marked_matching_inputs": len(surgeries),
                         "incidences": sum(incoming.values()),
                         "matching_degree": 10, "missing_degree": 7,
                         "fixed_pair_Hamming_vertices": 30,
                         "fixed_pair_Hamming_edges": 105,
                         "fixed_pair_bipartition": [15, 15],
                         "same_part_common_neighbours": 3},
           "projective_incidence_coordinates":{"left_points":points,"right_planes":planes},
           "mathieu": {"M12_order": len(m12), "M10_order": len(m10),
                       "M10_generators": gens10, "M10_code_orbits": orbit_census,
                       "W10_blocks": len(w10), "W10_binary_span_dimension": len(span(w10masks)).bit_length()-1,
                       "W10_ordered_distinct_block_intersections": dict(sorted(intersections.items())),
                       "W10_derived_matching_count": len(derivations),
                       "W10_selected_Hamming_orbit": len(selected_hamming),
                       "W10_derived_matching_fibres": dict(sorted(Counter(derivations.values()).items()))},
           "checks": CHECKS}
    encoded = json.dumps(out, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
