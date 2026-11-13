#!/usr/bin/env python3
"""Small exact controls; no producer imports, assertions or third-party packages.

Enumerates actual graceful completions, all C2 edge transports and vertex
gauges, and plus as a colimit over ALL covering sieves of finite discrete sites.
Prints deterministic JSON.  Mathematical generalizations remain in PROOF.md.
"""
from collections import Counter
from itertools import combinations, product
import hashlib
import json
from pathlib import Path


CHECKS = 0


def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


class UF:
    def __init__(self, n):
        self.p = list(range(n))

    def root(self, x):
        while self.p[x] != x:
            x = self.p[x]
        return x

    def join(self, x, y):
        x, y = self.root(x), self.root(y)
        if x != y:
            self.p[y] = x


def tree(n, edges):
    u = UF(n)
    for a, b in edges:
        if u.root(a) == u.root(b):
            return False
        u.join(a, b)
    return len(edges) == n - 1 and len({u.root(i) for i in range(n)}) == 1


def completion_control():
    fixed = ((0, 3),)
    choices = [tuple((i, i+d) for i in range(4-d)) for d in (1, 2)]
    states = [s for s in product(*choices) if tree(4, fixed + s)]
    expected = [((0, 1), (0, 2)), ((1, 2), (0, 2)),
                ((1, 2), (1, 3)), ((2, 3), (1, 3))]
    require(states == expected, "literal upper-tail completions")
    edges = [(i, j) for i, j in combinations(range(4), 2)
             if sum(a != b for a, b in zip(states[i], states[j])) <= 2]
    require(len(edges) == 6, "width-two graph is K4")
    all_bits = list(product((0, 1), repeat=4))
    orbit_keys = Counter()
    sections_histogram = Counter()
    cover_components_histogram = Counter()
    records = []
    for signs in product((0, 1), repeat=6):
        sections = [x for x in all_bits if all(
            x[i] ^ x[j] == s for (i, j), s in zip(edges, signs))]
        triangle_parities = [signs[edges.index((0, i))] ^
            signs[edges.index((0, j))] ^ signs[edges.index((i, j))]
            for i, j in combinations(range(1, 4), 2)]
        require(bool(sections) == (not any(triangle_parities)),
                "tree fundamental-cycle holonomy")
        require(len(sections) in (0, 2), "connected torsor section count")
        gauge_orbit = {tuple(s ^ x[i] ^ x[j]
            for (i, j), s in zip(edges, signs)) for x in all_bits}
        require(len(gauge_orbit) == 8, "global gauge kernel has size two")
        orbit_keys[min(gauge_orbit)] += 1
        sections_histogram[len(sections)] += 1
        cover = UF(8)
        for (i, j), s in zip(edges, signs):
            for b in (0, 1):
                cover.join(2*i+b, 2*j+(b ^ s))
        component_count = len({cover.root(i) for i in range(8)})
        require(component_count == (2 if sections else 1), "double-cover splitting")
        cover_components_histogram[component_count] += 1
        records.append({"edge_signs": signs, "fundamental_holonomies": triangle_parities,
                        "section_count": len(sections), "cover_components": component_count})
    require(len(orbit_keys) == 8, "eight gauge classes")
    require(set(orbit_keys.values()) == {8}, "equal gauge-orbit sizes")
    require(dict(sections_histogram) == {2: 8, 0: 56}, "positive and hostile totals")
    return {"fixed_edges": fixed, "free_colors": [1, 2], "completions": states,
            "exchange_edges": edges, "cycle_rank": 3, "gauge_classes": len(orbit_keys),
            "section_count_histogram": dict(sorted(sections_histogram.items())),
            "cover_components_histogram": dict(sorted(cover_components_histogram.items())),
            "all_transports": records}


def subset(v, u):
    return v & u == v


def sieves(n):
    result = {}
    for u in range(1 << n):
        objects = [v for v in range(1 << n) if subset(v, u)]
        covers = []
        for bits in product((0, 1), repeat=len(objects)):
            r = tuple(v for v, bit in zip(objects, bits) if bit)
            if any(w not in r for v in r for w in objects if subset(w, v)):
                continue
            union = 0
            for v in r:
                union |= v
            if union == u:
                covers.append(r)
        result[u] = covers
    return result


def constant(n, normalize_empty=False):
    values = {u: ((0,) if normalize_empty and not u else (0, 1))
              for u in range(1 << n)}
    restrictions = {(u, v): {s: (0 if normalize_empty and not v else s)
                           for s in values[u]}
                    for u in values for v in values if subset(v, u)}
    return values, restrictions


def matching(r, f):
    values, restrictions = f
    ordered = sorted(r, key=lambda u: (u.bit_count(), u))
    out = []
    def visit(chosen):
        if len(chosen) == len(ordered):
            out.append(tuple(chosen[v] for v in r))
            return
        v = ordered[len(chosen)]
        for s in values[v]:
            if all(restrictions[v, w][s] == chosen[w]
                   for w in chosen if subset(w, v)):
                chosen[v] = s
                visit(chosen)
                del chosen[v]
    visit({})
    return out


def plus(f, covers):
    """Literal finite colimit: glue matching families along all refinements."""
    tables, representatives, values = {}, {}, {}
    for u, cc in covers.items():
        nodes = [(r, a) for r in cc for a in matching(r, f)]
        index = {node: i for i, node in enumerate(nodes)}
        uf = UF(len(nodes))
        for i, (r, a) in enumerate(nodes):
            lookup = dict(zip(r, a))
            for s in cc:
                if set(s) <= set(r):
                    b = tuple(lookup[v] for v in s)
                    uf.join(i, index[s, b])
        groups = {}
        for i, node in enumerate(nodes):
            groups.setdefault(uf.root(i), []).append(node)
        ordered = sorted(groups.values(), key=lambda group: min(group))
        tables[u] = {node: k for k, group in enumerate(ordered) for node in group}
        representatives[u] = [min(group) for group in ordered]
        values[u] = tuple(range(len(ordered)))
    restrictions = {}
    for u in covers:
        for v in covers:
            if not subset(v, u):
                continue
            mapping = {}
            for (r, a), value in tables[u].items():
                pulled = tuple(w for w in r if subset(w, v))
                lookup = dict(zip(r, a))
                result = tables[v][pulled, tuple(lookup[w] for w in pulled)]
                if value in mapping:
                    require(mapping[value] == result, "restriction independent of representative")
                mapping[value] = result
            restrictions[u, v] = mapping
    return values, restrictions


def inspect_presheaf(f, covers):
    values, restrictions = f
    separated, sheaf = True, True
    for u in covers:
        for s in values[u]:
            require(restrictions[u, u][s] == s, "presheaf identity")
            for v in covers:
                if subset(v, u):
                    for w in covers:
                        if subset(w, v):
                            require(restrictions[v, w][restrictions[u, v][s]] ==
                                    restrictions[u, w][s], "presheaf composition")
        for r in covers[u]:
            restrictions_to_cover = [tuple(restrictions[u, v][s] for v in r)
                                     for s in values[u]]
            is_injective = len(set(restrictions_to_cover)) == len(values[u])
            is_surjective = set(restrictions_to_cover) == set(matching(r, f))
            separated &= is_injective
            sheaf &= is_injective and is_surjective
    return {"cardinalities": [len(values[u]) for u in covers],
            "separated": separated, "sheaf": sheaf}


def plus_control():
    results = []
    for n in (1, 2, 3):
        covers = sieves(n)
        for normalized in (False, True):
            f = constant(n, normalized)
            phases = []
            for k in range(4):
                phase = inspect_presheaf(f, covers)
                expected = [1 if u == 0 and (normalized or k > 0) else
                            2 ** u.bit_count() if k >= 2 or (normalized and k >= 1)
                            else 2 for u in covers]
                require(phase["cardinalities"] == expected, "plus cardinality formula")
                require(phase["separated"] == (normalized or k > 0), "separatedness stage")
                require(phase["sheaf"] == (k >= 2 or (k >= 1 and (normalized or n == 1))
                        or (normalized and n == 1)), "sheaf stage")
                phases.append(phase)
                if k < 3:
                    f = plus(f, covers)
            results.append({"points": n, "normalize_empty": normalized,
                            "cover_sieve_counts": [len(covers[u]) for u in covers],
                            "phases_0_through_3": phases})
    return results


def refinement_control():
    out = []
    for length in range(2, 10):
        # Observed state 0 is distinct; continuation moves one unit toward 0.
        colors = [0] + [1] * length
        counts = [len(set(colors))]
        while True:
            signatures = [(colors[i], colors[max(0, i-1)]) for i in range(length+1)]
            ids = {s: k for k, s in enumerate(sorted(set(signatures)))}
            refined = [ids[s] for s in signatures]
            if len(set(refined)) == len(set(colors)):
                break
            colors = refined
            counts.append(len(set(colors)))
        require(counts == list(range(2, length+2)), "arbitrary refinement depth")
        out.append({"path_length": length, "block_counts": counts,
                    "strict_refinement_rounds": len(counts)-1})
    return out


def main():
    result = {"status": "FINITE-EXACT", "completion_transport": completion_control(),
              "plus_all_covering_sieves": plus_control(),
              "deterministic_refinement": refinement_control()}
    result["checks"] = CHECKS
    result["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
