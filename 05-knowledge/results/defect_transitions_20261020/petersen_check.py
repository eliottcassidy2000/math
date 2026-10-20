"""Exact labelled boundary and subdivision checks; no production data imports.

Run: python 05-knowledge/results/defect_transitions_20261020/petersen_check.py
Dependency: networkx (only named-graph isomorphism and graph diagnostics).
All assertions use explicit exceptions and remain active under python -O.
"""
from collections import Counter
from itertools import combinations, product
import networkx as nx


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def edge(a, b):
    check(a != b, "loop")
    return tuple(sorted((a, b)))


P = ((0, 1), (0, 4), (0, 5), (1, 2), (1, 6), (2, 3), (2, 7),
     (3, 4), (3, 8), (4, 9), (5, 7), (5, 8), (6, 8), (6, 9), (7, 9))
B1 = ((0, 5), (1, 17), (2, 14), (3, 8), (4, 17), (6, 11), (7, 17),
      (9, 13), (10, 15), (12, 16),
      *((i, i + 1) for i in range(16)), (0, 16))
B2 = ((0, 1), (0, 2), (0, 14), (1, 5), (1, 11), (2, 3), (2, 6),
      (3, 4), (3, 9), (4, 5), (4, 7), (5, 6), (6, 8), (7, 8),
      (7, 17), (8, 9), (9, 15), (10, 11), (10, 14), (10, 16),
      (11, 12), (12, 13), (12, 17), (13, 14), (13, 15),
      (15, 16), (16, 17))


def pole_profile(vertices, edges, ports):
    """Ports are ordered semiedges, represented by their interior endpoint."""
    vertices = tuple(sorted(vertices))
    edges = tuple(edge(*e) for e in edges)
    ports = tuple(ports)
    check(len(set(edges)) == len(edges), "duplicate edge")
    incidence = Counter(x for e in edges for x in e) + Counter(ports)
    check(set(incidence) == set(vertices), "unlisted or isolated vertex")
    check(all(incidence[v] == 3 for v in vertices), "noncubic pole interior")
    variables = tuple((u, v) for u, v in edges) + tuple((v,) for v in ports)
    used = {v: 0 for v in vertices}
    assigned = [-1] * len(variables)
    result, witnesses = Counter(), {}

    def visit(left):
        if not left:
            word = tuple(assigned[len(edges):])
            result[word] += 1
            witnesses.setdefault(word, tuple(assigned))
            return
        best, bestmask = -1, 0
        for i in left:
            mask = 7
            for v in variables[i]:
                mask &= 7 ^ used[v]
            if not mask:
                return
            if best < 0 or mask.bit_count() < bestmask.bit_count():
                best, bestmask = i, mask
                if mask.bit_count() == 1:
                    break
        rest = tuple(i for i in left if i != best)
        for c in range(3):
            bit = 1 << c
            if not bestmask & bit:
                continue
            assigned[best] = c
            for v in variables[best]:
                used[v] |= bit
            visit(rest)
            for v in variables[best]:
                used[v] ^= bit
        assigned[best] = -1

    visit(tuple(range(len(variables))))
    return result, witnesses


def deleted_edge_pole(edges, deleted):
    vertices = set(x for e in edges for x in e)
    u, v = deleted
    check(edge(u, v) in set(map(lambda e: edge(*e), edges)), "missing deleted edge")
    adj = {x: [] for x in vertices}
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    ports = tuple(sorted(set(adj[u]) - {v}) + sorted(set(adj[v]) - {u}))
    kept = tuple(e for e in edges if u not in e and v not in e)
    return tuple(sorted(vertices - {u, v})), kept, ports


def kempe_wires(vertices, edges, ports, witness):
    """Extract and audit the two 0/1 paths for a boundary 0011 witness."""
    check(tuple(witness[len(edges):]) == (0, 0, 1, 1), "wrong witness word")
    labels = {v: ("v", v) for v in vertices}
    adj = {labels[v]: [] for v in vertices}
    for i in range(4):
        adj[("p", i)] = []
    for (u, v), c in zip(edges, witness):
        if c in (0, 1):
            adj[labels[u]].append(labels[v])
            adj[labels[v]].append(labels[u])
    for i, v in enumerate(ports):
        adj[("p", i)].append(labels[v])
        adj[labels[v]].append(("p", i))
    check(all(len(adj[labels[v]]) == 2 for v in vertices), "two-colour degree")
    paths = []
    for start, target in ((0, 1), (2, 3)):
        path = [("p", start)]
        previous = None
        current = path[0]
        while True:
            choices = [x for x in adj[current] if x != previous]
            check(len(choices) == 1, "path continuation")
            nxt = choices[0]
            path.append(nxt)
            if nxt[0] == "p":
                check(nxt == ("p", target), "forbidden cross-pair Kempe path")
                break
            previous, current = current, nxt
        paths.append(tuple(x[1] for x in path[1:-1]))
    check(not set(paths[0]) & set(paths[1]), "wire paths intersect")
    return tuple(paths)


def tait_count(edges):
    """Independent full-graph engine: match one colour, alternate residual cycles."""
    adj = {}
    for u, v in edges:
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)
    check(all(len(a) == 3 for a in adj.values()), "full graph noncubic")

    def match(left, chosen):
        if left:
            v = min(left, key=lambda x: len(adj[x] & left))
            return sum(match(left - {u, v}, chosen | {edge(u, v)})
                       for u in sorted(adj[v] & left))
        remaining = set(adj)
        cycles = 0
        while remaining:
            start = min(remaining)
            previous, current = None, start
            length = 0
            while True:
                remaining.discard(current)
                length += 1
                nxt = [v for v in adj[current]
                       if edge(current, v) not in chosen and v != previous]
                check(bool(nxt), "residual cycle missing")
                target = min(nxt)
                if target == start:
                    break
                previous, current = current, target
            if length % 2:
                return 0
            cycles += 1
        return 2 ** cycles

    return match(set(adj), set())


def insert_petersen(edges, first, second):
    """Fixed labelled pairing: first to 1,4; second to 7,8."""
    edges = {edge(*e) for e in edges}
    first, second = tuple(first), tuple(second)
    check(edge(*first) in edges and edge(*second) in edges, "missing source edge")
    check(edge(*first) != edge(*second), "identical source edges")
    fresh = max(x for e in edges for x in e) + 1
    ren = {v: fresh + i for i, v in enumerate(v for v in range(10) if v not in (0, 5))}
    result = edges - {edge(*first), edge(*second)}
    result.update(edge(ren[u], ren[v]) for u, v in P if 0 not in (u, v) and 5 not in (u, v))
    for x, port in zip(first + second, (1, 4, 7, 8)):
        result.add(edge(x, ren[port]))
    replacements = {}
    for source, local in ((first, (1, 2, 3, 4)), (second, (7, 9, 6, 8))):
        path = (source[0],) + tuple(ren[x] for x in local) + (source[1],)
        replacements[source] = path
        replacements[source[::-1]] = path[::-1]
    return tuple(sorted(result)), replacements


def lift_paths(paths, replacements):
    return tuple(tuple([path[0]] + [x for e in zip(path, path[1:])
                                   for x in replacements.get(e, e)[1:]])
                 for path in paths)


def audit_petersen_paths(edges, paths):
    check(len(paths) == 15, "wrong path bank length")
    available = {edge(*e) for e in edges}
    branch, interiors = {}, set()
    for (u, v), path in zip(P, paths):
        check(len(path) >= 2 and len(set(path)) == len(path), "not a simple path")
        check(all(edge(a, b) in available for a, b in zip(path, path[1:])), "absent path edge")
        for label, image in ((u, path[0]), (v, path[-1])):
            check(label not in branch or branch[label] == image, "inconsistent branch")
            branch[label] = image
        check(not interiors & set(path[1:-1]), "interiors overlap")
        interiors.update(path[1:-1])
    check(len(set(branch.values())) == 10, "branch images repeat")
    check(not set(branch.values()) & interiors, "branch/interior collision")
    return len(interiors)


def graph(edges):
    result = nx.Graph()
    result.add_edges_from(edges)
    return result


def flower(k):
    result = [(4 * i, 4 * i + j) for i in range(k) for j in (1, 2, 3)]
    result += [(4 * i + 1, 4 * ((i + 1) % k) + 1) for i in range(k)]
    ring = tuple(4 * i + j for j in (2, 3) for i in range(k))
    result += list(zip(ring, ring[1:] + ring[:1]))
    return tuple(map(lambda e: edge(*e), result))


def cycle_count(g, length):
    result = set()
    for start in g:
        def walk(path):
            if len(path) == length:
                if g.has_edge(path[-1], start):
                    result.add(frozenset(edge(a, b) for a, b in zip(path, path[1:] + path[:1])))
                return
            for v in g[path[-1]]:
                if v > start and v not in path:
                    walk(path + (v,))
        walk((start,))
    return len(result)


def replace_vertex(edges, vertex, pole_vertices, pole_edges, pole_ports):
    adjacent = sorted(v if u == vertex else u for u, v in edges if vertex in (u, v))
    check(len(adjacent) == len(pole_ports) == 3, "vertex substitution interface")
    fresh = max(x for e in edges for x in e) + 1
    ren = {v: fresh + i for i, v in enumerate(pole_vertices)}
    result = {edge(*e) for e in edges if vertex not in e}
    result.update(edge(ren[u], ren[v]) for u, v in pole_edges)
    result.update(edge(a, ren[p]) for a, p in zip(adjacent, pole_ports))
    return tuple(sorted(result))


def has_small_cyclic_cut(g):
    """Exhaust every cut of size at most 3, for the two named 18-vertex banks."""
    edges = tuple(g.edges())
    for size in range(1, 4):
        for cut in combinations(edges, size):
            h = g.copy()
            h.remove_edges_from(cut)
            cyclic = sum(h.subgraph(c).number_of_edges() >= len(c)
                         for c in nx.connected_components(h))
            if cyclic >= 2:
                return True
    return False


def main():
    B = {(a, a, b, b) for a, b in product(range(3), repeat=2)}
    D = {w for w in product(range(3), repeat=3) if len(set(w)) == 3}
    examples = [
        ("Petersen", P, (0, 5), 2),
        ("BlanusaI", B1, (7, 8), 4),
        ("BlanusaII", B2, (9, 15), 6),
        ("J5_spoke", flower(5), (0, 1), 12),
        ("J5_ring", flower(5), (1, 5), 4),
        ("J3_zero", flower(3), (1, 5), 0),
    ]
    for name, edges, cut, expected in examples:
        vertices, kept, ports = deleted_edge_pole(edges, cut)
        profile, witnesses = pole_profile(vertices, kept, ports)
        check(profile == Counter({w: expected for w in B}) if expected else not profile,
              name + " uniform profile")
        paths = kempe_wires(vertices, kept, ports, witnesses[(0, 0, 1, 1)]) if expected else ()
        print(f"pole {name} cut={cut} ports={ports} m={expected} kempe_paths={paths}")

    triangle = tuple(combinations(range(3), 2))
    k23 = tuple((u, v) for u in (0, 1) for v in (2, 3, 4))
    for name, vertices, edges, ports, expected in (
            ("triangle", range(3), triangle, (0, 1, 2), 1),
            ("K2_3", range(5), k23, (2, 3, 4), 2)):
        profile, _ = pole_profile(vertices, edges, ports)
        check(profile == Counter({w: expected for w in D}), name + " three-pole")
        print(f"three_pole {name} support={len(profile)} multiplicity={expected}")

    p3edges = tuple(e for e in P if 0 not in e)
    p3, _ = pole_profile(range(1, 10), p3edges, (1, 4, 5))
    check(not p3, "uncolourable Petersen three-pole hostile")
    hostile, _ = pole_profile((0, 1), ((0, 1),), (0, 0, 1, 1))
    check(len(hostile) == 12 and not set(hostile) & B, "wire-minor converse hostile")
    check((0, 0, 1, 1) in B and (0, 1, 0, 1) not in B, "pair-order hostile")
    print("hostiles empty Petersen three-pole; two-vertex wire-minor pole support=12 disjoint from B; crossed pairing rejected")

    # Product seven-pole: a triangle plus the ordered Petersen four-pole.
    vertices, kept, ports = deleted_edge_pole(P, (0, 5))
    shift = {v: v + 3 for v in vertices}
    seven_edges = triangle + tuple((shift[u], shift[v]) for u, v in kept)
    seven_ports = (0, 1, 2) + tuple(shift[v] for v in ports)
    seven_vertices = tuple(range(3)) + tuple(shift[v] for v in vertices)
    profile, _ = pole_profile(seven_vertices, seven_edges, seven_ports)
    J = {a + b for a in D for b in B}
    check(profile == Counter({w: 2 for w in J}), "seven-port repair profile")
    check(all(all(w.count(c) % 2 == 1 for c in range(3)) for w in J), "seven-port parity")
    print(f"seven_port triangle_x_Petersen support={len(J)} multiplicity=2 total={sum(profile.values())}")

    k4 = tuple(combinations(range(4), 2))
    check(tait_count(k4) == 6 and tait_count(P) == 0, "full graph controls")
    k4pairs = [(a, b) for a, b in combinations(k4, 2) if not set(a) & set(b)]
    for a, b in k4pairs:
        inserted, _ = insert_petersen(k4, a, b)
        check(tait_count(inserted) == 12, "positive factor-two control")
    print("colourable controls all 3 independent-edge K4 insertions: 6 -> 12")

    # Complete labelled universe: unordered independent edge pair, fixed endpoint order.
    named = (graph(B1), graph(B2))
    check(all(not has_small_cyclic_cut(g) for g in named), "named bank cyclic connectivity")
    buckets = Counter()
    named_examples = {}
    for a, b in combinations(P, 2):
        if set(a) & set(b):
            continue
        inserted, repl = insert_petersen(P, a, b)
        g = graph(inserted)
        audit_petersen_paths(inserted, lift_paths(P, repl))
        check(tait_count(inserted) == 0, "P insertion colouring")
        check(nx.is_connected(g) and not list(nx.bridges(g)), "snark connectivity")
        check(all(g.degree(v) == 3 for v in g), "snark cubicity")
        check(nx.girth(g) >= 5, "snark girth")
        hits = [i for i, h in enumerate(named, 1) if nx.is_isomorphic(g, h)]
        check(len(hits) == 1, "Blanusa named-bank isomorphism")
        buckets[hits[0]] += 1
        named_examples.setdefault(hits[0], (a, b, cycle_count(g, 5)))
    check(sum(buckets.values()) == 75 and len(buckets) == 2, "75 marked insertion census")
    print(f"Petersen independent-edge insertions=75 BlanusaI={buckets[1]} BlanusaII={buckets[2]}")
    for i in sorted(named_examples):
        print(f"named example Blanusa{i}: edges_and_five_cycles={named_examples[i]}")

    # The actual seven-boundary closure: replace vertex 0 and two distinct remote edges.
    for source_name, source, expected in (("K4", k4, 12), ("Petersen", P, 0)):
        inflated = replace_vertex(source, 0, range(3), triangle, (0, 1, 2))
        remote = [e for e in source if 0 not in e]
        inflated, _ = insert_petersen(inflated, remote[0], remote[-1])
        check(tait_count(inflated) == expected, "seven-boundary whole graph count")
        g = graph(inflated)
        check(nx.is_connected(g) and not list(nx.bridges(g)), "seven-boundary bridgeless")
        print(f"seven_port closure {source_name}: vertices={len(g)} Tait={expected}")

    # An explicit unbounded-depth mechanism, checked for 100 replacements.
    edges, paths = P, P
    for depth in range(1, 101):
        first = edges[(depth * 7) % len(edges)]
        second = next(e for e in reversed(edges) if not set(e) & set(first))
        edges, repl = insert_petersen(edges, first, second)
        paths = lift_paths(paths, repl)
        audit_petersen_paths(edges, paths)
    check(len(set(v for e in edges for v in e)) == 810, "iterated order")
    check(len(edges) == 1215, "iterated size")
    print("iterated Petersen substitution depth=100 vertices=810 edges=1215; all 100 subdivision certificates pass")
    print("PASS exact supports, independent whole-graph counts, labelled reductions, named-bank isomorphisms, and hostiles")


if __name__ == "__main__":
    main()
