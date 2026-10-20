#!/usr/bin/env python3
"""Labelled-support and rooted-minor audits beyond flower snarks.

Run from the repository root. Python standard library plus NetworkX 3.5.
All tests use explicit checks and therefore survive python -O.
"""

from collections import Counter
from itertools import combinations, product
from pathlib import Path
import runpy
import networkx as nx


CHECKS = 0
ROOT = Path(__file__).resolve().parents[1]
OLD = runpy.run_path(str(ROOT / "04-computation" /
                        "boundary_support_20261020_snarks.py"))
PETERSEN = set(OLD["PETERSEN"])


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def edge(a, b):
    return tuple(sorted((a, b)))


def graph(edges):
    edges = set(edge(*e) for e in edges)
    check(all(a != b for a, b in edges), "loop")
    return edges


def adjacency(edges):
    result = {}
    for a, b in edges:
        result.setdefault(a, set()).add(b)
        result.setdefault(b, set()).add(a)
    return result


def pole_check(edges, k):
    adj = adjacency(edges)
    check({v for v in adj if v < 0} == set(range(-k, 0)), "port labels")
    check(all(len(n) == (1 if v < 0 else 3) for v, n in adj.items()),
          "pole degrees")
    return adj


def pole_colourings(edges, k, on_colouring=None):
    """Enumerate all labelled extensions; also retain one witness per word."""
    adj = pole_check(edges, k)
    ordered = tuple(sorted(edges))
    used = {v: 0 for v in adj if v >= 0}
    assigned = [-1]*len(ordered)
    boundary = tuple(next(i for i, e in enumerate(ordered) if -p-1 in e)
                     for p in range(k))
    counts, witnesses = Counter(), {}

    def visit(left):
        if not left:
            word = tuple(assigned[i] for i in boundary)
            counts[word] += 1
            witnesses.setdefault(word, dict(zip(ordered, assigned)))
            if on_colouring is not None:
                on_colouring(word, dict(zip(ordered, assigned)))
            return
        best, options, size = None, 0, 4
        for i in left:
            mask = 0
            for v in ordered[i]:
                if v >= 0:
                    mask |= used[v]
            possible = 7 & ~mask
            if possible.bit_count() < size:
                best, options, size = i, possible, possible.bit_count()
                if size <= 1:
                    break
        if not options:
            return
        rest = left-{best}
        while options:
            bit = options & -options
            options ^= bit
            assigned[best] = bit.bit_length()-1
            for v in ordered[best]:
                if v >= 0:
                    used[v] |= bit
            visit(rest)
            for v in ordered[best]:
                if v >= 0:
                    used[v] ^= bit
        assigned[best] = -1

    visit(frozenset(range(len(ordered))))
    return counts, witnesses


def deleted_edge_pole(edges, removed):
    adj = adjacency(edges)
    u, v = removed
    endpoints = sorted(adj[u]-{v})+sorted(adj[v]-{u})
    result = {e for e in edges if u not in e and v not in e}
    result.update(edge(-i-1, w) for i, w in enumerate(endpoints))
    pole_check(result, 4)
    return result


def opened_edges_pole(edges, opened):
    check(len(set(sum((list(e) for e in opened), []))) == 4,
          "independent opened edges required")
    result = set(edges)-set(opened)
    for i, v in enumerate(sum((list(e) for e in opened), [])):
        result.add(edge(-i-1, v))
    pole_check(result, 4)
    return result


ISO = {w for w in product(range(3), repeat=4)
       if w[0] == w[1] and w[2] == w[3]}
HETERO = {w for w in product(range(3), repeat=4)
          if w[0] != w[1] and w[2] != w[3]
          and set(w[:2]) == set(w[2:])}
P8 = deleted_edge_pole(PETERSEN, (0, 1))  # Eight vertices, FOUR ports.


def kempe_pair_paths(pole, witness):
    selected = nx.Graph()
    selected.add_edges_from(e for e, c in witness.items() if c in (0, 1))
    paths = (tuple(nx.shortest_path(selected, -1, -2)),
             tuple(nx.shortest_path(selected, -3, -4)))
    check(not set(paths[0]) & set(paths[1]), "Kempe paths intersect")
    check(all(edge(a, b) in pole for p in paths for a, b in zip(p, p[1:])),
          "missing Kempe path edge")
    return paths


def suppress_interior(edges, keep=()):
    work = set(edges)
    keep = set(keep)
    while True:
        adj = adjacency(work)
        candidates = [v for v in adj if v >= 0 and v not in keep and len(adj[v]) == 2]
        if not candidates:
            return work
        v = min(candidates)
        a, b = sorted(adj[v])
        work -= {edge(v, a), edge(v, b)}
        check(a != b and edge(a, b) not in work, "invalid simple suppression")
        work.add(edge(a, b))


def iso_minor(pole, witness):
    paths = kempe_pair_paths(pole, witness)
    retained = {edge(a, b) for p in paths for a, b in zip(p, p[1:])}
    reduced = suppress_interior(retained)
    check(reduced == {(-2, -1), (-4, -3)}, "wrong two-wire minor")
    return paths


def hetero_minor(pole, witness):
    first, second = kempe_pair_paths(pole, witness)
    internal = nx.Graph()
    internal.add_edges_from(e for e in pole if min(e) >= 0)
    a_set, b_set = set(first)-{-1, -2}, set(second)-{-3, -4}
    candidates = []
    for a in a_set:
        for b in b_set:
            if nx.has_path(internal, a, b):
                candidates.append(tuple(nx.shortest_path(internal, a, b)))
    check(bool(candidates), "no connector between Kempe paths")
    connector = min(candidates, key=lambda p: (len(p), p))
    check(not set(connector[1:-1]) & (a_set | b_set), "connector intersects paths")
    retained = {edge(a, b) for p in (first, second, connector)
                for a, b in zip(p, p[1:])}
    reduced = suppress_interior(retained, (connector[0], connector[-1]))
    counts, _ = pole_colourings(reduced, 4)
    check(set(counts) == HETERO and set(counts.values()) == {1},
          "wrong adjacent-junction minor")
    check(sum(v >= 0 for v in adjacency(reduced)) == 2, "not two junctions")
    return first, second, connector


def insert_pole(base, opened, replacement=P8):
    """Substitute a four-pole for two labelled edges, including dangling wires."""
    check(len(opened) == 2 and set(opened) <= set(base), "missing insertion edges")
    base_vertices = set(adjacency(base))
    new_vertices = sorted(v for v in adjacency(replacement) if v >= 0)
    offset = max(base_vertices | {0})+1
    rename = {v: offset+i for i, v in enumerate(new_vertices)}
    attachment = {-(i+1): v for i, v in enumerate(
        [*opened[0], *opened[1]])}
    def renamed(v):
        return rename[v] if v >= 0 else attachment[v]
    new_edges = {edge(renamed(a), renamed(b)) for a, b in replacement}
    result = set(base)-set(opened)
    check(not result & new_edges, "insertion creates parallel edge")
    result |= new_edges
    return result, rename, attachment


def check_insertion_minor(base, opened, larger, replacement, paths, rename, attachment):
    def renamed(v):
        return rename[v] if v >= 0 else attachment[v]
    lifted = tuple(tuple(renamed(v) for v in p) for p in paths)
    check(all(edge(a, b) in larger for p in lifted for a, b in zip(p, p[1:])),
          "lifted minor path missing")
    check(not set(lifted[0][1:-1]) & set(lifted[1][1:-1]), "lifted path collision")
    retained = (set(base)-set(opened)) | {
        edge(a, b) for p in lifted for a, b in zip(p, p[1:])}
    reduced = suppress_interior(retained, keep=set(adjacency(base)))
    check(reduced == set(base), "insertion failed to reduce to exact base")


def girth(g):
    best = len(g)+1
    for start in g:
        distance, parent = {start: 0}, {start: None}
        todo = [start]
        for u in todo:
            for v in g[u]:
                if v not in distance:
                    distance[v] = distance[u]+1
                    parent[v] = u
                    todo.append(v)
                elif parent[u] != v and parent[v] != u:
                    best = min(best, distance[u]+distance[v]+1)
    return best


def cyclic_small_cut_audit(edges):
    """Every cut of size <=3 is covered by a bridge after <=2 deletions."""
    ordered = tuple(sorted(edges))
    found = set()
    for size in range(3):
        for removed in combinations(ordered, size):
            g = nx.Graph()
            g.add_edges_from(edges)
            g.remove_edges_from(removed)
            for e in list(nx.bridges(g)):
                h = g.copy()
                h.remove_edge(*e)
                cyclic = sum(h.subgraph(c).number_of_edges() >= len(c)
                             for c in nx.connected_components(h))
                if cyclic >= 2:
                    found.add(tuple(sorted((*removed, edge(*e)))))
    check(not found, f"cyclic cut of size <=3: {found}")


def family_audit(p8_paths):
    current = set(PETERSEN)
    history = []
    for step in range(1, 9):
        ordered = sorted(current)
        opened = next((a, b) for a, b in combinations(ordered, 2)
                      if not set(a) & set(b))
        larger, rename, attachment = insert_pole(current, opened)
        check_insertion_minor(current, opened, larger, P8, p8_paths, rename, attachment)
        g = nx.Graph()
        g.add_edges_from(larger)
        check(len(g) == 10+8*step and set(dict(g.degree()).values()) == {3},
              "dot family size or degree")
        check(nx.is_connected(g) and not list(nx.bridges(g)), "dot family bridge")
        check(girth(g) == 5, "dot family girth")
        count, perfect_matchings = OLD["matching_colour_count"](larger)
        check(count == 0, "independent dot-family colouring")
        cyclic_small_cut_audit(larger)
        print(f"dot family step={step}, n={len(g)}, m={len(larger)}, girth=5,"
              f" bridgefree, exact predecessor minor, PM={perfect_matchings}")
        history.append((current, opened, larger, rename, attachment))
        current = larger
    # Reverse every authenticated step all the way to the original Petersen.
    for base, opened, larger, rename, attachment in reversed(history):
        check(current == larger, "minor chain state")
        check_insertion_minor(base, opened, current, P8, p8_paths, rename, attachment)
        current = base
    check(current == PETERSEN, "minor chain Petersen endpoint")
    current = set(combinations(range(4), 2))
    for step in range(1, 4):
        opened = next((a, b) for a, b in combinations(sorted(current), 2)
                      if not set(a) & set(b))
        current, _, _ = insert_pole(current, opened)
        count, _ = OLD["matching_colour_count"](current)
        check(count == 6*2**step, "colourable count scaling control")
        print(f"K4 positive insertion control step={step}, C={count}")


def seven_pole_audit(p8_paths):
    junction = {(-1, 0), (-2, 0), (-3, 0), (-5, -4), (-7, -6)}
    counts, _ = pole_colourings(junction, 7)
    expected = {w for w in product(range(3), repeat=7)
                if len(set(w[:3])) == 3 and w[3] == w[4] and w[5] == w[6]}
    check(set(counts) == expected and len(expected) == 54, "junction kernel")
    first_opened = ((-1, 0), (-5, -4))
    first, r1, a1 = insert_pole(junction, first_opened)
    second_opened = ((-2, 0), (-7, -6))
    connected, r2, a2 = insert_pole(first, second_opened)
    counts2, _ = pole_colourings(connected, 7)
    check(counts2 == Counter({w: 4 for w in expected}),
          "connected seven-pole is not exactly 4J")
    g = nx.Graph()
    g.add_edges_from(connected)
    check(nx.is_connected(g) and sum(v >= 0 for v in g) == 17,
          "seven-pole not connected order17")
    check_insertion_minor(first, second_opened, connected, P8, p8_paths, r2, a2)
    check_insertion_minor(junction, first_opened, first, P8, p8_paths, r1, a1)
    print("connected seven-pole: internal_vertices=17, support=54, mass=216,"
          " exact vector=4J, two certified rooted reductions to junction+2 wires")
    tree_cases = tree_reductions = 0
    for tree in nx.nonisomorphic_trees(5):
        if max(dict(tree.degree()).values()) > 3:
            continue
        tree_cases += 1
        pole = set(edge(a, b) for a, b in tree.edges())
        port = 0
        for v in sorted(tree):
            for _ in range(3-tree.degree(v)):
                port += 1
                pole.add((-port, v))
        original, _ = pole_colourings(pole, 7)
        check(len(original) == 96 and set(original.values()) == {1},
              "tree seven-pole support")
        for a, b in combinations(sorted(edge(*e) for e in tree.edges()), 2):
            if set(a) & set(b):
                continue
            reduced = suppress_interior(pole-{a, b})
            response, _ = pole_colourings(reduced, 7)
            check(len(response) == 54, "tree smoothing not repaired junction type")
            check(not set(original) & set(response), "tree-charge hostile failed")
            tree_reductions += 1
    print(f"tree seven-pole hostile: shapes={tree_cases}, reductions={tree_reductions},"
          " original support96, every smoothed support54 disjoint")
    deleted_triples = Counter()
    for removed in combinations(range(10), 3):
        removed = set(removed)
        if sum(set(e) <= removed for e in PETERSEN) != 1:
            continue
        cut = [e for e in sorted(PETERSEN)
               if (e[0] in removed) != (e[1] in removed)]
        pole = {e for e in PETERSEN if not set(e) & removed}
        for i, e in enumerate(cut):
            pole.add(edge(-i-1, next(v for v in e if v not in removed)))
        response, _ = pole_colourings(pole, 7)
        deleted_triples[len(response), sum(response.values())] += 1
    check(deleted_triples == {(120, 120): 60}, "Petersen seven-pole census")
    print("Petersen seven-poles from triples inducing one edge: 60/60 support=mass=120,"
          " none equivalent to a 54-word junction kernel")
    # Pairing the ports of the two-vertex four-pole within its couples
    # creates a loop at each vertex. This is a genuine type failure for
    # ordinary simple bridgeless reduction, despite perfect support transport.
    coupling_endpoints = (0, 0, 1, 1)
    closure = {(0, 1),
               edge(coupling_endpoints[0], coupling_endpoints[1]),
               edge(coupling_endpoints[2], coupling_endpoints[3])}
    check(closure == {(0, 0), (0, 1), (1, 1)}, "dumbbell type hostile")
    print("type hostile: P8 -> two wires, closed by adjacent junctions,"
          " produces two loops plus bridge; reject as simple-snark reduction")


def kempe_defect_audit():
    """Independent physical enumeration of the parent's seven singleton defects."""
    junction = {(-1, 0), (-2, 0), (-3, 0), (-5, -4), (-7, -6)}
    first, _, _ = insert_pole(junction, ((-1, 0), (-5, -4)))
    connected, _, _ = insert_pole(first, ((-2, 0), (-7, -6)))
    tree = {(0, 1), (1, 2), (2, 3), (3, 4)}
    tree.update((-i-1, v) for i, v in enumerate((0, 0, 1, 2, 3, 4, 4)))
    removed = {0, 1, 3}
    cut = [e for e in sorted(PETERSEN) if (e[0] in removed) != (e[1] in removed)]
    p7 = {e for e in PETERSEN if not set(e) & removed}
    for i, e in enumerate(cut):
        p7.add(edge(-i-1, next(v for v in e if v not in removed)))
    for name, pole in (("connected_17_vertex_4J", connected),
                       ("five_vertex_path", tree),
                       ("Petersen_minus_0_1_3", p7)):
        actual_n = [0]*7
        def visit(word, coloring):
            for i, unique in enumerate(word):
                if word.count(unique) != 1:
                    continue
                h = nx.Graph()
                h.add_edges_from(e for e, color in coloring.items() if color != unique)
                components = [sorted(v for v in component if v < 0)
                              for component in nx.connected_components(h)
                              if any(v < 0 for v in component)]
                check(len(components) == 3 and all(len(c) == 2 for c in components),
                      "not three boundary Kempe paths")
                if all(word[-a-1] != word[-b-1] for a, b in components):
                    actual_n[i] += 1
        counts, _ = pole_colourings(pole, 7, visit)
        physical_dots = [sum(value * (1 if sorted(Counter(word).values()) == [1, 1, 5]
                                     else -1)
                             for word, value in counts.items()
                             if word.count(word[i]) == 1)
                         for i in range(7)]
        check(physical_dots == [-n for n in actual_n], "physical singleton/Kempe identity")
        check(all(d % 6 == 0 for d in physical_dots), "colour-orbit normalization")
        if name == "connected_17_vertex_4J":
            check(actual_n == [0]*7, "4J should be on Kempe equality boundary")
        else:
            check(any(actual_n), "negative-coordinate control absent")
        print(f"seven-port Kempe identity {name}: orbit_dots="
              f"{tuple(d//6 for d in physical_dots)}, N={tuple(actual_n)}")


def main():
    print(f"NetworkX {nx.__version__}; exact labelled supports, no random graph samples")
    counts, witnesses = pole_colourings(P8, 4)
    check(counts == Counter({w: 2 for w in ISO}), "Petersen eight-vertex four-pole profile")
    p8_paths = iso_minor(P8, witnesses[(0, 0, 1, 1)])
    print(f"Petersen eight-vertex four-pole: exact vector=2B, rooted wire paths={p8_paths}")
    inherited = runpy.run_path(str(ROOT / "04-computation" /
                                   "snark_boundary_fourpole_thm4116_independent.py"))
    for name, vertices, edges in inherited["graph_definitions"]():
        if name not in ("BlanusaI", "BlanusaII", "J5"):
            continue
        marked = min(edges)
        pole = deleted_edge_pole(set(edges), marked)
        response, witness = pole_colourings(pole, 4)
        check(set(response) == ISO and len(set(response.values())) == 1,
              "inherited isochromatic pole")
        paths = iso_minor(pole, witness[(0, 0, 1, 1)])
        print(f"{name} deleted {marked}: support=9, multiplicity={next(iter(response.values()))},"
              f" paired Kempe minor path lengths={tuple(len(p)-1 for p in paths)}")
    hetero_cases = 0
    for opened in combinations(sorted(PETERSEN), 2):
        if set(opened[0]) & set(opened[1]):
            continue
        pole = opened_edges_pole(PETERSEN, opened)
        response, witnesses = pole_colourings(pole, 4)
        check(set(response) == HETERO, "Petersen opened-edges profile")
        hetero_minor(pole, witnesses[(0, 1, 0, 1)])
        hetero_cases += 1
    check(hetero_cases == 75, "independent-edge census")
    print("Petersen opened-edge heterochromatic census: 75/75 full support and rooted junction minors")
    family_audit(p8_paths)
    seven_pole_audit(p8_paths)
    kempe_defect_audit()
    print(f"PASS: {CHECKS} exact checks")


if __name__ == "__main__":
    main()
