"""Exact cuff/seam holonomy and connected bipartite-fragment gluing.

Surface topology is proved in the companion note. The program checks its
finite F2 presentation and never infers plane colourability from surface genus.
"""
from dataclasses import dataclass, replace
from itertools import combinations, product
from math import comb
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def nat(x):
    need(type(x) is int and x >= 0, "exact nonnegative integer required")
    return x


def bits(value, length):
    need(type(value) is tuple and len(value) == length
         and all(type(x) is int and x in (0, 1) for x in value),
         "exact binary tuple required")
    return value


@dataclass(frozen=True)
class Graph:
    vertices: int
    edges: tuple


def audit_graph(g, cubic=False):
    need(type(g) is Graph and type(g.vertices) is int and g.vertices >= 1,
         "nonempty exact graph required")
    n = g.vertices
    need(type(g.edges) is tuple and
         all(type(e) is tuple and len(e) == 2
             and all(type(v) is int and 0 <= v < n for v in e)
             and e[0] <= e[1] for e in g.edges), "ordered indexed edges required")
    # Indexed edges permit parallel seams and loops; loops have two cuff ends.
    reached = {0}
    while True:
        out = reached | {v for u, v in g.edges if u in reached} | {
            u for u, v in g.edges if v in reached}
        if out == reached:
            break
        reached = out
    need(len(reached) == n, "connected graph required")
    if cubic:
        degree = [0] * n
        for u, v in g.edges:
            degree[u] += 1
            degree[v] += 1
        need(all(d == 3 for d in degree), "three cuff ends at each pants vertex required")
    return g


def tree_data(g):
    audit_graph(g)
    parent, edge_to_parent = [-1] * g.vertices, [-1] * g.vertices
    parent[0] = 0
    order, tree = [0], []
    for u in order:
        for i, (a, b) in enumerate(g.edges):
            if u not in (a, b):
                continue
            v = b if u == a else a
            if parent[v] == -1:
                parent[v], edge_to_parent[v] = u, i
                tree.append(i)
                order.append(v)
    chords = tuple(i for i in range(len(g.edges)) if i not in tree)
    return tuple(tree), chords, tuple(parent), tuple(edge_to_parent), tuple(order)


def gauge(g, transitions, vertex_flips):
    audit_graph(g)
    bits(transitions, len(g.edges))
    bits(vertex_flips, g.vertices)
    return tuple(t ^ vertex_flips[u] ^ vertex_flips[v]
                 for t, (u, v) in zip(transitions, g.edges))


def normalize(g, transitions):
    """Unique tree-zero representative, and its marked root-zero gauge."""
    bits(transitions, len(audit_graph(g).edges))
    tree, chords, parent, parent_edge, order = tree_data(g)
    flips = [0] * g.vertices
    for v in order[1:]:
        flips[v] = flips[parent[v]] ^ transitions[parent_edge[v]]
    normal = gauge(g, transitions, tuple(flips))
    need(all(normal[i] == 0 for i in tree), "tree normalization")
    return normal, tuple(flips), tuple(normal[i] for i in chords)


def incidence(g, values):
    audit_graph(g)
    bits(values, len(g.edges))
    out = [0] * g.vertices
    for value, (u, v) in zip(values, g.edges):
        out[u] ^= value
        out[v] ^= value
    return tuple(out)


def cycle_from_chords(g, values):
    tree, chords, parent, parent_edge, order = tree_data(g)
    bits(values, len(chords))
    result = [0] * len(g.edges)
    for i, value in zip(chords, values):
        result[i] = value
    residual = list(incidence(g, tuple(result)))
    for v in reversed(order[1:]):
        value = residual[v]
        result[parent_edge[v]] = value
        residual[v] ^= value
        residual[parent[v]] ^= value
    need(not any(residual), "the total parity residual vanishes")
    return tuple(result)


@dataclass(frozen=True)
class PantsCode:
    cuff_chords: tuple
    seam_chords: tuple


def encode_pants(g, cuffs, transitions):
    audit_graph(g, cubic=True)
    bits(cuffs, len(g.edges))
    need(not any(incidence(g, cuffs)), "each pair of pants has even cuff sum")
    _, chords, *_ = tree_data(g)
    seam = normalize(g, transitions)[2]
    return PantsCode(tuple(cuffs[i] for i in chords), seam)


def decode_pants(g, code):
    audit_graph(g, cubic=True)
    need(type(code) is PantsCode, "typed pants code required")
    _, chords, *_ = tree_data(g)
    bits(code.cuff_chords, len(chords))
    bits(code.seam_chords, len(chords))
    cuffs = cycle_from_chords(g, code.cuff_chords)
    transitions = [0] * len(g.edges)
    for i, value in zip(chords, code.seam_chords):
        transitions[i] = value
    return cuffs, tuple(transitions)


def seam_endpoints(multiplicities):
    need(type(multiplicities) is tuple and len(multiplicities) == 3,
         "three fixed seam multiplicities required")
    x, y, z = map(nat, multiplicities)
    return x + y, x + z, y + z


def recover_seams(endpoints):
    need(type(endpoints) is tuple and len(endpoints) == 3,
         "three endpoint counts required")
    a, b, c = map(nat, endpoints)
    numerators = a + b - c, a + c - b, b + c - a
    need(all(v >= 0 and v % 2 == 0 for v in numerators),
         "nonnegative triangle inequalities and even total required")
    return tuple(v // 2 for v in numerators)


def cuts(multiplicities):
    seam_endpoints(multiplicities)
    x, y, z = multiplicities
    return x + 1, x + y + 2


def recover_cuts(total, pair):
    nat(total)
    need(type(pair) is tuple and len(pair) == 2
         and all(type(v) is int for v in pair)
         and 1 <= pair[0] < pair[1] <= total + 2, "two marked cut positions required")
    i, j = pair
    return i - 1, j - i - 1, total + 2 - j


@dataclass(frozen=True)
class Fragment:
    vertices: tuple
    edges: frozenset
    colours: tuple


def audit_fragment(f):
    need(type(f) is Fragment and type(f.vertices) is tuple and f.vertices
         and all(type(v) is int and v >= 0 for v in f.vertices)
         and tuple(sorted(set(f.vertices))) == f.vertices, "marked fragment vertices required")
    need(type(f.edges) is frozenset and
         all(type(e) is tuple and len(e) == 2 and
             all(type(v) is int and v in f.vertices for v in e)
             and e[0] < e[1] for e in f.edges), "simple fragment edges required")
    bits(f.colours, len(f.vertices))
    values = dict(zip(f.vertices, f.colours))
    need(all(values[u] != values[v] for u, v in f.edges), "actual proper local colouring required")
    renaming = {v: i for i, v in enumerate(f.vertices)}
    audit_graph(Graph(len(f.vertices), tuple((renaming[u], renaming[v]) for u, v in sorted(f.edges))))
    return values


def glue_fragments(fragments):
    """One actual global lift, or None for a proved binary-colouring obstruction.

    Declared edges are exactly the union; unlisted geometric cross-edges are
    not silently added. Connected overlap graph is the stated input domain.
    """
    need(type(fragments) is tuple and fragments, "nonempty tuple of connected fragments required")
    local = tuple(audit_fragment(f) for f in fragments)
    seams, transitions = [], []
    for i, j in combinations(range(len(local)), 2):
        shared = sorted(local[i].keys() & local[j].keys())
        if not shared:
            continue
        values = {local[i][v] ^ local[j][v] for v in shared}
        if len(values) != 1:
            return None
        seams.append((i, j))
        transitions.append(values.pop())
    nerve = Graph(len(fragments), tuple(seams))
    normal, flips, _ = normalize(nerve, tuple(transitions))
    if any(normal):
        return None
    global_values = {}
    for values, flip in zip(local, flips):
        for v, colour in values.items():
            value = colour ^ flip
            need(v not in global_values or global_values[v] == value, "overlap agreement")
            global_values[v] = value
    edges = frozenset(e for f in fragments for e in f.edges)
    need(all(global_values[u] != global_values[v] for u, v in edges), "actual union colouring")
    return tuple(sorted(global_values.items()))


def all_bits(n):
    return product((0, 1), repeat=n)


def main():
    checks = 0
    def check(ok, label):
        nonlocal checks
        need(ok, label)
        checks += 1
    def rejects(fn):
        nonlocal checks
        try:
            fn()
        except (ValueError, TypeError):
            checks += 1
            return
        raise ValueError("hostile input accepted")

    examples = {
        "theta": Graph(2, ((0, 1),) * 3),
        "separating_cuff": Graph(2, ((0, 0), (0, 1), (1, 1))),
        "tetrahedral_dual": Graph(4, tuple(combinations(range(4), 2))),
        "K33_dual": Graph(6, tuple((i, j) for i in range(3) for j in range(3, 6))),
        "cube_dual": Graph(8, tuple((i, j) for i in range(8) for j in range(i + 1, 8)
                                  if (i ^ j).bit_count() == 1)),
    }
    summary = {}
    for name, graph in examples.items():
        audit_graph(graph, cubic=True)
        cycle_rank = len(graph.edges) - graph.vertices + 1
        cycles = [c for c in all_bits(len(graph.edges)) if not any(incidence(graph, c))]
        check(len(cycles) == 2 ** cycle_rank, "independent cuff enumeration")
        _, chords, *_ = tree_data(graph)
        for c in cycles:
            check(cycle_from_chords(graph, tuple(c[i] for i in chords)) == c,
                  "lossless cycle coordinates")
        normals = {}
        for transition in all_bits(len(graph.edges)):
            normal, flips, code = normalize(graph, transition)
            check(gauge(graph, transition, flips) == normal, "explicit gluing gauge")
            normals[normal] = normals.get(normal, 0) + 1
            check(normalize(graph, normal)[0] == normal, "canonical representative")
        check(len(normals) == 2 ** cycle_rank, "all seam classes")
        check(set(normals.values()) == {2 ** (graph.vertices - 1)}, "constant gauge orbit size")
        for c in cycles:
            for t in normals:
                encoded = encode_pants(graph, c, t)
                check(decode_pants(graph, encoded) == (c, t), "complete marked pants codec")
        # Independently compare gauge orbits on the three smallest examples.
        if graph.vertices <= 4:
            for t in all_bits(len(graph.edges)):
                orbit = {gauge(graph, t, v) for v in all_bits(graph.vertices)}
                same_code = {u for u in all_bits(len(graph.edges))
                             if normalize(graph, u)[2] == normalize(graph, t)[2]}
                check(orbit == same_code, "independent complete gauge orbit")
        summary[name] = {"pants": graph.vertices, "seams": len(graph.edges),
                         "genus": cycle_rank, "cuff_states": len(cycles),
                         "seam_states_per_cuff": len(normals),
                         "surface_characters": len(cycles) * len(normals)}
    theta = examples["theta"]
    trivial, twisted = (0, 0, 0), (0, 0, 1)
    check(encode_pants(theta, trivial, trivial) != encode_pants(theta, trivial, twisted),
          "trivial cuffs do not imply trivial global character")
    check(not any(normalize(theta, trivial)[0]) and any(normalize(theta, twisted)[0]),
          "smallest closed pants holonomy hostile")
    check(all(c[1] == 0 for c in all_bits(3)
              if not any(incidence(examples["separating_cuff"], c))),
          "a separating cuff has trivial binary character")

    for total in range(17):
        universe = [(x, y, total - x - y) for x in range(total + 1)
                    for y in range(total - x + 1)]
        check(len(universe) == comb(total + 2, 2), "triangular seam layer")
        for triple in universe:
            end = seam_endpoints(triple)
            check(recover_seams(end) == triple, "integer endpoint decoder")
            check(sum(end) == 2 * total, "endpoint budget")
            check(recover_cuts(total, cuts(triple)) == triple, "marked pair-index bijection")
        check({cuts(v) for v in universe} == set(combinations(range(1, total + 3), 2)),
              "complete stars and bars pair range")
    rejects(lambda: recover_seams((2, 0, 0)))  # Returning wave is outside the fixed seam bands.
    rejects(lambda: recover_seams((1, 1, 1)))  # Triangle inequalities alone omit index-two parity.

    # Independent signed-edge equations for every connected simple graph n<=4.
    signed_instances = 0
    for n in range(1, 5):
        possible = tuple(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            graph = Graph(n, tuple(e for i, e in enumerate(possible) if mask >> i & 1))
            try:
                audit_graph(graph)
            except ValueError:
                continue
            for transition in all_bits(len(graph.edges)):
                solutions = [s for s in all_bits(n)
                             if all((s[u] ^ s[v]) == a
                                    for (u, v), a in zip(graph.edges, transition))]
                normal, flips, _ = normalize(graph, transition)
                check(len(solutions) == (0 if any(normal) else 2),
                      "cycle holonomy iff literal binary solutions")
                if solutions:
                    check(flips in solutions, "decoded solution is literal")
                signed_instances += 1

    triangle = (Fragment((0, 1), frozenset({(0, 1)}), (0, 1)),
                Fragment((1, 2), frozenset({(1, 2)}), (0, 1)),
                Fragment((0, 2), frozenset({(0, 2)}), (0, 1)))
    check(glue_fragments(triangle) is None, "locally bipartite fragments, global triangle")
    square = triangle[:2] + (
        Fragment((2, 3), frozenset({(2, 3)}), (0, 1)),
        Fragment((0, 3), frozenset({(0, 3)}), (0, 1)))
    for flips in all_bits(4):
        shifted = tuple(replace(f, colours=tuple(c ^ s for c in f.colours))
                        for f, s in zip(square, flips))
        global_word = dict(glue_fragments(shifted))
        check(all(global_word[u] != global_word[v] for f in square for u, v in f.edges),
              "local gauge choices preserve a concrete union lift")
    unequal_ports = (
        Fragment((0, 1), frozenset({(0, 1)}), (0, 1)),
        Fragment((0, 1, 2), frozenset({(0, 2), (1, 2)}), (0, 0, 1)))
    check(glue_fragments(unequal_ports) is None, "multiple shared ports can already disagree")

    rejects(lambda: normalize(theta, (False, 0, 0)))
    rejects(lambda: gauge(theta, trivial, (0, 1.0)))
    rejects(lambda: audit_graph(Graph(2, ((False, 1),))))
    rejects(lambda: audit_graph(Graph(2, ())))
    rejects(lambda: encode_pants(Graph(3, ((0, 1), (1, 2))), (0, 0), (0, 0)))
    rejects(lambda: encode_pants(theta, (1, 0, 0), trivial))
    rejects(lambda: decode_pants(theta, PantsCode((0,), (0, 0))))
    rejects(lambda: recover_cuts(3, (True, 2)))
    rejects(lambda: seam_endpoints((0, -1, 1)))
    rejects(lambda: glue_fragments((replace(triangle[0], colours=(0, 0)),)))
    rejects(lambda: glue_fragments((replace(triangle[0], vertices=(False, 1)),)))
    rejects(lambda: glue_fragments((triangle[0], Fragment((3,), frozenset(), (0,)))))
    print(json.dumps({
        "scope": "marked F2 pants holonomy, fixed seam bands, and actual binary-fragment lifts",
        "checks": checks, "pants_controls": summary,
        "seam_layers": "N=0..16", "signed_connected_graph_instances": signed_instances,
        "hostiles": ["zero cuffs with nontrivial seam holonomy", "returning same-cuff wave",
                     "odd endpoint sum", "bipartite fragments forming a triangle",
                     "incompatible multiple ports", "typed malformed records"],
        "plane_colour_consequence": "none without actual colouring or noncolourability certificates",
    }, indent=2))


if __name__ == "__main__":
    main()
