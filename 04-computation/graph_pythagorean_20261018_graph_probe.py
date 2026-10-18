#!/usr/bin/env python3
"""Exact small graph/interface audit; no floating-point decision is used.

Run: python 04-computation/graph_pythagorean_20261018_graph_probe.py
The graphs are explicit NetworkX constructors, plus THM-4116's flower edges.
No filters beyond the stated named graph universe are applied.
"""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
import importlib.util
from pathlib import Path
import networkx as nx
import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def inherited_module(filename):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


legacy = inherited_module("snark_boundary_ap_shell_thm4116.py")
independent = inherited_module("snark_boundary_fourpole_thm4116_independent.py")


def edge_list(graph):
    return tuple(sorted(tuple(sorted(e)) for e in graph.edges()))


def matching_audit(graph):
    """Count edge colorings independently through perfect matching + 2-factor."""
    vertices = frozenset(graph)

    @lru_cache(None)
    def matchings(remaining):
        if not remaining:
            return ((),)
        v = min(remaining)
        answer = []
        for w in sorted(set(graph[v]) & remaining):
            for rest in matchings(remaining - {v, w}):
                answer.append(((v, w),) + rest)
        return tuple(answer)

    matchings_list = matchings(vertices)
    count = 0
    oddness = len(vertices)
    cycle_hist = Counter()
    for matching in matchings_list:
        residual = graph.copy()
        residual.remove_edges_from(matching)
        sizes = tuple(sorted(map(len, nx.connected_components(residual))))
        require(all(d == 2 for _, d in residual.degree()), "not a 2-factor")
        odd = sum(size % 2 for size in sizes)
        oddness = min(oddness, odd)
        cycle_hist[sizes] += 1
        if odd == 0:
            count += 2 ** len(sizes)
    return count, oddness, len(matchings_list), dict(sorted(cycle_hist.items()))


def interface(graph, side):
    edges = edge_list(graph)
    side = set(side)
    boundary = tuple(e for e in edges if (e[0] in side) != (e[1] in side))
    left = legacy.boundary_color_counts(graph, edges, side, boundary)
    right = legacy.boundary_color_counts(graph, edges, set(graph) - side, boundary)
    keys = set(left) | set(right)
    u = sum(v * v for v in left.values())
    v = sum(v * v for v in right.values())
    c = sum(left[k] * right[k] for k in keys)
    minus2 = sum((left[k] - right[k]) ** 2 for k in keys)
    plus2 = sum((left[k] + right[k]) ** 2 for k in keys)
    require(minus2 == u + v - 2 * c, "minus identity")
    require(plus2 == u + v + 2 * c, "plus identity")
    require(minus2 + plus2 == 2 * (u + v), "Pythagorean identity")
    require(plus2 - minus2 == 4 * c, "coloring gap")
    return {"side": tuple(sorted(side)), "boundary": boundary,
            "supports": (len(left), len(right)),
            "masses": (sum(left.values()), sum(right.values())),
            "norms_squared": (u, v), "dot": c,
            "triangle_squared": (minus2, plus2, 2 * (u + v)),
            "cosine_squared": str(Fraction(c * c, u * v)) if u * v else None}


def maximum_cut(graph):
    """Exhaust every cut, modulo global complement, using a Gray traversal."""
    n = len(graph)
    require(set(graph) == set(range(n)), "label convention")
    masks = [sum(1 << w for w in graph[v]) for v in range(n)]
    mask, cut, best, best_mask = 0, 0, 0, 0
    for i in range(1, 1 << (n - 1)):
        bit = (i & -i).bit_length() - 1
        neighbors_in = (masks[bit] & mask).bit_count()
        delta = graph.degree(bit) - 2 * neighbors_in
        cut += -delta if mask & (1 << bit) else delta
        mask ^= 1 << bit
        if cut > best:
            best, best_mask = cut, mask
    side = {i for i in range(n) if best_mask & (1 << i)}
    require(best == nx.cut_size(graph, side), "cut witness")
    return best, side


def flower_maxcut(k):
    """16-state exact transfer through flower modules; final c/d twist."""
    states = tuple(range(16))
    bits = lambda state: tuple((state >> i) & 1 for i in range(4))
    internal = {s: sum(bits(s)[0] != bits(s)[j] for j in (1, 2, 3)) for s in states}
    cross = lambda a, b: sum(bits(a)[j] != bits(b)[j] for j in (1, 2, 3))
    twist = lambda a, b: sum((bits(a)[i] != bits(b)[j]) for i, j in ((1, 1), (2, 3), (3, 2)))
    best = 0
    for first in states:
        score = {first: internal[first]}
        for _ in range(k - 1):
            score = {nxt: max(val + cross(last, nxt) + internal[nxt]
                              for last, val in score.items()) for nxt in states}
        best = max(best, max(val + twist(last, first) for last, val in score.items()))
    return best


def main():
    print("FINITE-EXACT named universe; Python integer/Fraction arithmetic")
    print("DEPENDENCIES: NetworkX", nx.__version__, "; SymPy", sp.__version__)
    graphs = {"K4": nx.complete_graph(4), "K3,3": nx.complete_bipartite_graph(3, 3),
              "triangular_prism": nx.circular_ladder_graph(3), "cube": nx.cubical_graph(),
              "pentagonal_prism": nx.circular_ladder_graph(5), "Petersen": nx.petersen_graph(),
              "Heawood": nx.heawood_graph(), "dodecahedron": nx.dodecahedral_graph(),
              "flower_J5": nx.Graph(legacy.flower_snark_edges(5))}
    for name, graph in graphs.items():
        count, oddness, nmatchings, factors = matching_audit(graph)
        maxcut, side = maximum_cut(graph)
        s = graph.number_of_edges() - maxcut
        q = Fraction(count * count * 3 ** s, 6 ** len(graph))
        require(all(graph.subgraph(part).degree(v) <= 1
                    for part in (side, set(graph) - side) for v in part), "maxcut side degree")
        e = edge_list(graph)[0]
        data = interface(graph, e)
        require(data["dot"] == count, "independent full coloring count")
        if len(graph) <= 10:
            direct = legacy.full_coloring_count(graph, edge_list(graph))
            require(direct == count, "direct edge-assignment count")
        signature = tuple(independent.matching_extension_count(tuple(graph), edge_list(graph), e, w)
                          for w in independent.REPRESENTATIVES)
        require(count == 6 * (signature[2] + signature[3]), "independent fourpole coloring")
        require(data["norms_squared"][1] == 3 * signature[0] ** 2 +
                6 * sum(x * x for x in signature[1:]), "independent fourpole norm")
        print(f"GRAPH {name}: n={len(graph)}, m={graph.number_of_edges()}, planar={nx.check_planarity(graph)[0]}, "
              f"bipartite={nx.is_bipartite(graph)}, 3edge_colors={count}, oddness={oddness}, "
              f"perfect_matchings={nmatchings}, maxcut={maxcut}, maxcut_fraction={Fraction(maxcut, graph.number_of_edges())}, s={s}, canonical_q={q}")
        print(f"  EDGE_INTERFACE {data}; signature={signature}")
        if name in {"K3,3", "Petersen", "pentagonal_prism", "cube"}:
            maxdata = interface(graph, side)
            require(min(maxdata["norms_squared"]) > 0, "nonvacuous maximum cut")
            require(maxdata["dot"] == count, "maximum cut gluing")
            require(Fraction(count * count, maxdata["norms_squared"][0] *
                             maxdata["norms_squared"][1]) == q, "canonical maxcut formula")
            print(f"  MAXCUT_INTERFACE {maxdata}")
        if name in {"Petersen", "pentagonal_prism"}:
            cycle_data = interface(graph, range(5))
            require(cycle_data["norms_squared"] == (30, 30), "C5 poles")
            require(cycle_data["dot"] == (0 if name == "Petersen" else 30), "C5 gluing")
            print(f"  C5_INTERFACE {cycle_data}")
    for name, graph in (("K5", nx.complete_graph(5)), ("K3,3", graphs["K3,3"]),
                        ("Petersen", graphs["Petersen"])):
        x = sp.Symbol("x")
        polynomial = sp.factor(sp.Matrix(nx.to_numpy_array(graph, dtype=int)).charpoly(x).as_expr())
        print(f"EXACT_SPECTRUM {name}: charpoly={polynomial}")
    for k in (5, 7, 9, 11, 21):
        maximum = flower_maxcut(k)
        if k == 5:
            require(maximum == maximum_cut(graphs["flower_J5"])[0], "flower transfer/bruteforce")
        print(f"FLOWER_MAXCUT J{k}: n={4*k}, m={6*k}, maxcut={maximum}, deficit={6*k-maximum}, ratio={Fraction(maximum,6*k)}")
    # Every labelled undirected Hamilton cycle on five fixed vertices, modulo
    # rotation/reversal, and the matching-glued two-cycle graph it specifies.
    base = nx.cycle_graph(5)
    base_edges = set(edge_list(base))
    cycle_overlap_histogram = Counter()
    for tail in permutations(range(1, 5)):
        cycle = (0,) + tail
        if cycle[1] > cycle[-1]:
            continue
        cycle_edges = {tuple(sorted((cycle[i], cycle[(i + 1) % 5]))) for i in range(5)}
        graph = nx.Graph(base_edges)
        graph.add_edges_from((a + 5, b + 5) for a, b in cycle_edges)
        graph.add_edges_from((i, i + 5) for i in range(5))
        overlap = len(base_edges & cycle_edges)
        count = matching_audit(graph)[0]
        cut_data = interface(graph, range(5))
        require(count == cut_data["dot"] == 6 * overlap, "five-pole intersection formula")
        cycle_overlap_histogram[overlap] += 1
    require(dict(sorted(cycle_overlap_histogram.items())) == {0: 1, 2: 5, 3: 5, 5: 1}, "cycle census")
    print("TWELVE_LABELLED_C5_GLUINGS overlap_histogram=", dict(sorted(cycle_overlap_histogram.items())))
    petersen = graphs["Petersen"]
    contracted_edges = {tuple(sorted((u % 5, v % 5))) for u, v in petersen.edges() if u % 5 != v % 5}
    require(contracted_edges == set(combinations(range(5), 2)), "spoke contraction K5")
    k33_paths = ((2, 3), (2, 1, 6), (2, 7), (8, 3), (8, 6),
                (8, 5, 7), (9, 4, 3), (9, 6), (9, 7))
    require({(path[0], path[-1]) for path in k33_paths} ==
            {(a, b) for a in (2, 8, 9) for b in (3, 6, 7)}, "K33 endpoints")
    used_internal = [v for path in k33_paths for v in path[1:-1]]
    require(len(used_internal) == len(set(used_internal)) == 3, "K33 disjoint path interiors")
    for path in k33_paths:
        require(all(petersen.has_edge(a, b) for a, b in zip(path, path[1:])), "K33 path edges")
    print("PETERSEN_K5_MINOR branch_sets=", [(i, i + 5) for i in range(5)])
    print("PETERSEN_K33_SUBDIVISION shores=((2,8,9),(3,6,7)), paths=", k33_paths)
    atlas_rows = []
    for graph in nx.graph_atlas_g():
        if not len(graph) or not all(d == 3 for _, d in graph.degree()):
            continue
        count = matching_audit(graph)[0]
        s = graph.number_of_edges() - maximum_cut(graph)[0]
        q = Fraction(count * count * 3 ** s, 6 ** len(graph))
        require(q <= Fraction(1, 4), "small atlas bound")
        atlas_rows.append((len(graph), count, s, str(q), nx.to_graph6_bytes(graph, header=False).decode().strip()))
    print("ALL_SIMPLE_CUBIC_GRAPHS_N_LE_7", atlas_rows)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
