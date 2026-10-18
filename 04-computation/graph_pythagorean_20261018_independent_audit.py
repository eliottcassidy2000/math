#!/usr/bin/env python3
"""Independent exact audit: cycle-space flows, cuts, and finite-field signs.

Universe: named cubic controls, all cubic graph-atlas members (orders <=7),
24 seeded regular-graph samples, six triangle expansions, and Paley squares
over F_(p^2) for p=3,5,7. No numerical eigensolver or edge-colouring recursion.
Run from repository root; stdout is the matching stored .out receipt.
"""
from fractions import Fraction
from itertools import combinations
import networkx as nx
import sympy as sp


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def canonical(g):
    return nx.convert_node_labels_to_integers(g, ordering="sorted")


def flow_count(g):
    """Nonzero F2^2 flows = ordered pairs of binary cycles covering E."""
    edges = sorted(tuple(sorted(e)) for e in g.edges())
    index = {e: i for i, e in enumerate(edges)}
    basis = []
    for cycle in nx.cycle_basis(g):
        mask = 0
        for a, b in zip(cycle, cycle[1:] + cycle[:1]):
            mask ^= 1 << index[tuple(sorted((a, b)))]
        basis.append(mask)
    check(len(basis) == len(edges) - len(g) + nx.number_connected_components(g),
          "wrong cycle-space rank")
    vectors = [0]
    for row in basis:
        vectors += [x ^ row for x in vectors]
    full = (1 << len(edges)) - 1
    return sum((x | y) == full for x in vectors for y in vectors)


def cut_audit(g):
    """Every cut, modulo complement; also certify every maximizing cut."""
    n = len(g)
    edges = sorted(g.edges())
    best, winners = -1, []
    for mask in range(1 << (n - 1)):
        size = sum(((mask >> a) ^ (mask >> b)) & 1 for a, b in edges)
        if size > best:
            best, winners = size, [mask]
        elif size == best:
            winners.append(mask)
    s = len(edges) - best
    products = set()
    for mask in winners:
        factors = []
        for bit in (0, 1):
            side = [v for v in g if ((mask >> v) & 1) == bit]
            h = g.subgraph(side)
            check(all(d <= 1 for _, d in h.degree()), "maximum-cut side not matching")
            j = h.number_of_edges()
            i = len(side) - 2 * j
            factors.append(6 ** i * 12 ** j)
        products.add(factors[0] * factors[1])
    check(products == {6 ** n // 3 ** s}, "cut-choice dependence")
    return best, s, len(winners), next(iter(products))


def expand_vertex(g):
    g = g.copy()
    v = min(g)
    neighbors = sorted(g[v])
    fresh = list(range(max(g) + 1, max(g) + 4))
    g.remove_node(v)
    g.add_edges_from(zip(neighbors, fresh))
    g.add_edges_from(combinations(fresh, 2))
    return canonical(g)


def seidel(g):
    g = canonical(g)
    n = len(g)
    return sp.ones(n) - sp.eye(n) - 2 * sp.Matrix(nx.to_numpy_array(g, dtype=int))


def switch(g, subset):
    h = g.copy()
    subset = set(subset)
    for x in subset:
        for y in set(g) - subset:
            if h.has_edge(x, y):
                h.remove_edge(x, y)
            else:
                h.add_edge(x, y)
    return h


def paley_square(p):
    d = next(a for a in range(2, p) if pow(a, (p - 1) // 2, p) == p - 1)
    points = [(a, b) for a in range(p) for b in range(p)]

    def multiply(x, y):
        a, b = x
        c, e = y
        return ((a * c + d * b * e) % p, (a * e + b * c) % p)

    def character(x):
        a, b = x
        norm = (a * a - d * b * b) % p
        if norm == 0:
            return 0
        return 1 if pow(norm, (p - 1) // 2, p) == 1 else -1

    g = nx.Graph()
    g.add_nodes_from(range(p * p + 1))
    idx = {x: i + 1 for i, x in enumerate(points)}
    for x, y in combinations(points, 2):
        diff = ((x[0] - y[0]) % p, (x[1] - y[1]) % p)
        if character(diff) == 1:
            g.add_edge(idx[x], idx[y])
    check(seidel(g) ** 2 == p * p * sp.eye(len(g)), "Paley conference identity")
    w = next(x for x in points if character(x) == -1)
    selected = {idx[multiply(w, (a, c))]
                for c in range((p - 1) // 2) for a in range(p)}
    h = switch(g, selected)
    degree = p * (p - 1) // 2
    check(set(dict(h.degree()).values()) == {degree}, "regular-switch degree")
    check(seidel(h) ** 2 == p * p * sp.eye(len(h)), "regular conference identity")
    c = (p * p + 1) // 2
    b = c - 1
    check(p * p + b * b == c * c, "spectral triple identity")
    if p == 3:
        check(nx.is_isomorphic(h, nx.petersen_graph()), "Paley root not Petersen")
    print("paley", p, "field_nonsquare", d, "switch_size", len(selected),
          "order_degree", (len(h), degree), "spectral_triple", (p, b, c))


def main():
    controls = [
        ("K4", nx.complete_graph(4), 6, 4),
        ("K33", nx.complete_bipartite_graph(3, 3), 12, 9),
        ("prism3", nx.circular_ladder_graph(3), 6, 7),
        ("cube", canonical(nx.cubical_graph()), 24, 12),
        ("prism5", nx.circular_ladder_graph(5), 30, 13),
        ("Petersen", nx.petersen_graph(), 0, 12),
        ("Heawood", nx.heawood_graph(), 48, 21),
    ]
    print("INDEPENDENT EXACT AUDIT: F2^2 cycle-space flow counts + exhaustive cuts")
    for name, raw, expected, expected_cut in controls:
        g = canonical(raw)
        c = flow_count(g)
        mc, s, cuts, m = cut_audit(g)
        check((c, mc) == (expected, expected_cut), "named control mismatch " + name)
        q = Fraction(c * c, m)
        check(q <= Fraction(1, 4), "range failure")
        print(name, "n", len(g), "flow_count", c, "maxcut", mc,
              "maximizing_cuts", cuts, "norm_product", m, "q", q)

    atlas = [g for g in nx.graph_atlas_g() if len(g) > 0
             and all(d == 3 for _, d in g.degree())]
    randoms = [nx.random_regular_graph(3, n, seed=seed)
               for n in (8, 10, 12) for seed in range(8)]
    for raw in atlas + randoms:
        g = canonical(raw)
        c = flow_count(g)
        mc, s, cuts, m = cut_audit(g)
        q = Fraction(c * c, m)
        check(0 <= q <= Fraction(1, 4), "sample range failure")
        check(q < Fraction(1, 4) or nx.is_isomorphic(g, nx.complete_graph(4)),
              "equality failure")
    print("all_cubic_atlas_orders_le7", len(atlas), "seeded_samples", len(randoms))

    g = nx.complete_graph(4)
    for k in range(6):
        c = flow_count(g)
        mc, s, cuts, m = cut_audit(g)
        check(c == 6 and nx.check_planarity(g)[0] and not list(nx.bridges(g)),
              "triangle-expansion invariant failure")
        print("triangle_expansion", k, "n", len(g), "C", c, "s", s,
              "q", Fraction(c * c, m))
        g = expand_vertex(g)

    for p in (3, 5, 7):
        paley_square(p)

    for name, raw, expected in [("K5", nx.complete_graph(5), 4),
                                ("K33", nx.complete_bipartite_graph(3, 3), 5),
                                ("Petersen", nx.petersen_graph(), 3)]:
        spectrum = seidel(raw).eigenvals()
        check(max(abs(x) for x in spectrum) == expected, "Seidel root")
        print("seidel_spectrum", name, dict(sorted(spectrum.items())))
    empty = switch(nx.complete_bipartite_graph(3, 3), range(3))
    check(empty.number_of_edges() == 0, "K33 switching hostile")
    check(nx.check_planarity(empty)[0] and not nx.check_planarity(nx.complete_bipartite_graph(3, 3))[0],
          "planarity hostile")
    print("same_Seidel_spectrum_planarity_hostile", "K33 -> empty6")
    print("PASS: all assertions are active under python -O too")


if __name__ == "__main__":
    main()
