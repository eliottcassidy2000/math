#!/usr/bin/env python3
"""Exact checks of the new lemmas of the Barnette preprint (openai/math #180) on duals of small
Barnette graphs G (cubic, bipartite, planar, 3-connected).  Dual T: faces of G = vertices of T;
vertices of G = triangles of T, black/white by the bipartition.  t0 = the triangle of a black
vertex g0; A = other black triangles; X = faces of G not incident to g0.
Checks: (i) J(r,s) = (1/3) sum delta(r(t),s(t)) is an integer for every pair (disk identity);
(ii) for every s with cyclic Q_s, sum_r i^J(r,s) = 0 (Lemma 4.2);
(iii) some pair has Q_s a forest (Prop 4.1), and the coefficient sum over forest s is nonzero;
(iv) the Hamiltonian cycle produced at the end (via Theorem 6.1's partition) -- we only check
that a forest pair exists and count them."""
import itertools, sys
import networkx as nx


def barnette_graphs():
    out = {}
    out["cube"] = nx.cubical_graph()
    for k in (3, 4, 5, 6):
        out[f"prism{2*k}"] = nx.circular_ladder_graph(2 * k)
    # permutohedron: Cayley(S4, adjacent transpositions) = truncated octahedron
    P = nx.Graph()
    for p in itertools.permutations(range(4)):
        for i in range(3):
            q = list(p); q[i], q[i + 1] = q[i + 1], q[i]
            P.add_edge(p, tuple(q))
    out["truncated_octahedron"] = P
    return out


def analyse(name, G):
    ok, emb = nx.check_planarity(G)
    assert ok and all(d == 3 for _, d in G.degree()) and nx.is_bipartite(G) and nx.node_connectivity(G) >= 3
    col = nx.bipartite.color(G)
    faceid = {}
    def face_of(v, w):
        f = frozenset(emb.traverse_face(v, w))
        if f not in faceid:
            faceid[f] = len(faceid)
        return faceid[f]
    around = {}
    for g in G.nodes():
        nb = list(emb.neighbors_cw_order(g))
        around[g] = [face_of(g, w) for w in nb]
        assert len(set(around[g])) == 3
    nF = len(faceid)
    black = [g for g in G if col[g] == 0]
    g0 = black[0]
    roots = set(around[g0])
    A = [g for g in black if g != g0]
    X = [f for f in range(nF) if f not in roots]
    assert len(A) == len(X)
    def delta(gb, v, w):
        c = around[gb]; i, j = c.index(v), c.index(w)
        return 1 if (i + 1) % 3 == j else -1
    # states: perfect matchings A -> X with r(t) in around[t]
    states = []
    def rec(i, used, cur):
        if i == len(A):
            states.append(tuple(cur)); return
        for f in around[A[i]]:
            if f in X and f not in used:
                used.add(f); cur.append(f); rec(i + 1, used, cur); cur.pop(); used.discard(f)
    rec(0, set(), [])
    def Q_cyclic(s):
        par = list(range(nF))
        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]; a = par[a]
            return a
        for t, f in zip(A, s):
            o = [h for h in around[t] if h != f]
            ra, rb = find(o[0]), find(o[1])
            if ra == rb:
                return True
            par[ra] = rb
        return False
    nonint = 0; cancel_fail = 0; forest_pairs = 0; forest_coeffs = []
    npairs = 0
    for s in states:
        c = 0
        # partners r: states with r(t) != s(t) for all t
        cyc = Q_cyclic(s)
        for r in states:
            if any(a == b for a, b in zip(r, s)):
                continue
            npairs += 1
            S3 = sum(delta(t, a, b) for t, a, b in zip(A, r, s))
            if S3 % 3:
                nonint += 1; continue
            c += 1j ** (S3 // 3)
            if not cyc:
                forest_pairs += 1
        if cyc and abs(c) > 1e-9:
            cancel_fail += 1
        if not cyc:
            forest_coeffs.append(c)
    nz = sum(1 for c in forest_coeffs if abs(c) > 1e-9)
    print(f"{name}: |V(G)|={G.number_of_nodes()}, k={len(black)}; states={len(states)}, pairs={npairs}; "
          f"non-integral J: {nonint}; cyclic-Q_s with nonzero partner sum: {cancel_fail}; "
          f"pairs with Q_s a forest: {forest_pairs}; forest s with nonzero coefficient: {nz}/{len(forest_coeffs)}",
          flush=True)


if __name__ == "__main__":
    for name, G in barnette_graphs().items():
        analyse(name, G)
