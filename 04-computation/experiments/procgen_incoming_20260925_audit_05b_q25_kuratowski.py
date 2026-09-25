#!/usr/bin/env python3
"""New observation check: the square-sum graph Q_N is planar for N<=24 and
nonplanar from N=25.  Extract a Kuratowski witness for Q_25, identify K5 vs K33
subdivision by branch-vertex degrees, and test whether the witness uses the ear
11-25-24 (the ear forced in every Hamiltonian path of Q_25)."""
import networkx as nx
from math import isqrt


def Q(N):
    g = nx.Graph()
    g.add_nodes_from(range(1, N + 1))
    for x in range(1, N + 1):
        for y in range(x + 1, N + 1):
            s = x + y
            if isqrt(s) ** 2 == s:
                g.add_edge(x, y)
    return g


def main():
    q24 = Q(24)
    ok, emb = nx.check_planarity(q24)
    print("Q24 planar:", ok, " V,E =", q24.number_of_nodes(), q24.number_of_edges())
    if ok:
        emb.check_structure()
        faces = 0
        seen = set()
        for u, v in emb.edges():
            if (u, v) not in seen:
                f = emb.traverse_face(u, v, mark_half_edges=seen)
                faces += 1
        comps = nx.number_connected_components(q24)
        print("  embedding faces F=%d, V-E+F=%d (components %d)" % (faces, q24.number_of_nodes() - q24.number_of_edges() + faces, comps))
    q24e = q24.copy()
    q24e.add_edge(11, 24)
    print("Q24 + marked edge {11,24} planar:", nx.check_planarity(q24e)[0])
    q25 = Q(25)
    ok25, cert = nx.check_planarity(q25, counterexample=True)
    print("Q25 planar:", ok25, " V,E =", q25.number_of_nodes(), q25.number_of_edges(), " N(25) =", sorted(q25.neighbors(25)))
    deg = dict(cert.degree())
    branch = sorted(v for v, d in deg.items() if d >= 3)
    print("Kuratowski witness: %d vertices, %d edges, branch vertices %s with degrees %s"
          % (cert.number_of_nodes(), cert.number_of_edges(), branch, [deg[v] for v in branch]))
    print("witness uses vertex 25:", 25 in cert, " edges at 25:", sorted(cert.edges(25)) if 25 in cert else None)
    # suppress degree-2 vertices to show the underlying K33/K5
    h = nx.Graph(cert)
    changed = True
    while changed:
        changed = False
        for v in list(h.nodes()):
            if h.degree(v) == 2:
                a, b = list(h.neighbors(v))
                if not h.has_edge(a, b):
                    h.remove_node(v)
                    h.add_edge(a, b)
                    changed = True
                    break
    print("suppressed witness: V=%d E=%d, is K33: %s, is K5: %s"
          % (h.number_of_nodes(), h.number_of_edges(),
             nx.is_isomorphic(h, nx.complete_bipartite_graph(3, 3)), nx.is_isomorphic(h, nx.complete_graph(5))))
    if nx.is_isomorphic(h, nx.complete_bipartite_graph(3, 3)):
        A, B = nx.bipartite.sets(h)
        print("  shores:", sorted(A), sorted(B))
    # list the subdivision paths
    # also check whether every Kuratowski subgraph must use 25: Q25 - 25 = Q24 planar, so yes.
    print("Q25 minus vertex 25 is Q24 (planar) -> every Kuratowski subgraph of Q25 uses vertex 25 (hence the ear 11-25-24).")
    # first nonplanar and Hamiltonicity thresholds side by side
    flags = []
    for N in range(10, 41):
        flags.append((N, nx.check_planarity(Q(N))[0]))
    print("planarity of Q_N, N=10..40:", flags)


if __name__ == "__main__":
    main()
