#!/usr/bin/env python3
"""Audit B, item 6: the 14-vertex Barnette graph remark ('15 of the 57 minimum blocking pairs are not stars').
Own construction: cube Q3 with vertex 0 replaced by a copy of Q3 minus a vertex, the three dangling edges
attached by each of the 6 bijections; keep the planar, bipartite, 3-connected ones; count edge pairs whose
deletion leaves no Hamiltonian cycle (exact DFS on 14 vertices)."""
import itertools
import networkx as nx

def hc_exists(n, adj, removed):
    adjl = [[w for w in adj[v] if frozenset((v, w)) not in removed] for v in range(n)]
    if any(len(a) < 2 for a in adjl):
        return False
    start = 0
    seen = [False] * n
    seen[start] = True
    def rec(v, k):
        if k == n:
            return start in adjl[v]
        for w in adjl[v]:
            if not seen[w]:
                seen[w] = True
                if rec(w, k + 1):
                    return True
                seen[w] = False
        return False
    return rec(start, 1)

results = set()
Q = [(a, b) for a in range(8) for b in range(8) if a < b and bin(a ^ b).count('1') == 1]
for perm in itertools.permutations([1, 2, 4]):
    G = nx.Graph()
    G.add_edges_from((('o', a), ('o', b)) for a, b in Q if a != 0 and b != 0)
    G.add_edges_from((('i', a), ('i', b)) for a, b in Q if a != 0 and b != 0)
    for x, y in zip([1, 2, 4], perm):
        G.add_edge(('o', x), ('i', y))
    planar = nx.check_planarity(G)[0]
    bip = nx.is_bipartite(G)
    conn = nx.node_connectivity(G)
    H = nx.convert_node_labels_to_integers(G)
    n = H.number_of_nodes()
    adj = [list(H[v]) for v in range(n)]
    E = [frozenset(e) for e in H.edges()]
    single = sum(1 for e in E if not hc_exists(n, adj, {e}))
    two = [p for p in itertools.combinations(E, 2) if not hc_exists(n, adj, set(p))]
    star = sum(1 for p in two if len(p[0] & p[1]) == 1)
    results.add((planar, bip, conn, n, len(E), single, len(two), star, len(two) - star))
    print(f"perm {perm}: planar={planar} bipartite={bip} 3-conn={conn>=3} V={n} E={len(E)} "
          f"single blockers={single} blocking pairs={len(two)} stars={star} non-stars={len(two)-star}")
