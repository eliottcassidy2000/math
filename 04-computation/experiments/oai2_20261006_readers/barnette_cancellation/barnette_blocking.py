# Barnette graph with a nontrivial 3-edge cut: cube with one vertex replaced by (cube minus a vertex).
import itertools, random
import networkx as nx
from ortools.sat.python import cp_model
from barnette_lemmas2 import substitute
def has_hc(G, removed=(), forced=()):
    m = cp_model.CpModel(); arcs = []; lit = {}
    nodes = list(G.nodes()); idx = {v: i for i, v in enumerate(nodes)}
    for (u, v) in G.edges():
        e = frozenset((u, v))
        if e in removed: continue
        a = m.NewBoolVar(""); b = m.NewBoolVar("")
        arcs.append((idx[u], idx[v], a)); arcs.append((idx[v], idx[u], b)); lit[e] = (a, b)
    for e in forced:
        if e not in lit: return False
        m.AddBoolOr(list(lit[e]))
    m.AddCircuit(arcs)
    s = cp_model.CpSolver(); s.parameters.num_workers = 2
    return s.Solve(m) in (cp_model.OPTIMAL, cp_model.FEASIBLE)
G = substitute(nx.cubical_graph(), 0, random.Random(1), 'g')
G = nx.convert_node_labels_to_integers(G)
print("V =", G.number_of_nodes(), "planar", nx.check_planarity(G)[0], "bipartite", nx.is_bipartite(G), "3-conn", nx.node_connectivity(G))
E = [frozenset(e) for e in G.edges()]
one = [e for e in E if not has_hc(G, removed={e})]
print("single edges whose deletion kills all HCs:", len(one))
two = [p for p in itertools.combinations(E, 2) if not has_hc(G, removed=set(p))]
star = [p for p in two if len(p[0] & p[1]) == 1]
print("blocking 2-sets:", len(two), "of which stars (share a vertex):", len(star), "non-star:", len(two) - len(star))
# P4 (3-edge paths): all contained in an HC?
bad = 0; tot = 0
for v in G.nodes():
    for a in G[v]:
        for b in G[a]:
            if b == v: continue
            for c in G[b]:
                if c in (v, a): continue
                tot += 1
                if not has_hc(G, forced=[frozenset((v, a)), frozenset((a, b)), frozenset((b, c))]): bad += 1
print("directed 3-edge paths:", tot, "not in any HC:", bad)
