"""Independent networkx check for HYP-9464 (A) (code independent of selfdual3_fewcycle_20261104.c).
(1) AT_n and TT_n are {-3}-self-dual for n = 9, 10, 11 (every T - X, |X| = 3, isomorphic to its converse);
(2) three random tournaments of order 10 (seed 1) are not;
(3) the two order-10 classes reported by `selfdual3 10` (rows in the record given as argv[1]) are TT_10 and AT_10.
Usage: python3 selfdual3_networkx_check_20261104.py selfdual3_fewcycle_20261104_records/sd3_10.out
kind-pasteur-2026-11-04-S10."""
import itertools, random, sys
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher

def tour_AT(n):
    G = nx.DiGraph(); G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) == (0, n - 1): G.add_edge(j, i)
            else: G.add_edge(i, j)
    return G

def tour_TT(n):
    G = nx.DiGraph(); G.add_nodes_from(range(n)); G.add_edges_from((i, j) for i in range(n) for j in range(i + 1, n))
    return G

def selfdual3(G):
    for X in itertools.combinations(G.nodes, 3):
        H = G.subgraph([v for v in G.nodes if v not in X]).copy()
        if not DiGraphMatcher(H, H.reverse(copy=True)).is_isomorphic():
            return False
    return True

for n in (9, 10, 11):
    print(f"n={n}: AT_n {{-3}}-self-dual: {selfdual3(tour_AT(n))}; TT_n {{-3}}-self-dual: {selfdual3(tour_TT(n))}")
random.seed(1)
for t in range(3):
    G = nx.DiGraph(); G.add_nodes_from(range(10))
    for i in range(10):
        for j in range(i + 1, 10):
            G.add_edge(i, j) if random.random() < 0.5 else G.add_edge(j, i)
    print(f"random order-10 tournament {t}: {{-3}}-self-dual: {selfdual3(G)}")
if len(sys.argv) > 1:
    for l in open(sys.argv[1]):
        if not l.startswith("SD3 n=10 "): continue
        rows = [int(x, 16) for x in l.split("rows:")[1].split()]
        n = len(rows)
        G = nx.DiGraph(); G.add_nodes_from(range(n))
        G.add_edges_from((v, w) for v in range(n) for w in range(n) if rows[v] >> w & 1)
        assert G.number_of_edges() == n * (n - 1) // 2
        print(l.split(" rows:")[0], "| iso AT_10:", DiGraphMatcher(G, tour_AT(n)).is_isomorphic(),
              "| iso TT_10:", DiGraphMatcher(G, tour_TT(n)).is_isomorphic())
