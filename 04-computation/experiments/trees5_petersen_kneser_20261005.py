#!/usr/bin/env python3
"""
trees5_petersen_kneser_20261005.py

The three unlabeled trees on 5 vertices against the Kuratowski graphs K_5, K_{3,3} and the Petersen
graph (= Kneser graph K(5,2): vertices = edges of K_5, adjacent iff disjoint; Aut = S_5).

Exact facts frozen here:
  A. spanning trees of K_5: 125 = 5^3 labelled (Cayley), 60 paths + 60 forks + 5 stars.
  B. 4-edge subgraphs of K_5 up to S_5 = 4-subsets of V(Petersen) up to Aut(Petersen): 6 orbits;
     induced Petersen subgraphs: star <-> 4K1 (the 5 maximum independent sets, Erdos--Ko--Rado),
     fork <-> P3 + K1, path <-> P4; the three non-trees: paw <-> K2 + 2K1, C4 <-> 2K2, C3+K2 <-> K_{1,3}.
     Characterization: a 4-subset of Petersen is a spanning tree of K_5 iff its induced subgraph is
     triangle-free-trivially (girth 5) AND ... -- we record the exact orbit table instead.
  C. containment of the three trees: K_5 contains all three; K_{3,3} and Petersen (cubic) contain the
     path and the fork but not the star; both are Erdos--Sos star-extremal graphs at t = 5
     (3n/2 edges, K_{1,4}-free): K_{3,3} at n = 6, Petersen at n = 10.
  D. Petersen has exactly 5 maximum independent sets (size 4), all stars of K_5; independence number 4.
  E. orientation sidecar: the 27 oriented 5-trees (A000238) and which embed in the Paley-type
     orientations?  (not here; see sumner_t5_20261005.py)
Reproduce: python3 trees5_petersen_kneser_20261005.py
"""
import itertools, math, sys
from collections import Counter
import networkx as nx

CHECKS = 0
def check(c, msg):
    global CHECKS
    CHECKS += 1
    if not c:
        print("CHECK FAILED:", msg); sys.exit(1)

def tree_name(T):
    degs = sorted((d for _, d in T.degree()), reverse=True)
    if degs[0] == T.number_of_nodes() - 1: return "star"
    if degs[0] == 2: return "path"
    return "fork"

trees5 = {tree_name(T): T for T in nx.nonisomorphic_trees(5)}
K5 = nx.complete_graph(5)
edges5 = list(itertools.combinations(range(5), 2))
assert len(edges5) == 10

print("=" * 78)
print("A. Spanning trees of K_5 by type (all 4-subsets of the 10 edges)")
print("=" * 78)
cnt = Counter(); nontree = Counter()
for S in itertools.combinations(edges5, 4):
    G = nx.Graph(); G.add_nodes_from(range(5)); G.add_edges_from(S)
    if nx.is_tree(G):
        cnt[tree_name(G)] += 1
    else:
        H = nx.Graph(S)  # drop isolated vertices for naming
        degs = tuple(sorted((d for _, d in H.degree()), reverse=True))
        cyc = nx.cycle_basis(H)
        name = {((2,2,2,2), 1): "C4", ((3,2,2,1), 1): "paw", ((2,2,2,1,1), 1): "C3+K2"}.get((degs, len(cyc)), str((degs, len(cyc))))
        nontree[name] += 1
check(sum(cnt.values()) == 125, "Cayley 5^3")
check(cnt == Counter({"path": 60, "fork": 60, "star": 5}), "60/60/5")
check(sum(cnt.values()) + sum(nontree.values()) == math.comb(10, 4), "210 four-subsets")
print("  spanning trees:", dict(cnt), " total 125 = 5^3")
print("  non-tree 4-edge subgraphs:", dict(nontree), " total", sum(nontree.values()))

print("\n" + "=" * 78)
print("B. Petersen = K(5,2): orbits of 4-subsets of vertices under Aut = S_5, induced subgraphs")
print("=" * 78)
P = nx.Graph()
P.add_nodes_from(edges5)
for e, f in itertools.combinations(edges5, 2):
    if not set(e) & set(f):
        P.add_edge(e, f)
check(P.number_of_nodes() == 10 and P.number_of_edges() == 15 and all(d == 3 for _, d in P.degree()), "Petersen cubic 10/15")
check(nx.is_isomorphic(P, nx.petersen_graph()), "is the Petersen graph")
check(nx.girth(P) == 5 if hasattr(nx, "girth") else True, "girth 5")
def induced_name(H):
    degs = tuple(sorted((d for _, d in H.degree()), reverse=True))
    return {(0,0,0,0): "4K1", (1,1,0,0): "K2+2K1", (1,1,1,1): "2K2", (2,1,1,0): "P3+K1", (2,2,1,1): "P4", (3,1,1,1): "K_{1,3}"}.get(degs, str(degs))
orbit = Counter(); witness = {}
for S in itertools.combinations(edges5, 4):
    G = nx.Graph(); G.add_nodes_from(range(5)); G.add_edges_from(S)
    k5type = tree_name(G) if nx.is_tree(G) else "non-tree"
    H = P.subgraph(S)
    orbit[(k5type, induced_name(H))] += 1
print("  (K_5 4-edge type, induced Petersen subgraph on the 4 vertices): count")
for k, v in sorted(orbit.items()):
    print(f"    {k}: {v}")
check(orbit[("path", "P4")] == 60 and orbit[("fork", "P3+K1")] == 60 and orbit[("star", "4K1")] == 5, "tree orbits")
check(set(k[1] for k in orbit if k[0] == "non-tree") == {"K2+2K1", "2K2", "K_{1,3}"}, "non-tree orbits")
print("  => star <-> 4K1 (independent 4-sets), fork <-> P3+K1, path <-> P4; the 4-subset is a spanning tree of K_5")
print("     iff its induced Petersen graph is one of 4K1, P3+K1, P4 (and not K2+2K1, 2K2, K_{1,3}).")
# independent sets
indep4 = [S for S in itertools.combinations(edges5, 4) if P.subgraph(S).number_of_edges() == 0]
indep5 = [S for S in itertools.combinations(edges5, 5) if P.subgraph(S).number_of_edges() == 0]
check(len(indep4) == 5 and len(indep5) == 0, "alpha(Petersen) = 4 with exactly 5 maximum independent sets")
check(all(len(set.intersection(*map(set, S))) == 1 for S in indep4), "each is a star of K_5 (EKR)")
print("  alpha(Petersen) = 4; the 5 maximum independent sets are the 5 stars of K_5 (Erdos--Ko--Rado for (5,2)).")

print("\n" + "=" * 78)
print("C. Containment of the three trees in K_5, K_{3,3}, Petersen; Erdos--Sos star-extremality")
print("=" * 78)
hosts = {"K_5": K5, "K_{3,3}": nx.complete_bipartite_graph(3, 3), "Petersen": nx.petersen_graph()}
def contains(G, T):
    return nx.algorithms.isomorphism.GraphMatcher(G, T).subgraph_is_monomorphic()
for hn, G in hosts.items():
    row = {tn: contains(G, T) for tn, T in trees5.items()}
    n, e = G.number_of_nodes(), G.number_of_edges()
    print(f"  {hn:8s} n={n:2d} e={e:2d} 3n/2={1.5*n:5.1f}  contains path:{row['path']} fork:{row['fork']} star:{row['star']}  maxdeg={max(d for _, d in G.degree())}")
check(all(contains(K5, T) for T in trees5.values()), "K5 contains all")
for hn in ("K_{3,3}", "Petersen"):
    G = hosts[hn]
    check(contains(G, trees5["path"]) and contains(G, trees5["fork"]) and not contains(G, trees5["star"]), f"{hn} star-free, path/fork present")
    check(G.number_of_edges() * 2 == 3 * G.number_of_nodes(), f"{hn} has exactly 3n/2 edges")
print("  K_{3,3} (n=6) and Petersen (n=10) are cubic, hence K_{1,4}-free with exactly 3n/2 = ex(n, K_{1,4}) edges:")
print("  both are Erdos--Sos star-extremal witnesses at t = 5; neither is path/fork-extremal (those need disjoint K_4's).")

print("\n" + "=" * 78)
print("D. Bipartition obstruction: a tree embeds in K_{3,3} iff both colour classes have <= 3 vertices")
print("=" * 78)
for tn, T in trees5.items():
    col = nx.bipartite.color(T)
    parts = sorted(Counter(col.values()).values())
    print(f"  {tn:5s} colour classes {parts}  embeds in K_(3,3): {parts[1] <= 3}")
    check((parts[1] <= 3) == contains(hosts["K_{3,3}"], T), "bipartition criterion")

print(f"\nALL {CHECKS} CHECKS PASSED")
