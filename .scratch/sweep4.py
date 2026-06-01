from lrc_coloring_core import *
import collections, itertools

# PLANARITY (part 3). G(t) is planar unless it contains K5 or K3,3 minor (Kuratowski/Wagner).
# Since chi=omega and max clique up to n, for n>=5 the t near 0 gives K_n >= K5 -> NON-PLANAR.
# Implement a planarity test from scratch (no networkx guarantee): use the
# Euler/Kuratowski via simple K5,K3,3 subgraph search is insufficient (need subdivisions),
# but for these tiny graphs (n<=7) we can test planarity by trying to find K5 or K3,3 MINOR
# through brute subgraph + the fact that for <=7 vertices nonplanarity always shows as a
# K5 or K3,3 subdivision; easier: use the boost-style left-right? Too much.
# Use: a graph on <=7 vertices, <= C(7,2)=21 edges. Planar iff |E|<=3|V|-6 is necessary not
# sufficient. We'll just detect K5 subgraph and K3,3 subgraph (not subdivision). For these
# unit-arc graphs the relevant obstruction IS a clique K5, so detecting omega>=5 suffices for
# our positive nonplanarity claims; we additionally do a real planarity check via networkx if available.
try:
    import networkx as nx
    HAVE_NX = True
except Exception:
    HAVE_NX = False
print("networkx available:", HAVE_NX)

def is_planar(N, edges):
    if HAVE_NX:
        G=nx.Graph(); G.add_nodes_from(range(N)); G.add_edges_from(tuple(e) for e in edges)
        ok,_=nx.check_planarity(G)
        return ok
    return None

# Report: for each set, fraction of cells nonplanar, and whether nonplanarity <=> omega>=5
for speeds in [(1,2,3,4),(1,2,3,5),(2,3,5,7),(1,2,3,4,5),(1,3,4,5,9),(1,2,3,4,5,6),(1,3,5,7,9,11)]:
    n=len(speeds)+1
    crit=critical_times(speeds,n); mids=cell_midpoints(crit)
    nonplanar=0; nonplanar_eq_clique5=0; total=len(mids)
    for t in mids:
        N,edges,pts=danger_graph(speeds,n,t)
        w=clique_number_arc(speeds,n,t)
        pl=is_planar(N,edges)
        if pl is False:
            nonplanar+=1
            if w>=5: nonplanar_eq_clique5+=1
    print(f"{speeds} n={n}: nonplanar {nonplanar}/{total}; of those omega>=5: {nonplanar_eq_clique5}")
