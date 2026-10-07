import itertools, random
import networkx as nx
import barnette_blocking as BB   # rebuilds G and the lists (prints its summary)
G, E, two = BB.G, BB.E, BB.two
nonstar = [p for p in two if len(p[0] & p[1]) == 0]
def in_3cut(p):
    for g in E:
        if g in p: continue
        H = G.copy(); H.remove_edges_from([tuple(e) for e in (p[0], p[1], g)])
        if not nx.is_connected(H): return True
    return False
print("non-star blocking pairs:", len(nonstar), "; contained in a 3-edge cut:", sum(in_3cut(p) for p in nonstar))
