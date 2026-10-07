import itertools, random, sys
import networkx as nx
from barnette_lemmas import analyse

def b3_cayley():
    G = nx.Graph()
    elems = []
    for p in itertools.permutations(range(3)):
        for sg in itertools.product((1, -1), repeat=3):
            elems.append(tuple(sg[i] * (p[i] + 1) for i in range(3)))
    def s1(x): return (x[1], x[0], x[2])
    def s2(x): return (x[0], x[2], x[1])
    def s3(x): return (x[0], x[1], -x[2])
    for x in elems:
        for s in (s1, s2, s3):
            G.add_edge(x, s(x))
    return G

def substitute(G, v, rng, tag):
    """replace vertex v by Q3 minus a vertex; keep if planar"""
    nb = list(G.neighbors(v))
    Q = nx.cubical_graph()   # nodes 0..7; remove node 0
    w = 0
    bnd = list(Q.neighbors(w))
    for perm in itertools.permutations(nb):
        H = G.copy(); H.remove_node(v)
        mp = {q: (tag, q) for q in Q.nodes() if q != w}
        for a, b in Q.edges():
            if w not in (a, b):
                H.add_edge(mp[a], mp[b])
        for bq, x in zip(bnd, perm):
            H.add_edge(mp[bq], x)
        if nx.check_planarity(H)[0] and nx.is_bipartite(H):
            return H
    return None

G = b3_cayley()
assert G.number_of_nodes() == 48
analyse("omnitruncated_cube(B3)", G)
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
for trial in range(6):
    H = nx.cubical_graph()
    steps = rng.randint(2, 4)
    for s in range(steps):
        v = rng.choice(list(H.nodes()))
        H2 = substitute(H, v, rng, (trial, s))
        if H2 is not None:
            H = H2
    H = nx.convert_node_labels_to_integers(H)
    analyse(f"random_substitution_{trial}", H)
