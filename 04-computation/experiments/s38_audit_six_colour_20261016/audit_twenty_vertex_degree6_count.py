#!/usr/bin/env python3
"""S38 audit (independent referee), exploratory: is a 20-vertex failed exterior possible at all?

A failed exterior B with |B| = 20 is vertex-minimal (no five-core has 17-19 vertices; the 16-vertex one is
degree-choosable).  With degrees 5 or 6 (degree 7 needs 22 vertices), h high vertices (degree 6, list 5), s low
degree-6 vertices (list 6) and 20 - k degree-5 vertices (k = h + s), the low-degree count of THM-4661 step (iv)
gives  5(20 - k) + 6 s <= (5/2)(20 - h - 1) + 6 h,  i.e.  8.5 k >= 52.5 + 9.5 s,  so k >= 7, and k is even
(degree sum 100 + k), so k >= 8.

SAT query: a graph on 20 vertices, triangle-free, all codegrees <= 2, degrees in {5,6}, at least K vertices of
degree 6.  Redundant constraints derived by hand (sound for every such graph): a degree-6 vertex z has at most two
degree-6 neighbours (|N_2(z)| >= (24 + t_z)/2 and |N_2(z)| <= 13).  Symmetry breaking: vertex 0 has degree 6 with
N(0) = {1..6}, and its (at most two) degree-6 neighbours are among {1, 2}.  Solver glucose4 (honours interrupt();
cadical195 does not in this pysat build -- a first version with cadical195 ran past its 480 s limit and was killed
by the outer timeout).  A TIMEOUT proves nothing.
"""
import itertools, json, os, sys, threading, time
import networkx as nx
from pysat.solvers import Solver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_nineteen_vertex_exhaustive import GraphSAT

HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'audit_twenty_vertex_degree6_count.json')


def log(*a):
    print(*a, flush=True)


def query(K, limit_s):
    n = 20
    gs = GraphSAT(n)
    ys = []
    for v in range(n):
        gs.card(gs.nbrs(v), 5, 'atleast')
        gs.card(gs.nbrs(v), 6, 'atmost')
        y = gs.pool.id(('y', v))
        gs.card(gs.nbrs(v), 6, 'atleast', cond=y)      # y_v -> degree 6
        ys.append(y)
    # not y_v -> degree exactly 5 (so y_v <-> degree 6)
    for v in range(n):
        nv = gs.pool.id(('ny', v))
        gs.cls.append([ys[v], nv])                       # not y_v -> nv
        gs.card(gs.nbrs(v), 5, 'atmost', cond=nv)
    # redundant: y_z -> at most two degree-6 neighbours
    for z in range(n):
        ds = []
        for u in range(n):
            if u == z:
                continue
            d = gs.pool.id(('d6', z, u))
            gs.cls.append([-gs.e(z, u), -ys[u], d])      # e(z,u) and y_u -> d
            ds.append(d)
        gs.card(ds, 2, 'atmost', cond=ys[z])
    gs.card(ys, K, 'atleast')
    extra = [[gs.e(0, u)] for u in range(1, 7)] + [[-gs.e(0, w)] for w in range(7, n)] + [[ys[0]]]
    extra += [[-ys[u]] for u in range(3, 7)]
    S = Solver(name='glucose4', bootstrap_with=gs.cls + extra)
    timer = threading.Timer(limit_s, lambda: S.interrupt())
    timer.start()
    t0 = time.time()
    r = S.solve_limited(expect_interrupt=True)
    timer.cancel()
    res = dict(K=K, status='TIMEOUT' if r is None else ('SAT' if r else 'UNSAT'), seconds=round(time.time() - t0, 1))
    if r:
        m = set(l for l in S.get_model() if l > 0)
        E = [(a, b) for a, b in itertools.combinations(range(n), 2) if gs.e(a, b) in m]
        G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
        assert sum(nx.triangles(G).values()) == 0
        assert all(len(set(G[a]) & set(G[b])) <= 2 for a, b in itertools.combinations(range(n), 2))
        assert all(5 <= d <= 6 for _, d in G.degree())
        res['degree6'] = sum(1 for _, d in G.degree() if d == 6)
        d5 = [v for v in G if G.degree(v) == 5]
        H = G.subgraph(d5)
        ok = True
        for comp in nx.connected_components(H):
            C = H.subgraph(comp)
            for blk in nx.biconnected_components(C):
                B = C.subgraph(blk)
                if not ((B.number_of_nodes() == 2) or (B.number_of_edges() == B.number_of_nodes()
                                                       and B.number_of_nodes() % 2 == 1)):
                    ok = False
        res['degree5_part_is_edge_oddcycle_forest'] = ok
        res['edges'] = E
    S.delete()
    log(json.dumps({k: v for k, v in res.items() if k != 'edges'}))
    return res


if __name__ == '__main__':
    out = [query(8, float(sys.argv[1]) if len(sys.argv) > 1 else 600), query(6, 240)]
    json.dump(out, open(OUTJ, 'w'), indent=1, default=str)
    log('wrote', OUTJ)
