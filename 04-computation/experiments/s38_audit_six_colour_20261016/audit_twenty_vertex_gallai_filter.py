#!/usr/bin/env python3
"""S38 audit (independent referee), exploratory: a cheap necessary filter for 20-vertex failed exteriors.

Since no five-core has 17-19 vertices and the 16-vertex one (Clebsch) is degree-choosable, a failed exterior on 20
points must itself be vertex-minimal.  Then (THM-4661 2(a)-(b), audited correct) every degree-5 vertex is low, and
the low vertices induce a Gallai forest whose blocks are single edges or odd cycles; equivalently, there must be a set
S of degree-6 vertices (the low degree-6 vertices, lists of size 6) such that the subgraph induced by
(degree-5 vertices) + S is such a forest, and at least one degree-6 vertex is outside S (a high vertex).
This script samples 20-vertex five-cores with a degree-6 vertex (SAT enumeration as in
audit_twenty_vertex_cegar_fast.py) and tests the weakest form of the filter: the degree-5 vertices alone must induce
a forest of edges and odd cycles (S empty is the most permissive choice, since adding vertices only adds blocks).
Graphs passing the filter are handed to the CEGAR of audit_twenty_vertex_cegar_fast.py with a round cap.
"""
import itertools, json, os, sys, threading, time
import networkx as nx
from pysat.solvers import Solver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_nineteen_vertex_exhaustive import GraphSAT, colourable

HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'audit_twenty_vertex_gallai_filter.json')
MAXM = int(sys.argv[1]) if len(sys.argv) > 1 else 3000


def log(*a):
    print(*a, flush=True)


def is_edge_oddcycle_forest(H):
    for comp in nx.connected_components(H):
        C = H.subgraph(comp)
        if C.number_of_nodes() == 1:
            continue
        for blk in nx.biconnected_components(C):
            B = C.subgraph(blk)
            nv, ne = B.number_of_nodes(), B.number_of_edges()
            if nv == 2 and ne == 1:
                continue
            if ne == nv and nv % 2 == 1 and all(d == 2 for _, d in B.degree()):
                continue
            return False
    return True


n = 20
gs = GraphSAT(n)
for v in range(n):
    gs.card(gs.nbrs(v), 5, 'atleast')
    gs.card(gs.nbrs(v), 6, 'atmost')
extra = [[gs.e(0, u)] for u in range(1, 7)] + [[-gs.e(0, w)] for w in range(7, n)]
S = Solver(name='cadical195', bootstrap_with=gs.cls + extra)
timer = threading.Timer(600, lambda: S.interrupt())
timer.start()
allpairs = list(itertools.combinations(range(n), 2))
t0 = time.time()
models = 0
passed = []
hist = {}
deg6hist = {}
while models < MAXM:
    r = S.solve_limited(expect_interrupt=True)
    if not r:
        break
    models += 1
    m = set(l for l in S.get_model() if l > 0)
    E = [(a, b) for a, b in allpairs if gs.e(a, b) in m]
    S.add_clause([(-gs.e(a, b) if (a, b) in E else gs.e(a, b)) for a, b in allpairs])
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
    d5 = [v for v in G if G.degree(v) == 5]
    k6 = n - len(d5)
    deg6hist[k6] = deg6hist.get(k6, 0) + 1
    ok = is_edge_oddcycle_forest(G.subgraph(d5))
    hist[ok] = hist.get(ok, 0) + 1
    if ok:
        passed.append(G)
timer.cancel(); S.delete()
log(f'{models} labelled 20-vertex five-cores sampled ({time.time()-t0:.1f}s); number of degree-6 vertices: '
    f'{dict(sorted(deg6hist.items()))}; degree-5 part is a forest of edges/odd cycles: {hist}')
# also: lower bound on the number of degree-6 vertices needed by the low-degree count, for reference
out = dict(models=models, degree6_histogram=deg6hist, filter_histogram={str(k): v for k, v in hist.items()},
           passed=len(passed))
classes = []
for G in passed:
    if not any(nx.is_isomorphic(G, H) for H in classes):
        classes.append(G)
out['passed_classes'] = len(classes)
log(f'{len(classes)} non-isomorphic sampled graphs pass the filter')
json.dump(out, open(OUTJ, 'w'), indent=1, default=str)
log('wrote', OUTJ)
