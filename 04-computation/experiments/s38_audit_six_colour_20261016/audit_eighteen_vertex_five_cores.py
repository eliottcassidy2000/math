#!/usr/bin/env python3
"""S38 audit (independent referee): do 18-vertex 'five-cores' of the THM-4661 type (triangle-free, all codegrees
<= 2, minimum degree >= 5) exist combinatorially?  (The note, section 6, says a five-core of 18 or 19 vertices
'exists combinatorially for all the counting cares'; the 19-vertex case is refuted by
audit_nineteen_vertex_exhaustive.py / audit_nineteen_vertex_backtrack.py.)

On 18 vertices a degree-6 vertex would need 1 + 6 + 12 = 19 vertices, so such a graph is 5-regular.  Root it at any
vertex z: N(z) = {u_1..u_5}; the other twelve vertices carry 2, 1 or 0 neighbours in N(z) with 2a + b = 20,
a + b + c = 12, and the a two-neighbour vertices sit on distinct pairs of {1..5}.  Up to S_5 the cases are
  A  (a,b,c) = (10,0,2);  B  (9,2,1), missing pair {1,2};  C1 (8,4,0), missing {1,2},{3,4};  C2 (8,4,0), missing
  {1,2},{1,3}.
Each case is a SAT instance on the 66 remaining pairs (exact degree 5, triangle-free, codegree <= 2).  For every
graph found, the script checks it independently (networkx), tests three-colourability and runs the CEGAR search for a
list assignment ([6] minus one colour per vertex) that kills every colouring.
"""
import itertools, json, os, sys, time, threading
import networkx as nx
from pysat.solvers import Solver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_nineteen_vertex_exhaustive import GraphSAT, colourable, blocking_list_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'audit_eighteen_vertex_five_cores.json')


def log(*a):
    print(*a, flush=True)


def case_structure(name):
    allp = list(itertools.combinations(range(1, 6), 2))
    if name == 'A':
        miss, singles, far = [], [], 2
    elif name == 'B':
        miss, singles, far = [(1, 2)], [1, 2], 1
    elif name == 'C1':
        miss, singles, far = [(1, 2), (3, 4)], [1, 2, 3, 4], 0
    else:
        miss, singles, far = [(1, 2), (1, 3)], [1, 1, 2, 3], 0
    pairs = [p for p in allp if p not in miss]
    att = [set(p) for p in pairs] + [{s} for s in singles] + [set() for _ in range(far)]
    assert len(att) == 12
    for u in range(1, 6):
        assert sum(1 for a in att if u in a) == 4
    return att


def solve_case(name, limit_s=600, max_models=50):
    n = 18
    att = case_structure(name)
    others = list(range(6, 18))
    gs = GraphSAT(n)
    fixed = []
    for i, j in itertools.combinations(range(n), 2):
        if i == 0:
            fixed.append([gs.e(i, j)] if 1 <= j <= 5 else [-gs.e(i, j)])
        elif 1 <= i <= 5 and 1 <= j <= 5:
            fixed.append([-gs.e(i, j)])
        elif 1 <= i <= 5 and j >= 6:
            fixed.append([gs.e(i, j)] if i in att[j - 6] else [-gs.e(i, j)])
    for v in others:
        gs.card(gs.nbrs(v), 5, 'equals')
    free = list(itertools.combinations(others, 2))
    S = Solver(name='glucose4', bootstrap_with=gs.cls + fixed)   # interruptible (cadical195 ignores interrupt here)
    timer = threading.Timer(limit_s, lambda: S.interrupt())
    timer.start()
    models = []
    status = None
    t0 = time.time()
    while len(models) < max_models:
        r = S.solve_limited(expect_interrupt=True)
        if r is None:
            status = 'TIMEOUT'
            break
        if not r:
            status = 'COMPLETE' if models else 'UNSAT'
            break
        m = set(l for l in S.get_model() if l > 0)
        es = [(a, b) for a, b in free if gs.e(a, b) in m]
        models.append(es)
        S.add_clause([(-gs.e(a, b) if (a, b) in es else gs.e(a, b)) for a, b in free])
    else:
        status = f'STOPPED_AFTER_{max_models}'
    timer.cancel()
    S.delete()
    graphs = []
    for es in models:
        E = [(0, u) for u in range(1, 6)] + [(u, v) for v in others for u in att[v - 6]] + es
        G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
        assert sum(nx.triangles(G).values()) == 0 and all(d == 5 for _, d in G.degree())
        assert all(len(set(G[a]) & set(G[b])) <= 2 for a, b in itertools.combinations(range(n), 2))
        graphs.append(G)
    classes = []
    for G in graphs:
        if not any(nx.is_isomorphic(G, H) for H in classes):
            classes.append(G)
    log(f'case {name}: status {status}, {len(models)} labelled completions found, {len(classes)} isomorphism classes '
        f'({time.time()-t0:.1f}s)')
    return status, classes


out = {}
allclasses = []
for name in ('A', 'B', 'C1', 'C2'):
    st, cl = solve_case(name)
    out[name] = dict(status=st, classes=len(cl))
    for G in cl:
        if not any(nx.is_isomorphic(G, H) for H in allclasses):
            allclasses.append(G)
log(f'total isomorphism classes found: {len(allclasses)}')
out['tests'] = []
for i, G in enumerate(allclasses[:10]):
    t0 = time.time()
    c3 = colourable(G, 3)
    f, it = blocking_list_assignment(G)
    rec = dict(index=i, three_colourable=c3 is not None, blocking_list_assignment_exists=f is not None,
               cegar_rounds=it, girth=nx.girth(G) if hasattr(nx, 'girth') else None,
               edges=sorted(tuple(sorted(e)) for e in G.edges()), seconds=round(time.time() - t0, 1))
    out['tests'].append(rec)
    log(f'18-vertex five-core #{i}: 3-colourable {rec["three_colourable"]}, blocking list assignment '
        f'{rec["blocking_list_assignment_exists"]} (CEGAR {it} rounds)')
json.dump(out, open(OUTJ, 'w'), indent=1)
log('wrote', OUTJ)
