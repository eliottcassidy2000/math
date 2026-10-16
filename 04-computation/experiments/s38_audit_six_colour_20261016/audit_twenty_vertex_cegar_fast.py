#!/usr/bin/env python3
"""S38 audit (independent referee), exploratory follow-up to audit_twenty_vertex_probe.py: a faster CEGAR for the
abstract question 'is there a list assignment L(v) = [6] minus one colour per vertex under which G has no proper
L-colouring?', applied to 20-vertex five-cores (triangle-free, codegree <= 2, degrees 5..6, a degree-6 vertex).

CEGAR: persistent colouring solver queried under assumptions (not x_{v, f(v)}); persistent f-solver with exactly one
forbidden colour per vertex, f(v_0) = 0 (colour symmetry), and for every colouring found all 720 colour permutations
added as hitting clauses.  The f-solver becoming UNSAT proves that no list assignment of this shape kills G
(G is then extendable for every old six-colouring pattern, in the abstract model); a SAT colouring query that fails
returns a killing assignment.  Graphs: up to MAX_GRAPHS distinct five-cores from a SAT enumeration with a degree-6
root; also the 3-colourability of each (THM-4661 2(c): a killing exterior is not 3-colourable).
Exploratory and FINITE-EXACT per graph; it samples, it does not enumerate all 20-vertex five-cores.
"""
import itertools, json, os, sys, threading, time
import networkx as nx
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_nineteen_vertex_exhaustive import GraphSAT, colourable

HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'audit_twenty_vertex_cegar_fast.json')
MAX_GRAPHS = int(sys.argv[1]) if len(sys.argv) > 1 else 40
CAP = float(sys.argv[2]) if len(sys.argv) > 2 else 150.0


def log(*a):
    print(*a, flush=True)


def blocking_fast(G, q=6, max_rounds=100000, time_cap=None):
    nodes = sorted(G.nodes())
    n = len(nodes)
    ix = {v: i for i, v in enumerate(nodes)}
    X = lambda i, c: i * q + c + 1
    F = lambda i, c: n * q + i * q + c + 1
    Sc = Solver(name='glucose4')   # interruptible
    for i in range(n):
        Sc.add_clause([X(i, c) for c in range(q)])
    for a, b in G.edges():
        for c in range(q):
            Sc.add_clause([-X(ix[a], c), -X(ix[b], c)])
    Sf = Solver(name='glucose4')   # interruptible
    for i in range(n):
        Sf.add_clause([F(i, c) for c in range(q)])
        for c1, c2 in itertools.combinations(range(q), 2):
            Sf.add_clause([-F(i, c1), -F(i, c2)])
    Sf.add_clause([F(0, 0)])
    perms = list(itertools.permutations(range(q)))
    t_start = time.time()
    timer = None
    if time_cap is not None:
        timer = threading.Timer(time_cap, lambda: (Sf.interrupt(), Sc.interrupt()))
        timer.start()
    for r in range(1, max_rounds + 1):
        rf = Sf.solve_limited(expect_interrupt=True)
        if rf is None or (time_cap is not None and time.time() - t_start > time_cap):
            if timer: timer.cancel()
            Sc.delete(); Sf.delete()
            return 'UNDECIDED', r
        if not rf:
            if timer: timer.cancel()
            Sc.delete(); Sf.delete()
            return None, r
        m = set(l for l in Sf.get_model() if l > 0)
        f = [next(c for c in range(q) if F(i, c) in m) for i in range(n)]
        rc = Sc.solve_limited(assumptions=[-X(i, f[i]) for i in range(n)], expect_interrupt=True)
        if rc is None:
            if timer: timer.cancel()
            Sc.delete(); Sf.delete()
            return 'UNDECIDED', r
        if not rc:
            if timer: timer.cancel()
            Sc.delete(); Sf.delete()
            return {nodes[i]: f[i] for i in range(n)}, r
        mc = set(l for l in Sc.get_model() if l > 0)
        col = [next(c for c in range(q) if X(i, c) in mc) for i in range(n)]
        for p in perms:
            Sf.add_clause([F(i, p[col[i]]) for i in range(n)])
    raise RuntimeError('no convergence')


def five_cores_20(max_graphs, limit_s=300):
    n = 20
    gs = GraphSAT(n)
    for v in range(n):
        gs.card(gs.nbrs(v), 5, 'atleast')
        gs.card(gs.nbrs(v), 6, 'atmost')
    extra = [[gs.e(0, u)] for u in range(1, 7)] + [[-gs.e(0, w)] for w in range(7, n)]
    S = Solver(name='cadical195', bootstrap_with=gs.cls + extra)
    timer = threading.Timer(limit_s, lambda: S.interrupt())
    timer.start()
    allpairs = list(itertools.combinations(range(n), 2))
    classes = []
    tried = 0
    while len(classes) < max_graphs and tried < 20 * max_graphs:
        r = S.solve_limited(expect_interrupt=True)
        if not r:
            break
        tried += 1
        m = set(l for l in S.get_model() if l > 0)
        E = [(a, b) for a, b in allpairs if gs.e(a, b) in m]
        S.add_clause([(-gs.e(a, b) if (a, b) in E else gs.e(a, b)) for a, b in allpairs])
        G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
        if not any(nx.is_isomorphic(G, H) for H in classes):
            classes.append(G)
    timer.cancel(); S.delete()
    return classes, tried


# positive control: K6 with every list [6] minus colour 0 must be killed
K6 = nx.complete_graph(6)
fK6, rK6 = blocking_fast(K6)
log(f'control K6: killing assignment found {fK6 is not None and fK6 != "UNDECIDED"} in {rK6} rounds')
assert fK6 is not None and fK6 != 'UNDECIDED'
t0 = time.time()
classes, tried = five_cores_20(MAX_GRAPHS)
log(f'{len(classes)} non-isomorphic 20-vertex five-cores with a degree-6 vertex from {tried} SAT models '
    f'({time.time()-t0:.1f}s)')
recs = []
for k, G in enumerate(classes):
    assert sum(nx.triangles(G).values()) == 0
    assert all(len(set(G[a]) & set(G[b])) <= 2 for a, b in itertools.combinations(G.nodes(), 2))
    t1 = time.time()
    c3 = colourable(G, 3)
    f, rounds = blocking_fast(G, time_cap=CAP)
    degs = sorted(d for _, d in G.degree())
    rec = dict(index=k, deg5=degs.count(5), deg6=degs.count(6), three_colourable=c3 is not None,
               killing_assignment=('UNDECIDED within cap' if f == 'UNDECIDED' else f is not None), cegar_rounds=rounds, seconds=round(time.time() - t1, 1),
               edges=sorted(G.edges()), assignment=(None if f == 'UNDECIDED' else f))
    recs.append(rec)
    log(f'#{k}: degrees 5^{rec["deg5"]} 6^{rec["deg6"]}, 3-colourable {rec["three_colourable"]}, '
        f'killing list assignment {rec["killing_assignment"]} (CEGAR {rounds} rounds, {rec["seconds"]}s)')
json.dump(dict(classes=len(classes), sat_models=tried, results=recs), open(OUTJ, 'w'), indent=1, default=str)
log('wrote', OUTJ)
