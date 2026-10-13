"""kappa(q) = chi(G_q) for the unit-circle Cayley graphs G_q = Cay(F_(q^2), mu_(q+1)): k-colourability by SAT with
symmetry breaking (colour(0) = 0, colour(1) = 1 since Aff(q) is arc-transitive, THM-4624) and lex-leader
constraints on the colour classes (in a fixed vertex order, colour c >= 2 may appear only after colour c-1).
The lex-leader constraints are what the mac-mini-2026-10-12-logic4 runs (kappa_sat.py, one hour undecided
at q = 13, 17, k = 5) lacked; with them G_13 is refuted 4-colourable in about a second.
Usage (from the repository root): python3 -B <script> q k [solver]      e.g. 13 4, 17 4, 13 5, 17 5, 17 6
Default solver: cadical195 (pysat).  Records: 05-knowledge/results/gq_kappa_sat_20261012.out
"""
import os, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
D = ROOT/'04-computation/experiments/logic4_20261012'
sys.path.insert(0, str(D)); os.chdir(D)
exec(open('gq_automorphisms.py').read().split('qs=[q for q')[0])  # defines gq(q)
from pysat.solvers import Solver

q = int(sys.argv[1]); k = int(sys.argv[2]); solver_name = sys.argv[3] if len(sys.argv) > 3 else 'cadical195'
adj, p, e, mu, elems, idx, add = gq(q); n = len(adj); one = idx[tuple([1]+[0]*(2*e-1))]
var = lambda v, col: v*k+col+1
cl = []
for v in range(n):
    cl.append([var(v, col) for col in range(k)])
    for col in range(k):
        for col2 in range(col+1, k):
            cl.append([-var(v, col), -var(v, col2)])
for v in range(n):
    for u_ in adj[v]:
        if u_ > v:
            for col in range(k):
                cl.append([-var(v, col), -var(u_, col)])
cl.append([var(0, 0)]); cl.append([var(one, 1)])
order = [0, one]+[v for v in range(n) if v not in (0, one)]
for i, v in enumerate(order):
    for col in range(2, k):
        cl.append([-var(v, col)]+[var(order[j], col-1) for j in range(i)])
t = time.time()
with Solver(name=solver_name, bootstrap_with=cl) as s:
    res = s.solve(); dt = time.time()-t
    print('q=%d n=%d k=%d solver=%s: %s (%.1fs)' % (q, n, k, solver_name, 'SAT (k-colourable)' if res else 'UNSAT (not k-colourable)', dt), flush=True)
    if res:
        mdl = s.get_model(); col = {v: [c_ for c_ in range(k) if mdl[var(v, c_)-1] > 0][0] for v in range(n)}
        assert all(col[v] != col[u_] for v in range(n) for u_ in adj[v]); print('colouring verified', flush=True)
        print('colouring:', [col[v] for v in range(n)], flush=True)
