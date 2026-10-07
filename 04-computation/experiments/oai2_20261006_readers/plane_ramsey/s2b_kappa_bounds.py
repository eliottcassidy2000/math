# Bounds for kappa(q) for larger q: upper bound via SAT (with a time budget), lower bound via clique & n/alpha (alpha by SAT/cardinality).
import sys, time, threading
from s2_kappa import kappa_graph, max_clique
from pysat.solvers import Cadical153, Glucose4
def try_colour(N, edges, k, clique, budget_s):
    v = lambda i, c: i*k + c + 1
    s = Glucose4()
    for i in range(N): s.add_clause([v(i, c) for c in range(k)])
    for a, b in edges:
        for c in range(k): s.add_clause([-v(a, c), -v(b, c)])
    for idx, x in enumerate(clique[:k]): s.add_clause([v(x, idx)])
    timer = threading.Timer(budget_s, lambda: s.interrupt()); timer.start()
    r = s.solve_limited(expect_interrupt=True); timer.cancel(); s.delete()
    return r   # True/False/None
for q in [int(a) for a in sys.argv[1:]]:
    t0 = time.time()
    N, edges, U = kappa_graph(q)
    cl = max_clique(N, edges)
    res = {}
    lo, hi = len(cl), None
    for k in range(len(cl), 12):
        r = try_colour(N, edges, k, cl, 90)
        res[k] = r
        if r is True: hi = k; break
        if r is False: lo = k + 1
    print(f"q={q}: |V|={N} deg={q+1} clique={len(cl)}  SAT results {res}  => kappa in [{lo},{hi}]  ({time.time()-t0:.0f}s)", flush=True)
