# Is kappa(p) <= 5 for p = 29, 41 (c-stable residue fields of the 'surviving' Polymath-ladder fields)?  Glucose with time budget.
import sys, threading, time
from pysat.solvers import Glucose4
from sympy import legendre_symbol
def graph(p):
    nu = next(a for a in range(2, p) if legendre_symbol(a, p) == -1)
    circ = [(a, b) for a in range(p) for b in range(p) if (a*a - nu*b*b) % p == 1]
    idx = lambda a, b: a*p + b
    E = set()
    for a in range(p):
        for b in range(p):
            for (c, d) in circ:
                i, j = idx(a, b), idx((a+c) % p, (b+d) % p)
                E.add((min(i, j), max(i, j)))
    return p*p, sorted(E)
for p in [int(x) for x in sys.argv[1:-1]]:
    budget = float(sys.argv[-1])
    N, E = graph(p); k = 5
    s = Glucose4(); v = lambda i, c: i*k + c + 1
    for i in range(N): s.add_clause([v(i, c) for c in range(k)])
    for a, b in E:
        for c in range(k): s.add_clause([-v(a, c), -v(b, c)])
    s.add_clause([v(0, 0)])
    t = threading.Timer(budget, lambda: s.interrupt()); t.start(); t0 = time.time()
    r = s.solve_limited(expect_interrupt=True); t.cancel()
    print(f"p={p}: 5-colourable? {r}  ({time.time()-t0:.0f}s)", flush=True)
