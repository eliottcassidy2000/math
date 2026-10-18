#!/usr/bin/env python3
"""Decide k-colourability of G_q = Cay(F_(q^2), mu_(q+1)) with CaDiCaL (pysat), symmetry broken by fixing the colours of
one triangle {0, 1, w} (w, w - 1 in mu). Writes a DRAT proof when UNSAT (proof file name printed).
Usage: python3 kappa_unsat.py q k"""
import sys, time
from pysat.solvers import Solver
q, k = int(sys.argv[1]), int(sys.argv[2])
r = next(x for x in range(2, q) if pow(x, (q - 1) // 2, q) == q - 1)
mu = [(a, b) for a in range(q) for b in range(q) if (a * a - r * b * b) % q == 1]
idx = lambda a, b: a * q + b
edges = set()
for a in range(q):
    for b in range(q):
        for (x, y) in mu:
            u, v = idx(a, b), idx((a + x) % q, (b + y) % q)
            if u < v: edges.add((u, v))
var = lambda v, c: v * k + c + 1
cl = [[var(v, c) for c in range(k)] for v in range(q * q)]
for (u, v) in edges:
    for c in range(k): cl.append([-var(u, c), -var(v, c)])
tri = next((x, y) for (x, y) in mu if ((x - 1) % q, y) in mu)
cl += [[var(idx(0, 0), 0)], [var(idx(1, 0), 1)], [var(idx(*tri), 2)]]
t0 = time.time()
with Solver(name='cadical153', bootstrap_with=cl, with_proof=False) as s:
    res = s.solve()
    print('q = %d, k = %d: %s in %.1f s (%d vertices, %d edges, %d clauses)' % (q, k, 'SAT (colourable)' if res else 'UNSAT (not colourable)', time.time() - t0, q * q, len(edges), len(cl)), flush=True)
    if res:
        model = set(l for l in s.get_model() if l > 0)
        col = [next(c for c in range(k) if var(v, c) in model) for v in range(q * q)]
        bad = sum(1 for (u, v) in edges if col[u] == col[v]); print('  verified proper:', bad == 0)
