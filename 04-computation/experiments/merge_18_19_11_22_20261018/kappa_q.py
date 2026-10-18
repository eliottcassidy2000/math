#!/usr/bin/env python3
"""kappa(q) = chi(Cay(F_(q^2), mu_(q+1))) for odd prime q (mac-mini-2026-10-18-sixcolour).
F_(q^2) = F_q[t]/(t^2 - r), r a non-residue; mu_(q+1) = {a + b t : a^2 - r b^2 = 1}.
Prints the Hoffman lower bound 1 - lmax/lmin (eigenvalues = additive-character sums over mu_(q+1)) and tries
k-colourings with a SAT solver (pysat; symmetry breaking: colour(0) = 0, colour(1) = 1, and the third vertex of a
triangle through 0 and 1 if one exists gets colour 2).
Usage: python3 kappa_q.py q k1 [k2 ...] [--time T]"""
import sys, time, cmath
from pysat.solvers import Solver
q = int(sys.argv[1]); ks = [int(a) for a in sys.argv[2:] if not a.startswith('--')]
T = float(sys.argv[sys.argv.index('--time') + 1]) if '--time' in sys.argv else 600.0
r = next(x for x in range(2, q) if pow(x, (q - 1) // 2, q) == q - 1)
mu = [(a, b) for a in range(q) for b in range(q) if (a * a - r * b * b) % q == 1]
assert len(mu) == q + 1
idx = lambda a, b: a * q + b
# Hoffman bound: eigenvalue for additive character psi_(c,d)(a + b t) = e((c a + d b)/q)
lam = []
for c in range(q):
    for d in range(q):
        s = sum(cmath.exp(2j * cmath.pi * ((c * a + d * b) % q) / q) for a, b in mu)
        lam.append(s.real)
lmax, lmin = max(lam), min(lam)
print('q = %d: |mu| = %d, eigenvalues in [%.4f, %.4f], Hoffman bound chi >= %.4f' % (q, len(mu), lmin, lmax, 1 - lmax / lmin), flush=True)
edges = set()
for a in range(q):
    for b in range(q):
        for (x, y) in mu:
            u, v = idx(a, b), idx((a + x) % q, (b + y) % q)
            if u < v: edges.add((u, v))
n = q * q
tri = next(((x, y) for (x, y) in mu if ((x - 1) % q, y) in mu), None)
for k in ks:
    var = lambda v, c: v * k + c + 1
    cl = [[var(v, c) for c in range(k)] for v in range(n)]
    for (u, v) in edges:
        for c in range(k): cl.append([-var(u, c), -var(v, c)])
    cl.append([var(idx(0, 0), 0)]); cl.append([var(idx(1, 0), 1)])
    if tri is not None and k >= 3: cl.append([var(idx(*tri), 2)])
    t0 = time.time()
    with Solver(name='cadical153', bootstrap_with=cl) as s:
        import threading
        timer = threading.Timer(T, lambda: s.interrupt()); timer.start()
        res = s.solve_limited(expect_interrupt=True); timer.cancel()
        dt = time.time() - t0
        if res is True:
            model = set(l for l in s.get_model() if l > 0)
            col = [next(c for c in range(k) if var(v, c) in model) for v in range(n)]
            bad = sum(1 for (u, v) in edges if col[u] == col[v])
            print('  k = %d: COLOURABLE (%.1f s); verified proper: %s' % (k, dt, bad == 0), flush=True)
            with open('kappa_%d_colouring_k%d.txt' % (q, k), 'w') as fo:
                fo.write('# q=%d r=%d vertex index a*q+b for a + b t, t^2 = r; colours 0..%d\n' % (q, r, k - 1))
                fo.write(' '.join(map(str, col)) + '\n')
        elif res is False:
            print('  k = %d: NOT colourable (UNSAT, %.1f s)' % (k, dt), flush=True)
        else:
            print('  k = %d: undecided after %.0f s' % (k, dt), flush=True)
