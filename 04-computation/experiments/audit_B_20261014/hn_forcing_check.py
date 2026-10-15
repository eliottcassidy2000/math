#!/usr/bin/env python3
"""Audit B: is radius 2 forced in a piece B (0 in B)?  FORCED means every proper 4-colouring gives some vertex with
|p|^2 = 4 the colour of 0.  We fix colour(0) = 0 and forbid colour 0 on the radius-2 orbit; SAT = not forced (and the
colouring is checked as a certificate), UNSAT = forced (solver verdict).  Edges and the orbit come from the independent
rebuild (hn_spindle_free_rebuild_zeta24.py, exact norms in Q(zeta_24)(sqrt5)).
Usage: python3 -u hn_forcing_check.py piece.npz SOLVER"""
import sys, time
import numpy as np
from pysat.solvers import Solver
sys.argv += [] 
path, name = sys.argv[1], sys.argv[2]
f = np.load(path)
X, D, I, J = f['X'], int(f['D'][0]), f['I'], f['J']
sys.path.insert(0, '.')
from importlib import util
spec = util.spec_from_file_location('rb', 'hn_spindle_free_rebuild_zeta24.py')
rb = util.module_from_spec(spec); sys.argv_saved = sys.argv; sys.argv = ['x', '/dev/null', 'none']; spec.loader.exec_module(rb); sys.argv = sys.argv_saved
n = len(X)
origin = int(np.nonzero(np.all(X == 0, axis=1))[0][0])
orbit = np.nonzero(rb.norm_is(X, D, 4, 1))[0]
k = 4
s = Solver(name=name)
x = lambda v, c: v * k + c + 1
for v in range(n):
    s.add_clause([x(v, c) for c in range(k)])
for a, b in zip(I.tolist(), J.tolist()):
    for c in range(k):
        s.add_clause([-x(a, c), -x(b, c)])
s.add_clause([x(origin, 0)])
for v in orbit.tolist():
    s.add_clause([-x(v, 0)])
t0 = time.time()
r = s.solve()
print(f"{path.split('/')[-1]}: n = {n}, |E| = {len(I)}, radius-2 orbit {len(orbit)}; colouring avoiding colour(0) on the orbit exists: {r} "
      f"-> {'NOT forced' if r else 'FORCED'} ({name}, {time.time() - t0:.1f}s)", flush=True)
if r:
    m = set(l for l in s.get_model() if l > 0)
    col = np.array([[c for c in range(k) if x(v, c) in m][0] for v in range(n)])
    print(f"   certificate: monochromatic edges {int(np.sum(col[I] == col[J]))}, orbit vertices with colour of 0: "
          f"{int(np.sum(col[orbit] == col[origin]))}", flush=True)
