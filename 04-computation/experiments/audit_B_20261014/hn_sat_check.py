#!/usr/bin/env python3
"""Audit B: 4-colourability of a graph from the independent rebuild (npz with edges I, J), with a chosen pysat solver.

Encoding: x(v,c) for 4 colours; one at-least-one clause per vertex, one conflict clause per edge and colour; symmetry
breaking by fixing the colours of the ORIGIN and of a triangle through it (origin -> 0, two adjacent unit neighbours
-> 1, 2), which is valid for any proper colouring up to a permutation of colours.  If SAT, the model is decoded and
checked to be a proper colouring of every edge (a certificate independent of the solver).

Usage: python3 -u hn_sat_check.py graph.npz SOLVER [timeout_s]
"""
import sys
import time
import threading

import numpy as np
from pysat.solvers import Solver

path, name = sys.argv[1], sys.argv[2]
tmo = float(sys.argv[3]) if len(sys.argv) > 3 else 0
f = np.load(path)
X, I, J = f['X'], f['I'].astype(np.int64), f['J'].astype(np.int64)
n = len(X)
k = 4
origin = int(np.nonzero(np.all(X == 0, axis=1))[0][0])
adj = [set() for _ in range(n)]
for a, b in zip(I.tolist(), J.tolist()):
    adj[a].add(b)
    adj[b].add(a)
tri = None
for u in sorted(adj[origin]):
    common = adj[origin] & adj[u]
    if common:
        tri = (origin, u, min(common))
        break
assert tri is not None
t0 = time.time()
s = Solver(name=name)
x = lambda v, c: v * k + c + 1
for v in range(n):
    s.add_clause([x(v, c) for c in range(k)])
for a, b in zip(I.tolist(), J.tolist()):
    for c in range(k):
        s.add_clause([-x(a, c), -x(b, c)])
for c, v in enumerate(tri):
    s.add_clause([x(v, c)])
print(f"{path}: n = {n}, |E| = {len(I)}, solver {name}, fixed triangle {tri}", flush=True)
if tmo > 0 and hasattr(s, 'interrupt'):
    timer = threading.Timer(tmo, lambda: s.interrupt())
    timer.start()
    try:
        r = s.solve_limited(expect_interrupt=True)
    except NotImplementedError:
        r = s.solve()
    timer.cancel()
else:
    r = s.solve()
dt = time.time() - t0
print(f"4-colourable: {r}  [{dt:.1f}s]", flush=True)
if r:
    model = set(l for l in s.get_model() if l > 0)
    col = np.full(n, -1)
    for v in range(n):
        cs = [c for c in range(k) if x(v, c) in model]
        col[v] = cs[0]
    bad = int(np.sum(col[I] == col[J]))
    print(f"certificate check: colouring of all {n} vertices, monochromatic edges = {bad}", flush=True)
    np.save(path.replace('.npz', f'_colouring_{name}.npy'), col)
s.delete()
