#!/usr/bin/env python3
"""Independent check of a candidate unit-distance graph saved by minimize_npz.py (JSON: den, basis, X, Y integer
numerators over den in the basis of Q(sqrt2, sqrt3, sqrt5), edges).

Checks, from scratch:
  1. every listed edge has exact length 1, and every exact unit pair among the vertices is listed (the graph is the
     full induced unit-distance graph) -- all pairs, integer arithmetic, no floating point;
  2. not 4-colourable, by two different solvers (CaDiCaL 1.9.5 and Glucose 4 or Kissat, via pysat);
  3. optionally vertex-critical (every single-vertex deletion is 4-colourable);
  4. no exact sqrt11 can occur: the coordinates live in Q(sqrt2, sqrt3, sqrt5), so no Moser spindle (rotation
     (5 + i sqrt11)/6) is possible -- recorded, not computed.

Usage: python3 -u verify_graph.py candidate.json [critical=1]
"""
import sys
import json
import time

import numpy as np
from pysat.solvers import Solver

PR = (2, 3, 5)
DIM = 8
MULC = [[1] * DIM for _ in range(DIM)]
for a in range(DIM):
    for b in range(DIM):
        c = 1
        for i in range(3):
            if (a & b) >> i & 1:
                c *= PR[i]
        MULC[a][b] = c


def unit_matrix(X, Y, den):
    n = len(X)
    unit = np.zeros((n, n), dtype=bool)
    target0 = den * den
    for i in range(n):
        dx = X[i + 1:] - X[i]
        dy = Y[i + 1:] - Y[i]
        acc = np.zeros((n - i - 1, DIM), dtype=np.int64)
        for b1 in range(DIM):
            for b2 in range(DIM):
                acc[:, b1 ^ b2] += MULC[b1][b2] * (dx[:, b1] * dx[:, b2] + dy[:, b1] * dy[:, b2])
        ok = (acc[:, 0] == target0) & np.all(acc[:, 1:] == 0, axis=1)
        unit[i, i + 1:] = ok
    return unit


def colourable(n, E, active, name):
    k = 4
    s = Solver(name=name)
    x = lambda v, c: v * k + c + 1
    act = set(active)
    for v in act:
        s.add_clause([x(v, c) for c in range(k)])
    for i, j in E:
        if i in act and j in act:
            for c in range(k):
                s.add_clause([-x(i, c), -x(j, c)])
    r = s.solve()
    s.delete()
    return r


def main():
    d = json.load(open(sys.argv[1]))
    crit = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    X = np.array(d['X'], dtype=np.int64)
    Y = np.array(d['Y'], dtype=np.int64)
    den = int(d['den'])
    n = len(X)
    E = [tuple(e) for e in d['edges']]
    t0 = time.time()
    U = unit_matrix(X, Y, den)
    listed = set((min(i, j), max(i, j)) for i, j in E)
    exact = set(zip(*np.nonzero(U)))
    exact = set((int(i), int(j)) for i, j in exact)
    print(f"{sys.argv[1]}: n = {n}, listed edges {len(listed)}, exact unit pairs {len(exact)}, equal: {listed == exact}  "
          f"[{time.time() - t0:.1f}s]", flush=True)
    assert listed == exact
    pts = set(zip(map(tuple, X.tolist()), map(tuple, Y.tolist())))
    assert len(pts) == n, "repeated vertex"
    for name in ('cadical195', 'glucose4'):
        t1 = time.time()
        r = colourable(n, E, range(n), name)
        print(f"   {name}: 4-colourable = {r}  [{time.time() - t1:.1f}s]", flush=True)
        assert r is False
    if crit:
        t1 = time.time()
        bad = [v for v in range(n) if not colourable(n, E, [u for u in range(n) if u != v], 'cadical195')]
        print(f"   vertex-critical: {not bad}" + (f" (non-critical vertices: {bad[:10]})" if bad else "") +
              f"  [{time.time() - t1:.1f}s]", flush=True)
    print("   coordinates in Q(sqrt2, sqrt3, sqrt5): no Moser spindle possible (its rotation needs sqrt-11)")


if __name__ == '__main__':
    main()
