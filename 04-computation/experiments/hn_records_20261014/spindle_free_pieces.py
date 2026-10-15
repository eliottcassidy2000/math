#!/usr/bin/env python3
"""Smaller spindle-free type-M candidates inside the Voronov-Neopryatnaya-Dergachev construction.

W = M3 u psi M3 (64,513 vertices) is a hinge at angular gap arccos(7/8): of its 120 new edges, 72 join z in M3 and psi w
with |z| = |w| = 2 and 48 join |z|^2 = 2/3 with |psi w|^2 = 8/3, 24 each way (audit B; corrected docstring).
Here we test smaller pieces A u psi B that keep radius-2 orbits, e.g. B = M1 + M1 (unclipped, radius <= 2), and report
size, exact edges and the 4-colourability verdict (CaDiCaL 1.9.5, a triangle fixed, conflict budget).

Usage: HN_PRIMES=2,3,5 python3 -u spindle_free_pieces.py
"""
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np
from pysat.solvers import Solver

assert os.environ.get('HN_PRIMES') == '2,3,5'
from hn_udg import el, cmul, cpow, cconj, union, minkowski, cfloat, cnorm, ONE, C0, is_unit, triangles_through
import spindle_free_build as sfb


def log(*a):
    print(*a, flush=True)


def gens():
    phi0 = (el(r6=Fr(1, 4), r2=Fr(1, 4)), el(r6=Fr(1, 4), r2=Fr(-1, 4)))
    phi1 = (el(r6=Fr(1, 3)), el(r3=Fr(1, 3)))
    psi = cmul(cconj(phi0), (el(r1=Fr(7, 8)), el(r15=Fr(1, 8))))
    clock = [cpow(phi0, a) for a in range(24)]
    M1 = frozenset([C0] + [cmul(c, cpow(phi1, b)) for c in clock for b in (-1, 0, 1)])
    M1clock = frozenset([C0] + clock)
    return phi0, phi1, psi, M1, M1clock


def clip(S, r2max):
    out = []
    for p in S:
        x, y = cfloat(p)
        r2 = x * x + y * y
        if abs(r2 - r2max) < 1e-9:
            out.append(p)
        elif r2 < r2max:
            out.append(p)
    return frozenset(out)


def test(name, A, B, psi, budget=20_000_000):
    t0 = time.time()
    G = union(A, frozenset(cmul(p, psi) for p in B))
    pts = sorted(G, key=cfloat)
    X, Y, den = sfb.to_int(pts)
    I, J = sfb.edges_of(pts, X, Y, den)
    E = list(zip(I.tolist(), J.tolist()))
    tris = triangles_through(len(pts), E, limit=1)
    k = 4
    s = Solver(name='cadical195')
    x = lambda v, c: v * k + c + 1
    for v in range(len(pts)):
        s.add_clause([x(v, c) for c in range(k)])
    for i, j in E:
        for c in range(k):
            s.add_clause([-x(i, c), -x(j, c)])
    for c, v in enumerate(tris[0]):
        s.add_clause([x(v, c)])
    s.conf_budget(budget)
    t1 = time.time()
    r = s.solve_limited()
    s.delete()
    verdict = {True: '4-colourable', False: 'NOT 4-colourable', None: '? (budget)'}[r]
    log(f"{name:40s} |V|={len(pts):6d} |E|={len(E):7d}  {verdict}  [build {t1 - t0:.1f}s, solve {time.time() - t1:.1f}s]")
    return r


def main():
    phi0, phi1, psi, M1, M1c = gens()
    M11 = minkowski(M1, M1)                      # radius <= 2, unclipped
    M2 = clip(M11, 1.0)
    M3 = minkowski(M2, M1)
    log(f"|M1| = {len(M1)}, |M1+M1| = {len(M11)}, |M2| = {len(M2)}, |M3| = {len(M3)}")
    test('(M1+M1) u psi (M1+M1)', M11, M11, psi)
    test('M3 u psi (M1+M1)', M3, M11, psi)
    M2c = clip(minkowski(M1c, M1c), 1.0)
    test('(M2+M1clock) u psi (M2+M1clock)', minkowski(M2, M1c), minkowski(M2, M1c), psi)
    test('(M2c+M1) u psi (M2c+M1)', minkowski(M2c, M1), minkowski(M2c, M1), psi)


if __name__ == '__main__':
    main()
