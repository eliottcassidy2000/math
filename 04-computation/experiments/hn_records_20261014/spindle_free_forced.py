#!/usr/bin/env python3
"""Is radius 2 forced in the spindle-free pieces?  For a piece B (containing 0), t = 4 is FORCED if every proper
4-colouring gives some vertex with |p|^2 = 4 the colour of 0 (SAT: forbid colour(0) on the radius-2 orbit).

Usage: HN_PRIMES=2,3,5 python3 -u spindle_free_forced.py
"""
import os
import time

from pysat.solvers import Solver

assert os.environ.get('HN_PRIMES') == '2,3,5'
from hn_udg import minkowski, cfloat, cnorm, C0, DIM
from fractions import Fraction as Fr
import spindle_free_build as sfb
from spindle_free_pieces import gens, clip


def forced(name, B, budget=50_000_000):
    t0 = time.time()
    pts = sorted(B, key=cfloat)
    X, Y, den = sfb.to_int(pts)
    I, J = sfb.edges_of(pts, X, Y, den)
    idx = {p: i for i, p in enumerate(pts)}
    c0 = idx[C0]
    four = (Fr(4),) + (Fr(0),) * (DIM - 1)
    orbit = [i for i, p in enumerate(pts) if cnorm(p) == four]
    k = 4
    s = Solver(name='cadical195')
    x = lambda v, c: v * k + c + 1
    for v in range(len(pts)):
        s.add_clause([x(v, c) for c in range(k)])
    for i, j in zip(I.tolist(), J.tolist()):
        for c in range(k):
            s.add_clause([-x(i, c), -x(j, c)])
    s.add_clause([x(c0, 0)])
    s.conf_budget(budget)
    col = s.solve_limited()
    if col is not True:
        print(f"{name:24s} |V|={len(pts):6d}: piece itself not shown 4-colourable ({col})", flush=True)
        return
    s.conf_budget(budget)
    r = s.solve_limited(assumptions=[-x(v, 0) for v in orbit])
    verdict = {False: 'FORCED', True: 'not forced', None: '? (budget)'}[r]
    print(f"{name:24s} |V|={len(pts):6d} |E|={len(I):7d} radius-2 orbit {len(orbit):4d}: {verdict}  "
          f"[{time.time() - t0:.1f}s]", flush=True)
    s.delete()


def main():
    phi0, phi1, psi, M1, M1c = gens()
    M11 = minkowski(M1, M1)
    M2 = clip(M11, 1.0)
    forced('M1+M1', M11)
    forced('M2+M1clock', minkowski(M2, M1c))
    forced('M3 = M2+M1', minkowski(M2, M1))


if __name__ == '__main__':
    main()
