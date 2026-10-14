#!/usr/bin/env python3
"""Hybrid type-M graphs: Parts' rigid piece L727 (Q(sqrt3, sqrt11)) hinged by rho = (7 + i sqrt15)/8 to small
15-degree 'clock' pieces (Q(sqrt2, sqrt3)) from the spindle-free family.  Field Q(sqrt2, sqrt3, sqrt5, sqrt11).

For each clock piece C: is radius 2 forced in C (every 4-colouring repeats the centre's colour on |p| = 2)?
Is L727 u rho C 4-colourable?  (CaDiCaL 1.9.5, triangle fixed.)

Usage: HN_PRIMES=2,3,5,11 python3 -u hybrid_probe.py
"""
import os
import time
from fractions import Fraction as Fr

from pysat.solvers import Solver

assert os.environ.get('HN_PRIMES') == '2,3,5,11'
from hn_udg import (el, cmul, cpow, union, minkowski, rotate, unit_graph, cfloat, cnorm, C0, DIM, RHO,
                    triangles_through, H, H_rot, I_OVER_SQ3)


def log(*a):
    print(*a, flush=True)


def clip(S, r2max):
    out = []
    for p in S:
        x, y = cfloat(p)
        r2 = x * x + y * y
        if r2 < r2max + 1e-9:
            out.append(p)
    return frozenset(out)


def colour_test(pts, E, extra_unit_clauses=None, budget=20_000_000):
    k = 4
    s = Solver(name='cadical195')
    x = lambda v, c: v * k + c + 1
    for v in range(len(pts)):
        s.add_clause([x(v, c) for c in range(k)])
    for i, j in E:
        for c in range(k):
            s.add_clause([-x(i, c), -x(j, c)])
    tri = triangles_through(len(pts), E, limit=1)
    if tri:
        for c, v in enumerate(tri[0]):
            s.add_clause([x(v, c)])
    s.conf_budget(budget)
    r = s.solve_limited()
    s.delete()
    return r


def forced_r2(B):
    pts, xy, E = unit_graph(B)
    idx = {p: i for i, p in enumerate(pts)}
    c0 = idx[C0]
    four = (Fr(4),) + (Fr(0),) * (DIM - 1)
    orbit = [i for i, p in enumerate(pts) if cnorm(p) == four]
    k = 4
    s = Solver(name='cadical195')
    x = lambda v, c: v * k + c + 1
    for v in range(len(pts)):
        s.add_clause([x(v, c) for c in range(k)])
    for i, j in E:
        for c in range(k):
            s.add_clause([-x(i, c), -x(j, c)])
    s.add_clause([x(c0, 0)])
    col = s.solve()
    r = s.solve(assumptions=[-x(v, 0) for v in orbit]) if col else None
    s.delete()
    return len(pts), len(E), len(orbit), col, (None if r is None else not r)


def main():
    phi0 = (el(r6=Fr(1, 4), r2=Fr(1, 4)), el(r6=Fr(1, 4), r2=Fr(-1, 4)))
    phi1 = (el(r6=Fr(1, 3)), el(r3=Fr(1, 3)))
    clock = frozenset([C0] + [cpow(phi0, a) for a in range(24)])
    clock6 = frozenset([C0] + [cpow(phi0, 4 * a) for a in range(6)])            # = H, the hexagonal wheel
    M1 = frozenset([C0] + [cmul(cpow(phi0, a), cpow(phi1, b)) for a in range(24) for b in (-1, 0, 1)])
    H1 = H_rot([-1, 0, 1]); H2 = H_rot([-2, -1, 0, 1, 2])
    V37 = union(H1, rotate(H1, I_OVER_SQ3)); V49 = union(H2, rotate(H1, I_OVER_SQ3))
    L = minkowski(V49, V37)
    assert len(L) == 727
    pieces = {
        'clock+clock': minkowski(clock, clock),
        'clip1(clock+clock)+clock': minkowski(clip(minkowski(clock, clock), 1.0), clock),
        'clock+clock+clock(r<=2)': clip(minkowski(minkowski(clock, clock), clock), 4.0),
        'clip1(M1+M1)+clock': minkowski(clip(minkowski(M1, M1), 1.0), clock),
    }
    for name, C in pieces.items():
        t0 = time.time()
        n, m, orb, col, forced = forced_r2(C)
        log(f"{name:28s} |V|={n:6d} |E|={m:6d} radius-2 orbit {orb:3d}: 4-colourable={col}, radius 2 forced={forced}  [{time.time() - t0:.1f}s]")
        if forced:
            t0 = time.time()
            G = union(L, rotate(C, RHO))
            pts, xy, E = unit_graph(G)
            r = colour_test(pts, E)
            verdict = {True: '4-colourable', False: 'NOT 4-colourable (5-chromatic)', None: '? (budget)'}[r]
            log(f"   L727 u rho [{name}]: |V|={len(pts)} |E|={len(E)}  {verdict}  [{time.time() - t0:.1f}s]")


if __name__ == '__main__':
    main()
