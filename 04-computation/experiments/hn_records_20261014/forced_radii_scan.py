#!/usr/bin/env python3
"""Forced monochromatic radii across many type-M base graphs (Minkowski sums of Parts' pieces V25, V31, V37, V49).

For a base graph B containing 0, a squared radius t is FORCED if every proper 4-colouring of the unit-distance graph on
B gives some vertex at squared distance t from 0 the colour of 0.  A type-M rotation w_t (cos = 1 - 1/(2t)) joins
p and w_t p at unit distance for |p|^2 = t, so radii forced in two pieces are where a rotation can hinge.
For each forced t we print the orbit size and the squarefree part of 4t - 1 (the rotation lives in Q(sqrt(-(4t-1)))).

Usage: python3 -u forced_radii_scan.py
"""
import time
from fractions import Fraction as Fr

import sympy
from pysat.solvers import Solver

from hn_udg import (H, H_rot, union, rotate, minkowski, unit_graph, cnorm, DIM, C0, I_OVER_SQ3, I_SQ3)


def pieces():
    H1 = H_rot([-1, 0, 1])
    H2 = H_rot([-2, -1, 0, 1, 2])
    return {
        'H1': H1, 'H2': H2,
        'V25': union(H1, rotate(H, I_OVER_SQ3)),
        'V31': union(H1, rotate(H_rot([-1, 1]), I_SQ3)),
        'V37': union(H1, rotate(H1, I_OVER_SQ3)),
        'V49': union(H2, rotate(H1, I_OVER_SQ3)),
    }


def forced_radii(B, k=4):
    pts, xy, E = unit_graph(B)
    idx = {p: i for i, p in enumerate(pts)}
    centre = idx[C0]
    radii = {}
    for i, p in enumerate(pts):
        t = cnorm(p)
        if all(t[b] == 0 for b in range(1, DIM)) and t[0] != 0:
            radii.setdefault(t[0], []).append(i)
    s = Solver(name='cadical195')
    x = lambda v, c: v * k + c + 1
    for v in range(len(pts)):
        s.add_clause([x(v, c) for c in range(k)])
    for i, j in E:
        for c in range(k):
            s.add_clause([-x(i, c), -x(j, c)])
    s.add_clause([x(centre, 0)])
    colourable = s.solve()
    out = []
    if colourable:
        for t in sorted(radii):
            if not s.solve(assumptions=[-x(v, 0) for v in radii[t]]):
                q = 4 * t - 1
                f = sympy.factorint(q.numerator * q.denominator)
                sf = 1
                for p_, e in f.items():
                    if e % 2:
                        sf *= p_
                out.append((t, len(radii[t]), sf))
    s.delete()
    return len(pts), len(E), colourable, out


def main():
    P = pieces()
    names = ['V25', 'V31', 'V37', 'V49']
    jobs = [(a,) for a in names] + [(a, b) for i, a in enumerate(names) for b in names[i:]]
    jobs += [('H1', 'H1', 'H1'), ('H1', 'H1', 'V25'), ('H2', 'H1', 'H1')]
    for job in jobs:
        t0 = time.time()
        B = P[job[0]]
        for nm in job[1:]:
            B = minkowski(B, P[nm])
        n, m, col, forced = forced_radii(B)
        desc = ', '.join(f"t={t} (|O|={o}, sqfree(4t-1)={sf})" for t, o, sf in forced)
        print(f"{'+'.join(job):14s} |V|={n:5d} |E|={m:6d} 4-colourable={col}  forced: {desc if forced else '-'}  "
              f"[{time.time() - t0:.1f}s]", flush=True)


if __name__ == '__main__':
    main()
