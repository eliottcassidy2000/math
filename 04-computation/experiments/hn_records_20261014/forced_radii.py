#!/usr/bin/env python3
"""Forced monochromatic radii of the type-M base graphs.

For a base graph B containing the centre 0 and a squared radius t, call t FORCED if every proper 4-colouring of the
unit-distance graph on B gives some vertex p with |p|^2 = t the colour of 0.  (SAT: forbid colour(0) on the whole
t-orbit; UNSAT means forced.)  A type-M rotation w_t joins p and w_t p at unit distance for |p|^2 = t, so a forced
radius in both pieces is the natural place to hinge a rotation.

Usage: python3 -u forced_radii.py
"""
import time
from fractions import Fraction as Fr

from pysat.solvers import Solver

from hn_udg import (H, H_rot, union, rotate, minkowski, unit_graph, cnorm, DIM, C0, I_OVER_SQ3, I_SQ3)


def base_graphs():
    H1 = H_rot([-1, 0, 1])
    H2 = H_rot([-2, -1, 0, 1, 2])
    V25 = union(H1, rotate(H, I_OVER_SQ3))
    V37 = union(H1, rotate(H1, I_OVER_SQ3))
    V49 = union(H2, rotate(H1, I_OVER_SQ3))
    V31 = union(H1, rotate(H_rot([-1, 1]), I_SQ3))
    return {'L727': minkowski(V49, V37), 'S361': minkowski(V31, V25)}


def forced(pts, E, centre, orbit, k=4):
    s = Solver(name='cadical195')
    x = lambda v, c: v * k + c + 1
    for v in range(len(pts)):
        s.add_clause([x(v, c) for c in range(k)])
    for i, j in E:
        for c in range(k):
            s.add_clause([-x(i, c), -x(j, c)])
    s.add_clause([x(centre, 0)])                      # colour of the centre is 0 (symmetry)
    for v in orbit:
        s.add_clause([-x(v, 0)])
    sat = s.solve()
    s.delete()
    return not sat


def main():
    for name, B in base_graphs().items():
        t0 = time.time()
        pts, xy, E = unit_graph(B)
        idx = {p: i for i, p in enumerate(pts)}
        centre = idx[C0]
        radii = {}
        for i, p in enumerate(pts):
            t = cnorm(p)
            if all(t[b] == 0 for b in range(1, DIM)) and t[0] != 0:
                radii.setdefault(t[0], []).append(i)
        print(f"{name}: {len(pts)} vertices, {len(E)} edges", flush=True)
        for t in sorted(radii):
            f = forced(pts, E, centre, radii[t])
            print(f"   t = {str(t):6s} |orbit| = {len(radii[t]):3d}  forced mono with centre: {f}", flush=True)
        print(f"   [{time.time() - t0:.1f}s]", flush=True)


if __name__ == '__main__':
    main()
