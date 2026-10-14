#!/usr/bin/env python3
"""Scan type-M graphs A u w_t B over base pairs (A, B) and reference radii t (Parts' open question: type-M
constructions with rotations other than rho = w_4).

w_t = exp(i theta_t), cos theta_t = 1 - 1/(2t), sin theta_t = sqrt(4t - 1)/(2t): |p - w_t p| = 1 whenever |p|^2 = t.
The base pieces live in Q(sqrt3, sqrt11); w_t needs sqrt(4t - 1), so HN_PRIMES must contain its squarefree primes.
For each graph: vertices, edges, and 4-colourability (CaDiCaL 1.9.5, a triangle fixed, conflict budget per call;
'?' if the budget runs out).

Usage: HN_PRIMES=3,5,11 python3 -u rotation_scan.py 4 [t2 ...]
"""
import sys
import time
from fractions import Fraction as Fr

from pysat.solvers import Solver

from hn_udg import union, rotate, unit_graph, triangles_through
from rotation_probe import w_of
from minimize import pieces, build

A_LIST = ['V49+V37', 'V49+V49', 'V37+V37', 'V31+V49', 'V25+V49']
B_LIST = ['V31+V25', 'V25+V25', 'V37+V25', 'V31+V31', 'V49+V37']
BUDGET = 2_000_000


def colourable(n, E, tri, k=4):
    s = Solver(name='cadical195')
    x = lambda v, c: v * k + c + 1
    for v in range(n):
        s.add_clause([x(v, c) for c in range(k)])
    for i, j in E:
        for c in range(k):
            s.add_clause([-x(i, c), -x(j, c)])
    for c, v in enumerate(tri):
        s.add_clause([x(v, c)])
    s.conf_budget(BUDGET)
    r = s.solve_limited()
    s.delete()
    return r


def main():
    ts = [Fr(a) for a in sys.argv[1:]] or [Fr(4)]
    P = pieces()
    built = {}
    for t in ts:
        w = w_of(t)
        for A in A_LIST:
            for B in B_LIST:
                t0 = time.time()
                GA = built.setdefault(A, build(A, P))
                GB = built.setdefault(B, build(B, P))
                G = union(GA, rotate(GB, w))
                pts, xy, E = unit_graph(G)
                tris = triangles_through(len(pts), E, limit=1)
                r = colourable(len(pts), E, tris[0])
                verdict = {True: '4-colourable', False: 'NOT 4-colourable (5-chromatic)', None: '? (budget)'}[r]
                print(f"t={str(t):5s} A={A:8s} B={B:8s} |V|={len(pts):5d} |E|={len(E):6d}  {verdict}  "
                      f"[{time.time() - t0:.1f}s]", flush=True)


if __name__ == '__main__':
    main()
