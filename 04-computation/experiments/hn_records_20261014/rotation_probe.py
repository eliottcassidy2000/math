#!/usr/bin/env python3
"""Type-M rotation probe at a forced radius: G = L727 u w_t X for X in {L727, S361}, where
w_t = exp(i theta_t), cos theta_t = 1 - 1/(2t), sin theta_t = sqrt(4t - 1)/(2t), so |p - w_t p| = 1 for |p|^2 = t.

The field must contain sqrt(4t - 1): run with HN_PRIMES covering its squarefree part, e.g.
  HN_PRIMES=3,5,7,11 python3 -u rotation_probe.py 4/9 LL

Reports the size, exact edges and whether the graph is 4-colourable (CaDiCaL 1.9.5, a triangle fixed).
"""
import sys
import time
from fractions import Fraction as Fr

import sympy

from hn_udg import (PRIMES, H, H_rot, union, rotate, minkowski, unit_graph, check_float_gap, triangles_through,
                    Colourer, is_unit, el, I_OVER_SQ3, I_SQ3)


def w_of(t):
    c = 1 - 1 / (2 * t)
    q = 4 * t - 1                       # sin = sqrt(q)/(2t)
    num, den = q.numerator, q.denominator
    f = sympy.factorint(num * den)
    s, k = 1, 1
    for p, e in f.items():
        k *= p ** (e // 2)
        if e % 2:
            s *= p
    assert all(p in PRIMES for p in sympy.factorint(s)), f"need primes of {s} in HN_PRIMES"
    coef = Fr(k, den) / (2 * t)         # sqrt(num/den) = k sqrt(s)/den
    w = (el(r1=c), el(**{f"r{s}": coef}))
    assert is_unit(w)
    return w


def main():
    t = Fr(sys.argv[1])
    pair = sys.argv[2] if len(sys.argv) > 2 else 'LS'
    t0 = time.time()
    H1 = H_rot([-1, 0, 1])
    H2 = H_rot([-2, -1, 0, 1, 2])
    V25 = union(H1, rotate(H, I_OVER_SQ3))
    V37 = union(H1, rotate(H1, I_OVER_SQ3))
    V49 = union(H2, rotate(H1, I_OVER_SQ3))
    V31 = union(H1, rotate(H_rot([-1, 1]), I_SQ3))
    L = minkowski(V49, V37)
    S = minkowski(V31, V25)
    X = {'LL': L, 'LS': S}[pair]
    w = w_of(t)
    G = union(L, rotate(X, w))
    pts, xy, E = unit_graph(G)
    n = len(pts)
    tris = triangles_through(n, E)
    print(f"t={t} pair={pair} w=({float(1 - 1/(2*t)):.6f}, ...) |G|={n} |E|={len(E)} triangles={len(tris)} "
          f"float-gap={check_float_gap(xy, E)}  [{time.time() - t0:.1f}s]", flush=True)
    col = Colourer(n, E, k=4, solver='cadical195', fixed_triangle=tris[0])
    t1 = time.time()
    ok, core = col.colourable(range(n))
    print(f"   4-colourable? {ok}" + ("" if ok else f"; core {len(core)}") + f"  [{time.time() - t1:.1f}s]", flush=True)


if __name__ == '__main__':
    main()
