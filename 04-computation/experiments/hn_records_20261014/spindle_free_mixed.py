#!/usr/bin/env python3
"""Mixed spindle-free type-M pieces, and the mono-set family of the clock piece at radius 2.

Pieces: P1 = M2 + M1clock (15,865 vertices, radius 2 forced on the 24-point clock orbit 2 phi0^a),
        M3 = M2 + M1 (32,257 vertices, radius 2 forced on a 72-point orbit).
Hinge psi = phi0^-1 (7/8 + i sqrt15/8) joins 2 phi0^(a-1) and 2 psi phi0^a (the rho angle, 28.955 degrees).

Usage: HN_PRIMES=2,3,5 python3 -u spindle_free_mixed.py
"""
import os
import time
import math
from fractions import Fraction as Fr

from pysat.solvers import Solver

assert os.environ.get('HN_PRIMES') == '2,3,5'
from hn_udg import minkowski, cfloat, cnorm, C0, DIM
import spindle_free_build as sfb
from spindle_free_pieces import gens, clip, test


def mono_sets(name, B, cap=5000, budget=50_000_000):
    t0 = time.time()
    pts = sorted(B, key=cfloat)
    X, Y, den = sfb.to_int(pts)
    I, J = sfb.edges_of(pts, X, Y, den)
    idx = {p: i for i, p in enumerate(pts)}
    c0 = idx[C0]
    four = (Fr(4),) + (Fr(0),) * (DIM - 1)
    orbit = [i for i, p in enumerate(pts) if cnorm(p) == four]
    xy = [cfloat(pts[i]) for i in orbit]
    ang = {i: round(math.degrees(math.atan2(y, x)) % 360, 3) for i, (x, y) in zip(orbit, xy)}
    k = 4
    s = Solver(name='cadical195')
    x = lambda v, c: v * k + c + 1
    for v in range(len(pts)):
        s.add_clause([x(v, c) for c in range(k)])
        for c in range(k):
            for d in range(c + 1, k):
                s.add_clause([-x(v, c), -x(v, d)])
    for i, j in zip(I.tolist(), J.tolist()):
        for c in range(k):
            s.add_clause([-x(i, c), -x(j, c)])
    s.add_clause([x(c0, 0)])
    fam = []
    while len(fam) < cap:
        s.conf_budget(budget)
        r = s.solve_limited()
        if r is not True:
            break
        model = set(l for l in s.get_model() if l > 0)
        M = tuple(i for i in orbit if x(i, 0) in model)
        fam.append(M)
        s.add_clause([(-x(i, 0) if i in M else x(i, 0)) for i in orbit])
    s.delete()
    from collections import Counter
    sizes = Counter(len(M) for M in fam)
    print(f"{name}: |V|={len(pts)}, radius-2 orbit {len(orbit)}; realizable mono-sets {len(fam)}"
          f"{' (cap)' if len(fam) >= cap else ''}; sizes {dict(sorted(sizes.items()))}  [{time.time() - t0:.1f}s]", flush=True)
    mn = min(sizes) if sizes else None
    for M in [M for M in fam if len(M) == mn][:8]:
        print("     ", sorted(ang[i] for i in M), flush=True)


def main():
    phi0, phi1, psi, M1, M1c = gens()
    M2 = clip(minkowski(M1, M1), 1.0)
    P1 = minkowski(M2, M1c)
    M3 = minkowski(M2, M1)
    test('M3 u psi P1', M3, P1, psi)
    test('P1 u psi M3', P1, M3, psi)
    mono_sets('P1 = M2+M1clock', P1, cap=300)


if __name__ == '__main__':
    main()
