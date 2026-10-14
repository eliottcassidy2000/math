#!/usr/bin/env python3
"""Mono-set families on a reference orbit.

For a base graph B (with 0) and a squared radius t, the mono-set of a proper 4-colouring is
M = {p in B : |p|^2 = t, colour(p) = colour(0)}.  We enumerate every realizable M (SAT with blocking clauses on the
orbit's colour-0 indicators), print the size distribution, the minimum size, and the angles (degrees) of the points
of each minimal M.  A rotation w_t makes p (in B1) and w_t q (in w_t B2) adjacent iff arg p - arg q is 0 or 2 theta_t,
so a type-M contradiction needs every pair (M1, M2) to contain such a pair.

Usage: HN_PRIMES=... python3 -u mono_sets.py SPEC t [max_sets]
"""
import sys
import math
import time
from collections import Counter
from fractions import Fraction as Fr

from pysat.solvers import Solver

from hn_udg import unit_graph, cnorm, cfloat, DIM, C0
from minimize import pieces, build


def main():
    spec, t = sys.argv[1], Fr(sys.argv[2])
    cap = int(sys.argv[3]) if len(sys.argv) > 3 else 100000
    P = pieces()
    Bset = build(spec, P)
    pts, xy, E = unit_graph(Bset)
    idx = {p: i for i, p in enumerate(pts)}
    c0 = idx[C0]
    orbit = [i for i, p in enumerate(pts) if cnorm(p) == (t,) + (Fr(0),) * (DIM - 1)]
    ang = {i: (math.degrees(math.atan2(xy[i][1], xy[i][0])) % 360) for i in orbit}
    orbit.sort(key=lambda i: ang[i])
    k = 4
    x = lambda v, c: v * k + c + 1
    s = Solver(name='cadical195')
    for v in range(len(pts)):
        s.add_clause([x(v, c) for c in range(k)])
        for c in range(k):
            for d in range(c + 1, k):
                s.add_clause([-x(v, c), -x(v, d)])
    for i, j in E:
        for c in range(k):
            s.add_clause([-x(i, c), -x(j, c)])
    s.add_clause([x(c0, 0)])
    t0 = time.time()
    fam = []
    while len(fam) < cap and s.solve():
        model = set(l for l in s.get_model() if l > 0)
        M = tuple(i for i in orbit if x(i, 0) in model)
        fam.append(M)
        s.add_clause([(-x(i, 0) if i in M else x(i, 0)) for i in orbit])   # block exactly this M
    s.delete()
    sizes = Counter(len(M) for M in fam)
    print(f"{spec} t={t}: |V|={len(pts)} |E|={len(E)} orbit {len(orbit)} points at angles "
          f"{[round(ang[i], 2) for i in orbit]}")
    print(f"   realizable mono-sets: {len(fam)}{' (cap reached)' if len(fam) >= cap else ''}; size distribution "
          f"{dict(sorted(sizes.items()))}  [{time.time() - t0:.1f}s]")
    mn = min(sizes)
    mins = [M for M in fam if len(M) == mn]
    print(f"   minimal size {mn}: {len(mins)} sets; first ones (angles):")
    for M in mins[:12]:
        print("     ", [round(ang[i], 2) for i in M])


if __name__ == '__main__':
    main()
