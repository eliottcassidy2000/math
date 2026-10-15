#!/usr/bin/env python3
"""Clip the rigid side: test M3' u psi P' for M3' = {p in M3 : |p|^2 <= r} u (radius-2 orbit of M3) and the clipped forcing
piece P' = {p in P1 : |p|^2 <= 3.0} u (radius-2 orbit), which already gives a 46,585-vertex 5-chromatic graph with full M3
(spindle_free_batch.log).  The cross edges of M3 u psi P1 join radius-2 points and |z|^2 = 8/3 points of M3 to |psi w|^2 = 2/3
points (audit B), so r >= 8/3 keeps every hinge edge.  CaDiCaL 1.9.5, triangle fixed, conflict budget as in spindle_free_pieces.

Usage: HN_PRIMES=2,3,5 python3 -u spindle_free_m3side.py >> spindle_free_m3side.log
"""
import os
import time
from fractions import Fraction as Fr

assert os.environ.get('HN_PRIMES') == '2,3,5'
from hn_udg import cmul, cfloat, cnorm, minkowski, DIM
from spindle_free_pieces import gens, clip, test


def log(*a):
    print(*a, flush=True)


def main():
    phi0, phi1, psi, M1, M1c = gens()
    M11 = minkowski(M1, M1)
    M2 = clip(M11, 1.0)
    M3 = minkowski(M2, M1)
    P1 = minkowski(M2, M1c)
    four = (Fr(4),) + (Fr(0),) * (DIM - 1)
    ring3 = frozenset(p for p in M3 if cnorm(p) == four)
    ringP = frozenset(p for p in P1 if cnorm(p) == four)
    r2 = lambda p: (lambda xy: xy[0] ** 2 + xy[1] ** 2)(cfloat(p))
    Pp = frozenset(p for p in P1 if r2(p) <= 3.0 + 1e-9) | ringP
    log(f"|M3| = {len(M3)}, radius-2 orbit of M3: {len(ring3)}; P' (|p|^2 <= 3.0) = {len(Pp)}")
    for r in (3.5, 3.0, 2.7):
        t0 = time.time()
        M3p = frozenset(p for p in M3 if r2(p) <= r + 1e-9) | ring3
        test(f"M3[|p|^2<={r}] u psi P'[3.0]", M3p, Pp, psi, budget=40_000_000)
        log(f"   [{time.time() - t0:.1f}s]")


if __name__ == '__main__':
    main()
