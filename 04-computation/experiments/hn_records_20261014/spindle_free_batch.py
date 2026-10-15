#!/usr/bin/env python3
"""Batch of piece-level reductions of the spindle-free 5-chromatic graph M3 u psi P1 (48,121 vertices).

For each candidate forcing piece P' inside P1 = M2 + M1clock: (1) test that radius 2 is still forced in P' (skip
otherwise); (2) test 4-colourability of M3 u psi P' (CaDiCaL 1.9.5, triangle fixed).  Every verdict is appended to
spindle_free_batch.log with the vertex/edge counts.  Candidates:
  P1(r2max): points of P1 with |p|^2 <= r2max, plus the radius-2 orbit;
  M2(s) + M1clock: M2 replaced by clip(M1 + M1, s) for s < 1.

Usage: HN_PRIMES=2,3,5 python3 -u spindle_free_batch.py
"""
import os
import time
from fractions import Fraction as Fr

assert os.environ.get('HN_PRIMES') == '2,3,5'
from hn_udg import cmul, cfloat, cnorm, union, minkowski, DIM
from spindle_free_pieces import gens, clip, test
from spindle_free_forced import forced


def log(*a):
    print(*a, flush=True)


def main():
    phi0, phi1, psi, M1, M1c = gens()
    M11 = minkowski(M1, M1)
    M2 = clip(M11, 1.0)
    M3 = minkowski(M2, M1)
    P1 = minkowski(M2, M1c)
    four = (Fr(4),) + (Fr(0),) * (DIM - 1)
    ring = frozenset(p for p in P1 if cnorm(p) == four)
    cands = []
    for r2 in (3.5, 3.0, 2.5):
        cands.append((f"P1 with |p|^2 <= {r2} plus the radius-2 orbit",
                      frozenset(p for p in P1 if (lambda xy: xy[0] ** 2 + xy[1] ** 2)(cfloat(p)) <= r2 + 1e-9) | ring))
    for s2 in (0.8, 0.6):
        cands.append((f"clip(M1+M1, |p|^2 <= {s2}) + M1clock", minkowski(clip(M11, s2), M1c)))
    for name, Pp in cands:
        t0 = time.time()
        log(f"== {name}: |P'| = {len(Pp)}")
        forced(name, Pp)
        test(f"M3 u psi [{name}]", M3, Pp, psi)
        log(f"   [{time.time() - t0:.1f}s]")


if __name__ == '__main__':
    main()
