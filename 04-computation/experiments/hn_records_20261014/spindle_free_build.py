#!/usr/bin/env python3
"""Rebuild the Moser-spindle-free 5-chromatic unit-distance graphs of Voronov, Neopryatnaya and Dergachev
(arXiv:2106.11824, series 2, case 10), exactly, and test 4-colourability.

Field: x, y in Q(sqrt2, sqrt3, sqrt5) (HN_PRIMES=2,3,5), so the plane field is Q(i, sqrt2, sqrt3, sqrt5); it has no
sqrt11, hence (their Proposition 6 / the rotation-field law) no Moser spindle.

  phi0 = e^(i pi/12) = (sqrt6 + sqrt2)/4 + i (sqrt6 - sqrt2)/4,   phi1 = sqrt6/3 + i sqrt3/3,
  M1 = {0} u {phi0^a phi1^b : 0 <= a <= 23, -1 <= b <= 1}       (73 points)
  M2 = clip(M1 + M1, 1),  M3 = M2 + M1                           (32,257 points expected)
  W  = M3 u psi M3,  psi = phi0^-1 (7/8 + i sqrt15/8)             (64,513 vertices expected)

Edges: chunked KD-tree ball queries give every pair at float distance within 1e-9 of 1; each candidate is then
decided exactly with integer arithmetic on numerators over a common denominator (vectorised over pairs).
Output: sizes, edge count, triangle used for symmetry breaking, and the CaDiCaL verdict.  The exact graph is saved
(npz: integer numerators X, Y of shape (n, 8), denominator, edges) to the path given as argv[1] (default: not saved).

Usage: HN_PRIMES=2,3,5 python3 -u spindle_free_build.py [out.npz] [solve=1]
"""
import os
import sys
import time
from fractions import Fraction as Fr
from math import lcm

import numpy as np
from scipy.spatial import cKDTree

assert os.environ.get('HN_PRIMES') == '2,3,5', "run with HN_PRIMES=2,3,5"
from hn_udg import (DIM, MULC, el, cmul, cpow, cconj, union, minkowski, cfloat, cnorm, ONE, C0, is_unit,
                    triangles_through)


def log(*a):
    print(*a, flush=True)


def build():
    phi0 = (el(r6=Fr(1, 4), r2=Fr(1, 4)), el(r6=Fr(1, 4), r2=Fr(-1, 4)))
    phi1 = (el(r6=Fr(1, 3)), el(r3=Fr(1, 3)))
    assert is_unit(phi0) and is_unit(phi1)
    assert cpow(phi0, 24) == (ONE, tuple(Fr(0) for _ in range(DIM)))
    gens = []
    for a in range(24):
        pa = cpow(phi0, a)
        for b in (-1, 0, 1):
            gens.append(cmul(pa, cpow(phi1, b)))
    M1 = frozenset([C0] + gens)
    assert len(M1) == 73, len(M1)
    S = minkowski(M1, M1)
    def inside(p):
        x, y = cfloat(p)
        r2 = x * x + y * y
        if abs(r2 - 1.0) < 1e-9:
            assert cnorm(p) == ONE, "ambiguous clip boundary"
            return True
        return r2 < 1.0
    M2 = frozenset(p for p in S if inside(p))
    M3 = minkowski(M2, M1)
    psi = cmul(cconj(phi0), (el(r1=Fr(7, 8)), el(r15=Fr(1, 8))))
    assert is_unit(psi)
    W = union(M3, frozenset(cmul(p, psi) for p in M3))
    return M1, M2, M3, W



def to_int(points):
    den = 1
    for p in points:
        for part in p:
            for c in part:
                den = lcm(den, c.denominator)
    X = np.array([[int(c * den) for c in p[0]] for p in points], dtype=np.int64)
    Y = np.array([[int(c * den) for c in p[1]] for p in points], dtype=np.int64)
    return X, Y, den


def exact_unit(X, Y, den, I, J):
    dx = X[I] - X[J]
    dy = Y[I] - Y[J]
    acc = np.zeros((len(I), DIM), dtype=object if den > 10 ** 4 else np.int64)
    for b1 in range(DIM):
        for b2 in range(DIM):
            m = MULC[b1][b2]
            acc[:, b1 ^ b2] += m * (dx[:, b1] * dx[:, b2] + dy[:, b1] * dy[:, b2])
    target = np.zeros(DIM, dtype=acc.dtype)
    target[0] = den * den
    return np.all(acc == target, axis=1)


def edges_of(points, X, Y, den, chunk=400):
    xy = np.array([cfloat(p) for p in points])
    tree = cKDTree(xy)
    I_all, J_all = [], []
    n = len(points)
    for s in range(0, n, chunk):
        idx = np.arange(s, min(n, s + chunk))
        nb = tree.query_ball_point(xy[idx], 1.0 + 1e-9, return_sorted=False)
        I, J = [], []
        for k, lst in enumerate(nb):
            i = idx[k]
            arr = np.asarray(lst, dtype=np.int64)
            arr = arr[arr > i]
            if len(arr) == 0:
                continue
            d = np.hypot(xy[arr, 0] - xy[i, 0], xy[arr, 1] - xy[i, 1])
            sel = arr[np.abs(d - 1.0) < 1e-9]
            I.extend([i] * len(sel))
            J.extend(sel.tolist())
        if I:
            I = np.array(I, dtype=np.int64)
            J = np.array(J, dtype=np.int64)
            ok = exact_unit(X, Y, den, I, J)
            I_all.append(I[ok])
            J_all.append(J[ok])
    I = np.concatenate(I_all) if I_all else np.zeros(0, dtype=np.int64)
    J = np.concatenate(J_all) if J_all else np.zeros(0, dtype=np.int64)
    return I, J


def main():
    out = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] != '-' else None
    solve = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    t0 = time.time()
    M1, M2, M3, W = build()
    log(f"|M1| = {len(M1)}, |M2| = {len(M2)}, |M3| = {len(M3)}, |W| = {len(W)}  [{time.time() - t0:.1f}s]")
    pts = sorted(W, key=cfloat)
    X, Y, den = to_int(pts)
    log(f"common denominator {den}; max |numerator| {int(max(np.abs(X).max(), np.abs(Y).max()))}  [{time.time() - t0:.1f}s]")
    I, J = edges_of(pts, X, Y, den)
    E = list(zip(I.tolist(), J.tolist()))
    log(f"|E| = {len(E)} exact unit edges  [{time.time() - t0:.1f}s]")
    if out:
        np.savez_compressed(out, X=X, Y=Y, den=np.array([den]), I=I, J=J)
        log(f"saved {out}")
    if not solve:
        return
    from pysat.solvers import Solver
    tri = triangles_through(len(pts), E, limit=1)[0]
    k = 4
    s = Solver(name='cadical195')
    x = lambda v, c: v * k + c + 1
    for v in range(len(pts)):
        s.add_clause([x(v, c) for c in range(k)])
    for i, j in E:
        for c in range(k):
            s.add_clause([-x(i, c), -x(j, c)])
    for c, v in enumerate(tri):
        s.add_clause([x(v, c)])
    t1 = time.time()
    r = s.solve()
    log(f"4-colourable? {r}  [{time.time() - t1:.1f}s solve, {time.time() - t0:.1f}s total]")


if __name__ == '__main__':
    main()
