#!/usr/bin/env python3
"""
procgen_tension_20260926_mip.py -- an independent, flow-based formulation of the Q3 queries (HiGHS MIP via
scipy.optimize.milp), used as a cross-check of the CP-SAT answers.

    exists a class-(i) level-k strategy sigma and a nonzero 0/1 circulation f on the edges of G_sigma with
    sum_e f_e (p0 [tail odd] - a0) >= 1 ?

A 0/1 circulation is an edge-disjoint union of cycles, so this holds iff some class-(i) sigma has a cycle of
density > a0/p0: the obstruction side is a FLOW, the class-(i) side a real-valued POTENTIAL (tension).
Run through the runner with --mip6 (level 6 may take long); levels 4, 5 are fast.
"""
import time
from fractions import Fraction

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

from procgen_tension_20260926_lib import best_lower, below_c, check, claim, karp


def flow_query(k, thr, tl=3600):
    M = 1 << k
    H = M >> 1
    F = best_lower(M)
    q, r = F.numerator, F.denominator
    a0, p0 = thr.numerator, thr.denominator
    odd = list(range(1, M, 2))
    xi = {s: i for i, s in enumerate(odd)}
    nx = len(odd)
    E = []
    for s in range(M):
        if s % 2 == 0:
            b = (s // 2) % H
            E += [(s, b, None), (s, b + H, None)]
        else:
            bp = ((3 * s + 1) // 2) % H
            bm = ((3 * s - 1) // 2) % H
            E += [(s, bp, 'p'), (s, bp + H, 'p'), (s, bm, 'm'), (s, bm + H, 'm')]
    nE = len(E)
    nv = nx + nE + M
    fo, po = nx, nx + nE
    U = (M // 2) * (r - q) + 1
    BIG = U + max(r - q, q) + 1
    rows, lb, ub = [], [], []

    def add(coefs, lo, hi):
        rows.append(coefs)
        lb.append(lo)
        ub.append(hi)
    for i, (s, t, lit) in enumerate(E):
        w = (r - q) if s % 2 else -q
        if lit == 'm':            # edge present iff sigma(s) = -  (x_s = 1)
            add({fo + i: 1, xi[s]: -1}, -np.inf, 0)
            add({po + t: 1, po + s: -1, xi[s]: BIG} if t != s else {xi[s]: BIG}, -np.inf, BIG - w)
        elif lit == 'p':          # edge present iff sigma(s) = +  (x_s = 0)
            add({fo + i: 1, xi[s]: 1}, -np.inf, 1)
            add({po + t: 1, po + s: -1, xi[s]: -BIG} if t != s else {xi[s]: -BIG}, -np.inf, -w)
        else:
            if t != s:
                add({po + t: 1, po + s: -1}, -np.inf, -w)
            else:
                check(-w >= 0, "even self-loop has non-positive weight")
    for v in range(M):
        c = {}
        for i, (s, t, lit) in enumerate(E):
            if t == v:
                c[fo + i] = c.get(fo + i, 0) + 1
            if s == v:
                c[fo + i] = c.get(fo + i, 0) - 1
        c = {j: x for j, x in c.items() if x != 0}
        if c:
            add(c, 0, 0)
    add({fo + i: ((p0 - a0) if s % 2 else -a0) for i, (s, t, lit) in enumerate(E)}, 1, np.inf)
    A = lil_matrix((len(rows), nv))
    for ri, c in enumerate(rows):
        for j, x in c.items():
            A[ri, j] = x
    integ = np.array([1] * (nx + nE) + [0] * M)
    lo = np.zeros(nv)
    hi = np.array([1.0] * (nx + nE) + [float(U)] * M)
    t0 = time.time()
    res = milp(np.zeros(nv), constraints=LinearConstraint(A.tocsr(), lb, ub), integrality=integ,
               bounds=Bounds(lo, hi), options={"time_limit": tl, "disp": False})
    dt = time.time() - t0
    if res.status == 0:
        mask = sum(1 << ((s - 1) // 2) for s in odd if res.x[xi[s]] > 0.5)
        return 'SAT', mask, dt
    if res.status == 2:
        return 'UNSAT', None, dt
    return 'UNKNOWN', None, dt


def run(levels=((4, Fraction(1, 2), 'SAT'), (4, Fraction(3, 5), 'UNSAT'), (5, Fraction(3, 5), 'SAT'),
                (5, Fraction(5, 8), 'UNSAT'), (6, Fraction(5, 8), 'UNSAT')), tl=3600):
    print("MIP cross-check (flow formulation, HiGHS)")
    for k, thr, expect in levels:
        st, mask, dt = flow_query(k, thr, tl)
        if st == 'SAT':
            rho = karp(k, mask)
            check(below_c(rho) and rho > thr, "MIP witness re-verified by exact Karp")
        check(st == expect, "flow MIP at level %d, threshold %s: %s (expected %s)" % (k, thr, st, expect))
        claim(True, "level %d: a class-(i) strategy with a cycle of density > %s %s (flow MIP, %.1f s)"
              % (k, thr, "exists (re-verified by exact Karp)" if st == 'SAT' else "does not exist", dt))


if __name__ == "__main__":
    run()
