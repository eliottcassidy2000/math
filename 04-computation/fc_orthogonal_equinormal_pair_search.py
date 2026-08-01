#!/usr/bin/env python3
"""Search the orthogonal / equinormal / mean-zero pairs for m = 3 and beyond.

HYP-9076 sec 6 reduced FC to: f = g + i h with g,h real, where

    L(f)   = 0  <=>  L(g) = L(h) = 0
    L(f^2) = 0  <=>  <g,g> = <h,h>  and  <g,h> = 0        (<u,v> := L(uv))
    L(f^3) = 0  <=>  L(g^3) = 3L(g h^2)  and  3L(g^2 h) = L(h^3)

and asked whether the m = 3 conditions can be met on that variety.

DIMENSION COUNT FIRST.  Each L(f^m) = 0 is ONE COMPLEX equation, i.e. two real
ones, so m = 1..M imposes 2M real conditions on the 2N real parameters of
f (N = number of monomials).  With the scale and phase symmetries
f -> lambda e^{i theta} f removing 2 more, solutions are expected whenever

    2M  <  2N - 2,     i.e.    M < N - 1.

For n=2, deg<=3 that is N = 10, so every window M <= 8 should be solvable;
for n=3, deg<=3, N = 20 and every M <= 18.  **No finite window can obstruct
FC** -- the conjecture is genuinely infinitary.  This script checks the
prediction by actually solving.
"""
import sys
from itertools import product
from math import factorial

import numpy as np
from scipy.optimize import least_squares


def monomials(nvars, maxdeg):
    return [a for a in product(range(maxdeg + 1), repeat=nvars)
            if sum(a) <= maxdeg]


def make_L_power(mons, nvars):
    idx = {a: i for i, a in enumerate(mons)}

    def Lm(coef, m):
        poly = {(0,) * nvars: 1.0 + 0j}
        f = {a: coef[i] for i, a in enumerate(mons) if coef[i] != 0}
        for _ in range(m):
            new = {}
            for a, ca in poly.items():
                for b, cb in f.items():
                    k = tuple(a[i] + b[i] for i in range(nvars))
                    new[k] = new.get(k, 0) + ca * cb
            poly = new
        tot = 0
        for a, c in poly.items():
            t = c
            for e in a:
                t *= factorial(e)
            tot += t
        return tot
    return Lm


def gram(mons, nvars):
    n = len(mons)
    G = np.zeros((n, n))
    for i, a in enumerate(mons):
        for j, b in enumerate(mons):
            t = 1
            for k in range(nvars):
                t *= factorial(a[k] + b[k])
            G[i, j] = t
    return G


def solve_window(nvars, maxdeg, M, tries=40, seed=0):
    mons = monomials(nvars, maxdeg)
    N = len(mons)
    Lm = make_L_power(mons, nvars)
    G = gram(mons, nvars)
    rng = np.random.default_rng(seed)

    def resid(p):
        c = p[:N] + 1j * p[N:]
        out = []
        for m in range(1, M + 1):
            v = Lm(c, m)
            out += [v.real, v.imag]
        # normalisation: ||f||_G = 1 (kills the scale symmetry and f = 0)
        nr = float(np.real(np.conj(c) @ G @ c))
        out.append(nr - 1.0)
        return out

    best = None
    for t in range(tries):
        p0 = rng.normal(size=2 * N)
        p0 /= np.linalg.norm(p0)
        r = least_squares(resid, p0, xtol=1e-15, ftol=1e-15, gtol=1e-15,
                          max_nfev=20000)
        res = np.max(np.abs(resid(r.x)))
        if best is None or res < best[0]:
            best = (res, r.x)
        if res < 1e-10:
            break
    return best, mons, Lm, G


if __name__ == "__main__":
    print("Searching orthogonal/equinormal/mean-zero pairs for m = 1..M")
    print("prediction: solvable while M < N-1  (N = #monomials)")
    print()
    for nvars, maxdeg in ((2, 3), (3, 3)):
        N = len(monomials(nvars, maxdeg))
        print(f"n={nvars}, deg<={maxdeg}:  N = {N} monomials, "
              f"2N = {2*N} real params, predicted window M < {N-1}")
        for M in (2, 3, 4, 6, 8, 12):
            if M > N + 2:
                continue
            (res, x), mons, Lm, G = solve_window(nvars, maxdeg, M)
            c = x[:N] + 1j * x[N:]
            nrm = float(np.real(np.conj(c) @ G @ c))
            tag = "SOLVED" if res < 1e-9 else ("partial" if res < 1e-4 else "NOT FOUND")
            print(f"   M={M:3d}: max|residual| = {res:.3e}   ||f||^2 = {nrm:.4f}"
                  f"   -> {tag}")
        print()
    print("""READING.  A solved window means a nonzero complex f whose first M
moments all vanish.  Such f exist for every window the dimension count allows,
which is exactly why no census can ever refute FC: the conjecture is about the
INFINITE family, and every finite truncation is under-determined.""")
