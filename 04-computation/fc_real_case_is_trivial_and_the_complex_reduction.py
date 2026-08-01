#!/usr/bin/env python3
"""The structural reduction for FC at n=3, degree 3 -- and why the censuses were vacuous.

L IS A POSITIVE FUNCTIONAL.  L(g) = int_{[0,inf)^n} g(x) e^{-(x_1+...+x_n)} dx,
so L(g) > 0 whenever g >= 0 and g is not identically 0.  Hence for REAL f,

    L(f^2) = int f(x)^2 e^{-|x|} dx > 0    unless f = 0.

**So FC is trivial over the reals: the m = 2 condition alone forces f = 0.**
Every census in HYP-9076 imposed m = 1..M with M >= 2 over INTEGER (hence
real) coefficients, so its emptiness was forced for this reason and carried NO
information about FC.  Those runs are hereby marked vacuous.

THE REAL REDUCTION.  Write f = g + i h with g, h real.  Then, with the
Gram/Bombieri-style inner product <u,v> := L(uv) -- positive definite by the
above --

    L(f)   = 0   <=>   L(g) = 0  and  L(h) = 0,
    L(f^2) = 0   <=>   <g,g> = <h,h>   and   <g,h> = 0,
    L(f^3) = 0   <=>   L(g^3) = 3L(g h^2)  and  3L(g^2 h) = L(h^3).

So FC becomes a question of EUCLIDEAN GEOMETRY in the space of real
polynomials with the inner product <u,v> = L(uv) = sum over monomial pairs of
(alpha+beta)!:  seek g PERPENDICULAR to h, of EQUAL NORM, both L-mean zero,
satisfying in addition the cubic pair of identities, and so on for every m.

That is the structural reduction: 19 complex unknowns collapse to a pair of
real vectors constrained by orthogonality and equal norm, both cheap linear /
quadratic conditions, with the content beginning at m = 3.
"""
import sys
from itertools import product
from math import factorial

import numpy as np


def monomials(nvars, maxdeg):
    return [a for a in product(range(maxdeg + 1), repeat=nvars)
            if sum(a) <= maxdeg]


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


def Lvec(mons, nvars):
    return np.array([np.prod([factorial(e) for e in a]) for a in mons], float)


def L_of_poly_power(coeffs, m, nvars):
    """exact L(f^m) for integer/complex dict f."""
    poly = {(0,) * nvars: 1}
    for _ in range(m):
        new = {}
        for a, ca in poly.items():
            for b, cb in coeffs.items():
                k = tuple(a[i] + b[i] for i in range(nvars))
                new[k] = new.get(k, 0) + ca * cb
        poly = {k: v for k, v in new.items() if v}
    tot = 0
    for a, c in poly.items():
        t = c
        for e in a:
            t *= factorial(e)
        tot += t
    return tot


if __name__ == "__main__":
    for nvars, maxdeg in ((2, 3), (3, 3)):
        mons = monomials(nvars, maxdeg)
        G = gram(mons, nvars)
        ev = np.linalg.eigvalsh(G)
        print(f"n={nvars}, deg<={maxdeg}: {len(mons)} monomials, "
              f"Gram <u,v>=L(uv) min eigenvalue {ev.min():.6g}  "
              f"-> positive definite: {ev.min() > 0}")
    print()
    print("CONSEQUENCE: for real f != 0, L(f^2) = <f,f> > 0, so no real f can")
    print("satisfy even the m=2 condition.  The HYP-9076 censuses were vacuous.")
    print()

    # An explicit complex witness at m = 1,2 -- impossible over the reals.
    nvars = 3
    mons = monomials(nvars, 3)
    Lv = Lvec(mons, nvars)
    G = gram(mons, nvars)
    # g, h with L(g)=L(h)=0, <g,h>=0, <g,g>=<h,h>: solve in a small subspace
    rng = np.random.default_rng(0)
    found = None
    for _ in range(20000):
        v = rng.integers(-3, 4, size=len(mons)).astype(float)
        w = rng.integers(-3, 4, size=len(mons)).astype(float)
        # project to L-mean zero
        v -= (Lv @ v) / (Lv @ Lv) * Lv
        w -= (Lv @ w) / (Lv @ Lv) * Lv
        if v @ G @ v < 1e-9 or w @ G @ w < 1e-9:
            continue
        # Gram-Schmidt w against v, then rescale to equal norm
        w = w - (v @ G @ w) / (v @ G @ v) * v
        if w @ G @ w < 1e-9:
            continue
        w *= np.sqrt((v @ G @ v) / (w @ G @ w))
        found = (v, w)
        break
    v, w = found
    print("explicit real pair (g,h) in n=3, deg<=3 with the m=1,2 conditions:")
    print(f"   L(g) = {Lv @ v:.3e}   L(h) = {Lv @ w:.3e}")
    print(f"   <g,h> = {v @ G @ w:.3e}   <g,g>-<h,h> = {v @ G @ v - w @ G @ w:.3e}")
    print("   so f = g + i h satisfies L(f) = L(f^2) = 0 while f != 0.")
    print()
    print("""=> The m <= 2 level is NON-EMPTY over C and empty over R for a
trivial reason.  Any meaningful FC census must search COMPLEX f, i.e. pairs
(g,h) already constrained to be L-mean-zero, orthogonal and equinormal, with
the first real content at m = 3.  That is the structural reduction the n = 3,
degree 3 case needs: not a bigger box, a different space.""")
