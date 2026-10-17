#!/usr/bin/env python3
"""S39 audit: reproduce THM-4664's 'Open' paragraph claim (uncommitted scratch run) for Gao's map G (p = w^3 - 6w^2 + 6w,
gamma = 2 - 4xy - x^2 z): the numbers of frequencies with |S| > p^2/2 are 340, 600, 384 at p = 11, 13, 17, with the line
(0, 0, gamma) among them.  The map is built from THM-4625's formula (as in a8), evaluated exactly mod p, FFT.
Usage: python -B a8b_gao_fourier.py
"""
import numpy as np
import sympy as sp

x, y, z, w, s = sp.symbols('x y z w s')
pw = w ** 3 - 6 * w ** 2 + 6 * w
qw = sp.integrate(s / 2 * sp.diff(pw.subs(w, s), s), (s, 0, w))
gam = 2 - 4 * x * y - x ** 2 * z
W = gam * (1 + x * y)
F = [sp.expand(sp.cancel((qw.subs(w, W) + gam * W) / (gam * x) ** 2)), sp.expand(sp.cancel((pw.subs(w, W) + 2 * gam) / (gam * x))),
     sp.expand(gam * x)]
dens = [int(sp.denom(sp.together(f))) for f in F]
polys = [sp.Poly(sp.expand(f * D), x, y, z) for f, D in zip(F, dens)]
print('det JF =', sp.simplify(sp.Matrix([[sp.diff(f, v) for v in (x, y, z)] for f in F]).det()))
for p in (11, 13, 17):
    a = np.arange(p, dtype=np.int64)
    X, Y, Z = np.meshgrid(a, a, a, indexing='ij')
    vals = []
    for P, D in zip(polys, dens):
        acc = np.zeros_like(X)
        for (i, j, k), cf in zip(P.monoms(), P.coeffs()):
            acc = (acc + (int(cf) % p) * (pow(1, 1) * (X ** i % p) * (Y ** j % p) % p * (Z ** k % p))) % p
        vals.append(acc * pow(D, -1, p) % p)
    N = np.bincount(((vals[0] * p + vals[1]) * p + vals[2]).ravel(), minlength=p ** 3).reshape(p, p, p).astype(float)
    S = np.abs(np.fft.fftn(N)); S[0, 0, 0] = 0
    big = S > p * p / 2
    line = all(big[0, 0, g] for g in range(1, p))
    print('p = %d: #{|S| > p^2/2} = %d; line (0,0,gamma) included: %s; fibre sizes: %s'
          % (p, int(big.sum()), line, np.bincount(N.astype(int).ravel()).tolist()))

# cross-check of the evaluation against exact sympy substitution at random points, and other thresholds
import random
random.seed(1)
for p in (11, 13):
    a = np.arange(p, dtype=np.int64)
    X, Y, Z = np.meshgrid(a, a, a, indexing='ij')
    vals = []
    for P, D in zip(polys, dens):
        acc = np.zeros_like(X)
        for (i, j, k), cf in zip(P.monoms(), P.coeffs()):
            acc = (acc + (int(cf) % p) * ((X ** i % p) * (Y ** j % p) % p * (Z ** k % p))) % p
        vals.append(acc * pow(D, -1, p) % p)
    for _ in range(200):
        xx, yy, zz = (random.randrange(p) for _ in range(3))
        exact = [int(sp.Rational(f.subs({x: xx, y: yy, z: zz})).p * pow(int(sp.Rational(f.subs({x: xx, y: yy, z: zz})).q), -1, p)) % p for f in F]
        assert exact == [int(v[xx, yy, zz]) for v in vals], (p, xx, yy, zz)
    N = np.bincount(((vals[0] * p + vals[1]) * p + vals[2]).ravel(), minlength=p ** 3).reshape(p, p, p).astype(float)
    S = np.abs(np.fft.fftn(N)); S[0, 0, 0] = 0
    srt = np.sort(S.ravel())[::-1]
    print('p = %d: evaluation matches exact substitution at 200 random points; counts above thresholds:' % p,
          {'p^2/2': int((S > p * p / 2).sum()), 'p^2/4': int((S > p * p / 4).sum()), 'p^1.5': int((S > p ** 1.5).sum()),
           '2p': int((S > 2 * p).sum()), 'p': int((S > p).sum())}, '; top |S|/p^2:', np.round(srt[:6] / p ** 2, 3).tolist())

# which frequencies are large (p = 11, 13)?  numpy's fftn uses e^(-2 pi i a.w/p): |S| is the same for a and -a.
for p in (11, 13):
    a = np.arange(p, dtype=np.int64)
    X, Y, Z = np.meshgrid(a, a, a, indexing='ij')
    vals = []
    for P, D in zip(polys, dens):
        acc = np.zeros_like(X)
        for (i, j, k), cf in zip(P.monoms(), P.coeffs()):
            acc = (acc + (int(cf) % p) * ((X ** i % p) * (Y ** j % p) % p * (Z ** k % p))) % p
        vals.append(acc * pow(D, -1, p) % p)
    N = np.bincount(((vals[0] * p + vals[1]) * p + vals[2]).ravel(), minlength=p ** 3).reshape(p, p, p).astype(float)
    S = np.abs(np.fft.fftn(N)); S[0, 0, 0] = 0
    big = sorted(tuple(int(t) for t in v) for v in np.argwhere(S > p * p / 2))
    print('p = %d: large frequencies (alpha, beta, gamma):' % p, big)
