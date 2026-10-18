#!/usr/bin/env python3
"""Referee helper (second S41 audit): exact arithmetic in K = Q(zeta15, j), j^2 = -11, written from scratch.

Elements are pairs (p0, p1) of length-8 Fraction tuples meaning p0(zeta) + p1(zeta) j, reduced modulo
Phi_15(x) = x^8 - x^7 + x^5 - x^4 + x^3 - x + 1.  No code is shared with escaped_ports_20261017d.py.
Also builds the golden source sets P, M, A, C, L from their definitions in the notes, and evaluates elements at the
intended complex embedding (zeta15 = exp(2 pi i/15), j = +i sqrt 11) or at a residue place (z0, j0) mod p.
"""
from fractions import Fraction as Fr
import cmath, math

PHI15 = [1, -1, 0, 1, -1, 1, 0, -1, 1]          # coefficients x^0..x^8


def _polmul(a, b):
    out = [Fr(0)] * 15
    for i, x in enumerate(a):
        if x:
            for k, y in enumerate(b):
                if y:
                    out[i + k] += x * y
    # reduce degree >= 8 using x^8 = -(PHI15[0..7]) . x^i
    for deg in range(14, 7, -1):
        cf = out[deg]
        if cf:
            out[deg] = Fr(0)
            for i in range(8):
                out[deg - 8 + i] -= cf * PHI15[i]
    return tuple(out[:8])


ZERO8 = tuple([Fr(0)] * 8)


def K(c0=0, c1=None):
    """Constant c0 (rational)."""
    p0 = [Fr(0)] * 8; p0[0] = Fr(c0)
    return (tuple(p0), ZERO8)


def zeta_pow(k):
    k %= 15
    p = [Fr(0)] * 15; p[k] = Fr(1)
    # reduce
    out = list(p)
    for deg in range(14, 7, -1):
        cf = out[deg]
        if cf:
            out[deg] = Fr(0)
            for i in range(8):
                out[deg - 8 + i] -= cf * PHI15[i]
    return (tuple(out[:8]), ZERO8)


J = (ZERO8, tuple([Fr(1)] + [Fr(0)] * 7))


def add(x, y):
    return (tuple(a + b for a, b in zip(x[0], y[0])), tuple(a + b for a, b in zip(x[1], y[1])))


def sub(x, y):
    return (tuple(a - b for a, b in zip(x[0], y[0])), tuple(a - b for a, b in zip(x[1], y[1])))


def neg(x):
    return (tuple(-a for a in x[0]), tuple(-a for a in x[1]))


def scale(x, q):
    q = Fr(q)
    return (tuple(a * q for a in x[0]), tuple(a * q for a in x[1]))


def mul(x, y):
    a0b0 = _polmul(x[0], y[0]); a1b1 = _polmul(x[1], y[1])
    a0b1 = _polmul(x[0], y[1]); a1b0 = _polmul(x[1], y[0])
    return (tuple(p - 11 * q for p, q in zip(a0b0, a1b1)), tuple(p + q for p, q in zip(a0b1, a1b0)))


def conj(x):
    """zeta -> zeta^-1 = zeta^14, j -> -j."""
    out = (ZERO8, ZERO8)
    zinv = [zeta_pow(14 * i) for i in range(8)]
    r0 = (ZERO8, ZERO8)
    for i, cf in enumerate(x[0]):
        if cf:
            r0 = add(r0, scale(zinv[i], cf))
    r1 = (ZERO8, ZERO8)
    for i, cf in enumerate(x[1]):
        if cf:
            r1 = add(r1, scale(zinv[i], cf))
    return add(r0, neg(mul(r1, J)))


def is_zero(x):
    return all(a == 0 for a in x[0]) and all(a == 0 for a in x[1])


def norm2(x):
    return mul(x, conj(x))


def is_unit(x):
    n = norm2(x)
    return n == K(1)


ZETA = cmath.exp(2j * math.pi / 15)
JN = 1j * math.sqrt(11)


def to_complex(x, zeta=ZETA, jn=JN):
    s0 = sum(float(c) * zeta ** i for i, c in enumerate(x[0]) if c)
    s1 = sum(float(c) * zeta ** i for i, c in enumerate(x[1]) if c)
    return s0 + s1 * jn


def to_mod(x, p, z0, j0):
    tot = 0
    for part, jp in ((x[0], 1), (x[1], j0)):
        for i, c in enumerate(part):
            if c:
                if c.denominator % p == 0:
                    raise ValueError('bad prime')
                tot += c.numerator * pow(c.denominator, -1, p) * pow(z0, i, p) * jp
    return tot % p


def uniq(lst):
    out, seen = [], set()
    for x in lst:
        if x not in seen:
            seen.add(x); out.append(x)
    return out


def golden_sets():
    """P (8 literal ports), M (spindle), A = P + M, C = M + pentagon, L = L+ u -L+ (13 points), d = 1 + zeta5."""
    one = K(1)
    alpha = zeta_pow(3)
    a2, a4 = zeta_pow(6), zeta_pow(12)
    phi = add(add(one, alpha), a4)
    root5 = sub(scale(phi, 2), one)
    jr5 = mul(J, root5)
    P = [K(0), one, add(one, alpha), add(add(one, alpha), a2), neg(a4),
         scale(sub(K(3), jr5), Fr(1, 8)), K(Fr(3, 4)), scale(add(K(3), jr5), Fr(1, 8))]
    eta = add(one, zeta_pow(5))
    rho = scale(add(K(5), J), Fr(1, 6))
    M = [K(0), one, eta, add(one, eta), rho, mul(rho, eta), mul(rho, add(one, eta))]
    A = uniq(P + [add(p, m) for p in P for m in M])
    C = uniq(M + P[:5])
    a5 = P[5]
    Lp = [K(0), one, sub(one, a5), K(Fr(1, 4)), K(Fr(5, 4)), sub(K(Fr(5, 4)), a5), K(Fr(1, 2))]
    L = uniq(Lp + [neg(x) for x in Lp])
    d = add(one, alpha)
    return dict(P=P, M=M, A=A, C=C, L=L, d=d, alpha=alpha, eta=eta, rho=rho)


def diffs(X):
    return uniq([sub(a, b) for a in X for b in X if a != b])


def golden_numbers():
    """Intended embedding: w the root of w^3+2w^2-4w-4 in (1, 2); u, Z, N, D, R, s, h, v, e as complex numbers."""
    import numpy as np
    w = [r.real for r in np.roots([1, 2, -4, -4]) if abs(r.imag) < 1e-12 and 1 < r.real < 2][0]
    # refine by Newton in float
    for _ in range(5):
        w -= (w ** 3 + 2 * w * w - 4 * w - 4) / (3 * w * w + 4 * w - 4)
    u = (w + 3 - JN) / (w + 3 + JN)
    d = 1 + cmath.exp(2j * math.pi / 5)
    Z = d + u
    N = abs(Z) ** 2
    D = 1 - 4 / N
    R = 4 - N
    s = 1j * math.sqrt(-D)
    h = math.sqrt(R)
    v = Z * (-1 + s) / 2
    e = Z * s / h
    return dict(w=w, u=u, d=d, Z=Z, N=N, D=D, R=R, s=s, h=h, v=v, e=e)
