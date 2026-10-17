#!/usr/bin/env python3
"""S39 audit, addendum (3) continued: identify the third moment of the four-bar family.

Eichler-Shimura/Deligne prediction: H^1(P^1, j_* Sym^3 R^1) for the universal curve over X_1(6) is
S_5(Gamma_1(6)) + conj = S_5(Gamma_0(6), (-3/.)) + conj (dimension 4), so
    t_p = -sum_c a_c^3 - (1 + (p/3))(4p + 2) = a_p(f) + conj(a_p(f))   for f the newform(s) of S_5(Gamma_0(6), (-3/.)).
For a form with real odd nebentypus (-3/.), conj(a_p) = (-3/p) a_p, so t_p = 0 at p = 2 mod 3 and 2 a_p at p = 1 mod 3.
Here: S_5(Gamma_0(6), (-3/.)) from eta quotients (weight 5, prod delta^(r_delta) in 3 Q^(x2)), its Hecke eigenforms
by T_7, and the comparison with t_p (recomputed exactly) at every prime 11 <= p < 700.
Usage: python -B a4d_weight5_level6.py
"""
import itertools
from fractions import Fraction
import numpy as np
import sympy
from sympy import primerange

NB = 800


def pentagonal(m, N):
    s = [0] * (N + 1)
    k = 0
    while True:
        hit = False
        for kk in ([k, -k] if k else [0]):
            e = m * kk * (3 * kk - 1) // 2
            if e <= N:
                s[e] += 1 - 2 * (abs(kk) % 2); hit = True
        if not hit and k > 0:
            break
        k += 1
    return s


def mul(a, b, N):
    out = [0] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(0, N + 1 - i):
                if b[j]:
                    out[i + j] += x * b[j]
    return out


def inv(a, N):
    r = [0] * (N + 1); r[0] = Fraction(1, a[0])
    for n in range(1, N + 1):
        r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1) if a[k]) / a[0]
    return [int(x) for x in r]


PEN = {m: pentagonal(m, NB) for m in (1, 2, 3, 6)}
PINV = {m: inv(PEN[m], NB) for m in (1, 2, 3, 6)}


def eta_quot(r, N=NB):
    ser = [1] + [0] * N
    for m, e in zip((1, 2, 3, 6), r):
        base = PEN[m] if e > 0 else PINV[m]
        for _ in range(abs(e)):
            ser = mul(ser, base, N)
    lead = sum(d * e for d, e in zip((1, 2, 3, 6), r)) // 24
    return ([0] * lead + ser)[:N + 1]


def orders(r):
    out = []
    for d in (1, 2, 3, 6):
        v = Fraction(0)
        for delta, e in zip((1, 2, 3, 6), r):
            g = int(sympy.gcd(d, delta))
            v += Fraction(g * g * e, d * delta)
        out.append(v / 4)
    return out


cands = []
for r in itertools.product(range(-8, 17), repeat=4):
    if sum(r) != 10:
        continue
    if sum(d * e for d, e in zip((1, 2, 3, 6), r)) % 24 or sum((6 // d) * e for d, e in zip((1, 2, 3, 6), r)) % 24:
        continue
    s = Fraction(1)
    for d, e in zip((1, 2, 3, 6), r):
        s *= Fraction(d) ** e
    q = s / 3
    if not sympy.sqrt(sympy.Rational(q.numerator, q.denominator)).is_rational:
        continue
    if all(o >= 1 for o in orders(r)):
        cands.append(r)
print('cuspidal eta quotients of weight 5, level 6, character (-3/.):', cands)
ser = [eta_quot(r) for r in cands]
Mx = sympy.Matrix([[s[n] for n in range(1, 60)] for s in ser])
print('dimension of their span:', Mx.rank())
basis = []
for s in ser:
    if sympy.Matrix([[b[n] for n in range(1, 60)] for b in basis + [s]]).rank() == len(basis) + 1:
        basis.append(s)
dim = len(basis)
chi = lambda n: 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)
# T_7 with character: (T_p f)_n = a_(pn) + chi(p) p^4 a_(n/p)
pT = 7
nn = 60
A = sympy.Matrix([[basis[j][n] for j in range(dim)] for n in range(1, nn + 1)])
cols = []
for j in range(dim):
    rhs = sympy.Matrix([basis[j][pT * n] + (chi(pT) * pT ** 4 * basis[j][n // pT] if n % pT == 0 else 0) for n in range(1, nn + 1)])
    sol = A.solve_least_squares(rhs)
    assert A * sol == rhs
    cols.append(sol)
T7 = sympy.Matrix.hstack(*cols)
print('T_7 characteristic polynomial on the span:', sympy.factor(T7.charpoly().as_expr()))
ev = T7.eigenvects()
forms = []
for val, m, vecs in ev:
    for vv in vecs:
        f = [sum(vv[j] * basis[j][n] for j in range(dim)) for n in range(NB + 1)]
        f = [sympy.nsimplify(sympy.simplify(x / f[1])) for x in f]
        forms.append((val, f))
for val, f in forms:
    print('eigenform (T_7 = %s): a_2..a_13 =' % val, [sympy.simplify(f[n]) for n in range(2, 14)])


# exact t_p
def qr_table(p):
    t = -np.ones(p, dtype=np.int64); t[0] = 0
    t[(np.arange(1, p) ** 2) % p] = 1
    return t


def tp_value(p):
    chi_ = qr_table(p)
    tt = np.arange(p, dtype=np.int64)
    S3 = 0
    for c in range(1, p):
        if c in (1, 9 % p):
            continue
        QE = ((tt * tt - 4 * c) % p) * ((tt - c - 1) % p) % p * ((tt - c + 3) % p) % p
        a = -1 - int(chi_[QE].sum())
        S3 += a ** 3
    e3 = 1 if p % 3 == 1 else -1
    return -S3 - (1 + e3) * (4 * p + 2)


print('NOTE: T_7 acts as the scalar 26 on this 2-dim space, so the vectors above are NOT Hecke eigenforms; the')
print('comparison below uses the trace of T_p on the space, which needs no eigenbasis.')


# ---- trace of T_p on the 2-dimensional space (no eigenform needed): echelon basis f1 = q + 0 q^2 + ..., f2 = q^2 + ...
NB2 = 1410
PEN = {m: pentagonal(m, NB2) for m in (1, 2, 3, 6)}
PINV = {m: inv(PEN[m], NB2) for m in (1, 2, 3, 6)}
g1 = eta_quot(cands[0], NB2); g2 = eta_quot(cands[1], NB2)
M2 = sympy.Matrix([[g1[1], g1[2]], [g2[1], g2[2]]])
Minv = M2.inv()
f1 = [Minv[0, 0] * x + Minv[0, 1] * y for x, y in zip(g1, g2)]
f2 = [Minv[1, 0] * x + Minv[1, 1] * y for x, y in zip(g1, g2)]
assert (f1[1], f1[2], f2[1], f2[2]) == (1, 0, 0, 1)
for g in (eta_quot(r, NB2) for r in cands[2:]):          # the other quotients lie in the span
    assert all(g[n] == g[1] * f1[n] + g[2] * f2[n] for n in range(NB2 + 1))
bad2 = []
for p in primerange(5, 700):
    trT = f1[p] + f2[2 * p]
    if p >= 11 and trT != tp_value(p):
        bad2.append(p)
print('\nt_p = trace(T_p | S_5(Gamma_0(6), (-3/.))) at every prime 11 <= p < 700: failures', bad2)
print('  trace(T_p) for p = 5, 7, 11, 13, 17, 19:', [f1[p] + f2[2 * p] for p in (5, 7, 11, 13, 17, 19)])
print('  so sum_c a_c^3 = -(1 + (p/3))(4p + 2) - tr(T_p | S_5(Gamma_1(6))) (FINITE-EXACT, 11 <= p < 700)' if not bad2 else '')
