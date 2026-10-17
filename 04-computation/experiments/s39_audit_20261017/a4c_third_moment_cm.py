#!/usr/bin/env python3
"""S39 audit, addendum (3): is the third moment of the four-bar family outside the modular span?

From a4b: t_p := -sum_c a_c^3 - (1 + (p/3))(4p + 2) (Lefschetz: the trace on the 4-dimensional weight-4 space
H^1(P^1, j_* Sym^3 R^1)) vanishes at every prime p = 2 mod 3 tested.  That suggests CM by Q(sqrt-3).  Test: for
p = 1 mod 3 write p = x^2 + 3y^2, pi = x + y sqrt(-3), A2 = pi^2 + pibar^2 = 2(x^2 - 3y^2),
A4 = pi^4 + pibar^4 = 2(x^4 - 18 x^2 y^2 + 9 y^4) (Hecke characters of Q(sqrt-3): weights 3 and 5), and search small
integer combinations t_p = alpha A4 + beta p A2 + gamma p^2 + delta p^2 (p/3)... over all primes 11 <= p < 700.
Same test for the twisted second moment residual T_p (6-dim weight 3) against p a_p(E) for a few curves and the
weight-4 forms is NOT attempted here (out of scope).
Usage: python -B a4c_third_moment_cm.py
"""
import itertools
import numpy as np
from sympy import primerange


def qr_table(p):
    t = -np.ones(p, dtype=np.int64); t[0] = 0
    t[(np.arange(1, p) ** 2) % p] = 1
    return t


def kr(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def tp_value(p):
    chi = qr_table(p)
    tt = np.arange(p, dtype=np.int64)
    S3 = 0
    for c in range(1, p):
        if c in (1, 9 % p):
            continue
        QE = ((tt * tt - 4 * c) % p) * ((tt - c - 1) % p) % p * ((tt - c + 3) % p) % p
        a = -(2 + int(chi[QE].sum())) + 1 + 0     # p + 1 - (p + 2 + sum) = -1 - sum
        a = -1 - int(chi[QE].sum())
        S3 += a ** 3
    e3 = kr(p, 3)
    return -S3 - (1 + e3) * (4 * p + 2)


def rep_x2_3y2(p):
    for y in range(1, int((p / 3) ** 0.5) + 2):
        x2 = p - 3 * y * y
        if x2 < 0:
            break
        x = int(round(x2 ** 0.5))
        if x * x == x2:
            return x, y
    return None


data = {}
for p in primerange(11, 700):
    data[p] = tp_value(p)
inert = [p for p in data if p % 3 == 2]
print('t_p at primes p = 2 mod 3 (inert in Q(sqrt-3)): all zero?', all(data[p] == 0 for p in inert), '(%d primes)' % len(inert))
split = [p for p in data if p % 3 == 1]
feat = {}
for p in split:
    x, y = rep_x2_3y2(p)
    A2 = 2 * (x * x - 3 * y * y)
    A4 = 2 * (x ** 4 - 18 * x * x * y * y + 9 * y ** 4)
    feat[p] = (A4, p * A2, p * p)
sols = []
R = range(-6, 7)
for al, be, ga in itertools.product(R, R, R):
    if all(data[p] == al * feat[p][0] + be * feat[p][1] + ga * feat[p][2] for p in split):
        sols.append((al, be, ga))
print('integer solutions (alpha, beta, gamma) of t_p = alpha A4 + beta p A2 + gamma p^2 on all %d split primes < 700:' % len(split), sols)
# least squares diagnostic
M = np.array([feat[p] for p in split], dtype=float); b = np.array([data[p] for p in split], dtype=float)
coef, res, rk, sv = np.linalg.lstsq(M, b, rcond=None)
print('least-squares coefficients:', np.round(coef, 6), ' max residual:', float(np.max(np.abs(M @ coef - b))))
print('samples (p, t_p, A4, p A2):', [(p, data[p], feat[p][0], feat[p][1]) for p in split[:6]])


# ---- second attempt: primary pi (pi = 2 mod 3 in Z[w]) and the cubic character of conductor 2 (lambda = pi mod 2)
def zmul(x, y):                      # (a + b w)(c + d w), w^2 = -1 - w
    a, b = x; c, d = y
    return (a * c - b * d, a * d + b * c - b * d)


def zpow(x, k):
    r = (1, 0)
    for _ in range(k):
        r = zmul(r, x)
    return r


def ztr(x):
    return 2 * x[0] - x[1]


def primary(p):
    R = int((4 * p / 3) ** 0.5) + 2
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            if a * a - a * b + b * b == p and a % 3 == 2 and b % 3 == 0:
                return (a, b)
    raise RuntimeError(p)


W = (0, 1)
feats2 = {}
for p in split:
    pi = primary(p)
    lam = {(1, 0): (1, 0), (0, 1): W, (1, 1): zmul(W, W)}[(pi[0] % 2, pi[1] % 2)]
    f = []
    for k in range(3):
        lk = zpow(lam, k)
        f.append(ztr(zmul(zpow(pi, 4), lk)))
    for k in range(3):
        lk = zpow(lam, k)
        f.append(p * ztr(zmul(zpow(pi, 2), lk)))
    f.append(p * p)
    feats2[p] = f
names = ['Tr pi^4', 'Tr pi^4 lam', 'Tr pi^4 lam^2', 'p Tr pi^2', 'p Tr pi^2 lam', 'p Tr pi^2 lam^2', 'p^2']
M2 = np.array([feats2[p] for p in split], dtype=float)
coef2, *_ = np.linalg.lstsq(M2, b, rcond=None)
print('\nprimary pi + conductor-2 cubic twists: least-squares coefficients', dict(zip(names, np.round(coef2, 6))),
      ' max residual %.3g' % float(np.max(np.abs(M2 @ coef2 - b))))
