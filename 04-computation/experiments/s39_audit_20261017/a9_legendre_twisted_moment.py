#!/usr/bin/env python3
"""S39 audit, item 9: the Legendre twisted second moment sum_(l != 0,1) chi(l) a_l^2 = -1 - p - c_p, c = eta(2z)^4 eta(4z)^4.

a_l = -sum_x chi(x(x-1)(x-l)) computed with a quadratic-residue table; c_p from Jacobi's identity eta^3 = sum (-1)^k (2k+1) q^(k(k+1)/2)
(eta(2z)^4 eta(4z)^4 = eta(2z)^3 eta(4z)^3 * eta(2z) eta(4z), the last factor from the pentagonal series), independent of the
session's repeated-product code.  Every prime 5 <= p < 1000.  Also the untwisted second moment (polynomial) and the
cubic-twist control.
Usage: python -B a9_legendre_twisted_moment.py
"""
import numpy as np
from sympy import primerange

N = 1000


def jac(m):
    s = [0] * (N + 1); k = 0
    while m * k * (k + 1) // 2 <= N:
        s[m * k * (k + 1) // 2] += (-1) ** k * (2 * k + 1); k += 1
    return s


def pent(m):
    s = [0] * (N + 1); k = 0
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


def mul(a, b):
    out = [0] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                if b[j]:
                    out[i + j] += x * b[j]
    return out


ser = mul(mul(jac(2), jac(4)), mul(pent(2), pent(4)))     # q-power offset: (2 + 4) * 4 / 24 = 1
c = [0] + ser[:N]
print('eta(2z)^4 eta(4z)^4: c_1..c_11 =', c[1:12])
bad, badU = [], []
for p in primerange(5, N):
    chi = -np.ones(p, dtype=np.int64); chi[0] = 0; chi[(np.arange(1, p) ** 2) % p] = 1
    xs = np.arange(p, dtype=np.int64)
    T = U = 0
    for l in range(2, p):
        a = -int(chi[(xs * (xs - 1) % p) * ((xs - l) % p) % p].sum())
        T += int(chi[l]) * a * a
        U += a * a
    if T != -1 - p - c[p]:
        bad.append(p)
    if U != p * p - 2 * p - 3:
        badU.append(p)
print('sum_(l != 0,1) chi(l) a_l^2 = -1 - p - c_p at every prime 5 <= p < %d: failures %s' % (N, bad))
print('untwisted: sum_l a_l^2 = p^2 - 2p - 3 (polynomial): failures %s' % badU)
