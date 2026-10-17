#!/usr/bin/env python3
"""Second moments of the two elliptic modular surfaces behind the Paley and the unit-distance ladders (S39).

(1) Tournaments: the Legendre family E_l : y^2 = x(x-1)(x-l). Its chi(l)-twisted second moment is the level-8 form of
    the Paley 4-wheel (S34; Ahlgren-Ono's 4F3(1) identity, CITED from memory):
        sum_(l != 0, 1) chi(l) a_l^2 = -1 - p - c_p,   sum_l c_p q^p ... = eta(2z)^4 eta(4z)^4;
    (the untwisted second moment is not of this form).
(2) Unit distances: the four-bar family E'_c : y^2 = (t^2 - 4c)(t - c - 1)(t - c + 3) (THM-4663 (9), the Gamma_1(6)
    surface). Its second moment is the hexagon form (THM-4663 (10)):
        sum_(c != 0, 1, 9) a_c^2 = (p + 1)(p - 4) - b_p,   b = eta(z)^2 eta(2z)^2 eta(3z)^2 eta(6z)^2.
Both checked exactly at every prime 5 <= p < PMAX (11 <= p for (2)).
Usage: python legendre_fourbar_second_moments_20261017.py [PMAX]
Session opus-2026-10-17-S39.
"""
import sys
import numpy as np
from sympy import jacobi_symbol, primerange


def kr(d, n):
    return int(jacobi_symbol(d % n, n)) if d % n else 0


def eta_product(factors, N):
    s = np.zeros(N + 1, dtype=object); s[0] = 1
    for m, e in factors:
        for n in range(1, N // m + 1):
            for _ in range(e):
                s[m * n:] = s[m * n:] - s[:N + 1 - m * n].copy()
    return s


if __name__ == '__main__':
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    s8 = eta_product([(2, 4), (4, 4)], PMAX)          # eta(2z)^4 eta(4z)^4 = q * s8
    s6 = eta_product([(1, 2), (2, 2), (3, 2), (6, 2)], PMAX)
    bad1, bad2 = [], []
    for p in primerange(5, PMAX):
        chi = [kr(x, p) for x in range(p)]
        T = 0
        for l in range(2, p):
            a = -sum(chi[(x * (x - 1) * (x - l)) % p] for x in range(p))
            T += chi[l] * a * a
        if T != -1 - p - int(s8[p - 1]):
            bad1.append(p)
        if p >= 11:
            S2 = 0
            for c in range(1, p):
                if c in (1, 9 % p):
                    continue
                a = p + 1 - (2 + sum(1 + chi[((t * t - 4 * c) * (t - c - 1) * (t - c + 3)) % p] for t in range(p)))
                S2 += a * a
            if S2 != (p + 1) * (p - 4) - int(s6[p - 1]):
                bad2.append(p)
    print('Legendre: sum chi(l) a_l^2 = -1 - p - c_p (eta(2z)^4 eta(4z)^4) for every prime 5 <= p < %d; failures %s' % (PMAX, bad1))
    print('four-bar: sum a_c^2 = (p+1)(p-4) - b_p (eta1^2 eta2^2 eta3^2 eta6^2) for every prime 11 <= p < %d; failures %s' % (PMAX, bad2))
    assert not bad1 and not bad2
    print('ALL CHECKS PASSED')
