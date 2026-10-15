#!/usr/bin/env python3
"""Audit (c2): the wheel traces are Greene's hypergeometric values, and the T_6 law is Koike's conjecture
(proved by Frechette-Ono-Papanikolas, IMRN 2004, Corollary 1.2).

Greene / FOP normalisation:  (A over B) = (1/p) sum_x A(x) Bbar(x - 1),
    nF_{n-1}(1) := p/(p-1) sum_chi (phi chi over chi)^n   (phi quadratic).
Substituting t = x/(x-1) gives (phi chi over chi) = (1/p) sum_t chi(t-1) (phi chi)(t), hence
    T_n = Tr((DS)^n) = sum_psi (sum_t chi(t-1) psi(t))^n = (p-1) p^(n-1) nF_{n-1}(1)      (exact, all n).
Known: p^3 4F3(1) = -a(p) - p (Ahlgren-Ono 2000);  b6(p) = -p^5 6F5(1) + p^4 4F3(1) + (1 - phi(-1)) p^2 (Koike, FOP).
Consequence: T_6 = (p-1)(-b6(p) - p a(p) - phi(-1) p^2), i.e. (p-1)(p^2 - p a - b6) for p = 3 mod 4 (S35's law).
Checked here with exact matrix traces at all primes 5 <= p <= 200 (both classes).
"""
import numpy as np
import sympy as sp
from audit_c_wheels import traces, level8_systems


def main():
    N = 402
    _, _, s4 = level8_systems(4, N)
    _, _, s6 = level8_systems(6, N)
    a = s4[-4]
    b6 = s6[20]
    bad4, bad6 = [], []
    n = 0
    for p in sp.primerange(5, 201):
        T = traces(p)
        phim1 = 1 if p % 4 == 1 else -1
        n += 1
        if T[4] != -(p - 1) * (p + int(a[p])):
            bad4.append(p)
        if T[6] != (p - 1) * (-int(b6[p]) - p * int(a[p]) - phim1 * p * p):
            bad6.append(p)
    print(f'{n} primes 5..199 (both classes): T_4 = -(p-1)(p + a_p) failures {bad4}; '
          f'T_6 = (p-1)(-b6_p - p a_p - phi(-1) p^2) [Koike/FOP] failures {bad6}')
    print('b6(13), b6(17) =', int(b6[13]), int(b6[17]))


if __name__ == '__main__':
    main()
