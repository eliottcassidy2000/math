#!/usr/bin/env python3
"""The wheel ladder: U(W_n) = q * Tr((D S)^n) for the Paley tournament, D = diag(chi), S_{xy} = chi(y - x).

Wheels W_n (hub + n-cycle rim, 2n edges, n + 1 vertices) satisfy |E| = 2|V| - 2: they are the classical 2D rigidity
circuits, and the double cover of P^{n-1} branched along their 2n graphic hyperplanes is Calabi-Yau.  The prediction
of the S35 rigidity/Calabi-Yau dictionary is motivic weight n - 1, i.e. modular contributions of weight n.
For q = 3 (mod 4), T_n = Tr((DS)^n) vanishes for odd n (S^T = -S).  Checks and fits:
  T_4 = -(q - 1)(q + a_q),  a_q = coefficient of eta(2z)^4 eta(4z)^4 (S34/Ahlgren-Ono);
  T_6, T_8: residual after subtracting the best polynomial-plus-a_q fit, against Hecke traces of weight-6 / weight-8
  eta-quotient eigenforms of level 8 and 16 (computed in paley_wheel_ladder_forms; here only the raw table).
Session opus-2026-10-14-S35.
"""
import sys
import numpy as np
from sympy import primerange


def chi_table(q):
    chi = np.zeros(q, dtype=np.int64)
    for x in range(1, q):
        chi[x] = 1 if pow(x, (q - 1) // 2, q) == 1 else -1
    return chi


def traces(q, nmax=10):
    chi = chi_table(q)
    idx = np.arange(q)
    S = chi[(idx[None, :] - idx[:, None]) % q].astype(object)
    DS = (chi[:, None].astype(object) * S)
    P = np.identity(q, dtype=object)
    out = {}
    for n in range(1, nmax + 1):
        P = P.dot(DS)
        out[n] = int(sum(P[i, i] for i in range(q)))
    return out


def beukers(N=400):
    P = [0] * (N + 1); P[0] = 1
    for n in range(1, N + 1):
        for step in (2 * n, 4 * n):
            if step > N:
                continue
            for _ in range(4):
                for i in range(N, step - 1, -1):
                    P[i] -= P[i - step]
    a = [0] * (N + 2)
    for i in range(N + 1):
        a[i + 1] = P[i]
    return a


if __name__ == '__main__':
    QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 131
    a = beukers()
    for q in primerange(5, QMAX + 1):
        T = traces(q, 8)
        chk = ''
        if q % 4 == 3:
            chk = '  T4 = -(q-1)(q+a_q): ' + str(T[4] == -(q - 1) * (q + a[q]))
        print(f'q = {q:4d} ({q % 4} mod 4): ' + '  '.join(f'T{n} = {T[n]}' for n in range(2, 9)) + chk, flush=True)
