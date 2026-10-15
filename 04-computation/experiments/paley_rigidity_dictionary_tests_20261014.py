#!/usr/bin/env python3
"""Tests of the rigidity / Calabi-Yau dictionary (HYP-9278) on six- and five-vertex graphs.

U6(G) = U(G)/(q(q-1)) = sum over x_2..x_{n-1} of prod_{(u,v) in E} chi(x_v - x_u), with x_0 = 0, x_1 = 1 ((0,1) in E).
Rows:
  - Laman-sparse, not series-parallel (|E| = 2|V| - 4, Fano double covers): K_{3,3} - e, triangular prism - e.
    Prediction: polynomial in q.  Test: interpolate on |V| + 1 primes, predict the rest exactly (each class mod 4).
  - circuit count but containing a smaller circuit: K_4 plus a degree-2 vertex (contains K_4), the octahedron minus the
    two edges 35, 45 (contains the wheel W_4).  Prediction: the smaller circuit's form (CM eta(4z)^6, i.e.
    b_q = 4a^2 - 2q for q = a^2 + b^2; and q + a_q with a = eta(2z)^4 eta(4z)^4).
Session opus-2026-10-14-S35.
"""
import numpy as np
import numba as nb
import sympy as sp
from sympy import primerange


@nb.njit(cache=True)
def U6(q, chi, eu, ev, nv):
    tot = 0
    m = nv - 2
    v = np.zeros(nv, dtype=np.int64); v[1] = 1
    total = q ** m
    for code in range(total):
        c = code
        for i in range(m):
            v[2 + i] = c % q; c //= q
        p = 1
        for k in range(eu.shape[0]):
            p *= chi[(v[ev[k]] - v[eu[k]]) % q]
            if p == 0:
                break
        tot += p
    return tot


def chi_table(q):
    chi = np.zeros(q, dtype=np.int64)
    for x in range(1, q):
        chi[x] = 1 if pow(x, (q - 1) // 2, q) == 1 else -1
    return chi


def beukers(N=200):
    P = [0] * (N + 1); P[0] = 1
    for n in range(1, N + 1):
        for step in (2 * n, 4 * n):
            if step > N:
                continue
            for _ in range(4):
                for i in range(N, step - 1, -1):
                    P[i] -= P[i - step]
    return [0] + P


def cm_b(q):
    for a in range(1, int(q ** 0.5) + 1, 2):
        b2 = q - a * a
        if b2 > 0 and int(round(b2 ** 0.5)) ** 2 == b2:
            return 4 * a * a - 2 * q
    return 0


GRAPHS = [
    ('K33 - e (6 v, 8 e; Laman-sparse, Fano)', 6, [(0, 1), (0, 3), (0, 5), (2, 1), (2, 3), (2, 5), (4, 1), (4, 3)], 'poly'),
    ('prism - e (6 v, 8 e; Laman-sparse, Fano)', 6, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4)], 'poly'),
    ('K4 + degree-2 vertex (5 v, 8 e; contains K4)', 5, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (4, 0), (4, 1)], 'K4'),
    ('octahedron - {35, 45} (6 v, 10 e; contains W4)', 6, [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 5), (2, 4), (2, 5), (3, 4)], 'W4'),
]

if __name__ == '__main__':
    a = beukers()
    x = sp.symbols('q')
    for name, nv, E, kind in GRAPHS:
        eu = np.array([e[0] for e in E], dtype=np.int64); ev = np.array([e[1] for e in E], dtype=np.int64)
        for cls in (3, 1):
            qs = [q for q in primerange(5, 80 if nv == 6 else 140) if q % 4 == cls]
            vals = {q: int(U6(q, chi_table(q), eu, ev, nv)) for q in qs}
            if kind == 'poly':
                base = qs[:nv + 1]
                P = sp.interpolate([(b, vals[b]) for b in base], x)
                ok = all(P.subs(x, q) == vals[q] for q in qs[nv + 1:])
                print(f'{name}, q = {cls} mod 4: polynomial through {nv + 1} primes predicts the other {len(qs) - nv - 1}: {ok};  U6 = {sp.factor(P)}')
            elif kind == 'K4':
                ok = all(vals[q] == (0 if q % 4 == 3 else -cm_b(q)) for q in qs)
                print(f'{name}, q = {cls} mod 4: U6 = -b_q (CM eta(4z)^6; 0 for q = 3 mod 4) at all {len(qs)} primes: {ok}')
            else:
                if cls == 3:
                    ok = all(vals[q] == q + a[q] for q in qs)
                    print(f'{name}, q = 3 mod 4: U6 = q + a_q (the W4 form) at all {len(qs)} primes: {ok}')
                else:
                    print(f'{name}, q = 1 mod 4: values {[vals[q] for q in qs]}')
