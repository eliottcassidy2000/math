#!/usr/bin/env python3
"""A_8 for the Paley tournament: sum over 9 distinct points v_0..v_8 of prod_{i=1..8} chi(v_i - v_(i-1)), v_0 = 0 fixed
(times q).  Tests the scoped conjecture of HYP-9274 (C): A_8 = q(q-1)(q-3) * quadratic with leading coefficient 14
(Catalan C_4, THM-438; A_(2k) ~ C_k q^(k+1)).  Session opus-2026-10-14-S34.
"""
import sys
import numpy as np
import numba as nb
from fractions import Fraction

QS = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [11, 19, 23, 31]


@nb.njit(cache=True)
def A8(q, chi):
    used = np.zeros(q, dtype=np.bool_)
    used[0] = True
    path = np.zeros(9, dtype=np.int64)
    sgn = np.zeros(9, dtype=np.int64)
    sgn[0] = 1
    nxt = np.zeros(10, dtype=np.int64)
    tot = 0
    depth = 1
    nxt[1] = 1
    while depth >= 1:
        v = nxt[depth]
        if v >= q:
            depth -= 1
            if depth >= 1:
                used[path[depth]] = False
                nxt[depth] += 1
            continue
        if used[v]:
            nxt[depth] += 1
            continue
        s = sgn[depth - 1] * chi[(v - path[depth - 1]) % q]
        if depth == 8:
            tot += s
            nxt[depth] += 1
            continue
        path[depth] = v
        sgn[depth] = s
        used[v] = True
        depth += 1
        nxt[depth] = 1
    return tot * q


vals = {}
for q in QS:
    chi = np.zeros(q, dtype=np.int64)
    sq = set((x * x) % q for x in range(1, q))
    for x in range(1, q):
        chi[x] = 1 if x in sq else -1
    a = int(A8(q, chi))
    vals[q] = a
    base = q * (q - 1) * (q - 3)
    print(f'   q = {q}: A_8 = {a},  A_8/(q(q-1)(q-3)) = {Fraction(a, base)},  A_8/q^5 = {a / q ** 5:.4f}', flush=True)
qs = sorted(vals)
# A_(2k) ~ C_k q^(k+1) (THM-438), so A_8/(q(q-1)(q-3)) has degree 2 with leading coefficient C_4 = 14:
# g(q) = 14 q^2 + b q + c through the first two primes, tested on the rest
if len(qs) >= 3:
    q1, q2 = qs[0], qs[1]
    y1 = Fraction(vals[q1], q1 * (q1 - 1) * (q1 - 3)) - 14 * q1 * q1
    y2 = Fraction(vals[q2], q2 * (q2 - 1) * (q2 - 3)) - 14 * q2 * q2
    b = (y1 - y2) / (q1 - q2)
    c = y1 - b * q1
    print(f'   quadratic through {qs[:2]}: A_8/(q(q-1)(q-3)) = 14 q^2 + ({b}) q + ({c})')
    for q in qs[2:]:
        g = 14 * q * q + b * q + c
        print(f'   check q = {q}: predicted A_8 = {g * q * (q - 1) * (q - 3)}, exact {vals[q]}, match: {g * q * (q - 1) * (q - 3) == vals[q]}')
    print('   factored: A_8 = 2 q (q-1)(q-3)(q-7)(7q-55)  [14q^2 - 208q + 770 = 2(q-7)(7q-55)]:',
          all(vals[q] == 2 * q * (q - 1) * (q - 3) * (q - 7) * (7 * q - 55) for q in qs))
