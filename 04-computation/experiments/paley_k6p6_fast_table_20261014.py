#!/usr/bin/env python3
"""Fast table of the K_6 - P_6 terminal Y_q = U(K_6 - P_6)/(q(q-1)) at odd primes q, in O(q^3).

Orientation and labelling as in paley_six_vertex_cores_20261014.py / the S35 note:
  E = (0,1),(0,3),(0,4),(0,5),(1,4),(1,5),(2,3),(2,4),(2,5),(3,4)   (edge (u,v) carries chi(x_v - x_u)),
  non-edges 0-2-1-3-5-4 (the removed P_6).
With x_0 = 0, x_1 = 1 (AGL(1,q) normalisation; (0,1) is an edge):
  Y_q = sum_{x2} g(x2) h(x2),
  g(x2) = sum_{x5} chi(x5) chi(x5 - 1) chi(x5 - x2)                    (a Legendre-curve trace),
  h(x2) = sum_{x3,x4} chi(x3) chi(x4) chi(x4 - 1) chi(x3 - x2) chi(x4 - x2) chi(x4 - x3).
Cross-checked against the brute-force values of the note (q <= 191).  Also prints v_2(Y_q - 1) and q mod 16.
Session opus-2026-10-14-S35.
"""
import sys, time
import numpy as np
import numba as nb
from sympy import primerange


@nb.njit(cache=True)
def Yq(q, chi):
    tot = 0
    for x2 in range(q):
        g = 0
        for x5 in range(q):
            g += chi[x5] * chi[(x5 - 1) % q] * chi[(x5 - x2) % q]
        if g == 0:
            continue
        h = 0
        for x4 in range(q):
            w = chi[x4] * chi[(x4 - 1) % q] * chi[(x4 - x2) % q]
            if w == 0:
                continue
            s = 0
            for x3 in range(q):
                s += chi[x3] * chi[(x3 - x2) % q] * chi[(x4 - x3) % q]
            h += w * s
        tot += g * h
    return tot


def chi_table(q):
    chi = np.zeros(q, dtype=np.int64)
    for x in range(1, q):
        chi[x] = 1 if pow(x, (q - 1) // 2, q) == 1 else -1
    return chi


def v2(n):
    if n == 0:
        return 99
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k


KNOWN = {3: 1, 7: 1, 11: -127, 19: -255, 23: 1025, 31: 1, 43: -1151, 47: -2047, 59: -895, 67: -4607, 71: -2047,
         79: -6143, 83: 12545, 103: 12289, 107: 12673, 127: -16383, 131: -12799, 139: 18817, 151: -15359, 191: 65537,
         5: -11, 13: 53, 17: 69, 29: -587, 37: -395, 41: -1403, 53: 1653, 61: 3125, 73: -827, 89: -6843, 97: -9019,
         101: -1163, 109: -11339}

if __name__ == '__main__':
    QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    t0 = time.time()
    for q in primerange(3, QMAX + 1):
        y = int(Yq(q, chi_table(q)))
        chk = '' if q not in KNOWN else ('  (matches brute force)' if KNOWN[q] == y else f'  MISMATCH vs {KNOWN[q]}')
        print(f'q = {q:4d}  q mod 16 = {q % 16:2d}  Y_q = {y:12d}  v2(Y_q - 1) = {v2(y - 1):2d}  Y/q^2 = {y / q**2:+.4f}{chk}', flush=True)
    print(f'[{time.time() - t0:.1f} s]')
