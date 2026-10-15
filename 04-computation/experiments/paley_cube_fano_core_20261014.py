#!/usr/bin/env python3
"""The cube Q_3: the first pure-chi core in the Fano regime (a pre-registered test of HYP-9278).

A pure core has minimum degree >= 3, so |E| >= 3|V|/2.  The Fano regime of HYP-9278 (double cover of P^{|V|-2} with
K = O(-1) or worse) needs |E| <= 2|V| - 4 with |E| even.  Both together force |V| >= 8.  The 3-cube Q_3 (8 vertices,
12 edges) is the first such core: Laman-sparse, bipartite, with no odd-reversal automorphism (every automorphism
reverses a multiple of 4 edges under the bit orientation), so THM-4652 does not kill it.
PREDICTION (HYP-9278, registered before this run): U(Q_3)/(q(q-1)) is a polynomial in q in each class mod 4.

Evaluation in O(q^4): top square 0-1-3-2 (x), bottom square 4-5-7-6 (y), rungs i -> i + 4, bit orientation.
  H(s, t) = sum_y chi(y) chi(y - s) chi(t - y)          (one q x q table),
  G(x) = sum_{a,b} chi(a - x_0) chi(b - x_3) H(a - x_1, b - x_1) H(a - x_2, b - x_2),
  U/(q(q-1)) = sum_{x_2, x_3} chi(x_2) chi(x_3 - 1) chi(x_3 - x_2) G(0, 1, x_2, x_3)    (x_0 = 0, x_1 = 1).
Here y_0 = a and y_3 = b, and the y_1 and y_2 sums are the two H factors.  Cross-checked by brute force at small q.
Session opus-2026-10-14-S35.
"""
import sys, time
import numpy as np
import numba as nb
import sympy as sp
from sympy import primerange


def chi_table(q):
    chi = np.zeros(q, dtype=np.int64)
    for x in range(1, q):
        chi[x] = 1 if pow(x, (q - 1) // 2, q) == 1 else -1
    return chi


@nb.njit(cache=True)
def cube_fast(q, chi):
    H = np.zeros((q, q), dtype=np.int64)
    for s in range(q):
        for t in range(q):
            acc = 0
            for y in range(q):
                acc += chi[y] * chi[(y - s) % q] * chi[(t - y) % q]
            H[s, t] = acc
    tot = 0
    for x2 in range(q):
        for x3 in range(q):
            w = chi[x2] * chi[(x3 - 1) % q] * chi[(x3 - x2) % q]
            if w == 0:
                continue
            G = 0
            for a in range(q):
                ca = chi[a]  # a - x0, x0 = 0
                if ca == 0:
                    continue
                for b in range(q):
                    cb = chi[(b - x3) % q]
                    if cb == 0:
                        continue
                    G += ca * cb * H[(a - 1) % q, (b - 1) % q] * H[(a - x2) % q, (b - x2) % q]
            tot += w * G
    return tot


@nb.njit(cache=True)
def cube_brute(q, chi, eu, ev):
    tot = 0
    v = np.zeros(8, dtype=np.int64); v[1] = 1
    for code in range(q ** 6):
        c = code
        for i in range(6):
            v[2 + i] = c % q; c //= q
        p = 1
        for k in range(eu.shape[0]):
            p *= chi[(v[ev[k]] - v[eu[k]]) % q]
            if p == 0:
                break
        tot += p
    return tot


E = [(u, u | (1 << k)) for u in range(8) for k in range(3) if not u >> k & 1]

if __name__ == '__main__':
    QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    eu = np.array([e[0] for e in E], dtype=np.int64); ev = np.array([e[1] for e in E], dtype=np.int64)
    for q in (5, 7, 11):
        assert cube_fast(q, chi_table(q)) == cube_brute(q, chi_table(q), eu, ev), q
    print('fast evaluator = brute force at q = 5, 7, 11')
    x = sp.symbols('q')
    for cls in (3, 1):
        t0 = time.time()
        qs = [q for q in primerange(5, QMAX + 1) if q % 4 == cls]
        vals = {q: int(cube_fast(q, chi_table(q))) for q in qs}
        print(f'q = {cls} mod 4: U(Q3)/(q(q-1)) at {len(qs)} primes: {[vals[q] for q in qs]}  [{time.time() - t0:.0f} s]', flush=True)
        for deg in range(0, 7):
            base = qs[:deg + 1]
            P = sp.interpolate([(b, vals[b]) for b in base], x) if deg > 0 else sp.Integer(vals[qs[0]])
            if all(P.subs(x, q) == vals[q] for q in qs):
                print(f'   polynomial of degree {deg} fits all {len(qs)} primes: {sp.factor(P)}')
                break
        else:
            print('   NOT a polynomial of degree <= 6 (prediction refuted in this class)')
