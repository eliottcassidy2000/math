#!/usr/bin/env python3
"""2-adic divisibility of the level-32 newform at Paley primes: v2(b(q)) >= min(v2(q + 1), K)?

c_q = |b(q)| = sqrt(M_q + (q + 1)^2) at q = 3 (mod 4) (HYP-9276 identity; square law re-verified here at every prime).
M_q by the THM-4653 FFT evaluator (as in paley_k6p6_two_adic_mersenne_20261014.py).  For every prime q = 3 (mod 4)
in [5, QMAX] this records (q, c_q), and reports, for each s = v2(q + 1), the minimum of v2(c_q), and the largest K
with v2(c_q) >= min(s, K) at all primes.  Writes paley_k6p6_cq_table_20261014.json (q -> c_q).
Session opus-2026-10-14-S35.
"""
import sys, json, math, time
from collections import defaultdict
import numpy as np
import sympy


def chi_table(q):
    c = np.full(q, -1, dtype=np.int64)
    c[np.unique((np.arange(1, q, dtype=np.int64) ** 2) % q)] = 1
    c[0] = 0
    return c


def M_fft(q):
    ci = chi_table(q); c = ci.astype(np.float64)
    x = np.arange(q)
    gi = ci * ci[(x - 1) % q]; g = gi.astype(np.float64)
    t = np.rint(np.fft.ifft(np.fft.fft(g) * np.conj(np.fft.fft(c))).real).astype(np.int64)
    ct = (ci * t).astype(np.float64)
    u = np.rint(np.fft.ifft(np.fft.fft(ct) * np.fft.fft(c)).real).astype(np.int64)
    return -int(np.sum(gi * t * u, dtype=np.int64))


def v2(n):
    return 99 if n == 0 else (n & -n).bit_length() - 1


if __name__ == '__main__':
    QMAX = int(float(sys.argv[1])) if len(sys.argv) > 1 else 200000
    OUT = sys.argv[2] if len(sys.argv) > 2 else '05-knowledge/results/paley_k6p6_cq_table_20261014.json'
    t0 = time.time()
    table = {}
    for q in sympy.primerange(7, QMAX + 1):
        if q % 4 != 3:
            continue
        M = M_fft(q)
        s = M + (q + 1) ** 2
        c = math.isqrt(s)
        assert c * c == s and c <= 2 * q, q
        table[q] = c
    json.dump({str(q): c for q, c in table.items()}, open(OUT, 'w'))
    by_s = defaultdict(list)
    for q, c in table.items():
        by_s[v2(q + 1)].append((v2(c), q, c))
    print(f'{len(table)} primes q = 3 mod 4 in [7, {QMAX}] (square law holds at all)  [{time.time() - t0:.0f} s]')
    Kmax = 99
    for s in sorted(by_s):
        lst = by_s[s]
        mn = min(v for v, _, _ in lst)
        worst = min(lst)
        print(f'  v2(q+1) = {s:2d}: {len(lst):5d} primes, min v2(c_q) = {mn if mn < 99 else "inf"} (e.g. q = {worst[1]}, c = {worst[2]}); '
              f'primes with c_q = 0: {[q for v, q, c in lst if c == 0]}')
        if mn < s:
            Kmax = min(Kmax, mn)
    print(f'largest K with v2(c_q) >= min(v2(q+1), K) at every prime: K = {Kmax}')
