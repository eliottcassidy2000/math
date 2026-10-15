#!/usr/bin/env python3
"""Deep 2-adic classes for the valuation law gcd(q + 1, 128) | b(q): primes q = -1 (mod 2^9) up to QMAX.

For such q, v2(q + 1) >= 9, so the law predicts v2(|b(q)|) >= 7 exactly at the saturation level.
|b(q)| = sqrt(M_q + (q + 1)^2) with M_q from the THM-4653 FFT evaluator (as in paley_k6p6_two_adic_mersenne_20261014.py).
Session opus-2026-10-14-S35.
"""
import sys, math, time
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
    QMAX = int(float(sys.argv[1])) if len(sys.argv) > 1 else 2_000_000
    t0 = time.time()
    rows = []
    q = 511
    while q <= QMAX:
        if sympy.isprime(q):
            M = M_fft(q)
            s = M + (q + 1) ** 2
            c = math.isqrt(s)
            assert c * c == s and c <= 2 * q, q
            rows.append((q, v2(q + 1), c, v2(c)))
        q += 512
    bad = [r for r in rows if r[3] < 7]
    print(f'{len(rows)} primes q = -1 mod 512 up to {QMAX}: square law holds at all; '
          f'min v2(|b(q)|) = {min(r[3] for r in rows)}; violations of v2 >= 7: {bad}  [{time.time() - t0:.0f} s]')
    from collections import Counter
    print('distribution of v2(|b(q)|):', sorted(Counter(min(r[3], 20) for r in rows).items()))
    print('zeros b(q) = 0:', [r[0] for r in rows if r[2] == 0])
