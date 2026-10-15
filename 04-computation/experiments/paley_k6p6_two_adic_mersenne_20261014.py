#!/usr/bin/env python3
"""The 2-adic law of the K_6 - P_6 terminal, its mechanism, and the Mersenne primes.

Uses the O(q log q) Legendre-bilinear evaluator of THM-4653 (mac-mini-2026-10-14-wilson; method credited there),
re-implemented with an exact int64 final sum:
    M_q = - sum_x g(x) t(x) u(x),  g(x) = chi(x) chi(x - 1),  t = g (correlated with) chi,  u = (chi t) * chi.
In this session's unrestricted normalisation, Y_q = U(K_6 - P_6)/(q(q - 1)) = 1 - chi(-1) M_q at Paley primes
(Y_q - 1 = M_q for q = 3 mod 4), and HYP-9276 (mac-mini) reads M_q = |b(q)|^2 - (q + 1)^2 there, b = the weight-3
level-32 newform.  Consequences tested here, far beyond the 30,000 of HYP-9276:
  (S) M_q + (q + 1)^2 is a perfect square c_q^2 with 0 <= c_q <= 2q   (q = 3 mod 4);
  (C) c_q = +-(q + 1) mod 2^m with m = 4 (q = 11 mod 16), 5 (q = 3 mod 16), 6 (q = 7 mod 8) -- the mechanism of the
      2^7 / 2^10 laws (v2(M_q) = v2(c - q - 1) + v2(c + q + 1));
  (Mersenne) q = 2^p - 1 for p = 5, 7, 13, 17, 19: c_q, v2(M_q).
Session opus-2026-10-14-S35.
"""
import sys, time, math
import numpy as np
import sympy


def chi_table(q):
    c = np.full(q, -1, dtype=np.int64)
    c[np.unique((np.arange(1, q, dtype=np.int64) ** 2) % q)] = 1
    c[0] = 0
    return c


def M_fft(q):
    ci = chi_table(q)
    c = ci.astype(np.float64)
    x = np.arange(q)
    gi = ci * ci[(x - 1) % q]
    g = gi.astype(np.float64)
    t = np.rint(np.fft.ifft(np.fft.fft(g) * np.conj(np.fft.fft(c))).real).astype(np.int64)
    ct = (ci * t).astype(np.float64)
    u = np.rint(np.fft.ifft(np.fft.fft(ct) * np.fft.fft(c)).real).astype(np.int64)
    return -int(np.sum(gi * t * u, dtype=np.int64))


def v2(n):
    return 99 if n == 0 else (n & -n).bit_length() - 1


def depth(c, q):
    return max(v2(c - q - 1), v2(c + q + 1))


def main():
    lo = int(float(sys.argv[1])) if len(sys.argv) > 1 else 30000
    hi = int(float(sys.argv[2])) if len(sys.argv) > 2 else 200000
    t0 = time.time()
    # sanity against this session's exact O(q^3) table (Y_q - 1 = M_q for q = 3 mod 4)
    exact = {11: -127, 19: -255, 23: 1025, 43: -1151, 127: -16383, 191: 65537, 499: 633601, 691: 1211137}
    for q, y in exact.items():
        assert M_fft(q) == y - 1, q
    print('FFT evaluator agrees with the exact O(q^3) table at q =', sorted(exact))
    from collections import defaultdict
    mind = defaultdict(lambda: 99); minv = defaultdict(lambda: 99); n = 0; bad = []
    for q in sympy.primerange(lo, hi + 1):
        if q % 4 != 3:
            continue
        M = M_fft(q)
        s = M + (q + 1) ** 2
        c = math.isqrt(s) if s >= 0 else -1
        if c < 0 or c * c != s or c > 2 * q:
            bad.append(q); continue
        n += 1
        r = q % 16
        mind[r] = min(mind[r], depth(c, q))
        minv[r] = min(minv[r], v2(M))
    print(f'primes q = 3 mod 4 in [{lo}, {hi}]: {n};  perfect-square law (S) failures: {bad[:10]} ({len(bad)})  [{time.time() - t0:.0f} s]')
    for r in (3, 7, 11, 15):
        print(f'  q = {r:2d} mod 16: min depth of c = +-(q+1) mod 2^m: {mind[r]};  min v2(M_q): {minv[r]}')
    print('Mersenne primes q = 2^p - 1:')
    for p in (3, 5, 7, 13, 17, 19):
        q = 2 ** p - 1
        M = M_fft(q)
        s = M + (q + 1) ** 2
        c = math.isqrt(s)
        print(f'  p = {p:2d}, q = {q}: M_q = {M}, perfect square: {c * c == s}, c_q = |b(q)| = {c} = {c // 2 ** v2(c) if c else 0} * 2^{v2(c)}, '
              f'v2(M_q) = {v2(M)}, depth = {depth(c, q)}', flush=True)


if __name__ == '__main__':
    main()
