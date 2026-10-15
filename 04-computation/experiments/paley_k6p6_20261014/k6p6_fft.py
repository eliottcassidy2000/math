#!/usr/bin/env python3
"""O(q log q) evaluation of the Paley six-vertex terminal V_q = I_{K6-P6}(q) through the Legendre family.

Identity (THM in this session; derivation in the results note):
  V_q = q(q-1) W_q,  W_q = -chi(-1) M_q + c(q),  c(q) = 2q - 6 (q = 3 mod 4), -2q + 20 (q = 1 mod 4),
  M_q = - sum_{a,x} chi(a) chi(x) chi(x-1) chi(x-a) t(a) t(x),   t(a) = sum_y chi(y) chi(y-1) chi(y-a)  (= -a_q(E_a)).
With g(y) = chi(y) chi(y-1):  t = correlation of g with chi,  u = (chi t) * chi (cyclic convolution on F_q),
M_q = - sum_x g(x) t(x) u(x).  FFTs in float64 with exact rounding (|entries| <= q^2, checked against the exact
O(q^3) values below 1500).

Usage: python3 k6p6_fft.py QMAX > k6p6_fft.out
"""
import sys
import json
import time
import numpy as np
import sympy


def chi_table(q):
    c = np.full(q, -1, dtype=np.int64)
    c[np.unique((np.arange(1, q, dtype=np.int64) ** 2) % q)] = 1
    c[0] = 0
    return c


def M_fft(q):
    c = chi_table(q).astype(np.float64)
    x = np.arange(q)
    g = c * c[(x - 1) % q]
    # t(a) = sum_y g(y) chi(y - a): correlation -> ifft(fft(g) * conj(fft(chi)))
    t = np.rint(np.fft.ifft(np.fft.fft(g) * np.conj(np.fft.fft(c))).real)
    ct = c * t
    u = np.rint(np.fft.ifft(np.fft.fft(ct) * np.fft.fft(c)).real)     # u(x) = sum_a chi(a) t(a) chi(x - a)
    return -int(np.rint((g * t * u).sum()))


def W_from_M(q, M):
    return (M + 2 * q - 6) if q % 4 == 3 else (-M - 2 * q + 20)


def main():
    QMAX = int(float(sys.argv[1])) if len(sys.argv) > 1 else 20000
    exact = {int(k): v for k, v in json.load(open('k6p6_W_allprimes.json')).items()}
    t0 = time.time()
    out = {}
    nchk = 0
    for q in sympy.primerange(5, QMAX + 1):
        M = M_fft(q)
        W = W_from_M(q, M)
        if q in exact:
            assert W == exact[q], q
            nchk += 1
        out[q] = (M, W)
    json.dump({str(q): [m, w] for q, (m, w) in out.items()}, open('k6p6_fft.json', 'w'))
    v73 = [q for q, (m, w) in out.items() if q % 8 == 3 and m % 128]
    v77 = [q for q, (m, w) in out.items() if q % 8 == 7 and m % 1024]
    print(f"{len(out)} primes 5 <= q <= {QMAX}; FFT values equal the exact O(q^3) values at all {nchk} primes below 1500  [{time.time() - t0:.1f}s]")
    print(f"2-adic divisibility: q = 3 mod 8 with 2^7 not dividing M_q: {v73[:10]} ({len(v73)}); "
          f"q = 7 mod 8 with 2^10 not dividing M_q: {v77[:10]} ({len(v77)})")
    for r in (1, 3, 5, 7):
        vals = [m for q, (m, w) in out.items() if q % 8 == r and m]
        if vals:
            print(f"  q = {r} mod 8: min v2(M_q) = {min((m & -m).bit_length() - 1 for m in vals)}")
    z = [q for q, (m, w) in out.items() if m == 0]
    print(f"M_q = 0 at {len(z)} primes: {z[:25]}")
    r = np.array([out[q][0] / q ** 2 for q in out])
    print(f"M_q/q^2: mean {r.mean():.4f}, rms {np.sqrt((r ** 2).mean()):.4f}, min {r.min():.3f}, max {r.max():.3f}")


if __name__ == '__main__':
    main()
