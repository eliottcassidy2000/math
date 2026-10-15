#!/usr/bin/env python3
"""Exact values of the Paley six-vertex terminal V_q = I_{K6-P6}(q) (paley_next_cores_20261014i.md, eq. (4)).

V_q = sum over injective x:{0..5} -> F_q of prod over the ten edges (a,b), a<b, of K6 minus the path 0-1-2-3-4-5
of chi(x_b - x_a), chi the quadratic character (chi(0) = 0), q = 3 mod 4 prime.

Affine normalisation (the summand is invariant under x -> lam x + c since the edge count is even): V_q = q(q-1) W_q,
with x0 = 0, x1 = 1 and

  W_q = sum_{x2 not in {0,1}, x3 not in {0,1}, x3 != x2} chi(x2) chi(x3) chi(x3-1) S(x2,x3),
  S   = F5(x2,x3) (F4(x2) - f4(x3; x2)) + chi(-x3) + chi(1-x3) + chi(x2-x3),
  f4(x; a)  = chi(x) chi(x-1) chi(x-a),            F4(a)   = sum_x f4(x; a),
  F5(a, b)  = sum_x chi(x) chi(x-1) chi(x-a) chi(x-b).

(Inclusion-exclusion over the four path non-edges 12, 23, 34, 45; the others are forced distinct by chi(0) = 0.)
The O(q^3) evaluation uses one float64 matrix product (exact: all entries are integers of size <= q).
A separate O(q^4) brute force over the normalised injections checks every q <= QBRUTE.

Usage: python3 k6p6_values.py [QMAX=1500] [QBRUTE=43] > k6p6_values.out
(also writes k6p6_values.json and k6p6_W_allprimes.json)
"""
import sys
import json
import time
import numpy as np
import sympy

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
QBRUTE = int(sys.argv[2]) if len(sys.argv) > 2 else 43
RECORDED = {7: 336, 11: -12320, 19: -76608, 31: 52080}   # paley_next_cores_20261014i.md section 5


def chi_table(q):
    c = np.zeros(q, dtype=np.int64)
    sq = (np.arange(1, q, dtype=np.int64) ** 2) % q
    c[np.unique(sq)] = 1
    c[(c == 0)] = -1
    c[0] = 0
    return c


def W_fast(q):
    c = chi_table(q)
    x = np.arange(q)
    base = c[x][None, :] * c[(x - 1) % q][None, :] * c[(x[None, :] - x[:, None]) % q]   # base[a, x] = f4(x; a)
    F4 = base.sum(axis=1)
    C = c[(x[None, :] - x[:, None]) % q]                    # C[b, x] = chi(x - b)
    F5 = np.rint(base.astype(np.float64) @ C.T.astype(np.float64)).astype(np.int64)   # F5[a, b]
    a = x[:, None]; b = x[None, :]
    mask = (a != 0) & (a != 1) & (b != 0) & (b != 1) & (a != b)
    w = c[a] * c[b] * c[(b - 1) % q]
    S = F5 * (F4[:, None] - base) + c[(-b) % q] + c[(1 - b) % q] + c[(a - b) % q]
    return int((mask * w * S).sum())


def W_brute(q):
    c = chi_table(q)
    total = 0
    r = range(q)
    for x2 in r:
        if x2 in (0, 1):
            continue
        for x3 in r:
            if x3 == x2:
                continue
            w3 = c[x2] * c[x3] * c[(x3 - 1) % q]
            if w3 == 0:
                continue
            for x4 in r:
                if x4 == x3:
                    continue
                w4 = c[x4] * c[(x4 - 1) % q] * c[(x4 - x2) % q]
                if w4 == 0:
                    continue
                for x5 in r:
                    if x5 == x4:
                        continue
                    total += w3 * w4 * c[x5] * c[(x5 - 1) % q] * c[(x5 - x2) % q] * c[(x5 - x3) % q]
    return int(total)


def main():
    t0 = time.time()
    out = {}
    for q in sympy.primerange(3, QMAX + 1):
        if q % 4 != 3:
            continue
        W = W_fast(q)
        V = q * (q - 1) * W
        if q <= QBRUTE:
            assert W == W_brute(q), q
        if q in RECORDED:
            assert V == RECORDED[q], (q, V, RECORDED[q])
        out[q] = V
    print(f"V_q for {len(out)} primes q = 3 mod 4 up to {QMAX}; brute force agrees for q <= {QBRUTE}; "
          f"recorded values {RECORDED} reproduced  [{time.time() - t0:.1f}s]")
    for q in list(out)[:12]:
        print(f"  q = {q:5d}: V_q = {out[q]},  W_q = V_q/(q(q-1)) = {out[q] // (q * (q - 1))}")
    json.dump({str(k): v for k, v in out.items()}, open('k6p6_values.json', 'w'))
    # W_q for every prime 5 <= q <= QMAX (both classes mod 4; the normalisation needs only the even edge count),
    # read by k6p6_fft.py as its exact reference
    allW = {str(q): W_fast(q) for q in sympy.primerange(5, QMAX + 1)}
    json.dump(allW, open('k6p6_W_allprimes.json', 'w'))
    print(f"W_q for all {len(allW)} primes 5 <= q <= {QMAX} written to k6p6_W_allprimes.json")


if __name__ == '__main__':
    main()
