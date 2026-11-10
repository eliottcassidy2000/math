#!/usr/bin/env python3
"""S91 (THM-5128): twisted transfer operators of the one-apex digit recursion modulo an odd m.

Notation.  psi(0) = 0, psi(2h) = psi(h) + h, psi(2h+1) = psi(h) + floor(h/2) + 1;  q(h) = h - psi(h), so
q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2).  The one-apex count is x(h) = 2 psi(h) + 1 = y(h) - 2 q(h), y(h) = 2h + 1.

Reading the binary digits of h from the top, the pair (q, h) moves by
        digit 0: (q, h) -> (q, 2h)            digit 1: (q, h) -> (q + ceil(h/2), 2h + 1).
Modulo an odd m, ceil(h/2) = (h + beta)/2 with beta = h mod 2, so (h mod m, h mod 2) is a Markov state with 2m values
and q is an additive functional of the path.  For a character k of Z/m the twisted transfer matrix is

        T_k[(eta, beta) -> (2 eta, 0)] = 1,      T_k[(eta, beta) -> (2 eta + 1, 1)] = zeta^(k (eta + beta) / 2),

zeta = exp(2 pi i / m), division by 2 modulo m, and
        S_n(k, l) = sum_{h < 2^n} zeta^(k q(h) + l h) = e_(0,0)^T  T_k^n  v_l,       v_l(eta, beta) = zeta^(l eta).

What this script does (all floating point: status VERIFIED, not exact):
  1. checks the matrix identity for S_n(k, l) against brute force (n <= NB) for every k, l and every m in the list;
  2. prints rho(T_k)/2 (the decay rate per binary digit of the exponential sum) for every k != 0;
  3. prints the proven strict bound check rho(T_k) < 2 and the worst rate over k for each m.
The strict inequality rho(T_k) < 2 for k != 0 mod m is PROVED in THM-5128 (two predecessors of a state differ by the
phase zeta^(k/2) != 1); the numbers printed here are its quantitative form.

Usage: python3 twisted_spectrum.py [NB]     (NB = number of bits for the brute-force check, default 16)
"""
import sys

import numpy as np

MODULI = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 73, 89, 127]


def q_table(n):
    size = 1 << n
    q = np.zeros(size, dtype=np.int64)
    for k in range(n):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    return q


def twisted_matrix(m, k):
    inv2 = pow(2, -1, m)
    zeta = np.exp(2j * np.pi / m)
    T = np.zeros((2 * m, 2 * m), dtype=np.complex128)
    for eta in range(m):
        for beta in range(2):
            i = 2 * eta + beta
            T[i, 2 * ((2 * eta) % m) + 0] += 1.0
            e = (k * ((eta + beta) * inv2 % m)) % m
            T[i, 2 * ((2 * eta + 1) % m) + 1] += zeta ** e
    return T


def main():
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    q = q_table(nb)
    h = np.arange(1 << nb, dtype=np.int64)
    ok = True
    print(f"brute-force check of S_n(k, l) = e^T T_k^n v_l for n = {nb} (all k, l), then rates rho(T_k)/2")
    print(f"{'m':>4} {'max_k rho/2':>12} {'argmax k':>9} {'min_k rho/2':>12}  {'max err':>9}   rates by k (k <= (m-1)/2)")
    for m in MODULI:
        zeta = np.exp(2j * np.pi / m)
        qm = q % m
        hm = h % m
        # brute force: counts of (q mod m, h mod m), then all character sums by a 2D DFT
        cnt = np.zeros((m, m), dtype=np.float64)
        np.add.at(cnt, (qm, hm), 1.0)
        brute = np.fft.fft2(cnt)            # brute[k, l] = sum cnt[a, b] exp(-2 pi i (k a + l b)/m)
        rates = {}
        maxerr = 0.0
        for k in range(m):
            T = twisted_matrix(m, k)
            vec = np.zeros(2 * m, dtype=np.complex128)
            vec[0] = 1.0
            for _ in range(nb):
                vec = vec @ T
            for l in range(m):
                v = np.array([zeta ** ((l * eta) % m) for eta in range(m) for _beta in range(2)])
                val = vec @ v
                ref = brute[(-k) % m, (-l) % m]
                maxerr = max(maxerr, abs(val - ref))
            if k:
                rates[k] = float(np.max(np.abs(np.linalg.eigvals(T)))) / 2.0
        if maxerr > 1e-6 * (1 << nb):
            ok = False
        worst = max(rates, key=rates.get)
        strict = all(r < 1.0 - 1e-9 for r in rates.values())
        ok &= strict
        half = [f"{rates[k]:.4f}" for k in range(1, (m - 1) // 2 + 1)]
        shown = " ".join(half[:12]) + (" ..." if len(half) > 12 else "")
        print(f"{m:4d} {rates[worst]:12.6f} {worst:9d} {min(rates.values()):12.6f}  {maxerr:9.2e}   {shown}")
        # symmetry k -> -k (complex conjugation)
        for k in range(1, m):
            if abs(rates[k] - rates[m - k]) > 1e-9:
                ok = False
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
