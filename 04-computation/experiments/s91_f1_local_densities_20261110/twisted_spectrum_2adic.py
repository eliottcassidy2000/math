#!/usr/bin/env python3
"""S91 (THM-5128): the prime 2.  Twisted transfer operators of the digit recursion of q modulo 2^j.

q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2).  Modulo 2^j the increment ceil(h/2) = (h + (h mod 2))/2 is a function of
w = h mod 2^(j+1) (the last j+1 binary digits), so the state is w, with 2^(j+1) values, and for a character k of Z/2^j

        T_k[w -> 2w mod 2^(j+1)] = 1,     T_k[w -> (2w + 1) mod 2^(j+1)] = exp(2 pi i k c(w) / 2^j),   c(w) = (w + (w & 1)) / 2.

The two predecessors of a state are w and w + 2^j; c differs by 2^(j-1) between them, so the phases differ by (-1)^k:
for odd k the argument of THM-5128 gives rho(T_k) < 2, and an even k is a character of a smaller modulus.

Checks (floating point, VERIFIED):
  1. S_n(k) = sum_{h < 2^n} exp(2 pi i k q(h) / 2^j) = e_0^T T_k^n 1 against brute force, n = NB, for every k and also
     jointly with the last J low digits of h (q mod 2^j is asymptotically independent of the low digits);
  2. rates rho(T_k)/2 for j = 1..6;
  3. consequence for F1: the proportion of each odd class modulo 8 and 16 among the F1 values over a box (brute force).
"""
import sys

import numpy as np


def q_table(n):
    size = 1 << n
    q = np.zeros(size, dtype=np.int64)
    for k in range(n):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    return q


def twisted_matrix(j, k):
    size = 1 << (j + 1)
    mod = 1 << j
    T = np.zeros((size, size), dtype=np.complex128)
    for w in range(size):
        T[w, (2 * w) % size] += 1.0
        c = ((w + (w & 1)) >> 1) % mod
        T[w, (2 * w + 1) % size] += np.exp(2j * np.pi * ((k * c) % mod) / mod)
    return T


def main():
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    q = q_table(nb)
    ok = True
    print(f"modulus 2^j: brute-force check at n = {nb}, then rates rho(T_k)/2 for odd k (even k belong to 2^(j-1))")
    for j in range(1, 7):
        mod = 1 << j
        rates = {}
        maxerr = 0.0
        for k in range(1, mod):
            T = twisted_matrix(j, k)
            vec = np.zeros(1 << (j + 1), dtype=np.complex128)
            vec[0] = 1.0
            for _ in range(nb):
                vec = vec @ T
            ref = np.exp(2j * np.pi * ((k * q) % mod) / mod).sum()
            maxerr = max(maxerr, abs(vec.sum() - ref))
            rates[k] = float(np.max(np.abs(np.linalg.eigvals(T)))) / 2.0
            if k % 2 == 1:
                # for odd k the rows w and w + 2^j are orthogonal: T_k / sqrt(2) is unitary (exact square-root law)
                ok &= bool(np.abs(T @ T.conj().T - 2 * np.eye(1 << (j + 1))).max() < 1e-12)
        ok &= maxerr < 1e-6 * (1 << nb)
        ok &= all(r < 1 - 1e-9 for r in rates.values())
        odd = {k: r for k, r in rates.items() if k % 2 == 1}
        worst = max(odd, key=odd.get)
        print(f"  j = {j}: max over odd k of rho/2 = {odd[worst]:.6f} (k = {worst}), min = {min(odd.values()):.6f}, "
              f"max over all k != 0 = {max(rates.values()):.6f}, brute-force error {maxerr:.1e}")
    # joint law of (q mod 8, h mod 8) over h < 2^nb
    h = np.arange(1 << nb, dtype=np.int64)
    # q mod 2 is the quadratic form sum_i h_i (h_(i+1) + h_(i+2)) in the binary digits (Rudin-Shapiro type);
    # so x = 2h + 1 - 2q satisfies (x - 1)/2 = h_0 + sum_i h_i (h_(i+1) + h_(i+2)) mod 2  (THM-5127 (C) on a chain)
    form = np.zeros(1 << nb, dtype=np.int64)
    for i in range(nb):
        form += ((h >> i) & 1) * (((h >> (i + 1)) & 1) + ((h >> (i + 2)) & 1))
    quad_ok = bool(((form - q) % 2 == 0).all())
    s_n = int(np.where(q % 2 == 0, 1, -1).sum())
    print(f"q mod 2 = sum_i h_i (h_(i+1) + h_(i+2)) for all h < 2^{nb}: {quad_ok};  "
          f"sum of (-1)^q over h < 2^{nb} = {s_n}  (bound 2^((n + 2)/2) = {2 ** ((nb + 2) / 2):.0f})")
    ok &= quad_ok and abs(s_n) <= 2 ** ((nb + 2) / 2)
    joint = np.zeros((8, 8))
    np.add.at(joint, (q % 8, h % 8), 1.0)
    joint /= (1 << nb)
    dev = float(np.abs(joint - 1 / 64).max() * 64)
    print(f"joint law of (q mod 8, h mod 8), h < 2^{nb}: max relative deviation from uniform {dev:.4f}")
    ok &= dev < 0.05
    # F1 over a box modulo 8 and 16
    n2 = 11
    y = 2 * h[: 1 << n2] + 1
    x = y - 2 * q[: 1 << n2]
    for s in (0, 1, 2, 3):
        N = np.multiply.outer(x, x) + (1 << s) * np.multiply.outer(y, y) - (1 if s == 0 else 0)
        for mod in (8, 16):
            cnt = np.bincount((N % mod).ravel(), minlength=mod)[1::2].astype(np.float64)
            rel = cnt / cnt.mean() - 1
            print(f"  F1 box [0, 2^{n2})^2, s = {s}, odd classes mod {mod:2d}: max relative deviation {np.abs(rel).max():.4f}")
            ok &= bool(np.abs(rel).max() < 0.08)
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
