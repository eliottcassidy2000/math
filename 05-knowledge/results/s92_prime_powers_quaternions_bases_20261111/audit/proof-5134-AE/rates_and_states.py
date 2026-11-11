# rates_and_states.py -- auditor proof-5134-AE.
# (1) Spectral rates max_{k != 0} rho(T_k)/2 of the twisted digit matrices of THM-5128 (B) at the composite
#     moduli used by THM-5134 (own matrices; numpy eigenvalues).  (B1) says the rate is < 1 for every odd m, which is
#     what (B4)(iv), hence THM-5134 (A), needs.  theta_m = ln(1/rho)/ln(2/rho^2).
# (2) Index-order equidistribution of the states modulo m (THM-5128 (B2)): on the blocks B_n the largest deviation
#     of a residue vector (x, y) mod m from 2^(n-1)/m^2, and the shares of the types gcd(g(h), m) = d against c_d.
import math, sys
import numpy as np
import numba as nb
from fractions import Fraction as Fr


def primes_of(n):
    ps = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        ps.append(n)
    return ps


def c_d(m, d):
    v = Fr(1, d * d)
    for p in primes_of(m // d):
        v *= (1 - Fr(1, p * p))
    return v


def rate(m):
    inv2 = (m + 1) // 2
    best = 0.0
    bestk = 0
    for k in range(1, m):
        T = np.zeros((2 * m, 2 * m), complex)
        for eta in range(m):
            for beta in range(2):
                i = 2 * eta + beta
                T[i, 2 * ((2 * eta) % m) + 0] += 1.0
                ph = (k * ((eta + beta) * inv2)) % m
                T[i, 2 * ((2 * eta + 1) % m) + 1] += np.exp(2j * np.pi * ph / m)
        r = max(abs(np.linalg.eigvals(T))) / 2
        if r > best:
            best, bestk = r, k
    return best, bestk


@nb.njit(cache=True)
def build_psi(H):
    psi = np.zeros(H + 1, np.int64)
    for h in range(1, H + 1):
        k = h >> 1
        if (h & 1) == 0:
            psi[h] = psi[k] + k
        else:
            psi[h] = psi[k] + (k >> 1) + 1
    return psi


@nb.njit(cache=True)
def block_hist(psi, n, m):
    lo = 1 << (n - 1)
    hi = 1 << n
    cnt = np.zeros(m * m, np.int64)
    for h in range(lo, hi):
        x = (2 * psi[h] + 1) % m
        y = (2 * h + 1) % m
        cnt[x * m + y] += 1
    return cnt


def main():
    mods = [9, 15, 21, 25, 27, 45, 49, 63, 75, 81, 125, 225]
    print("(1) rates max_k rho(T_k)/2 (own matrices), and theta_m/2 at rho_m = rate (supremum of the admissible exponents)")
    rates = {}
    for m in mods:
        r, k = rate(m)
        rates[m] = r
        th = math.log(1 / r) / math.log(2 / r ** 2)
        print("    m = %3d: rate %.6f at k = %3d (gcd(k, m) = %d)   theta_m < %.4f   theta_m/2 < %.4f   X^(-theta/2) at X = 2^28: %.3f" % (m, r, k, math.gcd(k, m), th, th / 2, 2 ** (-28 * th / 2)))
    NMAX = 26
    psi = build_psi(1 << NMAX)
    print()
    print("(2) blocks B_n: max over residue vectors of |#{(x, y) = w mod m} - 2^(n-1)/m^2| / (2^(n-1)/m^2), and the type shares")
    for m in (9, 27, 25, 45, 49):
        divs = [d for d in range(1, m + 1) if m % d == 0]
        print("  m = %d   c_d = %s   rate %.4f" % (m, ["%.6f" % float(c_d(m, d)) for d in divs], rates[m]))
        prev = None
        for n in (10, 12, 14, 16, 18, 20, 22, 24, 26):
            cnt = block_hist(psi, n, m).reshape(m, m)
            tot = 1 << (n - 1)
            dev = np.abs(cnt - tot / m ** 2).max() / (tot / m ** 2)
            sh = {d: 0 for d in divs}
            for a in range(m):
                for b in range(m):
                    sh[math.gcd(math.gcd(a, b), m)] += int(cnt[a, b])
            line = "  ".join("d=%d: %.6f" % (d, sh[d] / tot) for d in divs)
            emp = "" if prev is None else "  (decay per digit since n-2: %.4f)" % ((dev / prev) ** 0.5)
            prev = dev
            print("    n = %2d: rel. max dev %.5f%s   type shares %s" % (n, dev, emp, line))


if __name__ == "__main__":
    main()
