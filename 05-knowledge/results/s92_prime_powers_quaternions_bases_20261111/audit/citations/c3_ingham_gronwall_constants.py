#!/usr/bin/env python3
"""Auditor 'citations' - Ingham's asymptotic, Gronwall's constant for odd n, and the quoted constants of THM-5134 (E).

I1  D(n) = sum_{m=1}^{n-1} d(m) d(n-m) = #{(a1, a2, b1, b2) positive : a1 a2 + b1 b2 = n}  (checked by brute force for small n)
    against Ingham's main term (6/pi^2) sigma(n) ln^2 n; the ratio at n near 10^3 .. 10^6.
I2  The note's illustration: correlation of D(n)/n with sigma(n)/n on [20000, 20400).
G1  Gronwall for odd n: limsup sigma(n)/(n ln ln n) = e^gamma / 2, hence S = (8/pi^2) sigma(n)/n has maximal order
    (4 e^gamma / pi^2) ln ln n.  Evidence: odd 'primorial-type' n (products of the first odd primes), the Mertens
    product over odd primes, and the record ratios over all odd n <= 3 * 10^6.
G2  Mean of sigma(n)/n over odd n -> pi^2/8, and on the sheet 0 weighting (odd part of N + 1, N odd).
K1  The constants: 8/pi^2, w_0 S_0 + w_1 S_1 at N = 7 and N = 21 with the shares of THM-5128 (D1).
K2  Remark 1 of THM-5134: (1 - p^-2) sigma(p^v)/p^v against the squarefree-level factor 1 - p^-2 + 1/p.
"""
import math
import sys
from fractions import Fraction
import numpy as np

EULER_GAMMA = 0.57721566490153286060651209


def sieve_divisor_functions(L):
    d = np.zeros(L + 1, dtype=np.int64)
    s = np.zeros(L + 1, dtype=np.int64)
    for k in range(1, L + 1):
        d[k::k] += 1
        s[k::k] += k
    return d, s


def main():
    L = 1_000_100
    d, s = sieve_divisor_functions(L)

    # ---------- I1 ----------
    print("I1 Ingham: D(n) = sum_{m<n} d(m) d(n-m) against (6/pi^2) sigma(n) ln^2 n")
    # brute force identification for small n
    bad = 0
    for n in range(2, 61):
        brute = 0
        for a1 in range(1, n):
            for a2 in range(1, (n - 1) // a1 + 1):
                rest = n - a1 * a2
                if rest >= 1:
                    brute += int(d[rest])
        Dn = int(np.dot(d[1:n], d[n - 1:0:-1]))
        if brute != Dn:
            bad += 1
    print("   D(n) = number of positive solutions of a1 a2 + b1 b2 = n, 2 <= n <= 60: failures %d" % bad)

    def D(n):
        return int(np.dot(d[1:n], d[n - 1:0:-1]))

    print("   n, D(n), D(n) / ((6/pi^2) sigma(n) ln^2 n):")
    for n in (1000, 1001, 10007, 10080, 100003, 100800, 720720, 999983, 1000000):
        Dn = D(n)
        main_term = 6 / math.pi ** 2 * float(s[n]) * math.log(n) ** 2
        print("     %8d  %14d  %.4f" % (n, Dn, Dn / main_term))
    # ratio averaged over windows: converges like 1 + c/ln n
    for lo in (1000, 10000, 100000, 1000000 - 200):
        rs = []
        for n in range(lo, lo + 200):
            rs.append(D(n) / (6 / math.pi ** 2 * float(s[n]) * math.log(n) ** 2))
        print("   window [%d, %d): mean ratio %.4f, min %.4f, max %.4f" % (lo, lo + 200, np.mean(rs), min(rs), max(rs)))

    # ---------- I2 ----------
    xs = []
    ys = []
    for n in range(20000, 20400):
        xs.append(D(n) / n)
        ys.append(float(s[n]) / n)
    corr = float(np.corrcoef(xs, ys)[0, 1])
    print("I2 correlation of D(n)/n with sigma(n)/n on [20000, 20400): %.4f" % corr)
    # same with the odd part only / the series without the prime 2, for the record
    ys2 = []
    for n in range(20000, 20400):
        m = n
        while m % 2 == 0:
            m //= 2
        ys2.append(float(s[m]) / m)
    print("   (with sigma(n_odd)/n_odd instead: %.4f)" % float(np.corrcoef(xs, ys2)[0, 1]))

    # ---------- G1 ----------
    print("G1 Gronwall for odd n: limsup sigma(n)/(n ln ln n) = e^gamma/2 = %.6f ; all n: e^gamma = %.6f"
          % (math.exp(EULER_GAMMA) / 2, math.exp(EULER_GAMMA)))
    print("   (4 e^gamma / pi^2) = %.6f = (8/pi^2) * (e^gamma/2)" % (4 * math.exp(EULER_GAMMA) / math.pi ** 2))
    # Mertens over odd primes: prod_{2 < p <= x} (1 - 1/p)^-1 ~ (e^gamma/2) ln x
    P = 2_000_000
    sv = np.ones(P + 1, dtype=bool)
    sv[:2] = False
    for i in range(2, int(P ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    primes = np.nonzero(sv)[0]
    odd_primes = primes[1:]
    logprod = np.cumsum(-np.log1p(-1.0 / odd_primes))
    for x in (100, 1000, 10000, 100000, 1000000, 2000000):
        idx = np.searchsorted(odd_primes, x, side="right") - 1
        print("   x = %7d: prod_{2 < p <= x} (1 - 1/p)^-1 / ln x = %.6f" % (x, math.exp(logprod[idx]) / math.log(x)))
    # odd primorials n_k = 3 * 5 * ... * p_k : sigma(n)/n = prod (1 + 1/p); and with full prime powers p^a (a large):
    # sup over exponents is prod p/(p-1).  ln ln n_k ~ ln p_k.
    print("   odd primorial-type numbers n = prod_{2 < p <= y} p^a (a -> infinity gives prod p/(p-1)):")
    theta = np.cumsum(np.log(odd_primes.astype(float)))
    for y in (101, 1009, 10007, 100003, 1000003):
        idx = np.searchsorted(odd_primes, y, side="right") - 1
        # choose exponents a_p = floor(ln y / ln p) (so p^a <= y): n = lcm-like, ln n ~ psi(y) - ln-part of 2
        ps = odd_primes[: idx + 1].astype(float)
        a = np.floor(math.log(y) / np.log(ps))
        ln_n = float(np.sum(a * np.log(ps)))
        sig_over_n = float(np.exp(np.sum(np.log((1 - ps ** (-(a + 1))) / (1 - 1 / ps)))))
        print("     y = %8d: ln n = %.1f, sigma(n)/n = %.5f, ratio to ln ln n = %.5f, S/(ln ln n) = %.5f"
              % (y, ln_n, sig_over_n, sig_over_n / math.log(ln_n), 8 / math.pi ** 2 * sig_over_n / math.log(ln_n)))
    # records over all odd n <= 3e6 of sigma(n)/(n ln ln n)  (n >= 17 so that ln ln n > 1)
    L2 = 3_000_000
    s2 = np.zeros(L2 + 1, dtype=np.int64)
    for k in range(1, L2 + 1, 2):
        s2[k::k] += k
    n_odd = np.arange(17, L2 + 1, 2)
    rat = s2[n_odd] / n_odd / np.log(np.log(n_odd))
    for lo, hi in ((10 ** 3, 10 ** 4), (10 ** 4, 10 ** 5), (10 ** 5, 10 ** 6), (10 ** 6, 3 * 10 ** 6)):
        mask = (n_odd >= lo) & (n_odd < hi)
        i = int(np.argmax(rat[mask]))
        nn = int(n_odd[mask][i])
        print("   odd n in [%d, %d): max sigma(n)/(n ln ln n) = %.5f at n = %d" % (lo, hi, float(rat[mask][i]), nn))
    print("   (the maxima decrease towards e^gamma/2 = 0.89054 from above, as for Gronwall's theorem)")

    # ---------- G2 ----------
    for X in (10 ** 4, 10 ** 5, 10 ** 6, 3 * 10 ** 6):
        odd = np.arange(1, X + 1, 2)
        mean1 = float(np.mean(s2[odd] / odd))
        # sheet 0: N odd <= X, n = odd part of N + 1
        Np1 = np.arange(2, X + 2, 2)
        m = Np1.copy()
        while True:
            ev = (m % 2 == 0)
            if not ev.any():
                break
            m[ev] //= 2
        mask = m <= L2
        mean0 = float(np.mean(s2[m[mask]] / m[mask]))
        print("G2 X = %8d: mean of (8/pi^2) sigma(n)/n over odd N <= X: sheet 1 %.6f ; sheet 0 (n = odd part of N+1) %.6f"
              % (X, 8 / math.pi ** 2 * mean1, 8 / math.pi ** 2 * mean0))

    # ---------- K1 ----------
    w0, w1 = 0.43186124, 0.56813876
    c = 8 / math.pi ** 2
    print("K1 8/pi^2 = %.6f" % c)

    def sig(n):
        return sum(k for k in range(1, n + 1) if n % k == 0)

    def odd_part(n):
        while n % 2 == 0:
            n //= 2
        return n

    for N in (7, 21, 1, 3, 15, 63, 105, 127):
        n1 = odd_part(N)
        n0 = odd_part(N + 1)
        S1 = c * sig(n1) / n1
        S0 = c * sig(n0) / n0
        print("   N = %3d: S_1 = %.4f (n = %d), S_0 = %.4f (n = %d), w_0 S_0 + w_1 S_1 = %.4f" % (N, S1, n1, S0, n0, w0 * S0 + w1 * S1))
    # kappa values of THM-5128 (D1) from dilogarithms, and the shares
    def Li2(z, terms=4000):
        return sum(z ** k / k ** 2 for k in range(1, terms))

    k0 = math.pi ** 2 / 12 + 2 * Li2(-0.5) - Li2(-0.25)
    k1 = Li2(-0.25) - Li2(-0.5)
    print("   kappa_0 = %.10f, kappa_1 = %.10f, w_0 = %.8f, w_1 = %.8f" % (k0, k1, k0 / (k0 + k1), k1 / (k0 + k1)))

    # ---------- K2 ----------
    print("K2 Remark 1: factor at v = 1 minus squarefree-level factor = -p^-3; at v = 2: +p^-2 - p^-3 - p^-4")
    for p in (3, 5, 7, 11):
        sf = 1 - Fraction(1, p * p) + Fraction(1, p)
        f1 = (1 - Fraction(1, p * p)) * Fraction(1 + p, p)
        f2 = (1 - Fraction(1, p * p)) * Fraction(1 + p + p * p, p * p)
        # mean of the prime-power factor over the classes divisible by p (limit k -> infinity)
        mean = sum((1 - Fraction(1, p)) * Fraction(1, p ** (j - 1)) * (1 + Fraction(1, p)) * (1 - Fraction(1, p ** (j + 1)))
                   for j in range(1, 60))
        print("   p = %2d: f(v=1) - sf = %s (-p^-3 = %s); f(v=2) - sf = %s; mean over v >= 1 = %.12f, sf = %.12f"
              % (p, f1 - sf, Fraction(-1, p ** 3), f2 - sf, float(mean), float(sf)))
    # limits at 9
    print("   limits at 9: v = 0: %s, v = 1: %s, v >= 2: %s" % (Fraction(8, 9), Fraction(32, 27), Fraction(35, 27)),
          " = %.3f / %.3f / %.3f" % (8 / 9, 32 / 27, 35 / 27))


if __name__ == "__main__":
    main()
