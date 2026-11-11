# series_E.py -- auditor proof-5134-AE.  Numerical checks of THM-5134 (E).
import math
import numpy as np
import numba as nb
from fractions import Fraction as Fr

EULER_GAMMA = 0.5772156649015329


def li2(z, terms=4000):
    # dilogarithm for |z| <= 1/2 by the power series (fast convergence)
    s = 0.0
    zk = 1.0
    for k in range(1, terms):
        zk *= z
        s += zk / (k * k)
        if abs(zk) < 1e-18:
            break
    return s


def sigma_over_n(n):
    # exact sigma(n)/n
    s = Fr(0)
    for d in range(1, n + 1):
        if n % d == 0:
            s += Fr(1, d)
    return s


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


@nb.njit(cache=True)
def abundancy_odd(K):
    # ab[i] = sigma(n)/n for n = 2 i + 1 <= K   (sum of 1/d over odd divisors d)
    half = (K + 1) // 2
    ab = np.zeros(half, np.float64)
    for j in range(half):
        d = 2 * j + 1
        inv = 1.0 / d
        n = d
        while n <= K:
            ab[(n - 1) // 2] += inv
            n += 2 * d
    return ab


@nb.njit(cache=True)
def means(K, ab):
    # mean of sigma(n)/n over odd N <= K with n = N (sheet 1) and n = odd part of N + 1 (sheet 0)
    s1 = 0.0
    s0 = 0.0
    cnt = 0
    mn0 = 1e9
    mn1 = 1e9
    c0 = 0
    c1 = 0
    for N in range(1, K, 2):
        cnt += 1
        a1 = ab[(N - 1) // 2]
        t = N + 1
        while t % 2 == 0:
            t //= 2
        a0 = ab[(t - 1) // 2]
        s1 += a1
        s0 += a0
        if a0 == 1.0:
            c0 += 1
        if a1 == 1.0:
            c1 += 1
    return s0 / cnt, s1 / cnt, c0, c1, cnt


def main():
    c = 8 / math.pi ** 2
    print("8/pi^2 = %.6f   prod_{p odd}(1 - p^-2) = (6/pi^2)/(1 - 1/4) = %.6f" % (c, (6 / math.pi ** 2) / 0.75))
    # partial Euler product as a check
    pr = 1.0
    sieve = np.ones(2000001, bool)
    sieve[:2] = False
    for i in range(2, 1415):
        if sieve[i]:
            sieve[i * i::i] = False
    for p in np.nonzero(sieve)[0]:
        if p > 2:
            pr *= 1 - 1.0 / (p * p)
    print("   partial product over odd p < 2e6: %.7f" % pr)
    k0 = math.pi ** 2 / 12 + 2 * li2(-0.5) - li2(-0.25)
    k1 = li2(-0.25) - li2(-0.5)
    w0 = k0 / (k0 + k1)
    w1 = k1 / (k0 + k1)
    print("kappa_0 = %.10f  kappa_1 = %.10f  (THM-5128 (D1): 0.1615389173, 0.2125139092)   w_0 = %.8f  w_1 = %.8f" % (k0, k1, w0, w1))
    for N in (7, 21, 1, 3, 15, 31, 63, 105, 119, 945):
        n0 = odd_part(N + 1)
        n1 = N
        S0 = c * float(sigma_over_n(n0))
        S1 = c * float(sigma_over_n(n1))
        print("N = %4d:  sheet 0: n = %4d  S_0 = %.5f    sheet 1: n = %4d  S_1 = %.5f    w_0 S_0 + w_1 S_1 = %.5f" % (N, n0, S0, n1, S1, w0 * S0 + w1 * S1))
    for L in (16, 20, 24):
        K = 1 << L
        ab = abundancy_odd(K + 2)
        m0, m1, c0, c1, cnt = means(K, ab)
        print("odd N < 2^%d: mean of S_0 = %.6f  mean of S_1 = %.6f  (pi^2/8 * 8/pi^2 = 1);  #{S_0 minimal} = %d (N + 1 a power of 2: %d values)  #{S_1 minimal} = %d"
              % (L, c * m0, c * m1, c0, L, c1))
    # maximal order for odd n: odd primorials
    print("maximal order: sigma(n)/(n ln ln n) at the odd primorials n = 3*5*...*p  against e^gamma/2 = %.5f ; S/(ln ln n) against 4 e^gamma/pi^2 = %.5f"
          % (math.exp(EULER_GAMMA) / 2, 4 * math.exp(EULER_GAMMA) / math.pi ** 2))
    logn = 0.0
    ab_ = 1.0
    k = 0
    for p in np.nonzero(sieve)[0]:
        if p == 2:
            continue
        logn += math.log(p)
        ab_ *= 1 + 1.0 / p
        k += 1
        if k in (5, 10, 100, 1000, 10000, 100000, 148932):
            print("   %6d odd primes (p <= %7d): ln n = %.1f   sigma(n)/n = %.4f   ratio to ln ln n = %.5f   [squarefree n; prime powers add the factor zeta-like prod 1/(1-p^-2) -> pi^2/8: %.5f]"
                  % (k, p, logn, ab_, ab_ / math.log(logn), ab_ * (math.pi ** 2 / 8) / math.log(logn)))
    # the loose reading of "a sequence of m whose prime factors exhaust the odd primes"
    print("limit along m_j = p_j^2 (single primes; union of prime factors = all odd primes) for n = 1: L_(p^2)(unit) = 1 - p^-2 -> 1, not 8/pi^2 = %.4f" % c)
    print("limit along m_j = (3*5*...*p_j)^2 for n = 1: prod (1 - p^-2) -> %.4f" % c)


if __name__ == "__main__":
    main()
