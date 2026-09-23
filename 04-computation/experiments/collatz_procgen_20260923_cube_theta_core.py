#!/usr/bin/env python3
"""
collatz_procgen_20260923_cube_theta_core.py

Shared exact helpers for the cube-theta lane (HYP-9127).

Conventions
-----------
* rho = 2^L / M0 with M0 odd; for the cube-swap word Y3 under 3x+1, L = 10, M0 = 3^9.
* A 2-adic integer is stored as a Python/gmpy2 integer in [0, 2^N) ("mod 2^N").
* theta_P(rho) = sum_{k >= kmin} rho^{P(k)} for an increasing, nonnegative, integer-valued P.
  Every term rho^e with e >= 0 is a 2-adic integer of valuation L*e, so only the finitely
  many k with L*P(k) < N matter modulo 2^N.
* X := theta_{k^3}(rho) = sum_{k>=0} rho^(k^3) is the number of HYP-9127.

Everything here is exact integer arithmetic (gmpy2); no floating point enters any certificate.
"""
import math
import gmpy2
from gmpy2 import mpz

LOG2_3 = math.log2(3.0)
PHI = (1 + 5 ** 0.5) / 2

L_Y3 = 10
M0_Y3 = 3 ** 9          # 19683
A_Y3 = 9                # number of 1s per block (3x+1: M0 = 3^A)


def mubar(L, M0):
    """mu_bar = max(1, log2(M0)/L): the height exponent of rho = 2^L/M0."""
    return max(1.0, math.log2(M0) / L)


def v2(x):
    """2-adic valuation of a nonzero integer (inf for 0)."""
    x = mpz(x)
    if x == 0:
        return math.inf
    return int(gmpy2.bit_scan1(abs(x)))


def inv_mod2N(a, N):
    return gmpy2.invert(mpz(a), mpz(1) << N)


def rho_pow_mod(e, N, L=L_Y3, M0=M0_Y3, invM0=None):
    """rho^e mod 2^N for an integer e >= 0 (0 if L*e >= N)."""
    if L * e >= N:
        return mpz(0)
    mod = mpz(1) << N
    if invM0 is None:
        invM0 = inv_mod2N(M0, N)
    return ((mpz(1) << (L * e)) * gmpy2.powmod(invM0, e, mod)) % mod


def theta_mod(N, P, L=L_Y3, M0=M0_Y3, kmin=0, weight=None):
    """sum_{k>=kmin} w(k) rho^{P(k)} mod 2^N, P increasing with P(k) >= 0 for k >= kmin."""
    mod = mpz(1) << N
    invM0 = inv_mod2N(M0, N)
    s = mpz(0)
    k = kmin
    while True:
        e = P(k)
        if e < 0:
            raise ValueError("P(k) < 0")
        if L * e >= N:
            # P increasing: all later terms vanish mod 2^N
            break
        t = rho_pow_mod(e, N, L, M0, invM0)
        if weight is not None:
            t = t * weight(k)
        s += t
        k += 1
    return s % mod


def X_mod(N, L=L_Y3, M0=M0_Y3):
    """X = sum_{k>=0} rho^(k^3) mod 2^N."""
    return theta_mod(N, lambda k: k ** 3, L, M0)


def aux_mod(N, j, L=L_Y3, M0=M0_Y3):
    """Y_j = sum_{k>=0} k^j rho^(k^3) mod 2^N  (values of u-derivatives of G at (1,1))."""
    return theta_mod(N, lambda k: k ** 3, L, M0, weight=lambda k: mpz(k) ** j)


def theta2_mod(N, L=L_Y3, M0=M0_Y3):
    """sum_{k>=0} rho^(k^2) mod 2^N (square-swap number of Theorem Y)."""
    return theta_mod(N, lambda k: k ** 2, L, M0)


def partial_sum_fraction(K, P=lambda k: k ** 3, L=L_Y3, M0=M0_Y3, kmin=0):
    """S_K = sum_{k=kmin..K} rho^{P(k)} = A/M0^{P(K)} as exact integers (A, D)."""
    EK = P(K)
    A = mpz(0)
    for k in range(kmin, K + 1):
        e = P(k)
        A += (mpz(2) ** (L * e)) * (mpz(M0) ** (EK - e))
    return A, mpz(M0) ** EK


def log2abs(x):
    x = abs(mpz(x))
    if x == 0:
        return -math.inf
    # exact enough for reporting: bit_length - 1 + log2(mantissa)
    b = int(x.bit_length())
    if b <= 1000:
        return math.log2(int(x))
    sh = b - 60
    return math.log2(int(x >> sh)) + sh


def cube_word(nblocks, swap=lambda n: False, normal="1111111110", swapped="1111111101"):
    """Concatenate 10-letter blocks: 'swapped' at indices n with swap(n) True, else 'normal'."""
    return "".join(swapped if swap(n) else normal for n in range(nblocks))


def is_cube(n):
    if n <= 0:
        return False
    r = int(round(n ** (1.0 / 3)))
    for c in (r - 1, r, r + 1):
        if c > 0 and c * c * c == n:
            return True
    return False


def is_square(n):
    if n <= 0:
        return False
    r = math.isqrt(n)
    return r * r == n


def bernstein_3x1_mod(word, N):
    """Phi(w) mod 2^N for T(x) = x/2, (3x+1)/2: Phi(w) = -sum_l 2^{d_l} / 3^{l+1}
    (d_l = position of the l-th 1). Requires len(word) >= N (terms with d_l >= N vanish)."""
    mod = mpz(1) << N
    inv3 = inv_mod2N(3, N)
    p = inv3            # 3^{-(l+1)} for l = 0
    s = mpz(0)
    for d, ch in enumerate(word[:N]):
        if ch == "1":
            s += (mpz(1) << d) * p
            p = (p * inv3) % mod
    return (-s) % mod


def parity_check_3x1(x, word, N):
    """Iterate T on x mod 2^N; each step loses one bit of precision. Return the number of
    initial letters of 'word' that agree with the parity vector of x (at most N)."""
    x = mpz(x)
    prec = N
    agree = 0
    for ch in word:
        if prec <= 0:
            break
        b = int(x & 1)
        if str(b) != ch:
            return agree
        agree += 1
        if b:
            x = 3 * x + 1
        x >>= 1
        prec -= 1
        x &= (mpz(1) << prec) - 1
    return agree
