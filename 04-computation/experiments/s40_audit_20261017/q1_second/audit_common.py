"""Independent helpers for the S40 audit (written from scratch; no imports from the session's scripts).

- prime sieve
- Legendre symbol by Euler's criterion (numba modpow)
- Kronecker character of a fundamental discriminant via its factorisation into prime discriminants
- class numbers h(D) by Cohen's reduced-form count (loop over b, then divisors a of (b^2 - D)/4)
- L(1, chi_D) by the theta-function (erfc) series, an analytic third method
"""
import math
import numpy as np
from numba import njit


def sieve(n):
    s = np.ones(n + 1, dtype=np.bool_)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


@njit(cache=True)
def modpow(b, e, m):
    r = 1
    b %= m
    while e > 0:
        if e & 1:
            r = (r * b) % m
        b = (b * b) % m
        e >>= 1
    return r


@njit(cache=True)
def legendre(a, p):
    """(a/p) for an odd prime p by Euler's criterion."""
    a %= p
    if a == 0:
        return 0
    r = modpow(a, (p - 1) // 2, p)
    return 1 if r == 1 else -1


@njit(cache=True)
def h_forms(D):
    """Class number of primitive positive definite forms of discriminant D < 0 (Cohen, Alg. 5.3.5 style)."""
    n = -D
    h = 0
    b = n % 2
    B = int(math.sqrt(n / 3.0)) + 1
    while 3 * b * b <= n:
        q = (b * b + n) // 4
        a = b if b > 1 else 1
        while a * a <= q:
            if q % a == 0:
                c = q // a
                # gcd(a, b, c)
                g = a
                x = b
                while x:
                    g, x = x, g % x
                x = c
                while x:
                    g, x = x, g % x
                if g == 1:
                    if a == b or a == c or b == 0:
                        h += 1
                    else:
                        h += 2
            a += 1
        b += 2
    return h


def odd_prime_factors(m):
    f = []
    p = 3
    while p * p <= m:
        if m % p == 0:
            f.append(p)
            while m % p == 0:
                m //= p
        p += 2
    if m > 1:
        f.append(m)
    return f


def fundamental_neg(D):
    """True iff D < 0 is a fundamental discriminant (independent implementation)."""
    n = -D
    if n % 4 == 3:          # D = 1 mod 4
        m = n
        e2 = 0
    elif n % 4 == 0:
        m4 = n // 4          # D = 4 * (-m4), need -m4 = 2,3 mod 4  <=>  m4 = 2, 1 mod 4
        if m4 % 4 not in (1, 2):
            return False
        m = m4
    else:
        return False
    # squarefree check of m
    k = 2
    while k * k <= m:
        if m % (k * k) == 0:
            return False
        k += 1
    return True


def char_table(D):
    """chi_D(a) for 0 <= a < |D| as int8 array, D fundamental < 0, via prime discriminants."""
    n = -D
    odd = n
    while odd % 2 == 0:
        odd //= 2
    ps = odd_prime_factors(odd)
    prod = 1
    for p in ps:
        prod *= (p if p % 4 == 1 else -p)
    d2 = D // prod        # in {1, -4, 8, -8}
    assert d2 * prod == D and d2 in (1, -4, 8, -8), (D, d2)
    a = np.arange(n, dtype=np.int64)
    chi = np.ones(n, dtype=np.int64)
    for p in ps:
        chi *= leg_vec(a, p)
    if d2 != 1:
        odd_mask = (a % 2 == 1)
        c4 = np.where(a % 4 == 1, 1, -1)
        c8 = np.where((a % 8 == 1) | (a % 8 == 7), 1, -1)
        if d2 == -4:
            c = c4
        elif d2 == 8:
            c = c8
        else:
            c = c4 * c8
        chi = np.where(odd_mask, chi * c, 0)
    return chi.astype(np.int8)


@njit(cache=True)
def leg_vec(a, p):
    out = np.empty(a.shape[0], dtype=np.int64)
    for i in range(a.shape[0]):
        out[i] = legendre(a[i], p)
    return out


@njit(cache=True)
def L1_theta(chi, q):
    """L(1, chi) for a real odd primitive character chi mod q (chi given as a table mod q):
    L(1,chi) = sum_n chi(n) [exp(-pi n^2/q)/n + (pi/sqrt q) erfc(n sqrt(pi/q))]."""
    s = 0.0
    n = 1
    lim = int(6.5 * math.sqrt(q)) + 10
    c = math.sqrt(math.pi / q)
    while n <= lim:
        x = chi[n % q]
        if x != 0:
            s += x * (math.exp(-math.pi * n * n / q) / n + (math.pi / math.sqrt(q)) * math.erfc(n * c))
        n += 1
    return s


@njit(cache=True)
def paley_chi(p):
    """(a/p) for 0 <= a < p via squares (used only after cross-checking against Euler's criterion)."""
    chi = np.full(p, -1, dtype=np.int64)
    chi[0] = 0
    for x in range(1, (p - 1) // 2 + 1):
        chi[(x * x) % p] = 1
    return chi
