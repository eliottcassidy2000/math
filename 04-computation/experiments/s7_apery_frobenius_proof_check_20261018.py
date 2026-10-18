#!/usr/bin/env python3
"""Check that the proof of THM-4669 extends to every Frobenius point n = (p^r - 1)/2 (S43, HYP-9291 -> THM).

p-adic checks (rationals reduced in Z_(p) modulo p^2; a value with p in the denominator is reported):
 (C1') c_k(n) = C(n,k)C(n+k,k) == (-1)^k a_k^2 (mod p^2), 0 <= k <= n   [ratio prod_j (1 - p^(2r)/(2j-1)^2)];
 (C2') a_m == (-1)^m C(n,m) (1 - eps S_m) (mod p^2), eps = p^r/2, S_m = sum_{i<m} 1/(n-i)  [eps/(n-i) in p Z_p];
 (L')  sum_k a_k^2 a_(n-k)^2 == sum_k a_k^4 (mod p^2);
 (T')  s_7(n) == (-1)^n A(n) (mod p^2).
Session opus-2026-10-17-S43.
"""
from math import comb
from fractions import Fraction as Fr

def vp(x, p):
    if x == 0:
        return 10 ** 9
    v = 0
    while x % p == 0:
        x //= p; v += 1
    return v

def red(x, p):
    """reduce a rational with p-free denominator mod p^2; raise if p divides the denominator."""
    m = p * p
    if x.denominator % p == 0:
        raise ValueError('p in denominator')
    return x.numerator * pow(x.denominator, -1, m) % m

def a(k): return Fr(comb(2 * k, k), 4 ** k)

for p, r in ((5, 2), (7, 2), (11, 2), (13, 2), (5, 3), (7, 3)):
    n = (p ** r - 1) // 2
    eps = Fr(p ** r, 2)
    c1 = c2 = True
    for k in range(n + 1):
        ck = comb(n, k) * comb(n + k, k)
        if (ck - red((-1) ** k * a(k) ** 2, p)) % (p * p):
            c1 = False
    S = Fr(0)
    for m in range(n + 1):
        if m > 0:
            S += Fr(1, n - (m - 1))
        val = (-1) ** m * comb(n, m) * (1 - eps * S)
        try:
            if (red(a(m), p) - red(val, p)) % (p * p):
                c2 = False
        except ValueError:
            c2 = False
    L = (red(sum(a(k) ** 2 * a(n - k) ** 2 for k in range(n + 1)), p) - red(sum(a(k) ** 4 for k in range(n + 1)), p)) % (p * p) == 0
    s7 = sum(comb(n, k) ** 2 * comb(n + k, k) * comb(2 * k, n) for k in range((n + 1) // 2, n + 1))
    A = sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1))
    T = (s7 - (-1) ** n * A) % (p * p) == 0
    print('p = %d, r = %d, n = %d: C1 %s, C2 %s, L %s, T %s' % (p, r, n, c1, c2, L, T), flush=True)
