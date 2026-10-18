#!/usr/bin/env python3
"""The zeta(2)-class Clausen-twin found by clausen_twin_search_20261018.py, and its reduction (S43).

Claim (FINITE-EXACT test): for n = (p-1)/2,  sum_{k=0}^n C(n+k,k)^3 == b(n) (mod p^2),  b(n) = sum_k C(n,k)^2 C(n+k,k)
(Apery zeta(2) numbers, M24 class 4C at the half-period). Reduction: C(n+k,k) == a_k (1 + p O_k), C(n,k) == (-1)^k a_k
(1 - p O_k) (mod p^2), O_k = sum_{j<=k} 1/(2j-1), so the claim is equivalent to  sum_{k<=n} a_k^3 O_k == 0 (mod p),
and, with a_k == (-1)^k C(n,k), O_k == -(H_n - H_(n-k))/2 (mod p), to  H_n D_n == S_n (mod p), D_n = sum (-1)^k C(n,k)^3
(Dixon), S_n = sum (-1)^k C(n,k)^3 H_(n-k).
Session opus-2026-10-17-S43.
"""
from math import comb
from fractions import Fraction as Fr
primes = [p for p in range(5, 1000) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
def red(x, m): return x.numerator * pow(x.denominator, -1, m) % m
t1 = t2 = t3 = True
for p in primes:
    n = (p - 1) // 2
    u = sum(comb(n + k, k) ** 3 for k in range(n + 1))
    b = sum(comb(n, k) ** 2 * comb(n + k, k) for k in range(n + 1))
    if (u - b) % (p * p): t1 = False
    if p < 400:
        O = Fr(0); s = Fr(0)
        for k in range(n + 1):
            if k > 0: O += Fr(1, 2 * k - 1)
            s += Fr(comb(2 * k, k), 4 ** k) ** 3 * O
        if red(s, p) % p: t2 = False
        H = lambda m: sum(Fr(1, j) for j in range(1, m + 1))
        Dn = sum((-1) ** k * comb(n, k) ** 3 for k in range(n + 1))
        Sn = sum((-1) ** k * comb(n, k) ** 3 * H(n - k) for k in range(n + 1))
        if (red(H(n) * Dn, p) - red(Sn, p)) % p: t3 = False
print('sum C(n+k,k)^3 == b(n) mod p^2, all primes 5..997:', t1)
print('sum a_k^3 O_k == 0 mod p, primes < 400:', t2)
print('H_n D_n == S_n mod p, primes < 400:', t3)
