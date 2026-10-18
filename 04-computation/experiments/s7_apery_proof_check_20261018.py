#!/usr/bin/env python3
"""Step-by-step check of the proof of s7((p-1)/2) == (-1)^((p-1)/2) A((p-1)/2) (mod p^2) (S43, THM for HYP-9291 r = 1).

a_k = C(2k,k)/4^k. Steps:
 (E1) exact: C(n,k) C(2k,n) = C(2k,k) C(k, n-k);
 (E2) exact (Clausen + quadratic transformation): sum_k (-1)^(n-k) 4^k a_k^3 C(k, n-k) = sum_k a_k^2 a_(n-k)^2, all n;
 (C1) n = (p-1)/2: C(n,k) C(n+k,k) == (-1)^k a_k^2 (mod p^2), 0 <= k <= n;
 (C2) n = (p-1)/2: a_m == (-1)^m C(n,m) (1 - (p/2)(H_n - H_(n-m))) (mod p^2), 0 <= m <= n;
 (S)  exact: sum_k C(n,k)^4 (H_k - H_(n-k)) = 0, all n;
 (L)  n = (p-1)/2: sum_k a_k^2 a_(n-k)^2 == sum_k a_k^4 (mod p^2);
 (T)  the theorem.
Session opus-2026-10-17-S43.
"""
from math import comb
from fractions import Fraction as Fr

def a(k): return Fr(comb(2 * k, k), 4 ** k)
ok = {}
ok['E1'] = all(comb(n, k) * comb(2 * k, n) == comb(2 * k, k) * comb(k, n - k) for n in range(40) for k in range(n + 1))
ok['E2'] = all(sum((-1) ** (n - k) * 4 ** k * a(k) ** 3 * comb(k, n - k) for k in range(n + 1))
               == sum(a(k) ** 2 * a(n - k) ** 2 for k in range(n + 1)) for n in range(41))
H = lambda m: sum(Fr(1, j) for j in range(1, m + 1))
ok['S'] = all(sum(comb(n, k) ** 4 * (H(k) - H(n - k)) for k in range(n + 1)) == 0 for n in range(40))
primes = [p for p in range(3, 400) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
def modp2(x, p):
    m = p * p
    return x.numerator * pow(x.denominator, -1, m) % m
c1 = c2 = L = T = True
for p in primes:
    n = (p - 1) // 2
    for k in range(n + 1):
        if (comb(n, k) * comb(n + k, k) - modp2((-1) ** k * a(k) ** 2, p)) % (p * p):
            c1 = False
        if (modp2(a(k), p) - modp2((-1) ** k * comb(n, k) * (1 - Fr(p, 2) * (H(n) - H(n - k))), p)) % (p * p):
            c2 = False
    if (modp2(sum(a(k) ** 2 * a(n - k) ** 2 for k in range(n + 1)), p) - modp2(sum(a(k) ** 4 for k in range(n + 1)), p)) % (p * p):
        L = False
    s7 = sum(comb(n, k) ** 2 * comb(n + k, k) * comb(2 * k, n) for k in range(n + 1))
    A = sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1))
    if (s7 - (-1) ** n * A) % (p * p):
        T = False
ok.update(C1=c1, C2=c2, L=L, T=T)
print('E1 (exact, n < 40):', ok['E1'])
print('E2 (exact Clausen identity, n <= 40):', ok['E2'])
print('S  (exact symmetry, n < 40):', ok['S'])
print('C1, C2, L, T (all odd primes 3 <= p < 400):', ok['C1'], ok['C2'], ok['L'], ok['T'])
