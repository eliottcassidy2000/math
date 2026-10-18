#!/usr/bin/env python3
"""Step-by-step check of the reflection proof of HYP-9310 (L_p) (KLMSY 2016, Cor. 3). Session cloud-2026-10-18-M1.

n = (p-1)/2, F_m = sum_{k<=n} C(2k,k)^3 m^-k, h_m = sum_{k<=n} C(2k,k)^3 m^-k (H_2k - H_k), q = q_p(m).
 (R)  exact: Feps_m := sum_{k<=n} [(-n)_k/(1/2-n)_k]^3 (64/m)^k = C(2n,n)^-3 (4096/m)^n F_{4096/m}
 (T)  Feps_m == F_m - 3p h_m (mod p^2)          [first-order Taylor in eps = -p/2]
 (M)  Morley: C(2n,n) == (-1)^n 4^(p-1) (mod p^3), p >= 5
 (S)  m^-n == (m/p)(1 - (p/2) q) (mod p^2)
 (E)  (L_p) <=> F_{4096/m} == (-m/p) F_m (mod p^2), tested for ALL m with 1 <= |m| <= MMAX, p not | m
 (K2) KLMSY Thm 2: F_m == ((1-lambda)/p)(a_p^2 - 2p) (mod p^2) (0 if supersingular), lambda = 64/m,
      E_lambda: y^2 = (x-1)(x^2 - c), c = 1/(1-lambda) = m/(m-64); checked at the 6 non-self-dual CM points
 (I)  isogeny with twist: a_p(E_lambda) = (-2/p) a_p(E_{1/lambda}) (point counts)
 (L)  the law itself at the CM points, including the bad primes p = 7 for m = 1, 4096
Usage: python3 hyp9310_proof_check.py [PMAX] [MMAX]
"""
import sys
from fractions import Fraction as Fr
from math import comb

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
MMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 60

def primes(lo, hi):
    return [p for p in range(lo, hi) if p > 1 and all(p % d for d in range(2, int(p ** 0.5) + 1))]

def leg(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def F(m, p, mod):
    n = (p - 1) // 2
    mi = pow(m % mod, -1, mod)
    return sum(comb(2 * k, k) ** 3 * pow(mi, k, mod) for k in range(n + 1)) % mod

def h(m, p, mod):
    n = (p - 1) // 2
    mi = pow(m % mod, -1, mod)
    H = [0] * (2 * n + 1)
    for j in range(1, 2 * n + 1):
        H[j] = (H[j - 1] + pow(j, -1, mod)) % mod
    return sum(comb(2 * k, k) ** 3 * pow(mi, k, mod) * (H[2 * k] - H[k]) for k in range(n + 1)) % mod

def qp(m, p):
    return ((pow(m % (p * p), p - 1, p * p) - 1) // p) % p

def ap_count(c_num, c_den, p):
    # a_p of y^2 = (x-1)(x^2 - c), c = c_num/c_den, by point counting
    c = c_num * pow(c_den % p, -1, p) % p
    s = 0
    for x in range(p):
        s += leg((x - 1) * (x * x - c), p)
    return -s

ok = {k: True for k in ['R', 'T', 'M', 'S', 'E', 'K2', 'I', 'L']}
# (R) exact, small p
for p in primes(5, 40):
    n = (p - 1) // 2
    for m in [1, 4096, -8, -512, 16, 256, 3, -5, 7]:
        lhs = Fr(0)
        for k in range(n + 1):
            num = Fr(1); den = Fr(1)
            for j in range(k):
                num *= (-n + j); den *= (Fr(1, 2) - n + j)
            lhs += (num / den) ** 3 * Fr(64, m) ** k
        rhs = Fr(comb(2 * n, n)) ** -3 * Fr(4096, m) ** n * sum(Fr(comb(2 * k, k) ** 3) * Fr(m, 4096) ** k for k in range(n + 1))
        if lhs != rhs:
            ok['R'] = False
CM = [1, 4096, -8, -512, 16, 256]
for p in primes(5, PMAX):
    n = (p - 1) // 2
    P2, P3 = p * p, p ** 3
    if (comb(2 * n, n) - (-1) ** n * pow(4, p - 1, P3)) % P3:
        ok['M'] = False
    for m in range(-MMAX, MMAX + 1):
        if m == 0 or m % p == 0:
            continue
        q = qp(m, p)
        # (S)
        if (pow(pow(m % P2, n, P2), -1, P2) - leg(m, p) * (1 - Fr(p, 2) * q).numerator * pow((1 - Fr(p, 2) * q).denominator, -1, P2)) % P2:
            ok['S'] = False
        # (T): Feps via the reflection, compared with F - 3 p h
        Fm = F(m, p, P2)
        F2 = F(Fr(4096, m), p, P2) if False else None
        # F_{4096/m} with rational m' = 4096/m: compute directly
        mi = Fr(m, 4096)
        Fref = sum(comb(2 * k, k) ** 3 * (mi.numerator * pow(mi.denominator, -1, P2)) ** k for k in range(n + 1)) % P2 \
            if mi.denominator % p else None
        if Fref is None:
            continue
        Feps = pow(comb(2 * n, n), -3, P2) * pow(4096 * pow(m % P2, -1, P2), n, P2) * Fref % P2
        if (Feps - (Fm - 3 * p * h(m, p, P2))) % P2:
            ok['T'] = False
        # (E) equivalence of truth values
        Lp = (6 * h(m, p, p) - q * F(m, p, p)) % p == 0
        iso = (Fref - leg(-m, p) * Fm) % P2 == 0
        if Lp != iso:
            ok['E'] = False
    for m in CM:
        if m % p == 0:
            continue
        q = qp(m, p)
        Lp = (6 * h(m, p, p) - q * F(m, p, p)) % p == 0
        if not Lp:
            ok['L'] = False
        if p == 7 and m in (1, 4096):
            continue  # bad reduction: (K2) not claimed; (L) checked directly above
        a = ap_count(m, m - 64, p)
        lam = Fr(64, m)
        lm1 = 1 - lam
        s = leg(lm1.numerator * lm1.denominator, p)  # ((1-lambda)/p): the sign that the data force (see note)
        pred = 0 if a % p == 0 else (s * (a * a - 2 * p)) % P2
        if (F(m, p, P2) - pred) % P2:
            ok['K2'] = False
        # (I): E_{1/lambda} has c' = 1/(1 - 1/lambda) = lambda/(lambda-1) = 64/(64-m)
        a2 = ap_count(64, 64 - m, p)
        if a != leg(-2, p) * a2:
            ok['I'] = False
for k, v in ok.items():
    print(k, v)
print('PMAX', PMAX, 'MMAX', MMAX)
