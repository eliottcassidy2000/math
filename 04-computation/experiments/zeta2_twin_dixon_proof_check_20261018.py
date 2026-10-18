#!/usr/bin/env python3
"""Check of the proof of HYP-9292 (S43): closed forms from the derivative of Dixon's summation.

Dixon: 3F2(a,b,c;1+a-b,1+a-c;1) = G(1+a/2)G(1+a-b)G(1+a-c)G(1+a/2-b-c)/(G(1+a)G(1+a/2-b)G(1+a/2-c)G(1+a-b-c)).
With a = -n + d, b = c = -n and d/dd at d = 0, with D_n = sum (-1)^k C(n,k)^3, S_n = sum (-1)^k C(n,k)^3 H_(n-k),
E_n = S_n - H_n D_n:
  (odd n)   S_n = G(1-n/2) G(1+3n/2) / (3n G(1+n/2)^2)          [exact; p = 2n+1 divides the numerator once]
  (n = 2m)  E_n = (D_n/2) (H_m - H_(3m)),  D_n = (-1)^m (3m)!/m!^3   [exact; H_(3m) - H_m == 0 mod p = 4m+1]
and the theorem: sum_{k<=n} C(n+k,k)^3 == sum_k C(n,k)^2 C(n+k,k) (mod p^2), n = (p-1)/2, p >= 5.
Session opus-2026-10-17-S43.
"""
from math import comb, factorial
from fractions import Fraction as Fr

def H(m): return sum(Fr(1, j) for j in range(1, m + 1))
def dfact(x):   # odd double factorial x!! for odd x >= -1
    r = 1
    while x > 1:
        r *= x; x -= 2
    return r
def gamma_half_ratio(n):
    """G(1-n/2) G(1+3n/2) / G(1+n/2)^2 for odd n, as an exact rational (pi factors cancel)."""
    m = (n - 1) // 2
    g1 = Fr((-1) ** m * 2 ** m, dfact(n - 2))           # G(1/2 - m)/sqrt(pi)
    g2 = Fr(dfact(3 * n), 2 ** ((3 * n + 1) // 2))       # G(3n/2 + 1)/sqrt(pi)
    g3 = Fr(dfact(n), 2 ** ((n + 1) // 2))               # G(n/2 + 1)/sqrt(pi)
    return g1 * g2 / (g3 * g3)
ok_odd = ok_even = True
for n in range(1, 41):
    D = sum((-1) ** k * comb(n, k) ** 3 for k in range(n + 1))
    S = sum((-1) ** k * comb(n, k) ** 3 * H(n - k) for k in range(n + 1))
    if n % 2:
        if S != gamma_half_ratio(n) / (3 * n):
            ok_odd = False
    else:
        m = n // 2
        if D != (-1) ** m * Fr(factorial(3 * m), factorial(m) ** 3) or S - H(n) * D != Fr(D, 2) * (H(m) - H(3 * m)):
            ok_even = False
print('odd-n closed form for S_n (n <= 40):', ok_odd)
print('even-n closed form for E_n and Dixon value (n <= 40):', ok_even)
primes = [p for p in range(5, 1200) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
ok_thm = all((sum(comb((p - 1) // 2 + k, k) ** 3 for k in range((p + 1) // 2))
              - sum(comb((p - 1) // 2, k) ** 2 * comb((p - 1) // 2 + k, k) for k in range((p + 1) // 2))) % (p * p) == 0
             for p in primes)
print('theorem for all primes 5 <= p < 1200:', ok_thm)
ok_sym = all(sum(Fr(1, j) for j in range(m + 1, 3 * m + 1)).numerator % (4 * m + 1) == 0
             for m in range(1, 300) if all((4 * m + 1) % d for d in range(2, int((4 * m + 1) ** 0.5) + 1)))
print('H_(3m) - H_m == 0 mod p for all primes p = 4m + 1 < 1200:', ok_sym)
