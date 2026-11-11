#!/usr/bin/env python3
"""Auditor proof-5128-CF, check A7: THM-5128 (F), divisible states (own code, exact integers).

(1) x(h) = y(h) iff h = 0 or h is a power of two, for all h < 2^24; then x = y = 2h + 1.
(2) first p-divisible state (p | gcd(x, y)) for the primes quoted in (F), and, for every odd prime p < 300 that
    divides some 2^m + 1, the least such m, the state 2^(m-1), the first p-divisible state and p^2 (the generic scale):
    is "the prime factors of 2^m + 1 have early p-divisible states" true as a general statement?
(3) proportion of h < 2^24 with gcd(x, y) = 1 (file: 0.81176), density of p-divisible states times p^2,
    density of gcd(g, m) = 1 for m = 15, 105 against prod (1 - 1/p^2).
(4) sign of A_p (partial sums to 2^24) against "p-divisible state early" for the primes < 100.
"""
import sys
import numpy as np
import sympy
from a_lib import psi_table

NB = 24
H = 1 << NB
psi = psi_table(H)
h = np.arange(H, dtype=np.int64)
x = 2 * psi + 1
y = 2 * h + 1
q = h - psi
g = np.gcd(x, y)
ok = True

eq = np.nonzero(x == y)[0]
pw = [0] + [1 << j for j in range(NB)]
good1 = list(eq) == pw
ok &= good1
print("(1) {h < 2^%d : x = y} = {0} + powers of two: %s  (%d values); x = y = 2^(j+1) + 1 there: %s"
      % (NB, good1, len(eq), all(int(x[1 << j]) == (1 << (j + 1)) + 1 for j in range(NB))))
# the mechanism: q(h) = sum over one-bits i of ceil(floor(h / 2^(i+1)) / 2)
qq = np.zeros(H, dtype=np.int64)
for i in range(NB):
    Hi = h >> (i + 1)
    qq += ((h >> i) & 1) * ((Hi + 1) >> 1)
print("    q(h) = sum over one-bits i of ceil(H_i / 2) for all h < 2^%d: %s" % (NB, bool(np.all(qq == q))))


def first_div(p):
    idx = np.nonzero(g % p == 0)[0]
    return int(idx[0]) if len(idx) else None


quoted = {7: 38, 19: 218, 23: 103, 29: 1435, 31: 8664, 37: 1350, 41: 266, 127: 109283}
got = {p: first_div(p) for p in quoted}
print("(2) first p-divisible states, file:", quoted)
print("                              mine:", got)
ok &= got == quoted
fermat = {3: 1, 5: 2, 17: 8, 11: 16, 13: 32, 43: 64, 257: 128, 65537: 2 ** 15}
print("    states quoted for prime factors of 2^m + 1 (p: h): p | g(h) holds:",
      all(int(g[hv]) % p == 0 for p, hv in fermat.items()),
      "; is it the FIRST p-divisible state:", {p: first_div(p) == hv for p, hv in fermat.items()})
print("    p | least m with p | 2^m + 1 | state 2^(m-1) | first p-divisible state | p^2 | 2^(m-1) < p^2 ? | first < p^2 ?")
for p in sympy.primerange(3, 300):
    o = sympy.n_order(2, p)
    fd = first_div(p)
    if o % 2 == 0:
        m = o // 2
        st = 1 << (m - 1)
        print("   %3d | %3d | %s | %s | %6d | %s | %s" % (p, m, st if m <= 40 else "2^%d" % (m - 1), fd, p * p,
                                                      st < p * p, (fd is not None and fd < p * p)))
    else:
        print("   %3d | none (ord_p(2) = %d odd) | - | %s | %6d | - | %s" % (p, o, fd, p * p,
                                                                     (fd is not None and fd < p * p)))

print("(3) proportion of h < 2^%d with gcd(x, y) = 1: %.5f   (file: 0.81176; 8/pi^2 = %.5f)"
      % (NB, float(np.mean(g == 1)), 8 / np.pi ** 2))
for p in (3, 5, 7, 11, 13, 17, 31, 127):
    print("    p = %3d: #p-divisible states below 2^%d = %d, times p^2 / 2^%d = %.4f"
          % (p, NB, int(np.sum(g % p == 0)), NB, float(np.sum(g % p == 0)) * p * p / H))
for m, primes in ((15, (3, 5)), (105, (3, 5, 7))):
    print("    density of gcd(g, %d) = 1 below 2^%d: %.5f, product %.5f"
          % (m, NB, float(np.mean(np.gcd(g, m) == 1)), float(np.prod([1 - 1 / p ** 2 for p in primes]))))

print("(4) sign of A_p^0 (partial sum over h < 2^%d) against the position of the first p-divisible state" % NB)
xf = x.astype(np.float64)
yf = y.astype(np.float64)
w0 = np.log(2 * (xf + yf) / (xf + 2 * yf)) / xf
print("    p | first p-div state | first / p^2 | number of p-div states below p^2 | A_p^0 (2^%d) | divides some 2^m + 1" % NB)
for p in sympy.primerange(3, 100):
    div = (g % p == 0)
    A0 = float(np.sum(div * w0) - np.sum(w0) / p ** 2)
    fd = first_div(p)
    print("   %3d | %6s | %7.3f | %2d | %+.5f | %s" % (p, fd, fd / p ** 2, int(div[: p * p].sum()), A0,
                                                   sympy.n_order(2, p) % 2 == 0))
print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
sys.exit(0 if ok else 1)
