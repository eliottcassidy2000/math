#!/usr/bin/env python3
"""(1) McCarthy-Papanikolas Thm 1.8 at prime powers: H_Q := 3F2^*(chi4, phi, phi; eps, eps | 1)_Q equals alpha^r + beta^r
(trace of Frob_Q on rho_g) for Q = p^r = 1 mod 4, including the inert squares p^2 with p = 3 mod 4.
(2) The descent form of HYP-9276: M_q = -H_(q^2) - chi_(-4)(q) q^2 - 2q - 1 for primes q (M_q from k6p6_fft.json).
Complex-float Gauss sums; values are compared with exact integers to 1e-6 relative.  Usage: python3 k6p6_descent_check.py"""
import json
from mp_3f2_field import F32_star, B

ap = {}
for line in open('g_coeffs.out'):
    parts = line.split()
    if len(parts) == 3 and parts[0].isdigit():
        ap[int(parts[0])] = complex(int(parts[1]), int(parts[2]))


def sat(p, r):
    b = ap[p]
    chi = 1 if p % 4 == 1 else -1
    s = [2, b]
    for _ in range(2, r + 1):
        s.append(b * s[-1] - chi * p * p * s[-2])
    return s[r]


ok1 = 0
cases = [(5, 1), (13, 1), (17, 1), (29, 1), (37, 1), (41, 1), (3, 2), (3, 4), (3, 6), (5, 2), (5, 3), (5, 4), (7, 2),
         (11, 2), (13, 2), (17, 2), (19, 2), (23, 2), (29, 2), (31, 2)]
for p, r in cases:
    v, t = F32_star(p, r), sat(p, r)
    assert abs(v - t) < 1e-6 * max(1, abs(t)), (p, r, v, t)
    ok1 += 1
print(f"(1) H_Q = alpha^r + beta^r at all {ok1} prime powers Q = 1 mod 4: {[p ** r for p, r in cases]}")

M = {int(k): v[0] for k, v in json.load(open('k6p6_fft.json')).items()}
ok2 = []
for q in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
    H = F32_star(q, 2)
    chi = 1 if q % 4 == 1 else -1
    rhs = -H - chi * q * q - 2 * q - 1
    assert abs(M[q] - rhs) < 1e-6 * max(1, abs(M[q])), (q, M[q], rhs)
    ok2.append(q)
print(f"(2) descent identity M_q = -H_(q^2) - chi(q) q^2 - 2q - 1 holds at q = {ok2}")
