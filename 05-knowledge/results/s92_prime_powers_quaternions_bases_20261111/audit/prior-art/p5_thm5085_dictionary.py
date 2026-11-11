#!/usr/bin/env python3
"""Auditor prior-art (S92).  THM-5134 (H) against THM-5085 (2) (kind-pasteur S19, fixed-port tails in every base).

THM-5085: port state (h, R, Delta); letters S: (h, 2R, 0), B: (h + x (R + Delta)/2, 2R - h_u, -h_u), h_u fixed;
closed form h = h0 + (x/2) [R0 N(I) + Delta0 [1 in I] - h_u Q(I)], Q(2m) = Q(m), Q(2m+1) = Q(m) + m + (m mod 2).
Claim checked: for a fixed smaller state (h1, s) the base-b count of THM-5134 (H),
    N_b(h1, m, s) = (1 + b psi(h1)) (1 + b psi(m)) + b (2^s y1 y(m) - [s = 0])/2,
is the THM-5085 tail with (h0, R0, Delta0, h_u) = (X1 + (b/2)(2^s y1 - [s = 0]), 2 (X1 + 2^s y1), 0, X1), X1 = 1 + b psi(h1),
both by the closed form and by running the letters on the binary digits of m (least significant first).
Exact rational arithmetic (b odd makes b/2 a half-integer; the total is an integer)."""
import sys
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def psi(h):
    if h == 0:
        return 0
    if h % 2 == 0:
        return psi(h // 2) + h // 2
    k = h // 2
    return psi(k) + k // 2 + 1


def Q(m):
    if m == 0:
        return 0
    if m % 2 == 0:
        return Q(m // 2)
    k = m // 2
    return Q(k) + k + (k % 2)


ok_closed = ok_run = True
n = 0
for b in (2, 3, 4, 5, 6, 7, 10):
    for h1 in range(0, 24):
        for s in range(0, 4):
            X1 = 1 + b * psi(h1)
            y1 = 2 * h1 + 1
            h0 = Fraction(X1) + Fraction(b, 2) * ((1 << s) * y1 - (1 if s == 0 else 0))
            R0 = Fraction(2 * (X1 + (1 << s) * y1))
            hu = Fraction(X1)
            for m in range(0, 200):
                Nb = (1 + b * psi(h1)) * (1 + b * psi(m)) + Fraction(b * ((1 << s) * y1 * (2 * m + 1) - (1 if s == 0 else 0)), 2)
                closed = h0 + Fraction(b, 2) * (R0 * m - hu * Q(m))
                if closed != Nb:
                    ok_closed = False
                h, R, D = h0, R0, Fraction(0)
                mm = m
                while mm:
                    if mm & 1:
                        h, R, D = h + Fraction(b, 2) * (R + D), 2 * R - hu, -hu
                    else:
                        h, R, D = h, 2 * R, Fraction(0)
                    mm >>= 1
                if h != Nb or Nb.denominator != 1:
                    ok_run = False
                n += 1
print("cases:", n)
print("closed form of THM-5085 (2) = N_b of THM-5134 (H):", "PASS" if ok_closed else "FAIL")
print("letter recursion of THM-5085 (1) on the digits of m = N_b:", "PASS" if ok_run else "FAIL")
