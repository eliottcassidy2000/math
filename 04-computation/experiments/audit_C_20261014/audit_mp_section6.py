#!/usr/bin/env python3
"""Audit C: what McCarthy-Papanikolas (arXiv:1205.1006) Theorem 1.6 and their Section 6 observation say numerically.

Greene's functions with all parameters quadratic (phi) / trivial (e, e(0) = 0) over F_q:
  S1(x) = sum_y phi(y)phi(1-y)phi(1-xy) (x != 0), S1(0) = 0;  2F1(x) = phi(-1) S1(x)/q
  S2(x) = sum_y S1(xy) phi(y)phi(1-y);                         3F2(x) = S2(x)/q^2
  S3(x) = sum_y S2(xy) phi(y)phi(1-y);                         4F3(x) = phi(-1) S3(x)/q^3
McCarthy's normalisation for these parameters is (-q)^3 x Greene (binom(phi, e) = -1/q), so
  4F3^McC(phi,phi,phi,phi; e,e,e | -1)_q = -phi(-1) S3(-1).
Tests:
  (T1) q = p:   MP Theorem 1.6 + (1.5) [conditional on van Geemen-van Straten]: value = a(p) b(p) (p = 1 mod 4), 0 (p = 3 mod 4)
  (T2) q = p^2: MP Section 6 observation (p < 20): value = a_p^2 - 2 a_{p^2}; if L(F_7) = L(f1 x f2) this is the trace of
       Frob_p^2 on rho_f1 (x) rho_f2 = (a(p)^2 - 2p)(b(p)^2 - 2 chi_-4(p) p^2).  For p = 3 mod 4 (a(p) = 0) that is
       2p(|b(p)|^2 - 2p^2): an F_{p^2} hypergeometric expression for |b(p)|^2 (prior to HYP-9276's F_p expression).
a(p) of f1 = eta(4z)^2 eta(8z)^2 from y^2 = x^3 - x; b(p) from audit_newform_bp.json.
"""
import sys
import json
import numpy as np
sys.path.insert(0, '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_C_20261014')
from audit_common import GF, chi_prime

D = '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_C_20261014/'
bp = {int(k): tuple(v) for k, v in json.load(open(D + 'audit_newform_bp.json')).items()}


def a_f1(p):
    c = chi_prime(p)
    x = np.arange(p)
    return -int(c[(x ** 3 - x) % p].sum())


def S3_at_minus1(F):
    q = F.q
    c = F.chi
    R = np.arange(q)
    one_minus = np.array([F.sub_scalar(1, int(y)) for y in R])
    k = c * c[one_minus]                                   # phi(y) phi(1-y)
    # multiplication table rows via logs: xy
    lg, ex = F.log, F.exp
    def mulrow(x):
        if x == 0:
            return np.zeros(q, dtype=np.int64)
        out = np.zeros(q, dtype=np.int64)
        nz = R[1:]
        out[1:] = ex[(lg[x] + lg[nz]) % (q - 1)]
        return out
    S1 = np.zeros(q, dtype=np.int64)
    for x in range(1, q):
        S1[x] = int(np.dot(k, c[one_minus[mulrow(x)]]))
    S2 = np.zeros(q, dtype=np.int64)
    for x in range(q):
        S2[x] = int(np.dot(k, S1[mulrow(x)]))
    m1 = F.neg_one
    return int(np.dot(k, S2[mulrow(m1)]))


def main():
    print("T1 (q = p): MP Thm 1.6 value -phi(-1) S3(-1) vs a(p)b(p) [p = 1 mod 4] or 0 [p = 3 mod 4]")
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 37, 41):
        F = GF(p, 1)
        val = -F.eps * S3_at_minus1(F)
        a = a_f1(p)
        b = bp[p]
        pred = a * b[0] if p % 4 == 1 else 0
        print(f"  p = {p:2d}: value {val:8d}   a(p) = {a:3d}, b(p) = {b[0]}{b[1]:+d}i   predicted {pred:8d}   "
              f"{'ok' if val == pred else 'differs'}")
    print("T2 (q = p^2): MP Sec. 6 value -S3(-1) over F_{p^2} vs (a(p)^2 - 2p)(b(p)^2 - 2 chi(p) p^2)")
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        F = GF(p, 2)
        val = -F.eps * S3_at_minus1(F)
        a = a_f1(p)
        b = bp[p]
        chi4 = 1 if p % 4 == 1 else -1
        b2 = b[0] ** 2 - b[1] ** 2                         # b(p)^2 (real)
        pred = (a * a - 2 * p) * (b2 - 2 * chi4 * p * p)
        extra = ''
        if p % 4 == 3:
            extra = f"  -> |b(p)|^2 = 2p^2 + value/(2p) = {2 * p * p + val / (2 * p):.1f} (true {b[1] ** 2})"
        print(f"  p = {p:2d}: value {val:12d}   predicted {pred:12d}   {'ok' if val == pred else 'differs'}{extra}")


if __name__ == '__main__':
    main()
