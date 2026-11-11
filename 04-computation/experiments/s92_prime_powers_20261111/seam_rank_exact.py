#!/usr/bin/env python3
"""S92 (THM-5134 (F), THM-5135 (1), (4)): exact class laws of bilinear seams with uniform states modulo a prime power.

(1)  [THM-5135 (1)] For u, v independent and uniform in (Z/p^k)^r and a bilinear form of unit determinant, a class c with
     v = v_p(c) < k has   p^k P(u . v = c) = (1 - p^-r) sigma_(r-1)(p^v) / p^((r-1) v).
(F)  [THM-5134 (F)] For odd p the quaternary forms a1 a2 + u b1 b2 (u a unit) and a^2 + b^2 + c^2 + d^2 take every class modulo
     p^k equally often (they are equivalent over Z_p): the F1 seam has the local densities of the four-square form.

Exact rational arithmetic (class counts by convolution of the distributions of the summands).  Checked:
  * THM-5135 (1) for r = 1, 2, 3, 4 at p^k = 9, 27, 25, 49 and for r = 1, 2, 3 at 2^k = 8, 16 (the statement holds at p = 2
    for a seam with unit determinant; the F1 seam is not of that kind at 2);
  * a non-standard kernel of unit determinant ([[1, 2], [3, 5]], det = -1) gives the same law;
  * (F): class counts of a1 a2 + 2 b1 b2 and of the sum of four squares agree modulo 9, 27, 25, 49, 121;
  * the degenerate kernel [[149, 72], [72, 37]] of Session261 (det = 329 = 7 * 47) is 2 (h + c)(g + d) modulo 7:
    its class law modulo 7 is the rank-one law P(0) = 2/p - 1/p^2, not the rank-two law.
Exit status 1 on failure.
"""
import sys
from fractions import Fraction as F

import numpy as np


def conv(a, b, q):
    out = np.zeros(q, dtype=object)
    for i in range(q):
        if a[i]:
            for j in range(q):
                if b[j]:
                    out[(i + j) % q] += a[i] * b[j]
    return out


def dist_product(q, coeff=1):
    """distribution (counts) of coeff * a * b over a, b in Z/q."""
    d = np.zeros(q, dtype=object)
    for a in range(q):
        for b in range(q):
            d[(coeff * a * b) % q] += 1
    return d


def dist_square(q):
    d = np.zeros(q, dtype=object)
    for a in range(q):
        d[(a * a) % q] += 1
    return d


def vp(c, p, k):
    if c == 0:
        return k
    v = 0
    while c % p == 0:
        c //= p
        v += 1
    return v


def predicted(p, k, r, v):
    return (1 - F(1, p ** r)) * sum(F(1, p ** ((r - 1) * i)) for i in range(v + 1))


def main():
    ok = True
    print("THM-5135 (1): p^k * P(u . v = c) by valuation of c, uniform states, seam of rank r: exact = formula ?")
    for p, k in ((3, 2), (3, 3), (5, 2), (7, 2), (2, 3), (2, 4)):
        q = p ** k
        one = dist_product(q)
        d = None
        for r in (1, 2, 3, 4):
            if p == 2 and r == 4:
                continue
            d = one if r == 1 else conv(d, one, q)
            tot = q ** (2 * r)
            good = True
            vals = {}
            for c in range(q):
                v = vp(c, p, k)
                f = F(int(d[c]) * q, tot)
                vals.setdefault(v, set()).add(f)
                if v < k and f != predicted(p, k, r, v):
                    good = False
            ok &= good and all(len(s) == 1 for s in vals.values())
            print(f"  p^k = {q:3d}, r = {r}: " + ", ".join(f"v={v}: {sorted(s)[0]}" for v, s in sorted(vals.items()))
                  + f"   {'ok' if good else 'FAIL'}")
    # a non-standard unimodular kernel
    q = 27
    d = np.zeros(q, dtype=object)
    K = ((1, 2), (3, 5))
    for u1 in range(q):
        for u2 in range(q):
            a = (K[0][0] * u1 + K[1][0] * u2) % q
            b = (K[0][1] * u1 + K[1][1] * u2) % q
            for v1 in range(q):
                base = a * v1
                for v2 in range(q):
                    d[(base + b * v2) % q] += 1
    good = all(F(int(d[c]) * q, q ** 4) == predicted(3, 3, 2, vp(c, 3, 3)) for c in range(q) if vp(c, 3, 3) < 3)
    ok &= good
    print(f"kernel [[1, 2], [3, 5]] (det -1) modulo 27: the law of rank 2: {good}")
    # (F): the F1 seam against the sum of four squares
    print("THM-5134 (F): class counts of a1 a2 + 2 b1 b2 and of a^2 + b^2 + c^2 + d^2 modulo p^k:")
    for q in (9, 27, 25, 49, 121):
        f1 = conv(dist_product(q), dist_product(q, 2), q)
        sq = dist_square(q)
        s2 = conv(sq, sq, q)
        s4 = conv(s2, s2, q)
        same = all(int(f1[c]) == int(s4[c]) for c in range(q))
        ok &= same
        print(f"  modulo {q:3d}: equal for every class: {same}")
    # the kernel of Session261 modulo 7 and 47 (degenerate) and modulo 5 (unit determinant)
    print("kernel K = [[149, 72], [72, 37]] (det 329 = 7 * 47): p * P(K(u, v) = c) for uniform u, v modulo p")
    for p in (5, 7, 47):
        d = np.zeros(p, dtype=object)
        for h in range(p):
            for c in range(p):
                a = (149 * h + 72 * c) % p
                b = (72 * h + 37 * c) % p
                for g in range(p):
                    for dd in range(p):
                        d[(a * g + b * dd) % p] += 1
        f0 = F(int(d[0]) * p, p ** 4)
        f1 = F(int(d[1]) * p, p ** 4)
        rank2 = (1 + F(1, p) - F(1, p * p), 1 - F(1, p * p))
        rank1 = (2 - F(1, p), 1 - F(1, p))
        kind = "rank 2" if (f0, f1) == rank2 else ("rank 1" if (f0, f1) == rank1 else "other")
        print(f"  p = {p:2d}: class 0: {f0}, class 1: {f1}  -> {kind} law")
        ok &= kind == ("rank 2" if p == 5 else "rank 1")
    fac = all((149 * h * g + 72 * (h * d + c * g) + 37 * c * d - 2 * (h + c) * (g + d)) % 7 == 0
              for h in range(7) for c in range(7) for g in range(7) for d in range(7))
    print(f"  K(h, c; g, d) = 2 (h + c)(g + d) modulo 7: {fac}")
    ok &= fac
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
