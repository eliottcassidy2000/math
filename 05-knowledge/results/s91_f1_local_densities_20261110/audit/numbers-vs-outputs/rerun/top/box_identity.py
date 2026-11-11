#!/usr/bin/env python3
"""S91 (THM-5128): the exact class law of a bilinear join over two boxes, and its first-order form.

For an odd prime p let D_n be the distribution of v(h) = (x(h), y(h)) mod p over h in [0, 2^n) (computed EXACTLY by the
digit recursion on the state (h mod p, h mod 2, q mod p); no enumeration of h), pi_n = D_n(0, 0) the proportion of
p-divisible states.  For independent h1 in [0, 2^n1), h2 in [0, 2^n2) and a unit u = 2^s:

    P( x1 x2 + u y1 y2 = c' )  =  1/p + (pi_n1 + pi_n2 - 1/p^2) ( [c' = 0] - 1/p ) + Rhat(c'),             (*)

    Rhat(c') = (1/p) sum_{t != 0} e(-t c'/p) sum_{v != 0} delta_1(v) eps_2(t v_u),

delta_1 = D_n1 - 1/p^2, eps_2 = the non-trivial Fourier coefficients of D_n2, v_u = (a, u b).  The F1 class is
c = c' - [s = 0], so the distinguished class is N = 0 for s >= 1 and N = -1 for s = 0.

This script
  1. checks the exact digit recursion for D_n against brute force (n <= 14);
  2. checks (*) as an identity in exact rational arithmetic (Fractions) for small p, and in floating point for the
     other primes, and checks that sum_c' Rhat(c') = 0;
  3. checks the box law against a brute-force count of the F1 values over a box (n1 = n2 = 9);
  4. prints, for n1 = n2 = n, the first-order coefficient pi_n1 + pi_n2 - 1/p^2, the size of Rhat and the proven
     bound |Rhat| <= (p - 1) max|eps_1| max|eps_2| ... (see THM-5128), for n = 10, 20, 40, 80: the limit law is
     the generic one, P(c' = 0) -> 1/p + (1/p^2)(1 - 1/p), P(c' != 0) -> 1/p - 1/p^3.

Usage: python3 box_identity.py
"""
import sys
from fractions import Fraction

import numpy as np

PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 61]


def q_table(n):
    size = 1 << n
    q = np.zeros(size, dtype=np.int64)
    for k in range(n):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    return q


def exact_counts(p, nmax):
    """counts[n][a][b] = #{h < 2^n : x(h) = a, y(h) = b mod p}, exact Python integers, for n = 0..nmax."""
    inv2 = pow(2, -1, p)
    # state: (eta, beta, r) with eta = h mod p, beta = h mod 2, r = q mod p
    cur = {(0, 0, 0): 1}
    out = []
    for n in range(nmax + 1):
        tab = [[0] * p for _ in range(p)]
        for (eta, beta, r), c in cur.items():
            b = (2 * eta + 1) % p
            a = (b - 2 * r) % p
            tab[a][b] += c
        out.append(tab)
        nxt = {}
        for (eta, beta, r), c in cur.items():
            k0 = ((2 * eta) % p, 0, r)
            nxt[k0] = nxt.get(k0, 0) + c
            k1 = ((2 * eta + 1) % p, 1, (r + (eta + beta) * inv2) % p)
            nxt[k1] = nxt.get(k1, 0) + c
        cur = nxt
    return out


def class_law_direct(D1, D2, p, u, exact=False):
    """distribution of a1 a2 + u b1 b2 mod p for independent (a1,b1) ~ D1, (a2,b2) ~ D2."""
    zero = Fraction(0) if exact else 0.0
    out = [zero] * p
    for a1 in range(p):
        for b1 in range(p):
            w1 = D1[a1][b1]
            if w1 == 0:
                continue
            for a2 in range(p):
                for b2 in range(p):
                    w2 = D2[a2][b2]
                    if w2:
                        out[(a1 * a2 + u * b1 * b2) % p] += w1 * w2
    return out


def fourier(D, p):
    """Dhat[alpha][beta] = sum D[a][b] e((alpha a + beta b)/p) as a numpy array."""
    arr = np.array(D, dtype=np.float64)
    return np.conj(np.fft.fft2(arr))


def main():
    ok = True
    # 1. digit recursion against brute force
    nb = 14
    q = q_table(nb)
    h = np.arange(1 << nb, dtype=np.int64)
    y = 2 * h + 1
    x = y - 2 * q
    for p in PRIMES:
        cnts = exact_counts(p, nb)
        for n in (0, 1, 5, 9, nb):
            ref = np.zeros((p, p), dtype=np.int64)
            np.add.at(ref, (x[: 1 << n] % p, y[: 1 << n] % p), 1)
            ok &= bool((np.array(cnts[n], dtype=np.int64) == ref).all())
    print("1. exact digit recursion for the distribution of (x, y) mod p agrees with brute force (n <= 14):", ok)

    # 2. identity (*) in exact arithmetic for p = 3, 5, 7 and several (n1, n2, u)
    exact_ok = True
    for p in (3, 5, 7):
        cnts = exact_counts(p, 12)
        for (n1, n2) in ((3, 4), (6, 9), (12, 12)):
            D1 = [[Fraction(c, 1 << n1) for c in row] for row in cnts[n1]]
            D2 = [[Fraction(c, 1 << n2) for c in row] for row in cnts[n2]]
            for s in (0, 1, 2, 5):
                u = pow(2, s, p)
                law = class_law_direct(D1, D2, p, u, exact=True)
                exact_ok &= sum(law) == 1
                pi1, pi2 = D1[0][0], D2[0][0]
                first = [Fraction(1, p) + (pi1 + pi2 - Fraction(1, p * p)) * ((1 if c == 0 else 0) - Fraction(1, p))
                         for c in range(p)]
                rhat = [law[c] - first[c] for c in range(p)]
                exact_ok &= sum(rhat) == 0
                # Rhat must vanish identically when block 2 is exactly uniform: test with the uniform law
                U = [[Fraction(1, p * p)] * p for _ in range(p)]
                lawu = class_law_direct(D1, U, p, u, exact=True)
                firstu = [Fraction(1, p) + (pi1 + Fraction(1, p * p) - Fraction(1, p * p)) * ((1 if c == 0 else 0) - Fraction(1, p))
                          for c in range(p)]
                exact_ok &= lawu == firstu
    ok &= exact_ok
    print("2. (*) with Rhat summing to 0, and Rhat = 0 exactly when one block is uniform (exact rationals, p = 3, 5, 7):", exact_ok)

    # 2b. Fourier form of Rhat and the bound, floating point, all primes
    fl_ok = True
    worst_ratio = 0.0
    for p in PRIMES:
        cnts = exact_counts(p, 24)
        for (n1, n2) in ((8, 8), (12, 20), (24, 24)):
            D1 = np.array(cnts[n1], dtype=np.float64) / float(1 << n1)
            D2 = np.array(cnts[n2], dtype=np.float64) / float(1 << n2)
            F1, F2 = fourier(D1, p), fourier(D2, p)
            e1 = np.abs(F1).copy(); e1[0, 0] = 0.0
            e2 = np.abs(F2).copy(); e2[0, 0] = 0.0
            for s in (0, 1, 3):
                u = pow(2, s, p)
                law = np.array(class_law_direct(D1.tolist(), D2.tolist(), p, u))
                pi1, pi2 = D1[0, 0], D2[0, 0]
                first = np.array([1 / p + (pi1 + pi2 - 1 / p ** 2) * ((1.0 if c == 0 else 0.0) - 1 / p) for c in range(p)])
                rhat = law - first
                # Fourier evaluation of Rhat
                delta1 = D1 - 1 / p ** 2
                rf = np.zeros(p, dtype=np.complex128)
                for t in range(1, p):
                    tot = 0j
                    for a in range(p):
                        for b in range(p):
                            if a == 0 and b == 0:
                                continue
                            tot += delta1[a, b] * F2[(t * a) % p, (t * u * b) % p]
                    for c in range(p):
                        rf[c] += np.exp(-2j * np.pi * t * c / p) * tot / p
                fl_ok &= bool(np.max(np.abs(rf - rhat)) < 1e-12)
                bound = (p - 1) * e1.max() * e2.max()
                if bound > 0:
                    worst_ratio = max(worst_ratio, float(np.max(np.abs(rhat)) / bound))
                fl_ok &= bool(np.max(np.abs(rhat)) <= bound + 1e-15)
    ok &= fl_ok
    print("2b. Fourier formula for Rhat and the bound |Rhat| <= (p - 1) max|eps_1| max|eps_2| (float):", fl_ok,
          " worst |Rhat| / bound = %.3f" % worst_ratio)

    # 3. brute-force F1 box
    nb2 = 9
    xs, ys = x[: 1 << nb2], y[: 1 << nb2]
    box_ok = True
    for p in (3, 5, 7, 11, 17, 43):
        cnts = exact_counts(p, nb2)
        D = [[Fraction(c, 1 << nb2) for c in row] for row in cnts[nb2]]
        for s in (0, 1, 4):
            N = np.multiply.outer(xs, xs) + (1 << s) * np.multiply.outer(ys, ys) - (1 if s == 0 else 0)
            ref = np.bincount((N % p).ravel(), minlength=p)
            law = class_law_direct(D, D, p, pow(2, s, p), exact=True)
            shift = 1 if s == 0 else 0
            pred = [law[(c + shift) % p] * (1 << (2 * nb2)) for c in range(p)]
            box_ok &= all(int(pred[c]) == int(ref[c]) and pred[c].denominator == 1 for c in range(p))
    ok &= box_ok
    print("3. box law = brute-force count of F1 residues over [0, 2^9)^2 for s = 0, 1, 4 (p = 3, 5, 7, 11, 17, 43):", box_ok)

    # 4. approach to the generic law
    print("4. n1 = n2 = n, unit u = 2:  A = pi_n1 + pi_n2 - 1/p^2 (limit 1/p^2),  max|Rhat|,  E = (P(0) - P(c))/P(c)"
          " for c' != 0 averaged (limit p/(p^2-1))")
    print(f"   {'p':>3} {'limit E':>8} | " + " | ".join(f"n={n:<3d} A*p^2   max|R|*p   E     " for n in (10, 20, 40, 80)))
    for p in PRIMES:
        cnts = exact_counts(p, 80)
        cells = []
        for n in (10, 20, 40, 80):
            D = [[Fraction(c, 1 << n) for c in row] for row in cnts[n]]
            Df = [[float(v) for v in row] for row in D]
            law = class_law_direct(Df, Df, p, 2 % p)
            A = float(2 * D[0][0] - Fraction(1, p * p))
            first = [1 / p + A * ((1.0 if c == 0 else 0.0) - 1 / p) for c in range(p)]
            R = max(abs(law[c] - first[c]) for c in range(p))
            other = sum(law[1:]) / (p - 1)
            E = (law[0] - other) / other
            cells.append(f"      {A * p * p:6.3f}  {R * p:9.2e}  {E:7.4f}")
        print(f"   {p:3d} {p / (p * p - 1):8.4f} | " + " | ".join(cells))
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
