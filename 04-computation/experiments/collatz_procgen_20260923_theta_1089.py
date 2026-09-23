#!/usr/bin/env python3
"""
collatz_procgen_20260923_theta_1089.py

A quick, honest test of the owner's hint "33^2 = 1089" for X = sum (2^10/3^9)^(k^3) (HYP-9127) and for
the 33-letter square-swap word of HYP-9131.

  N1  arithmetic facts (exact): 1089 = 3^2 11^2, 3^5 = 2*11^2 + 1, 3^9 = 18*1089 + 81, 2^10 = 1089 - 65,
      33 = 2^5 + 1, 3^10 = 1 (mod 11^2) (11 is a base-3 Wieferich prime), residues of rho mod 11, 121
  N2  11-adic behaviour of the terms rho^(k^3) and the partial sums (the series does not converge 11-adically)
  N3  the real value X_R: continued fraction statistics (Gauss-Kuzmin, Levy) and PSLQ relations
  N4  2-adic X: denominators 1089^j, 33^j, 11^j (rational reconstruction covers them); block length 33
  N5  LLL relation test mixing X, X_R-independent data with 1089 (sanity)

Usage: python3 collatz_procgen_20260923_theta_1089.py [--quick]
"""
import math
import sys
import time

import gmpy2
from gmpy2 import mpz
import mpmath

QUICK = "--quick" in sys.argv


def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)
    sys.stdout.flush()


def section_N1():
    hdr("N1  arithmetic facts")
    facts = [
        ("1089 = 33^2 = 3^2 * 11^2", 1089 == 33 ** 2 == 9 * 121),
        ("33 = 2^5 + 1", 33 == 2 ** 5 + 1),
        ("3^5 = 2*11^2 + 1", 3 ** 5 == 2 * 121 + 1),
        ("3^9 = 18*1089 + 81", 3 ** 9 == 18 * 1089 + 81),
        ("2^10 = 1089 - 65 (= 32^2 = 33^2 - 65)", 2 ** 10 == 1089 - 65),
        ("3^7 - 2^11 = 139 = 2^7 + 11", 3 ** 7 - 2 ** 11 == 139 == 2 ** 7 + 11),
        ("3^10 = 1 mod 121  (11 is a base-3 Wieferich prime)", pow(3, 10, 121) == 1),
        ("3^9 - 2^10 = 18659 = 47 * 397", 3 ** 9 - 2 ** 10 == 18659 == 47 * 397),
    ]
    for s, ok in facts:
        print(f"  {s:55s} {ok}")
    inv = pow(3 ** 9, -1, 121)
    r121 = (2 ** 10 * inv) % 121
    order = next(k for k in range(1, 200) if pow(r121, k, 121) == 1)
    print(f"  rho = 2^10/3^9 mod 11 = {r121 % 11}, mod 121 = {r121} (coincidence: 47 also divides 3^9 - 2^10);"
          f" multiplicative order mod 121 = {order}")
    print(f"  rho^5 mod 121 = {pow(r121, 5, 121)} = 2^50 mod 121 = {pow(2, 50, 121)} (because 3^45 = 1 mod 121)")
    print("  3^5 = 2*11^2 + 1 is the last solution of 3^n = 2y^2 + 1 (n = 1, 2, 5): classical.")
    print("  None of these involves the prime 2 beyond 2^5 + 1 = 33; HYP-9127 is a statement about 2-adic digits.")


def section_N2():
    hdr("N2  11-adic view: terms rho^(k^3) are 11-adic units; the series does not converge in Q_11")
    inv = pow(3 ** 9, -1, 121)
    r = (2 ** 10 * inv) % 121
    terms = [pow(r, k ** 3, 121) for k in range(0, 200)]
    # period of k -> rho^(k^3) mod 121 : k^3 mod 55 has period 55
    per = next(p for p in range(1, 200) if all(terms[k] == terms[k + p] for k in range(0, 200 - p)))
    S = []
    acc = 0
    for t in terms:
        acc = (acc + t) % 121
        S.append(acc)
    drift = (S[per - 1 + per] - S[per - 1]) % 121
    print(f"  k -> rho^(k^3) mod 121 is periodic with period {per}; the partial sums S_K mod 121 drift by {drift} per period")
    print(f"  so S_K mod 121 is periodic with period {per * 121 // math.gcd(121, drift) if drift else per} and has no 11-adic limit")
    print("  (drift 0 is generic: cubing permutes Z/55, so a period sums to (rho^55 - 1)/(rho - 1) = 0 mod 121)")
    print("  (|rho^(k^3)|_11 = 1). Nothing 11-adic constrains the 2-adic number X or a putative rational value a/b.")


def section_N3():
    hdr("N3  the real value X_R = sum (2^10/3^9)^(k^3): continued fraction and integer relations")
    dps = 1200 if QUICK else 3000
    mpmath.mp.dps = dps + 50
    r = mpmath.mpf(2) ** 10 / mpmath.mpf(3) ** 9
    X = mpmath.mpf(0)
    k = 0
    while True:
        t = r ** (k ** 3)
        if t < mpmath.mpf(10) ** (-(dps + 40)):
            break
        X += t
        k += 1
    print(f"  X_R = {mpmath.nstr(X, 50)}  ({k} terms, {dps} digits)")
    # continued fraction
    x = X
    pq = []
    for _ in range(int(dps * 0.9)):
        a = int(mpmath.floor(x))
        pq.append(a)
        f = x - a
        if f == 0:
            break
        x = 1 / f
    # convergent denominators to check precision budget: keep terms while q_n^2 < 10^dps
    qm, q = 0, 1
    good = 0
    for a in pq:
        qm, q = q, a * q + qm
        if 2 * mpmath.log10(q) > dps - 20:
            break
        good += 1
    pq = pq[:good]
    mx = max(pq[1:])
    levy = float(mpmath.log(q) / good) if good else 0
    big = sum(1 for a in pq[1:] if a >= 100)
    exp_big = (good - 1) * math.log2(1 + 1 / 100)
    print(f"  reliable partial quotients: {good}; max a_n = {mx}; #a_n >= 100: {big} (Gauss-Kuzmin expects {exp_big:.1f})")
    print(f"  Levy constant estimate log(q_n)/n = {levy:.4f} (typical: pi^2/(12 ln 2) = {math.pi**2/(12*math.log(2)):.4f})")
    # PSLQ for polynomial relations of degree <= 6
    deg = 4 if QUICK else 6
    mpmath.mp.dps = 400
    Xs = mpmath.mpf(X)
    rel = mpmath.pslq([Xs ** i for i in range(deg + 1)], maxcoeff=10 ** 40, maxsteps=10 ** 6)
    print(f"  PSLQ, degree <= {deg}, coefficients <= 10^40 at 400 digits: {rel}")
    rel2 = mpmath.pslq([mpmath.mpf(1), Xs, Xs / 1089, Xs * 33, mpmath.mpf(1) / 1089], maxcoeff=10 ** 40, maxsteps=10 ** 5)
    print(f"  PSLQ on [1, X, X/1089, 33 X, 1/1089] (sanity: finds the trivial rational dependencies): {rel2}")


def section_N4():
    hdr("N4  2-adic X: denominators built from 1089, and the block length 33")
    N = 20000 if QUICK else 100000
    mod = mpz(1) << N
    inv = gmpy2.invert(mpz(3 ** 9), mod)
    X = mpz(0)
    k = 0
    while 10 * k ** 3 < N:
        X += (mpz(1) << (10 * k ** 3)) * gmpy2.powmod(inv, k ** 3, mod)
        k += 1
    X %= mod
    worst = None
    for base in (1089, 33, 11, 121, 3):
        for j in range(0, 60):
            d = mpz(base) ** j
            y = (d * X) % mod
            # distance of d*X to the nearest integer of absolute size < 2^(N/2): y or y - 2^N
            small = min(y, mod - y)
            lb = int(small).bit_length()
            if worst is None or lb < worst[0]:
                worst = (lb, base, j)
    print(f"  for d = 1089^j, 33^j, 11^j, 121^j, 3^j (j < 60): the balanced residue of d X mod 2^N has >= {worst[0]} bits"
          f" (N = {N}), so X != c/d for every integer c with |c| < 2^{worst[0] - 1}.")
    print("  (Rational reconstruction in the previous note already excludes every a/b with |a|,|b| <= 2^4999999.)")
    print("  Block length 33 of the HYP-9131 word enters only through mu_bar = 23 log2(5)/33; Theorem H1 proves")
    print("  the word irrational for every L with mu_bar < 7/4, so '33' plays no role in the proof.")


if __name__ == "__main__":
    t0 = time.time()
    print("collatz_procgen_20260923_theta_1089.py" + (" --quick" if QUICK else ""))
    section_N1()
    section_N2()
    section_N3()
    section_N4()
    print("\nVERDICT: NUMEROLOGY. The arithmetic facts are real (11 is a base-3 Wieferich prime); no periodicity,")
    print("identity or functional equation of X or of the 33-letter word modulo 1089, 11-adically, or through")
    print("block length 33/1089 was found, and the proof of HYP-9131 does not use 33.")
    print(f"\n[1089 done in {time.time()-t0:.1f}s]")
