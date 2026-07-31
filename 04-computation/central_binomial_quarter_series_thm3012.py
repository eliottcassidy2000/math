#!/usr/bin/env python3
"""Referee for THM-3012: the series S(k) = sum C(2n,n)C(4n,2n)/((kn+1)64^n).

Checks (exact where possible, else high-precision mpmath):
  C1  C(2n,n)C(4n,2n)/64^n = (1/4)_n(3/4)_n/(n!)^2  as exact rationals, n<60.
  C2  S(k) = 3F2(1/4,3/4,1/k; 1,1+1/k; 1) and the three known closed forms
      S(1)=8sqrt2/(3pi), S(2)=(4/pi)log(1+sqrt2),
      S(3)=-(1/pi)(sqrt3 log(5-2sqrt6)+2arctan(sqrt2/5)), to 200 digits.
  C3  the quadratic transformation (b-a=1/2) and its elementary kernel
      F(z) = (1/pi) int_0^pi dphi/sqrt(1+sqrt z cos phi).
  C4  the K-moment form.
  C5  THOMAE NORMAL FORM  S(k) = (8 sqrt2/(3 pi k)) 3F2(1-1/k,1,1; 5/4,7/4; 1)
      -- only ONE parameter carries k.
  C6  the uniform integral representation
      S(k) = (2/pi) int_0^{pi/2} 2F1(1, 4/k; 1+2/k; (1-sqrt2 sin th)/2) dth.
  C7  SCHWARZ ANGLE SUM of 2F1(1,4/k;1+2/k;.) equals 8/k-1 for k<=3 (>1,
      spherical/finite monodromy) and EXACTLY 1 for every k>=4 (Euclidean):
      k=4 is the wall.
  C8  arctan(sqrt2/5) = pi - 3 arctan(sqrt2), giving the reduced form
      pi S(3) = 2 sqrt3 log(sqrt3+sqrt2) - 2 pi + 6 arctan(sqrt2).
  C9  PSLQ exclusion battery for k = 4,5,6,12 with a k=3 positive control.
"""

from fractions import Fraction as Fr
from math import comb, factorial

from mpmath import (mp, mpf, hyp2f1, hyp3f2, quad, pi, sqrt, log, atan, sin,
                    ellipk, nstr, pslq, catalan, gamma)


def require(c, m):
    if not c:
        raise RuntimeError(m)


def poch(a, n):
    r = Fr(1)
    for i in range(n):
        r *= (a + i)
    return r


def c1():
    for n in range(60):
        lhs = Fr(comb(2 * n, n) * comb(4 * n, 2 * n), 64 ** n)
        rhs = poch(Fr(1, 4), n) * poch(Fr(3, 4), n) / Fr(factorial(n)) ** 2
        require(lhs == rhs, f"C1 fails at n={n}")
    print("C1 C(2n,n)C(4n,2n)/64^n = (1/4)_n(3/4)_n/(n!)^2 exactly, n<60: OK")


def S(k):
    return hyp3f2(mpf(1) / 4, mpf(3) / 4, mpf(1) / k, 1, 1 + mpf(1) / k, 1)


def c2():
    mp.dps = 220
    forms = {
        1: 8 * sqrt(2) / (3 * pi),
        2: 4 / pi * log(1 + sqrt(2)),
        3: -(sqrt(3) * log(5 - 2 * sqrt(6)) + 2 * atan(sqrt(2) / 5)) / pi,
    }
    for k, v in forms.items():
        require(abs(S(k) - v) < mpf(10) ** -195, f"C2 closed form fails k={k}")
    print("C2 S(k)=3F2(1/4,3/4,1/k;1,1+1/k;1); k=1,2,3 closed forms to 195 digits: OK")
    mp.dps = 60


def c3():
    for z in [mpf("0.1"), mpf("0.37"), mpf("0.8")]:
        F = hyp2f1(mpf(1) / 4, mpf(3) / 4, 1, z)
        q = (1 + sqrt(z)) ** mpf(-0.5) * hyp2f1(mpf(1) / 2, mpf(1) / 2, 1,
                                                2 * sqrt(z) / (1 + sqrt(z)))
        ker = 1 / pi * quad(lambda ph: 1 / sqrt(1 + sqrt(z) * __import__(
            "mpmath").cos(ph)), [0, pi])
        require(abs(F - q) < mpf(10) ** -40, "C3 quadratic transformation")
        require(abs(F - ker) < mpf(10) ** -30, "C3 elementary kernel")
    print("C3 quadratic transformation + elementary kernel "
          "F(z)=(1/pi)int_0^pi dphi/sqrt(1+sqrt z cos phi): OK")


def c4():
    for k in [1, 2, 3, 4, 5, 6]:
        v = 16 / (k * pi * sqrt(2)) * quad(
            lambda m: m ** (mpf(4) / k - 1) * (2 - m ** 2) ** (-mpf(2) / k - mpf(1) / 2)
            * ellipk(m ** 2), [0, 1])
        require(abs(v - S(k)) < mpf(10) ** -20, f"C4 K-moment fails k={k}")
    print("C4 K-moment S(k)=(16/(k pi sqrt2))int mu^{4/k-1}(2-mu^2)^{-2/k-1/2}K: OK")


def c5():
    for k in [1, 2, 3, 4, 5, 6, 8, 12]:
        v = 8 * sqrt(2) / (3 * pi * k) * hyp3f2(1 - mpf(1) / k, 1, 1,
                                                mpf(5) / 4, mpf(7) / 4, 1)
        require(abs(v - S(k)) < mpf(10) ** -22, f"C5 Thomae form fails k={k}")
    print("C5 Thomae normal form S(k)=(8sqrt2/(3 pi k))3F2(1-1/k,1,1;5/4,7/4;1): OK")


def c6():
    for k in [1, 2, 3, 4, 5, 6]:
        v = 2 / pi * quad(lambda th: hyp2f1(1, mpf(4) / k, 1 + mpf(2) / k,
                                            (1 - sqrt(2) * sin(th)) / 2), [0, pi / 2])
        require(abs(v - S(k)) < mpf(10) ** -15, f"C6 uniform rep fails k={k}")
    print("C6 uniform rep S(k)=(2/pi)int_0^{pi/2}2F1(1,4/k;1+2/k;(1-sqrt2 sin th)/2): OK")


def c7():
    for k in range(1, 40):
        s = Fr(2, k) + Fr(2, k) + abs(Fr(1) - Fr(4, k))
        if k <= 3:
            require(s == Fr(8, k) - 1 and s > 1, f"C7 spherical fails k={k}")
        else:
            require(s == 1, f"C7 Euclidean fails k={k}")
    print("C7 Schwarz angle sum of 2F1(1,4/k;1+2/k;.) = 8/k-1 > 1 (spherical) for "
          "k<=3, EXACTLY 1 (Euclidean) for all k>=4: OK  <- k=4 is the wall")


def c8():
    require(abs(atan(sqrt(2) / 5) - (pi - 3 * atan(sqrt(2)))) < mpf(10) ** -40,
            "C8 arctan identity")
    require(abs(pi * S(3) - (2 * sqrt(3) * log(sqrt(3) + sqrt(2)) - 2 * pi
                             + 6 * atan(sqrt(2)))) < mpf(10) ** -40, "C8 reduced S(3)")
    print("C8 arctan(sqrt2/5)=pi-3arctan(sqrt2); "
          "pi S(3)=2sqrt3 log(sqrt3+sqrt2)-2pi+6arctan(sqrt2): OK")


def c9():
    import itertools
    mp.dps = 150
    s2, s3, s6 = sqrt(2), sqrt(3), sqrt(6)
    C = {"pi": pi, "log(5+2r6)": log(5 + 2 * s6), "r3log(5+2r6)": s3 * log(5 + 2 * s6),
         "atan(r2/5)": atan(s2 / 5), "log(1+r2)": log(1 + s2),
         "r3log(1+r2)": s3 * log(1 + s2), "log(2+r3)": log(2 + s3),
         "atan(r2)": atan(s2), "atan(r2/3)": atan(s2 / 3), "log2": log(2),
         "r3log2": s3 * log(2), "1": mpf(1), "G": catalan,
         "varpi": gamma(mpf(1) / 4) ** 2 / (2 * sqrt(2 * pi))}
    names = list(C)

    def sweep(T):
        hits = []
        for size in (2, 3):
            for combo in itertools.combinations(names, size):
                r = pslq([T] + [C[n] for n in combo], tol=mpf(10) ** -120,
                         maxcoeff=800, maxsteps=50000)
                if r and r[0] != 0:
                    hits.append((combo, r))
        return hits
    ctrl = sweep(pi * S(3))
    require(ctrl, "C9 POSITIVE CONTROL FAILED: sweep must recover k=3")
    print(f"C9 positive control: k=3 recovered ({len(ctrl)} relations), e.g. "
          f"{ctrl[0][1]} on {ctrl[0][0]}")
    for k in (4, 5, 6, 12):
        require(not sweep(pi * S(k)), f"C9 unexpected relation at k={k}")
    print("C9 no relation for k=4,5,6,12 in the same pool (search validated): OK")
    mp.dps = 60


if __name__ == "__main__":
    mp.dps = 60
    c1(); c2(); c3(); c4(); c5(); c6(); c7(); c8(); c9()
    print("ALL THM-3012 REFEREE CHECKS PASSED")
