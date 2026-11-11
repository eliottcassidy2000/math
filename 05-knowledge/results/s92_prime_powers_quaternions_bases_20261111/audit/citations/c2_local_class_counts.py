#!/usr/bin/env python3
"""Auditor 'citations' - local class counts behind the classical statements of THM-5134 (F) and THM-5135 (1)-(3).

L1  Multiplication table modulo m: N_m(t) = #{(a, b) : ab = t} = sum_{d | gcd(t, m)} d phi(m/d)   (all m <= 200).
L2  Rank-r law at prime powers, ANY prime (2 included): r-fold additive power of the table against
        p^k P = (1 - p^-r) sigma_(r-1)(p^v) / p^((r-1) v)          (v = v_p(c) < k),
    and the value at v >= k.
L3  (N_m * N_m)(c) = m^3 L_m(c) for every odd m <= 201; also with one table multiplied by a unit.
L4  Sum of 2r squares against the hyperbolic form of rank r modulo p^k: equal class counts iff r even or p = 1 (mod 4);
    zeros modulo p: p^(2r-1) + eta (p^r - p^(r-1)).
L5  The prime 2: four squares against a1 a2 + b1 b2 modulo 2, 4, 8, 16, 32; 2-adic density of four squares at odd n;
    eight squares at odd n; the seam factors (1 - 2^-r) sigma_(r-1)(2^v)/2^((r-1)v).
L6  E8 form x^T C x / 2 modulo 2, 4, 8 (full enumeration, numba), modulo 16 (full enumeration, 2^32 vectors, parallel),
    modulo 3, 9, 5, 7: against the rank-4 law.
Exact integer arithmetic.
"""
import math
import sys
import time
from fractions import Fraction
import numpy as np
from numba import njit, prange


def phi(n):
    r = n
    m = n
    p = 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def table_counts(m):
    """brute-force multiplication table counts modulo m (python ints)."""
    t = [0] * m
    for a in range(m):
        for b in range(m):
            t[(a * b) % m] += 1
    return t


def cconv(a, b, m):
    out = [0] * m
    for i in range(m):
        ai = a[i]
        if ai == 0:
            continue
        for j in range(m):
            out[(i + j) % m] += ai * b[j]
    return out


def vp(c, p, k):
    """valuation of the class c modulo p^k, capped at k."""
    if c % (p ** k) == 0:
        return k
    v = 0
    while c % p == 0:
        c //= p
        v += 1
    return v


def sigma_pow(p, v, e):
    return sum(p ** (e * i) for i in range(v + 1))


def rank_law(p, k, r, c):
    """exact p^k * P(B = c) for the unimodular seam of rank r modulo p^k (Fraction), from the theorem's formula for
    v < k, and from my own derivation for v >= k:  sum_{i<k} p^(-(r-1)i)(1 - p^-r) + p^(k - rk)."""
    v = vp(c, p, k)
    if v < k:
        return (1 - Fraction(1, p ** r)) * Fraction(sigma_pow(p, v, r - 1), p ** ((r - 1) * v))
    return (1 - Fraction(1, p ** r)) * sum(Fraction(1, p ** ((r - 1) * i)) for i in range(k)) + Fraction(p ** k, p ** (r * k))


def L_m(m, e):
    """L_m(e) = sum_{d | gcd(m, e)} d^-1 prod_{p | m/d} (1 - p^-2)."""
    g = math.gcd(m, e % m) if e % m else m
    tot = Fraction(0)
    for d in divisors(g):
        q = m // d
        f = Fraction(1, d)
        pp = 2
        qq = q
        while qq > 1:
            if qq % pp == 0:
                f *= (1 - Fraction(1, pp * pp))
                while qq % pp == 0:
                    qq //= pp
            pp += 1
        tot += f
    return tot


def squares_counts(m, nsq):
    s1 = [0] * m
    for a in range(m):
        s1[(a * a) % m] += 1
    out = [0] * m
    out[0] = 1
    for _ in range(nsq):
        out = cconv(out, s1, m)
    return out


def hyperbolic_counts(m, r, unit=1):
    t = table_counts(m)
    out = [0] * m
    out[0] = 1
    for i in range(r):
        if i == 1 and unit != 1:
            tu = [0] * m
            for c in range(m):
                tu[(unit * c) % m] += t[c]
            out = cconv(out, tu, m)
        else:
            out = cconv(out, t, m)
    return out


@njit(parallel=True, cache=False)
def e8_class_counts(C, m):
    """class counts of x^T C x / 2 modulo m over all x in (Z/m)^8, in lattice coordinates (C the even Cartan matrix).
    x^T C x / 2 = sum_i x_i^2 + sum_{i<j} C_ij x_i x_j  (C_ii = 2)."""
    res = np.zeros((m, m), dtype=np.int64)  # first index: x_0 (for prange), second: class
    for x0 in prange(m):
        loc = np.zeros(m, dtype=np.int64)
        x = np.zeros(8, dtype=np.int64)
        x[0] = x0
        for x1 in range(m):
            x[1] = x1
            for x2 in range(m):
                x[2] = x2
                for x3 in range(m):
                    x[3] = x3
                    # partial form in x0..x3
                    q4 = 0
                    for i in range(4):
                        q4 += x[i] * x[i]
                        for j in range(i + 1, 4):
                            q4 += C[i, j] * x[i] * x[j]
                    # linear coefficients of x4..x7 from x0..x3
                    l4 = 0
                    l5 = 0
                    l6 = 0
                    l7 = 0
                    for i in range(4):
                        l4 += C[i, 4] * x[i]
                        l5 += C[i, 5] * x[i]
                        l6 += C[i, 6] * x[i]
                        l7 += C[i, 7] * x[i]
                    for x4 in range(m):
                        a4 = q4 + x4 * x4 + l4 * x4
                        for x5 in range(m):
                            a5 = a4 + x5 * x5 + (l5 + C[4, 5] * x4) * x5
                            for x6 in range(m):
                                a6 = a5 + x6 * x6 + (l6 + C[4, 6] * x4 + C[5, 6] * x5) * x6
                                c7 = l7 + C[4, 7] * x4 + C[5, 7] * x5 + C[6, 7] * x6
                                for x7 in range(m):
                                    v = a6 + x7 * x7 + c7 * x7
                                    loc[v % m] += 1
        for c in range(m):
            res[x0, c] = loc[c]
    out = np.zeros(m, dtype=np.int64)
    for x0 in range(m):
        for c in range(m):
            out[c] += res[x0, c]
    return out


def main():
    ok = True
    args = sys.argv[1:]
    do16 = "--mod16" in args

    # ---------- L1 ----------
    bad = 0
    for m in range(1, 201):
        t = table_counts(m)
        for c in range(m):
            g = math.gcd(c, m) if c else m
            f = sum(d * phi(m // d) for d in divisors(g))
            if t[c] != f:
                bad += 1
    print("L1 N_m(t) = sum_{d | gcd(t, m)} d phi(m/d), all m <= 200, all classes: failures %d" % bad)
    ok &= bad == 0

    # ---------- L2 ----------
    print("L2 rank-r law at prime powers (r-fold additive power of the table), exact:")
    cases = [(2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (3, 1), (3, 2), (3, 3), (3, 4), (5, 1), (5, 2), (5, 3), (7, 1), (7, 2),
             (11, 1), (11, 2), (13, 1)]
    bad = 0
    for p, k in cases:
        m = p ** k
        t = table_counts(m)
        cur = [0] * m
        cur[0] = 1
        for r in range(1, 6):
            cur = cconv(cur, t, m)
            tot = m ** (2 * r)
            for c in range(m):
                lhs = Fraction(cur[c] * m, tot)
                if lhs != rank_law(p, k, r, c):
                    bad += 1
    print("   primes 2, 3, 5, 7, 11, 13, ranks r = 1..5, %d (p, k) cases: failures %d" % (len(cases), bad))
    ok &= bad == 0
    for p, k in [(2, 3), (3, 2)]:
        m = p ** k
        for r in (2, 4):
            vals = sorted({(vp(c, p, k), str(rank_law(p, k, r, c))) for c in range(m)})
            print("   p^k = %d, r = %d: (valuation, p^k P) = %s" % (m, r, vals))

    # ---------- L3 ----------
    bad = 0
    nm = 0
    for m in range(3, 202, 2):
        t = table_counts(m)
        sq = cconv(t, t, m)
        nm += 1
        for c in range(m):
            if Fraction(sq[c], m ** 3) != L_m(m, c):
                bad += 1
    print("L3 (N_m * N_m)(c) = m^3 L_m(c), every odd 3 <= m <= 201 (%d moduli): failures %d" % (nm, bad))
    ok &= bad == 0
    bad = 0
    for m in (9, 15, 21, 25, 27, 35, 45, 49, 63, 75, 81, 105):
        t = table_counts(m)
        for u in range(2, m):
            if math.gcd(u, m) != 1:
                continue
            tu = [0] * m
            for c in range(m):
                tu[(u * c) % m] += t[c]
            sq = cconv(t, tu, m)
            for c in range(m):
                if Fraction(sq[c], m ** 3) != L_m(m, c):
                    bad += 1
    print("   with one table multiplied by any unit u (12 odd moduli, all units): failures %d" % bad)
    ok &= bad == 0
    # L_m for an EVEN modulus is not the table square (the theorem states odd m only): show the 2-part
    for m in (2, 4, 8, 6, 12):
        t = table_counts(m)
        sq = cconv(t, t, m)
        print("   m = %2d (even, outside the statement): (N*N)(c)/m^3 = %s ; L_m(c) = %s"
              % (m, [str(Fraction(sq[c], m ** 3)) for c in range(m)], [str(L_m(m, c)) for c in range(m)]))

    # ---------- L4 ----------
    print("L4 sum of 2r squares against the hyperbolic form of rank r:")
    bad_iff = 0
    bad_zero = 0
    rows = []
    for p in (3, 5, 7, 11, 13, 17, 19):
        for r in (1, 2, 3, 4):
            sq = squares_counts(p, 2 * r)
            hy = hyperbolic_counts(p, r)
            equal = sq == hy
            expect = (r % 2 == 0) or (p % 4 == 1)
            if equal != expect:
                bad_iff += 1
            eta_sq = 1 if expect else -1
            if sq[0] != p ** (2 * r - 1) + eta_sq * (p ** r - p ** (r - 1)):
                bad_zero += 1
            if hy[0] != p ** (2 * r - 1) + (p ** r - p ** (r - 1)):
                bad_zero += 1
            rows.append((p, r, equal))
    print("   modulo p in {3,5,7,11,13,17,19}, r = 1..4: 'equal iff r even or p = 1 (4)' failures %d; zero-count formula "
          "failures %d" % (bad_iff, bad_zero))
    ok &= (bad_iff == 0 and bad_zero == 0)
    bad = 0
    for m, p in ((9, 3), (27, 3), (81, 3), (25, 5), (125, 5), (49, 7), (121, 11), (169, 13)):
        for r in (1, 2, 3, 4):
            if m > 81 and r > 2:
                continue
            sq = squares_counts(m, 2 * r)
            hy = hyperbolic_counts(m, r)
            equal = sq == hy
            expect = (r % 2 == 0) or (p % 4 == 1)
            if equal != expect:
                bad += 1
    print("   modulo 9, 27, 81, 25, 125, 49, 121, 169 (r <= 4, r <= 2 above 81): failures %d" % bad)
    ok &= bad == 0
    bad = 0
    for m in (9, 27, 25, 49, 121):
        four = squares_counts(m, 4)
        for u in range(1, m):
            if math.gcd(u, m) != 1:
                continue
            if hyperbolic_counts(m, 2, unit=u) != four:
                bad += 1
    print("   four squares = a1 a2 + u b1 b2 for every unit u, modulo 9, 27, 25, 49, 121: failures %d" % bad)
    ok &= bad == 0
    # the odd-p local density of four squares at n: (1 - p^-2) sigma(p^v)/p^v
    bad = 0
    for p, k in ((3, 4), (5, 3), (7, 2), (11, 2)):
        m = p ** k
        four = squares_counts(m, 4)
        for c in range(m):
            if Fraction(four[c] * m, m ** 4) != rank_law(p, k, 2, c):
                bad += 1
    print("   four squares modulo 81, 125, 49, 121: m * share = (1 - p^-2) sigma(p^v)/p^v (v < k): failures %d" % bad)
    ok &= bad == 0

    # ---------- L5 ----------
    print("L5 the prime 2:")
    for k in (1, 2, 3, 4, 5):
        m = 2 ** k
        four = squares_counts(m, 4)
        seam = hyperbolic_counts(m, 2)
        dens4 = [str(Fraction(four[c] * m, m ** 4)) for c in range(m)]
        denss = [str(Fraction(seam[c] * m, m ** 4)) for c in range(m)]
        odd_ok = all(Fraction(four[c] * m, m ** 4) == 1 for c in range(1, m, 2))
        print("   modulo %2d: four squares differ from a1 a2 + b1 b2: %s ; four-square factor at every odd class = 1: %s"
              % (m, four != seam, odd_ok))
        if k <= 3:
            print("      four squares, m * share by class:", dens4)
            print("      a1a2 + b1b2 , m * share by class:", denss)
        ok &= (four != seam) and odd_ok
    # even classes of the four-square form: 2-adic density 3/2^v? (Jacobi: r_4 = 24 sigma(odd part) for even n)
    m = 64
    four = squares_counts(m, 4)
    row = {}
    for c in range(1, m):
        v = vp(c, 2, 6)
        row.setdefault(v, set()).add(Fraction(four[c] * m, m ** 4))
    print("   four squares modulo 64, m * share by 2-adic valuation of the class:",
          {v: sorted(str(x) for x in s) for v, s in sorted(row.items())})
    for k in (1, 2, 3, 4):
        m = 2 ** k
        eight = squares_counts(m, 8)
        seam4 = hyperbolic_counts(m, 4)
        odd_ok = all(Fraction(eight[c] * m, m ** 8) == 1 for c in range(1, m, 2))
        print("   modulo %2d: eight squares differ from the rank-4 seam: %s ; eight-square factor at odd classes = 1: %s"
              % (m, eight != seam4, odd_ok))

    # ---------- L6 ----------
    print("L6 E8 form x^T C x / 2 against the rank-4 seam:")
    C = np.zeros((8, 8), dtype=np.int64)
    for i in range(8):
        C[i, i] = 2
    for a, b in [(1, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (2, 4)]:
        C[a - 1, b - 1] = -1
        C[b - 1, a - 1] = -1
    mods = [2, 4, 8, 3, 9, 5, 7]
    if do16:
        mods.append(16)
    for m in mods:
        t0 = time.time()
        cnt = e8_class_counts(C, m)
        assert int(cnt.sum()) == m ** 8
        seam4 = hyperbolic_counts(m, 4)
        same = [int(x) for x in cnt] == seam4
        p = 2 if m % 2 == 0 else (3 if m % 3 == 0 else m)
        k = round(math.log(m) / math.log(p))
        law_ok = all(Fraction(int(cnt[c]) * m, m ** 8) == rank_law(p, k, 4, c) for c in range(m))
        print("   modulo %2d: E8 class counts = rank-4 seam counts: %s ; = formula (1 - p^-4) sigma_3(p^v)/p^(3v): %s  [%.1fs]"
              % (m, same, law_ok, time.time() - t0))
        if m in (2, 4, 8):
            print("      counts:", [int(x) for x in cnt])
        ok &= same and law_ok
        sys.stdout.flush()
    # Z^8 at 2 for contrast
    for m in (2, 4, 8):
        eight = squares_counts(m, 8)
        print("   modulo %d: Z^8 (eight squares) counts: %s" % (m, eight))

    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")


if __name__ == "__main__":
    main()
