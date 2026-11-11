#!/usr/bin/env python3
"""Auditor proof-5135, script a1: THM-5135 (1).

Independent recomputation (written before reading the session scripts) of
  - the table multiplicity N_m(t) = #{(a,b) mod m : ab = t} against sum_{d | gcd(t,m)} d phi(m/d);
  - the rank-r class law modulo p^k, r = 1..4 (and 5, 6), p^k in a list including 4, 8, 9, 27, 25,
    as the r-fold cyclic additive convolution of N_m (exact integers), against
        p^k P(c) = (1 - p^-r) sigma_{r-1}(p^v) / p^((r-1) v)          (v = v_p(c) < k)
    and the value at the class 0 (v >= k), which the theorem does not state;
  - brute force over all (u, v) for random matrices G of unit determinant, and of non-unit determinant
    (to see that the hypothesis is needed), r = 2, 3 (and r = 4 at 4 and 8... small);
  - (N_m * N_m)(c) = m^3 L_m(c), L_m(e) = sum_{d | gcd(m,e)} d^-1 prod_{p | m/d} (1 - p^-2), for all m <= 210
    (odd AND even), with one table multiplied by a unit;
  - the Euler product sigma_{r-1}(n) / (zeta(r) n^(r-1)) (numerically, partial products).
All counts are exact integers / Fractions.
"""
import sys
import itertools
import random
from fractions import Fraction
from math import gcd, pi

import numpy as np


def phi(n):
    r = n
    x = n
    p = 2
    while p * p <= x:
        if x % p == 0:
            while x % p == 0:
                x //= p
            r -= r // p
        p += 1
    if x > 1:
        r -= r // x
    return r


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def prime_factors(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


def table_brute(m):
    """N_m(t) by direct enumeration of the m x m multiplication table."""
    cnt = [0] * m
    for a in range(m):
        for b in range(m):
            cnt[(a * b) % m] += 1
    return cnt


def table_formula(m):
    out = []
    for t in range(m):
        g = gcd(t, m)  # gcd(0, m) = m
        out.append(sum(d * phi(m // d) for d in divisors(g)))
    return out


def cyc_conv(a, b):
    m = len(a)
    out = [0] * m
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[(i + j) % m] += ai * bj
    return out


def vp(c, p, k):
    """valuation of the class c modulo p^k, capped at k (class 0 -> k)."""
    c %= p ** k
    if c == 0:
        return k
    v = 0
    while c % p == 0:
        c //= p
        v += 1
    return v


def law_factor(p, r, v):
    """(1 - p^-r) sigma_{r-1}(p^v) / p^((r-1)v) as a Fraction."""
    s = sum(Fraction(1, p ** ((r - 1) * i)) for i in range(v + 1))
    return (1 - Fraction(1, p ** r)) * s


def zero_class_factor(p, r, k):
    """my own derivation for the class 0 mod p^k (v >= k):
       p^k P(0) = (1 - p^-r) sum_{i<k} p^-(r-1)i + p^-(r-1)k."""
    return law_factor(p, r, k - 1) + Fraction(1, p ** ((r - 1) * k))


def main():
    out = []
    P = out.append

    # ---------------------------------------------------------------- table multiplicity
    bad = 0
    for m in range(1, 201):
        if table_brute(m) != table_formula(m):
            bad += 1
            P("  TABLE MISMATCH at m = %d" % m)
    P("[T] N_m(t) = sum_{d | gcd(t,m)} d phi(m/d): brute force for all m in 1..200: mismatches = %d" % bad)

    # ---------------------------------------------------------------- rank-r law by convolution
    P("")
    P("[R] rank-r law modulo p^k as the r-fold additive convolution of the table (exact)")
    moduli = [(2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (3, 1), (3, 2), (3, 3), (3, 4), (5, 1), (5, 2), (5, 3),
              (7, 1), (7, 2), (11, 1), (11, 2), (13, 1)]
    allok = True
    for (p, k) in moduli:
        m = p ** k
        N1 = table_formula(m)
        cur = N1[:]
        for r in range(1, 7):
            if r > 1:
                cur = cyc_conv(cur, N1)
            tot = m ** (2 * r)
            assert sum(cur) == tot
            ok = True
            for c in range(m):
                v = vp(c, p, k)
                fac = Fraction(cur[c] * m, tot)  # p^k P(c)
                if v < k:
                    pred = law_factor(p, r, v)
                else:
                    pred = zero_class_factor(p, r, k)
                if fac != pred:
                    ok = False
            allok = allok and ok
            if r <= 4 and m in (4, 8, 9, 27, 25):
                row = []
                for v in range(k + 1):
                    c = (p ** v) % m
                    row.append("v=%d: count %d, factor %s" % (v, cur[c], Fraction(cur[c] * m, tot)))
                P("  p^k = %d, r = %d: law %s | %s" % (m, r, "OK" if ok else "FAIL", "; ".join(row)))
    P("  all (p^k, r) pairs, r = 1..6, p^k in %s: %s" % ([p ** k for p, k in moduli], "ALL OK" if allok else "SOME FAIL"))
    P("  (class 0, v >= k, compared with my formula (1-p^-r) sigma_{r-1}(p^(k-1))/p^((r-1)(k-1)) + p^-(r-1)k)")

    # ---------------------------------------------------------------- brute force with matrices G
    P("")
    P("[G] brute force of u^T G v over all (u, v), G random integer matrices")
    rng = random.Random(5135)

    def det_int(M):
        n = len(M)
        if n == 1:
            return M[0][0]
        if n == 2:
            return M[0][0] * M[1][1] - M[0][1] * M[1][0]
        s = 0
        for j in range(n):
            minor = [row[:j] + row[j + 1:] for row in M[1:]]
            s += (-1) ** j * M[0][j] * det_int(minor)
        return s

    def counts_brute(G, m):
        r = len(G)
        G = np.array(G, dtype=np.int64) % m
        grids = np.array(list(itertools.product(range(m), repeat=r)), dtype=np.int64)  # m^r x r
        Gv = (grids @ G.T) % m  # rows: G v
        cnt = np.zeros(m, dtype=np.int64)
        # u . (G v) over all u, v
        for u in grids:
            vals = (Gv @ u) % m
            cnt += np.bincount(vals, minlength=m)
        return [int(x) for x in cnt]

    def law_counts(p, k, r):
        m = p ** k
        N1 = table_formula(m)
        cur = N1[:]
        for _ in range(r - 1):
            cur = cyc_conv(cur, N1)
        return cur

    tests = [(2, 2, 2), (2, 3, 2), (3, 2, 2), (3, 3, 2), (5, 2, 2), (2, 2, 3), (2, 3, 3), (3, 2, 3), (2, 2, 4), (2, 3, 4), (3, 1, 4)]
    for (p, k, r) in tests:
        m = p ** k
        ref = law_counts(p, k, r)
        n_unit = 0
        n_unit_ok = 0
        n_non = 0
        n_non_diff = 0
        trials = 0
        while (n_unit < 6 or n_non < 6) and trials < 4000:
            trials += 1
            G = [[rng.randrange(-6, 7) for _ in range(r)] for _ in range(r)]
            d = det_int(G)
            if d % p != 0:
                if n_unit >= 6:
                    continue
                n_unit += 1
                if counts_brute(G, m) == ref:
                    n_unit_ok += 1
            else:
                if n_non >= 6:
                    continue
                n_non += 1
                if counts_brute(G, m) != ref:
                    n_non_diff += 1
        P("  p^k = %d, r = %d: unit determinant: %d/%d matrices have the law of (1); p | det: %d/%d differ from it"
          % (m, r, n_unit_ok, n_unit, n_non_diff, n_non))

    # an explicit family with p | det G where the law of (1) fails
    P("  explicit: G = diag(1, p) modulo p^2, classes by valuation (counts) vs the unimodular rank-2 law")
    for p in (2, 3, 5):
        m = p * p
        cb = counts_brute([[1, 0], [0, p]], m)
        ref = law_counts(p, 2, 2)
        P("    p = %d: diag(1,p): %s ; law (1): %s" % (p, [cb[(p ** v) % m] for v in range(3)], [ref[(p ** v) % m] for v in range(3)]))

    # ---------------------------------------------------------------- additive square = m^3 L_m
    P("")
    P("[L] (N_m * N_m)(c) = m^3 L_m(c), L_m(e) = sum_{d | gcd(m,e)} d^-1 prod_{p | m/d} (1 - p^-2)")

    def L(m, e):
        g = gcd(m, e % m) if e % m else m
        s = Fraction(0)
        for d in divisors(g):
            t = Fraction(1, d)
            for p in prime_factors(m // d):
                t *= 1 - Fraction(1, p * p)
            s += t
        return s

    bad_odd = bad_even = 0
    n_odd = n_even = 0
    bad_unit = 0
    for m in range(1, 211):
        N1 = table_formula(m)
        sq = cyc_conv(N1, N1)
        okm = all(Fraction(sq[c], m ** 3) == L(m, c) for c in range(m))
        if m % 2:
            n_odd += 1
            bad_odd += (not okm)
        else:
            n_even += 1
            bad_even += (not okm)
        # one table multiplied by a unit u: #{a1 a2 + u b1 b2 = c}
        for u in range(2, m):
            if gcd(u, m) != 1:
                continue
            Nu = [0] * m
            for t in range(m):
                Nu[(u * t) % m] += N1[t]
            if cyc_conv(N1, Nu) != sq:
                bad_unit += 1
            if u > 6:
                break
    P("  odd m in 1..209: %d moduli, mismatches %d; even m in 2..210: %d moduli, mismatches %d" % (n_odd, bad_odd, n_even, bad_even))
    P("  with one table multiplied by a unit u (u = 2..7 prime to m): mismatches %d" % bad_unit)
    # multiplicativity of L_m and the prime-power values
    badL = 0
    for (p, k) in [(3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (7, 1), (7, 2), (2, 1), (2, 2), (2, 3), (2, 4)]:
        m = p ** k
        for e in range(m):
            v = vp(e, p, k)
            if v < k:
                pred = (1 - Fraction(1, p * p)) * sum(Fraction(1, p ** i) for i in range(v + 1))
            else:
                pred = 1 + Fraction(1, p) - Fraction(1, p ** (k + 1))
            if L(m, e) != pred:
                badL += 1
    P("  L_(p^k)(e) = (1-p^-2) sigma(p^v)/p^v (v<k), 1 + 1/p - p^-(k+1) (v>=k): mismatches %d" % badL)

    # ---------------------------------------------------------------- composite moduli, general r (CRT)
    P("")
    P("[C] composite moduli: r-fold power of the table is the product of the prime-power laws (CRT)")
    badc = 0
    for m in (6, 12, 15, 18, 20, 36, 45, 60, 63, 75, 105):
        N1 = table_formula(m)
        cur = N1[:]
        pk = []
        x = m
        for p in prime_factors(m):
            k = 0
            while x % p == 0:
                x //= p
                k += 1
            pk.append((p, k))
        for r in range(1, 5):
            if r > 1:
                cur = cyc_conv(cur, N1)
            for c in range(m):
                pred = Fraction(1)
                for (p, k) in pk:
                    v = vp(c % (p ** k), p, k)
                    pred *= law_factor(p, r, v) if v < k else zero_class_factor(p, r, k)
                if Fraction(cur[c] * m, m ** (2 * r)) != pred:
                    badc += 1
    P("  m in (6, 12, 15, 18, 20, 36, 45, 60, 63, 75, 105), r = 1..4: mismatches %d" % badc)

    # ---------------------------------------------------------------- Euler product
    P("")
    P("[E] Euler product prod_p (1 - p^-r) sigma_{r-1}(p^v)/p^((r-1)v) vs sigma_{r-1}(n)/(zeta(r) n^(r-1))")
    zeta = {2: pi ** 2 / 6, 3: 1.2020569031595942, 4: pi ** 4 / 90}
    LIM = 200000
    sieve = np.ones(LIM + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(LIM ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    primes = [int(x) for x in np.nonzero(sieve)[0]]
    for r in (2, 3, 4):
        for n in (1, 2, 12, 45, 360, 1001):
            prod = 1.0
            for p in primes:
                v = 0
                x = n
                while x % p == 0:
                    x //= p
                    v += 1
                prod *= (1 - p ** (-r)) * sum(p ** (-(r - 1) * i) for i in range(v + 1))
            sig = sum(d ** (r - 1) for d in divisors(n))
            closed = sig / (zeta[r] * n ** (r - 1))
            P("  r = %d, n = %4d: partial product (p < %d) %.9f, closed form %.9f, diff %.2e" % (r, n, LIM, prod, closed, prod - closed))
    # r = 1: product over p <= x of (1 - 1/p)(v+1) tends to 0
    for n in (1, 12):
        for cut in (100, 10000, 200000):
            prod = 1.0
            for p in primes:
                if p > cut:
                    break
                v = 0
                x = n
                while x % p == 0:
                    x //= p
                    v += 1
                prod *= (1 - 1.0 / p) * (v + 1)
            P("  r = 1, n = %d: product over p <= %d = %.6f (tends to 0 like d(n) e^-gamma / ln x)" % (n, cut, prod))

    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
