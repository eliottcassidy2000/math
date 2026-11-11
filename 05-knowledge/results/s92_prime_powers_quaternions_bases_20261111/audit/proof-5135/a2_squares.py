#!/usr/bin/env python3
"""Auditor proof-5135, script a2: THM-5135 (2), (3)(a), (3)(c), (3)(d) and the counting formulas of (5).

Independent recomputation:
  [S] class counts of x_1^2 + ... + x_2r^2 modulo p^k (cyclic convolution of the distribution of squares)
      against the hyperbolic form sum_{i<=r} a_i b_i (convolution of the multiplication table), r = 1..6,
      p odd, k = 1..3 (or as far as cheap): equal iff r even or p = 1 (mod 4);
      zero counts modulo p against p^(2r-1) + eta (p^r - p^(r-1));
      a direct brute force (no convolution) for small cases.
  [2] the prime 2: four squares vs a1 a2 + b1 b2 modulo 2, 4, 8, 16, 32 (which classes differ);
      2-adic density of four / eight squares at odd classes; eight squares vs rank 4 modulo 2^k.
  [J] r_4(n) = 8 sigma(n) (n odd), general Jacobi; Hurwitz order 24 sigma(odd part);
      r_8(n) = 16 sigma_3(n) (n odd), r_8(2) = 112, general formula; additive convolution of r_2.
  [B] bicomplex counts: sum_{ab=n} r_2(a) r_2(b) on Z[i]^2, zero for n = 3 (mod 4); the count on the order
      Z[1, i, j, ij] spanned by the four units (index 4 in Z[i]^2) for comparison.
  [Q] all quadratic forms in 4 variables over F_2 and F_3: which have the class counts of a1 a2 + b1 b2
      (supports (3)(d): a form with the seam law modulo p is non-degenerate of hyperbolic type).
"""
import sys
import itertools
from fractions import Fraction
from math import gcd, isqrt

import numpy as np


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def table(m):
    cnt = [0] * m
    for a in range(m):
        for b in range(m):
            cnt[(a * b) % m] += 1
    return cnt


def squares(m):
    cnt = [0] * m
    for a in range(m):
        cnt[(a * a) % m] += 1
    return cnt


def cyc_conv(a, b):
    m = len(a)
    A = np.array(a, dtype=object)
    out = [0] * m
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[(i + j) % m] += ai * bj
    return out


def conv_power(a, r):
    cur = a[:]
    for _ in range(r - 1):
        cur = cyc_conv(cur, a)
    return cur


def main():
    out = []
    P = out.append

    # ------------------------------------------------------------------ [S]
    P("[S] sum of 2r squares vs hyperbolic form of rank r modulo p^k, p odd")
    P("    criterion claimed: equal class counts iff r even or p = 1 (mod 4)")
    viol = 0
    ncase = 0
    zero_bad = 0
    detail = []
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        for k in (1, 2, 3):
            m = p ** k
            if m > 1400:
                continue
            sq = squares(m)
            tb = table(m)
            S = None
            H = None
            for r in range(1, 7):
                S = conv_power(sq, 2) if r == 1 else cyc_conv(S, conv_power(sq, 2))
                H = tb[:] if r == 1 else cyc_conv(H, tb)
                same = (S == H)
                pred = (r % 2 == 0) or (p % 4 == 1)
                ncase += 1
                if same != pred:
                    viol += 1
                    P("    VIOLATION p^k = %d r = %d: same = %s predicted %s" % (m, r, same, pred))
                if k == 1:
                    # number of zeros modulo p
                    eta_h = 1
                    eta_s = 1 if pred else -1
                    zh = p ** (2 * r - 1) + eta_h * (p ** r - p ** (r - 1))
                    zs = p ** (2 * r - 1) + eta_s * (p ** r - p ** (r - 1))
                    if H[0] != zh or S[0] != zs:
                        zero_bad += 1
                if (not same) and r <= 3 and m in (3, 7, 9, 27, 49, 11):
                    ndiff = sum(1 for c in range(m) if S[c] != H[c])
                    detail.append("    p^k = %d, r = %d: %d of %d classes differ (e.g. class 0: squares %d, hyperbolic %d; class 1: %d, %d)"
                                  % (m, r, ndiff, m, S[0], H[0], S[1], H[1]))
    P("    cases (p in 3..41, k <= 3 with p^k <= 1400, r = 1..6): %d, violations of the criterion: %d" % (ncase, viol))
    P("    zero counts modulo p vs p^(2r-1) + eta (p^r - p^(r-1)): mismatches %d" % zero_bad)
    out.extend(detail)

    # direct brute force, no convolution: 4 squares vs a1a2+b1b2 mod 9 and 25, 2 squares vs a1 b1 mod 9, 25, 49, 13, 169
    def brute_counts(fn, m, nvar):
        cnt = np.zeros(m, dtype=np.int64)
        rng = np.arange(m, dtype=np.int64)
        grids = np.meshgrid(*([rng] * nvar), indexing="ij")
        vals = fn(*grids) % m
        cnt = np.bincount(vals.ravel(), minlength=m)
        return [int(x) for x in cnt]

    P("    direct enumeration (no convolution):")
    for m in (9, 25, 27, 49):
        a = brute_counts(lambda x, y, z, w: x * x + y * y + z * z + w * w, m, 4)
        b = brute_counts(lambda x, y, z, w: x * y + z * w, m, 4)
        P("      four squares vs a1a2+b1b2 modulo %d: %s" % (m, "EQUAL" if a == b else "DIFFERENT"))
    for m in (5, 25, 13, 169, 3, 9, 7, 49):
        a = brute_counts(lambda x, y: x * x + y * y, m, 2)
        b = brute_counts(lambda x, y: x * y, m, 2)
        P("      two squares vs a b modulo %d: %s" % (m, "EQUAL" if a == b else "DIFFERENT"))
    for m in (3, 5, 7, 9):
        a = brute_counts(lambda a1, a2, a3, a4, a5, a6: a1 * a1 + a2 * a2 + a3 * a3 + a4 * a4 + a5 * a5 + a6 * a6, m, 6)
        b = brute_counts(lambda a1, a2, a3, a4, a5, a6: a1 * a2 + a3 * a4 + a5 * a6, m, 6)
        P("      six squares vs rank-3 hyperbolic modulo %d: %s" % (m, "EQUAL" if a == b else "DIFFERENT"))

    # ------------------------------------------------------------------ [2]
    P("")
    P("[2] the prime 2")
    for k in (1, 2, 3, 4, 5, 6):
        m = 2 ** k
        sq = squares(m)
        tb = table(m)
        S4 = conv_power(sq, 4)
        H2 = conv_power(tb, 2)
        S8 = cyc_conv(S4, S4)
        H4 = cyc_conv(H2, H2)
        diff4 = [c for c in range(m) if S4[c] != H2[c]]
        diff8 = [c for c in range(m) if S8[c] != H4[c]]
        dens4 = set(Fraction(S4[c], m ** 3) for c in range(1, m, 2))
        dens8 = set(Fraction(S8[c], m ** 7) for c in range(1, m, 2))
        P("  modulo %d: four squares vs rank 2: %d of %d classes differ; eight squares vs rank 4: %d of %d differ"
          % (m, len(diff4), m, len(diff8), m))
        if m <= 8:
            P("     four squares counts %s" % S4)
            P("     a1a2+b1b2 counts   %s" % H2)
        P("     2-adic density of four squares on odd classes: %s; of eight squares on odd classes: %s"
          % (sorted(dens4), sorted(dens8)))

    # ------------------------------------------------------------------ [J]
    P("")
    P("[J] Jacobi, Hurwitz, eight squares")
    NMAX = 3000
    # r_2 by enumeration
    r2 = [0] * (NMAX + 1)
    lim = isqrt(NMAX)
    for x in range(-lim, lim + 1):
        for y in range(-lim, lim + 1):
            n = x * x + y * y
            if n <= NMAX:
                r2[n] += 1
    r2a = np.array(r2, dtype=np.int64)
    r4 = np.convolve(r2a, r2a)[:NMAX + 1]
    r8 = np.convolve(r4, r4)[:NMAX + 1]
    sig = lambda n, k=1: sum(d ** k for d in divisors(n))
    bad4 = sum(1 for n in range(1, NMAX + 1, 2) if r4[n] != 8 * sig(n))
    bad4g = sum(1 for n in range(1, NMAX + 1) if r4[n] != 8 * sum(d for d in divisors(n) if d % 4))
    bad8 = sum(1 for n in range(1, NMAX + 1, 2) if r8[n] != 16 * sig(n, 3))
    bad8g = sum(1 for n in range(1, NMAX + 1) if r8[n] != 16 * sum((-1) ** (n + d) * d ** 3 for d in divisors(n)))
    even_eq = [n for n in range(2, NMAX + 1, 2) if r8[n] == 16 * sig(n, 3)]
    P("  r_4(n) = sum_{a+b=n} r_2(a) r_2(b) (a, b >= 0, r_2(0) = 1) = 8 sigma(n), odd n <= %d: mismatches %d" % (NMAX, bad4))
    P("  r_4(n) = 8 sum_{d | n, 4 !| d} d, all n <= %d: mismatches %d" % (NMAX, bad4g))
    P("  r_8(n) = 16 sigma_3(n), odd n <= %d: mismatches %d; r_8(2) = %d, 16 sigma_3(2) = %d" % (NMAX, bad8, r8[2], 16 * sig(2, 3)))
    P("  r_8(n) = 16 sum_{d|n} (-1)^(n+d) d^3, all n <= %d: mismatches %d; even n with r_8 = 16 sigma_3: %s" % (NMAX, bad8g, even_eq))
    # Hurwitz order: Lipschitz plus the coset (Z+1/2)^4 : norm n <=> sum of four odd squares = 4n
    HN = 400
    odd_sq = [0] * (4 * HN + 1)
    o = 1
    while o * o <= 4 * HN:
        odd_sq[o * o] += 2  # +-o
        o += 2
    osq = np.array(odd_sq, dtype=np.int64)
    o2 = np.convolve(osq, osq)[:4 * HN + 1]
    o4 = np.convolve(o2, o2)[:4 * HN + 1]
    badh = 0
    for n in range(1, HN + 1):
        hur = int(r4[n]) + int(o4[4 * n])
        oddpart_sigma = sum(d for d in divisors(n) if d % 2)
        if hur != 24 * oddpart_sigma:
            badh += 1
    P("  Hurwitz order: #{norm n} = 24 * (sum of odd divisors of n), n <= %d: mismatches %d (= 24 sigma(n) for odd n)" % (HN, badh))

    # ------------------------------------------------------------------ [B]
    P("")
    P("[B] bicomplex counts")
    BN = 600
    mult = [0] * (BN + 1)
    for a in range(1, BN + 1):
        if r2[a] == 0:
            continue
        for b in range(1, BN // a + 1):
            mult[a * b] += r2[a] * r2[b]
    nz3 = [n for n in range(3, BN + 1, 4) if mult[n] != 0]
    P("  M(n) = sum_{ab=n} r_2(a) r_2(b), n <= %d: nonzero values at n = 3 (mod 4): %s" % (BN, nz3))
    P("  first values M(1..20) = %s" % mult[1:21])
    # the order Z[1, i, j, ij] (i^2 = j^2 = -1, ij = ji): w = a + b i + c j + d ij -> (z1, z2) = ((a+d) + (b-c) i, (a-d) + (b+c) i)
    # norm = |z1|^2 |z2|^2. Count by enumeration of (a, b, c, d) in a box: |z1|^2 |z2|^2 = n, both factors >= 1,
    # so |z1|^2, |z2|^2 <= n and |a +- d|, |b -+ c| <= sqrt(n).
    BN2 = 120
    R = isqrt(BN2) + 1
    basis = [0] * (BN2 + 1)
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            for c in range(-R, R + 1):
                for d in range(-R, R + 1):
                    n1 = (a + d) ** 2 + (b - c) ** 2
                    n2 = (a - d) ** 2 + (b + c) ** 2
                    n = n1 * n2
                    if 1 <= n <= BN2:
                        basis[n] += 1
    P("  order Z[1,i,j,ij] (integer coordinates in the four units; index 4 in Z[i]^2): B(n) = #{norm n}")
    P("    B(1..20) = %s" % basis[1:21])
    odd_half = all(2 * basis[n] == mult[n] for n in range(1, BN2 + 1, 2))
    P("    for odd n <= %d: B(n) = M(n)/2: %s; B(1) = %d (the 8 units +-1, +-i, +-j, +-ij), M(1) = %d (the 16 units of Z[i]^2)"
      % (BN2, odd_half, basis[1], mult[1]))
    ev = [(n, basis[n], mult[n]) for n in range(2, 34, 2)]
    P("    even n: (n, B(n), M(n)) = %s" % ev)
    # the caveat on the box: both |a+d|,|a-d| <= sqrt(n) so |a|,|d| <= sqrt(n); same for b, c. R suffices.

    # ------------------------------------------------------------------ [Q]
    P("")
    P("[Q] quaternary quadratic forms over F_p with the class counts of a1 a2 + b1 b2 (p = 2, 3)")
    for p in (2, 3):
        vecs = np.array(list(itertools.product(range(p), repeat=4)), dtype=np.int64)
        monos = [(i, j) for i in range(4) for j in range(i, 4)]
        mono_vals = np.stack([vecs[:, i] * vecs[:, j] for (i, j) in monos], axis=1) % p  # p^4 x 10
        ref = np.bincount((vecs[:, 0] * vecs[:, 1] + vecs[:, 2] * vecs[:, 3]) % p, minlength=p)
        n_same = 0
        n_same_nondeg = 0
        n_total = 0
        zero_only_same = 0
        for coeffs in itertools.product(range(p), repeat=10):
            n_total += 1
            cf = np.array(coeffs, dtype=np.int64)
            vals = (mono_vals @ cf) % p
            cnt = np.bincount(vals, minlength=p)
            if cnt[0] == ref[0]:
                zero_only_same += 1
            if (cnt == ref).all():
                n_same += 1
                # polar form b(x,y) = Q(x+y) - Q(x) - Q(y): Gram matrix with 2 c_ii on the diagonal, c_ij off it
                G = np.zeros((4, 4), dtype=np.int64)
                for (i, j), c in zip(monos, coeffs):
                    if i == j:
                        G[i, i] = 2 * c
                    else:
                        G[i, j] = c
                        G[j, i] = c
                det = int(round(np.linalg.det(G.astype(float)))) % p
                if det != 0:
                    n_same_nondeg += 1
                    if p == 3:
                        assert det % 3 == 1  # square discriminant: det of the polar form of a1a2+b1b2 is 1
        P("  p = %d: %d forms; same number of zeros as the seam: %d; same class counts: %d, of which with non-degenerate polar form: %d"
          % (p, n_total, zero_only_same, n_same, n_same_nondeg))
    P("  (so over F_2 and F_3 the seam counts force a non-degenerate polar form, i.e. an even unimodular form at p)")

    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
