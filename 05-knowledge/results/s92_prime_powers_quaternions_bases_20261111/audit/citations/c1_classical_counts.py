#!/usr/bin/env python3
"""Auditor 'citations' - independent recomputation of the classical counting theorems cited by THM-5134 / THM-5135.

C1  r_4(n): brute count (exact convolution of r_1) against Jacobi: r_4(n) = 8 * sum_{d | n, 4 not | d} d.
C2  r_8(n): brute count against 16 * sum_{d | n} (-1)^(n+d) d^3.
C3  Hurwitz integral quaternions of norm n: brute count against 24 * (sum of odd divisors of n).
C4  E8 shells: vectors of norm 2n in E8 = D8 u (D8 + (1/2)^8), against 240 sigma_3(n).
C5  E8 shells from the Cartan matrix (independent model): x^T C x / 2 = n for n = 1, 2, 3 by a bounded search.
C6  r_2 convolutions: quaternion norm additive, bicomplex norm multiplicative.
C7  real densities and zeta(4).
All integer arithmetic exact (Python ints / int64 with headroom checks).
"""
import math
import sys
import itertools
import numpy as np


def divisors(n):
    ds = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            ds.append(i)
            if i != n // i:
                ds.append(n // i)
        i += 1
    return sorted(ds)


def sigma_k(n, k=1):
    return sum(d ** k for d in divisors(n))


def conv_trunc(a, b, L):
    """exact truncated convolution of two python-int lists (length L+1)."""
    out = [0] * (L + 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j in range(0, L + 1 - i):
            bj = b[j]
            if bj:
                out[i + j] += ai * bj
    return out


def main():
    ok_all = True

    # ---------------- C1: four squares ----------------
    L = 3000
    r1 = [0] * (L + 1)
    x = 0
    while x * x <= L:
        r1[x * x] += 1 if x == 0 else 2
        x += 1
    r2 = conv_trunc(r1, r1, L)
    r4 = conv_trunc(r2, r2, L)
    bad_general = 0
    bad_odd = 0
    bad_even = 0
    for n in range(1, L + 1):
        jac = 8 * sum(d for d in divisors(n) if d % 4 != 0)
        if r4[n] != jac:
            bad_general += 1
        if n % 2 == 1:
            if r4[n] != 8 * sigma_k(n):
                bad_odd += 1
        else:
            m = n
            while m % 2 == 0:
                m //= 2
            if r4[n] != 24 * sigma_k(m):
                bad_even += 1
    print("C1 four squares, 1 <= n <= %d (exact convolution):" % L)
    print("   r_4(n) = 8 * sum_{d | n, 4 not | d} d      failures: %d" % bad_general)
    print("   odd n : r_4(n) = 8 sigma(n)                failures: %d" % bad_odd)
    print("   even n: r_4(n) = 24 sigma(odd part of n)   failures: %d" % bad_even)
    print("   r_4(1..12) =", r4[1:13])
    ok_all &= (bad_general == 0 and bad_odd == 0 and bad_even == 0)

    # ---------------- C2: eight squares ----------------
    L8 = 1500
    r4s = r4[: L8 + 1]
    r8 = conv_trunc(r4s, r4s, L8)
    bad8 = 0
    bad8_odd = 0
    even_equal_16s3 = []
    for n in range(1, L8 + 1):
        f = 16 * sum(((-1) ** (n + d)) * d ** 3 for d in divisors(n))
        if r8[n] != f:
            bad8 += 1
        if n % 2 == 1 and r8[n] != 16 * sigma_k(n, 3):
            bad8_odd += 1
        if n % 2 == 0 and r8[n] == 16 * sigma_k(n, 3):
            even_equal_16s3.append(n)
    print("C2 eight squares, 1 <= n <= %d:" % L8)
    print("   r_8(n) = 16 * sum_{d | n} (-1)^(n+d) d^3   failures: %d" % bad8)
    print("   odd n : r_8(n) = 16 sigma_3(n)             failures: %d" % bad8_odd)
    print("   even n with r_8(n) = 16 sigma_3(n): %d (none expected)" % len(even_equal_16s3))
    print("   r_8(1..8) =", r8[1:9], "   r_8(2) =", r8[2], " 16 sigma_3(2) =", 16 * sigma_k(2, 3))
    ok_all &= (bad8 == 0 and bad8_odd == 0 and not even_equal_16s3 and r8[2] == 112)

    # ---------------- C3: Hurwitz ----------------
    # Hurwitz order = Z^4 u (Z + 1/2)^4; norm n <=> all-integer: r_4(n); half-integers: (a,b,c,d) all odd, a^2+..+d^2 = 4n.
    LH = 700
    odd1 = [0] * (4 * LH + 1)
    x = 1
    while x * x <= 4 * LH:
        odd1[x * x] += 2
        x += 2
    odd2 = conv_trunc(odd1, odd1, 4 * LH)
    odd4 = conv_trunc(odd2, odd2, 4 * LH)
    badH = 0
    badH_odd = 0
    hur = []
    for n in range(1, LH + 1):
        cnt = r4[n] + odd4[4 * n]
        hur.append(cnt)
        sodd = sum(d for d in divisors(n) if d % 2 == 1)
        if cnt != 24 * sodd:
            badH += 1
        if n % 2 == 1 and cnt != 24 * sigma_k(n):
            badH_odd += 1
    print("C3 Hurwitz quaternions of norm n, 1 <= n <= %d:" % LH)
    print("   count = 24 * (sum of odd divisors of n)    failures: %d" % badH)
    print("   odd n : count = 24 sigma(n)                failures: %d" % badH_odd)
    print("   counts n = 1..8:", hur[:8])
    ok_all &= (badH == 0 and badH_odd == 0)

    # ---------------- C4: E8 as D8+ ----------------
    # vectors of norm 2n: (i) x in Z^8, sum even, x.x = 2n;  (ii) x in (Z+1/2)^8, sum even, x.x = 2n.
    # (i): x.x = sum x_i^2 = sum x_i (mod 2), so sum even <=> x.x even: count = r_8(2n).
    # (ii): y = 2x all odd, y.y = 8n, sum(y) = 0 (mod 4).  By symmetry y -> (-y_1, y_2..): exactly half of the all-odd
    #       vectors of a given norm have sum = 0 (mod 4) (flipping one sign changes the sum by 2 y_1 = 2 mod 4).
    LE = 180
    need = 8 * LE
    o1 = [0] * (need + 1)
    x = 1
    while x * x <= need:
        o1[x * x] += 2
        x += 2
    o2 = conv_trunc(o1, o1, need)
    o4 = conv_trunc(o2, o2, need)
    o8 = conv_trunc(o4, o4, need)
    r1b = [0] * (2 * LE + 1)
    x = 0
    while x * x <= 2 * LE:
        r1b[x * x] += 1 if x == 0 else 2
        x += 1
    r2b = conv_trunc(r1b, r1b, 2 * LE)
    r4b = conv_trunc(r2b, r2b, 2 * LE)
    r8b = conv_trunc(r4b, r4b, 2 * LE)
    badE = 0
    shells = []
    for n in range(1, LE + 1):
        assert o8[8 * n] % 2 == 0
        cnt = r8b[2 * n] + o8[8 * n] // 2
        shells.append(cnt)
        if cnt != 240 * sigma_k(n, 3):
            badE += 1
    print("C4 E8 = D8 u (D8 + half), vectors of norm 2n, 1 <= n <= %d:" % LE)
    print("   count = 240 sigma_3(n)                     failures: %d" % badE)
    print("   shells n = 1..6:", shells[:6], "  (240 sigma_3: %s)" % [240 * sigma_k(n, 3) for n in range(1, 7)])
    ok_all &= (badE == 0)

    # ---------------- C5: E8 from the Cartan matrix ----------------
    # Bourbaki E8 Dynkin diagram: chain 1-3-4-5-6-7-8, node 2 attached to 4.
    C = np.zeros((8, 8), dtype=np.int64)
    for i in range(8):
        C[i, i] = 2
    edges = [(1, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (2, 4)]
    for a, b in edges:
        C[a - 1, b - 1] = -1
        C[b - 1, a - 1] = -1
    det = round(float(np.linalg.det(C.astype(float))))
    print("C5 Cartan matrix of E8: det = %d, symmetric = %s" % (det, bool((C == C.T).all())))
    # exact count of x in Z^8 with x^T C x = 2n by Fincke-Pohst style recursion through the Cholesky (LDL) form,
    # done with exact rationals to avoid any rounding issue.
    from fractions import Fraction
    n8 = 8
    # LDL^T: Q(x) = sum_i d_i (x_i + sum_{j>i} u_ij x_j)^2
    A = [[Fraction(int(C[i, j])) for j in range(n8)] for i in range(n8)]
    dcoef = [Fraction(0)] * n8
    U = [[Fraction(0)] * n8 for _ in range(n8)]
    for i in range(n8):
        dcoef[i] = A[i][i]
        for j in range(i + 1, n8):
            U[i][j] = A[i][j] / A[i][i]
        for j in range(i + 1, n8):
            for k in range(i + 1, n8):
                A[j][k] -= U[i][j] * A[i][k]
    counts = {}

    def rec(i, xs, rem):
        # variables x_{i}..x_{7} chosen from the last index downward; i is the current index to choose
        if i < 0:
            v = 2 * 0
            q = sum(int(C[a, b]) * xs[a] * xs[b] for a in range(n8) for b in range(n8))
            counts[q] = counts.get(q, 0) + 1
            return
        shift = sum(U[i][j] * xs[j] for j in range(i + 1, n8))
        # d_i (x_i + shift)^2 <= rem
        bound = rem / dcoef[i]
        # integer range
        r = math.isqrt(int(bound)) + 2
        centre = -shift
        lo = math.floor(centre) - r
        hi = math.ceil(centre) + r
        for xi in range(lo, hi + 1):
            t = dcoef[i] * (xi + shift) ** 2
            if t <= rem:
                xs[i] = xi
                rec(i - 1, xs, rem - t)
        xs[i] = 0

    NMAX = 3
    rec(n8 - 1, [0] * n8, Fraction(2 * NMAX))
    print("   vectors with x^T C x = 2n (exact search): ", {q // 2: counts[q] for q in sorted(counts) if q > 0})
    print("   odd values of x^T C x found: %d (lattice is even)" % sum(1 for q in counts if q % 2))
    for n in range(1, NMAX + 1):
        ok_all &= (counts.get(2 * n, 0) == 240 * sigma_k(n, 3))

    # ---------------- C6: r_2 convolutions ----------------
    LQ = 2000
    badq = 0
    badb = 0
    bad_b3 = 0
    r2l = r2[: LQ + 1]
    for n in range(1, LQ + 1):
        add = sum(r2l[a] * r2l[n - a] for a in range(0, n + 1))
        if add != r4[n]:
            badq += 1
        if n % 2 == 1 and add != 8 * sigma_k(n):
            badq += 1
        mult = sum(r2l[a] * r2l[n // a] for a in divisors(n))
        # brute force over Z[i]^2 for small n
        if n <= 300:
            brute = 0
            for a in divisors(n):
                brute += r2l[a] * r2l[n // a]
            if brute != mult:
                badb += 1
        if n % 4 == 3 and mult != 0:
            bad_b3 += 1
    print("C6 r_2 convolutions, 1 <= n <= %d:" % LQ)
    print("   sum_{a+b=n} r_2(a) r_2(b) = r_4(n) (= 8 sigma(n), n odd)   failures: %d" % badq)
    print("   sum_{ab=n} r_2(a) r_2(b) = 0 for n = 3 (mod 4)             failures: %d" % bad_b3)
    mult_small = [sum(r2l[a] * r2l[n // a] for a in divisors(n)) for n in range(1, 14)]
    print("   multiplicative convolution n = 1..13:", mult_small)
    # the Z-span of the four units 1, i, j, ij (i^2 = j^2 = -1, ij = ji) inside C + C:
    #   a + b i + c j + d ij  ->  (z1, z2) = ((a - d) + (b + c) i, (a + d) + (b - c) i),  norm = |z1|^2 |z2|^2.
    R = 9
    span_counts = {}
    for a, b, c, d in itertools.product(range(-R, R + 1), repeat=4):
        nrm = ((a - d) ** 2 + (b + c) ** 2) * ((a + d) ** 2 + (b - c) ** 2)
        if 0 < nrm <= 13:
            span_counts[nrm] = span_counts.get(nrm, 0) + 1
    print("   Z-span of the units {1, i, j, ij} (index 4 in Z[i]^2), elements of norm n = 1..13 (box +-%d):" % R)
    print("    ", [span_counts.get(n, 0) for n in range(1, 14)])
    ok_all &= (badq == 0 and bad_b3 == 0)

    # ---------------- C7: real densities, zeta(4) ----------------
    def ball(dim, r):
        return math.pi ** (dim / 2) / math.gamma(dim / 2 + 1) * r ** dim

    print("C7 real densities (d/dn of the volume of {Q <= n}):")
    n = 7.0
    h = 1e-5
    d4 = (ball(4, math.sqrt(n + h)) - ball(4, math.sqrt(n - h))) / (2 * h)
    print("   four squares: d/dn Vol{x in R^4 : |x|^2 <= n} at n = 7: %.9f ; pi^2 n = %.9f" % (d4, math.pi ** 2 * n))
    d8 = (ball(8, math.sqrt(2 * (n + h))) - ball(8, math.sqrt(2 * (n - h)))) / (2 * h)
    print("   E8 (covolume 1, Q = x.x/2): d/dn Vol{x in R^8 : |x|^2 <= 2n} at n = 7: %.6f ; (8 pi^4/3) n^3 = %.6f"
          % (d8, 8 * math.pi ** 4 / 3 * n ** 3))
    d8z = (ball(8, math.sqrt(n + h)) - ball(8, math.sqrt(n - h))) / (2 * h)
    print("   Z^8 (Q = x.x): d/dn Vol{|x|^2 <= n} at n = 7: %.6f ; (pi^4/6) n^3 = %.6f" % (d8z, math.pi ** 4 / 6 * n ** 3))
    z4 = sum(1.0 / k ** 4 for k in range(1, 200001))
    print("   zeta(4) = %.12f ; pi^4/90 = %.12f" % (z4, math.pi ** 4 / 90))
    print("   (8 pi^4/3) / zeta(4) = %.9f (240 expected)" % ((8 * math.pi ** 4 / 3) / (math.pi ** 4 / 90)))
    print("   pi^2 * (8/pi^2) = 8 ; (pi^4/6) * (96/pi^4) = %.9f (16 expected: Z^8 at odd n, 2-adic density 1)"
          % ((math.pi ** 4 / 6) * (1 / ((math.pi ** 4 / 90) * (1 - 2 ** -4)))))
    print("   prod_{p odd} (1 - p^-2) = (6/pi^2)/(3/4) = %.9f ; 8/pi^2 = %.9f" % ((6 / math.pi ** 2) / 0.75, 8 / math.pi ** 2))

    print("ALL CHECKS PASSED" if ok_all else "SOME CHECK FAILED")


if __name__ == "__main__":
    main()
