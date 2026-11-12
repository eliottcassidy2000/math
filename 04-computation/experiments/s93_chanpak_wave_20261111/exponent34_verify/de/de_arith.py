#!/usr/bin/env python3
"""Independent verifier (part D/E) - the counting and exponent arithmetic of (D), (E), exactly.

1. |W| >= 2^(n-2) with R = ceil(log2 n) + 1: exact count of run-free words by a transfer recursion (cross-checked by
   enumeration for n <= 18), the union bound (n - R) 2^(n-1-R), and 2n <= 2^R < 4n.
2. The cells (t', s): t <= t' <= n-1, 0 <= s <= n-1-t', |s - t'| <= R+3.  Count <= n(2R+7); w >= 1;
   w <= n - 2t' + R + 3; max(t'-1, s) <= t' + R + 3; 4^w <= 4^(R+3) 2^n; exponent of the gcd term.
3. The elementary variant (beta0) and the trivial variant (beta = 2, Remark 1).
4. The polynomial factors.
5. The passage from n to X.
"""
import math
from fractions import Fraction
from itertools import product


def clog2(n):
    return (n - 1).bit_length()


def runfree_count(n, R):
    """binary strings of length n with first bit fixed and every run of equal bits of length <= R"""
    a = [0] * (n + 1)
    a[0] = 1
    for m in range(1, n + 1):
        a[m] = sum(a[m - j] for j in range(1, min(R, m) + 1))
    return a[n]


def runfree_brute_all(n):
    """cum[R] = number of n-bit words (top bit 1) whose longest run of equal bits is <= R, by enumeration"""
    hist = [0] * (n + 1)
    for h in range(1 << (n - 1), 1 << n):
        best = cur = 1
        prev = h & 1
        for i in range(1, n):
            b = (h >> i) & 1
            if b == prev:
                cur += 1
                if cur > best:
                    best = cur
            else:
                cur = 1
                prev = b
        hist[best] += 1
    cum = [0] * (n + 1)
    for R in range(1, n + 1):
        cum[R] = cum[R - 1] + hist[R]
    return cum


print("== 1. |W| >= 2^(n-2) ==")
for n in range(2, 19):
    cum = runfree_brute_all(n)
    for R in range(1, n + 1):
        assert runfree_count(n, R) == cum[R], (n, R)
    # exact count of words with a run of >= R+1 equal bits starting at a given position: 2^(n-1-R) each
    for R in range(1, n):
        for pos in range(0, n - R):
            c = 0
            for h in range(1 << (n - 1), 1 << n):
                blk = (h >> pos) & ((1 << (R + 1)) - 1)
                if blk == 0 or blk == (1 << (R + 1)) - 1:
                    c += 1
            assert c == 2 ** (n - 1 - R), (n, R, pos, c)
        if n > 12:
            break
print("transfer recursion = enumeration for n <= 18, all R; each of the n - R positions carries exactly 2^(n-1-R) words (n <= 12 all R)")
worst = (10, None)
for n in list(range(8, 3001)) + [4096, 4097, 8192, 8193, 20000]:
    R = clog2(n) + 1
    assert 2 * n <= 2 ** R < 4 * n, n
    Wn = runfree_count(n, R)
    longrun = 2 ** (n - 1) - Wn
    ub = (n - R) * 2 ** (n - 1 - R)
    assert longrun <= ub <= 2 ** (n - 2), n
    assert Wn >= 2 ** (n - 2)
    frac = Fraction(Wn, 2 ** (n - 1))
    if frac < worst[0]:
        worst = (frac, n)
    # exactness of the per-position count: a run of >= R+1 equal bits at a given position
print("n = 8..3000 and 4096, 4097, 8192, 8193, 20000: 2n <= 2^R < 4n, long-run words <= (n-R)2^(n-1-R) <= 2^(n-2), |W| >= 2^(n-2): OK")
print("   smallest |W|/2^(n-1) = %.4f at n = %d (the bound needs >= 0.5)" % (float(worst[0]), worst[1]))
for n in (8, 9, 16, 17, 22, 26, 32, 33, 64, 65, 1024, 1025):
    R = clog2(n) + 1
    Wn = runfree_count(n, R)
    print("   n=%5d R=%2d  |W|/2^(n-1) = %.5f   union bound for the complement: %.5f" %
          (n, R, Wn / 2 ** (n - 1), (n - R) * 2.0 ** (-R)))
# would R' = R - 1 = ceil(log2 n) also do?  (information only)
fails = [n for n in range(8, 600) if runfree_count(n, clog2(n)) < 2 ** (n - 2)]
print("   (information) with R = ceil(log2 n) instead: |W| >= 2^(n-2) fails for %d values of n in 8..599, first few %s" % (len(fails), fails[:8]))

print("\n== 2. cells and exponents, t = ceil(n/4) ==")
beta0 = math.log2(6) / 2 + 0.25
rho0 = (beta0 - 1) / (2 * beta0 - 1)


def cells(n, t, R):
    for tp in range(t, n):
        for s in range(max(0, tp - R - 3), min(n - 1 - tp, tp + R + 3) + 1):
            yield tp, s


max_lit = Fraction(-100)
max_sharp = Fraction(-100)
max_first = Fraction(-100)
arg = None
for n in list(range(8, 601)) + [1000, 1023, 1024, 1025, 2048, 2049, 4095, 4096, 4097]:
    R = clog2(n) + 1
    t = -(-n // 4)
    assert n <= 4 * t and 4 * t <= n + 3
    ncell = 0
    for tp, s in cells(n, t, R):
        ncell += 1
        w = n - tp - s
        assert w >= 1
        assert w <= n - 2 * tp + R + 3
        mx = max(tp - 1, s)
        assert mx <= tp + R + 3
        assert 2 * w <= n + 2 * (R + 3)                       # 4^w <= 4^(R+3) 2^n
        assert 2 * mx + 3 * w <= 2 * n + 5 * (R + 3)          # literal: 2^mx 2^(3w/2) <= 2^((5/2)(R+3)) 2^n
        e = Fraction(2 * mx + 3 * w - 2 * n, 2 * (R + 3))     # (mx + 1.5 w - n)/(R+3)
        if e > max_lit:
            max_lit = e
            arg = (n, tp, s)
        f = Fraction(2 * w - n, 2 * (R + 3))                  # (w - n/2)/(R+3) <= 1
        if f > max_first:
            max_first = f
    assert ncell <= n * (2 * R + 7), (n, ncell)
print("n = 8..600 and 9 larger n: w >= 1, w <= n-2t'+R+3, max(t'-1,s) <= t'+R+3, 4^w <= 4^(R+3) 2^n, "
      "2^max 2^(3w/2) <= 2^((5/2)(R+3)) 2^n, #cells <= n(2R+7): all OK")
print("   max over cells of (max(t'-1,s) + 1.5 w - n)/(R+3) = %s = %.4f at (n,t',s)=%s   [the proof uses 2.5; 1.5 suffices]"
      % (max_lit, float(max_lit), arg))
print("   max over cells of (w - n/2)/(R+3) = %.4f   [the proof uses 1]" % float(max_first))

print("\n== 3. variants ==")
print("beta0 = %.12f ; (beta0-1)/(2 beta0-1) = %.12f ; beta0/(2 beta0-1) = %.12f" % (beta0, rho0, beta0 / (2 * beta0 - 1)))
for beta, name in ((beta0, 'beta0 (B1)'), (2.0, 'beta=2 (trivial S_w <= 4^w)'), (1.5, 'beta=3/2 (B2)')):
    lit = -100.0
    sharp_ok = True
    worst_sharp = -100.0
    n4 = True
    crit = True
    for n in list(range(8, 401)) + [1000, 1024, 1025, 4096, 4097]:
        R = clog2(n) + 1
        t = math.ceil(n * (beta - 1) / (2 * beta - 1) - 1e-12)
        if 4 * t < n:
            n4 = False
        for tp, s in cells(n, t, R):
            w = n - tp - s
            mx = max(tp - 1, s)
            if tp + beta * (n - 2 * tp) > n + 1e-9:
                crit = False
            v = (mx + beta * w - n) / (R + 3)
            lit = max(lit, v)
            worst_sharp = max(worst_sharp, v)
            if 2 * w > n + 2 * (R + 3):
                n4 = False
    print("  %-28s t = ceil(n(beta-1)/(2beta-1)): n <= 4t and 4^w <= 4^(R+3) 2^n: %s ; t'+beta(n-2t') <= n: %s ; "
          "max (max(t'-1,s) + beta w - n)/(R+3) = %.4f  (literal argument gives 1+beta = %.4f, sharp split gives beta = %.4f) ; "
          "exponent beta/(2beta-1) = %.6f" % (name, n4, crit, lit, 1 + beta, beta, beta / (2 * beta - 1)))

print("\n== 4. polynomial factors ==")
print("   n, R, n(2R+7)[4^(R+4) + n(32n)^(5/2)] / (n^(9/2) ln n), and the beta0-variant with (32n)^(1+beta0)")
for n in (8, 16, 17, 64, 1024, 10 ** 4, 10 ** 6, 10 ** 9, 10 ** 12, 10 ** 18):
    R = clog2(n) + 1
    a = n * (2 * R + 7) * (4.0 ** (R + 4) + n * (32.0 * n) ** 2.5) / (n ** 4.5 * math.log(n))
    b = n * (2 * R + 7) * (4.0 ** (R + 4) + n * (32.0 * n) ** (1 + beta0)) / (n ** 4.5 * math.log(n))
    c = n * (2 * R + 7) * (4.0 ** (R + 4) + n * (32.0 * n) ** (beta0)) / (n ** 4.5 * math.log(n))
    print("   n=%-20d R=%2d  beta=3/2 literal: %.4e   beta0 literal: %.4e (unbounded, ~ n^%.4f)   beta0 sharp split: %.4e"
          % (n, R, a, b, beta0 - 1.5, c))

print("\n== 5. from n to X ==")
# for X in [2^(n+1), 2^(n+2)): values <= 2^(n+1) - 1 < X; ratio of 2^(3n/4)/(n^(9/2) ln n) to X^(3/4)/(ln X)^5
mn = (1e99, None)
for n in range(8, 20001):
    for X in (2.0 ** 0 * 1.0,):
        pass
    lnXlo = (n + 1) * math.log(2)
    lnXhi = (n + 2) * math.log(2)
    # worst case in the dyadic block: X -> 2^(n+2) for the power, X = 2^(n+1) for the log
    r = math.exp(0.75 * n * math.log(2) - 0.75 * lnXhi) * lnXlo ** 5 / (n ** 4.5 * math.log(n))
    if r < mn[0]:
        mn = (r, n)
print("   min over n in 8..20000 of [2^(3n/4)/(n^(9/2) ln n)] / [X^(3/4)/(ln X)^5] (worst X in [2^(n+1), 2^(n+2))) = %.5f at n = %d"
      % mn)
print("   (the ratio grows like sqrt(ln X)/ln ln X, so (log X)^5 is a valid, slightly lossy, exponent; (log X)^(9/2) log log X is exact)")
print("   x(0) = 1, x(1) = 3: for 3 <= X < 2^9 the count #{x(h) <= X} >= 2 > 0, so c can be adjusted")
