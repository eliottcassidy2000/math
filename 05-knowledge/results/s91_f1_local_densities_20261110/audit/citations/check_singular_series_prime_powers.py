#!/usr/bin/env python3
"""Auditor 'citations' (S91): is the Euler product of THM-5128 (D3) a singular series in the classical sense
(product of p-adic densities), or only the law modulo squarefree moduli?

Part 1 (exact).  For independent uniform states (a1, b1), (a2, b2) in (Z/p^k)^2 and a unit u, the law of
a1 a2 + u b1 b2 modulo p^k (the generic limit of THM-5128 extended to prime-power moduli; (B2) gives the
equidistribution of the states for every odd modulus).  Compared with
   (i)  'uniform inside the class modulo p':  p^-k (1 - 1/p^2 + [p | c]/p)      (what S_sigma(N) of (D3) would give)
   (ii) the classical determinant-form density  p^-k (1 - p^-2) sum_{j <= v_p(c)} p^-j   for c != 0 mod p^k.

Part 2 (data).  All F1 triples (h1, h2, s) with N <= X, X = 2^22: among the triples of the sheet s >= 1 with
3 | N, the shares of N = 0, 3, 6 (mod 9); among those of the sheet s = 0 with 3 | N + 1, the shares of
N + 1 = 0, 3, 6 (mod 9).  'Uniform inside the class' predicts 1/3 each.
No code of the audited session is imported.
"""
import sys
from fractions import Fraction

import numpy as np
from numba import njit

print("Part 1: exact law of a1 a2 + u b1 b2 modulo p^k for uniform states")
for p, k, u in ((3, 2, 1), (3, 2, 2), (3, 3, 1), (5, 2, 1), (5, 2, 2), (7, 2, 3)):
    m = p ** k
    cnt = [0] * m
    # law of a1*a2 mod m, then convolve
    prod = [0] * m
    for a in range(m):
        for b in range(m):
            prod[a * b % m] += 1
    for c1 in range(m):
        if prod[c1]:
            for c2 in range(m):
                if prod[c2]:
                    cnt[(c1 + u * c2) % m] += prod[c1] * prod[c2]
    tot = m ** 4
    ok_i = True
    ok_ii = True
    for c in range(m):
        P = Fraction(cnt[c], tot)
        naive = Fraction(1, m) * (1 - Fraction(1, p * p) + (Fraction(1, p) if c % p == 0 else 0))
        if P != naive:
            ok_i = False
        if c % m != 0:
            v = 0
            cc = c
            while cc % p == 0:
                cc //= p
                v += 1
            det = Fraction(1, m) * (1 - Fraction(1, p * p)) * sum(Fraction(1, p ** j) for j in range(v + 1))
            if P != det:
                ok_ii = False
    by_v = {}
    for c in range(m):
        v = 0
        cc = c
        while cc and cc % p == 0:
            cc //= p
            v += 1
        key = v if c else "0"
        by_v.setdefault(key, set()).add(Fraction(cnt[c], tot) * m)
    print("  p^k = %d^%d, u = %d: m * P(c) by v_p(c): %s" % (p, k, u, {kk: [str(x) for x in sorted(vv)] for kk, vv in by_v.items()}))
    print("     (i) uniform inside the class mod p: %s   (ii) determinant-form density for c != 0: %s" % (ok_i, ok_ii))
    avg = sum(Fraction(cnt[c], tot) for c in range(0, m, p)) * p
    print("     mass of the class 0 mod p, times p: %s  (= 1 - 1/p^2 + 1/p = %s)" % (avg, 1 - Fraction(1, p * p) + Fraction(1, p)))


@njit(cache=False)
def psi_table(n):
    psi = np.zeros(n, dtype=np.int64)
    for h in range(1, n):
        k = h >> 1
        if h & 1:
            psi[h] = psi[k] + (k >> 1) + 1
        else:
            psi[h] = psi[k] + k
    return psi


@njit(cache=False)
def census(X, x, y, hmax):
    # counts[sheet][residue mod 9 of M], M = x1 x2 + 2^s y1 y2  (N = M - [s = 0])
    counts = np.zeros((2, 9), dtype=np.int64)
    total = np.zeros(2, dtype=np.int64)
    s = 0
    while (1 << s) + 1 <= X + 1:
        pw = 1 << s
        lim = X + 1 if s == 0 else X
        sheet = 0 if s == 0 else 1
        for h1 in range(hmax):
            x1 = x[h1]
            y1 = pw * y[h1]
            # x1 > y(h1)/2, so every M with this or a later h1 exceeds y(h1)//2 + pw*y(h1)
            if (y[h1] // 2) + y1 > lim:
                break
            if x1 + y1 > lim:
                continue
            for h2 in range(hmax):
                M = x1 * x[h2] + y1 * y[h2]
                if M > lim:
                    # x, y are not monotone in h jointly, but M >= (x1/2 + y1) * y2 grows with h2: stop when the
                    # lower bound exceeds the limit
                    if (x1 // 2 + y1) * y[h2] > lim:
                        break
                    continue
                counts[sheet, M % 9] += 1
                total[sheet] += 1
        s += 1
    return counts, total


print()
print("Part 2: F1 triples with N <= 2^22 (ordered triples; M = N + [s = 0])")
X = 1 << 22
hmax = X // 2 + 2
psi = psi_table(hmax)
hh = np.arange(hmax, dtype=np.int64)
x = 2 * psi + 1
y = 2 * hh + 1
counts, total = census(X, x, y, hmax)
for sheet in (0, 1):
    c = counts[sheet]
    T = total[sheet]
    cls = c[0] + c[3] + c[6]
    print("  sheet %d: %d triples; share of 3 | M: %.5f (generic limit (1 - 1/9 + 1/3)/3 = %.5f)" % (sheet, T, cls / T, (1 - 1 / 9 + 1 / 3) / 3))
    sh = [c[0] / cls, c[3] / cls, c[6] / cls]
    se = (1 / 3 * 2 / 3 / cls) ** 0.5
    print("     inside the class 3 | M: shares of M = 0, 3, 6 (mod 9): %.5f %.5f %.5f   (uniform: 0.33333; binomial s.e. %.5f)" % (sh[0], sh[1], sh[2], se))
    print("     determinant-form generic limit: %.5f %.5f %.5f" % (0.14403292 / 0.40740741, 0.13168724 / 0.40740741, 0.13168724 / 0.40740741))
    other = [c[r] / T for r in (1, 2, 4, 5, 7, 8)]
    print("     classes prime to 3 (mod 9), shares: " + " ".join("%.5f" % v for v in other))
sys.exit(0)
