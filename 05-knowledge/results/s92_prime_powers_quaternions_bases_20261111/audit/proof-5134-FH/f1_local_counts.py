#!/usr/bin/env python3
"""Auditor proof-5134-FH.  THM-5134 (F): local class counts of the two quaternary forms.

Independent code (written before reading the session scripts).

For a modulus m:
  D_u(c) = #{(a1,a2,b1,b2) mod m : a1 a2 + u b1 b2 = c}          (u a unit)
  S(c)   = #{(a,b,c,d) mod m    : a^2+b^2+c^2+d^2 = c}
  m^3 L_m(c), L_m(e) = sum_{d | gcd(m,e)} d^-1 prod_{p | m/d} (1 - p^-2)   (THM-5134 definition)
and the closed prime-power form
  L_{p^k}(e) = (1+1/p)(1-p^(-v-1))   (v = v_p(e) < k),     1 + 1/p - p^(-k-1)   (v >= k).

Two methods for the counts: (1) cyclic convolution of exact integer histograms,
(2) direct enumeration of all m^4 tuples with numpy (for m <= 49).
"""
import sys
from fractions import Fraction
from math import gcd
import numpy as np


def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def divisors(n):
    ds = [1]
    for p, k in factor(n).items():
        ds = [d * p ** i for d in ds for i in range(k + 1)]
    return sorted(ds)


def L_def(m, e):
    """L_m(e) from the definition in the Setting of THM-5134."""
    g = gcd(m, e % m) if e % m else m
    tot = Fraction(0)
    for d in divisors(g):
        t = Fraction(1, d)
        for p in factor(m // d):
            t *= 1 - Fraction(1, p * p)
        tot += t
    return tot


def L_closed(m, e):
    """Product of the closed prime-power values of statement (D)."""
    tot = Fraction(1)
    for p, k in factor(m).items():
        ee = e % (p ** k)
        if ee == 0:
            v = k
        else:
            v = 0
            while ee % p == 0:
                ee //= p
                v += 1
        if v < k:
            tot *= (1 + Fraction(1, p)) * (1 - Fraction(1, p ** (v + 1)))
        else:
            tot *= 1 + Fraction(1, p) - Fraction(1, p ** (k + 1))
    return tot


def conv(f, g, m):
    h = [0] * m
    for i, fi in enumerate(f):
        if fi == 0:
            continue
        for j, gj in enumerate(g):
            if gj:
                h[(i + j) % m] += fi * gj
    return h


def hist_products(m, u):
    f = [0] * m
    for a in range(m):
        for b in range(m):
            f[(u * a * b) % m] += 1
    return f


def hist_squares(m):
    f = [0] * m
    for a in range(m):
        f[(a * a) % m] += 1
    return f


def det_counts(m, u):
    return conv(hist_products(m, 1), hist_products(m, u), m)


def foursq_counts(m):
    s = hist_squares(m)
    s2 = conv(s, s, m)
    return conv(s2, s2, m)


def direct_det(m, u):
    a = np.arange(m, dtype=np.int64)
    p1 = (a[:, None] * a[None, :]) % m           # a1 a2
    p2 = (u * p1) % m                           # u b1 b2
    tot = (p1.reshape(-1)[:, None] + p2.reshape(-1)[None, :]) % m
    return np.bincount(tot.reshape(-1), minlength=m).tolist()


def direct_foursq(m):
    a = np.arange(m, dtype=np.int64)
    s2 = ((a * a)[:, None] + (a * a)[None, :]) % m
    tot = (s2.reshape(-1)[:, None] + s2.reshape(-1)[None, :]) % m
    return np.bincount(tot.reshape(-1), minlength=m).tolist()


def main():
    out = []
    P = out.append
    P("THM-5134 (F): class counts of a1 a2 + u b1 b2 and of a^2+b^2+c^2+d^2 modulo m")
    P("=" * 78)
    mods = [3, 5, 7, 9, 11, 13, 15, 21, 25, 27, 33, 35, 45, 49, 63, 75, 81, 105, 121, 125]
    allok = True
    for m in mods:
        units = [u for u in range(1, m) if gcd(u, m) == 1]
        fs = foursq_counts(m)
        target = [L_def(m, c) * m ** 3 for c in range(m)]
        closed_ok = all(L_def(m, c) == L_closed(m, c) for c in range(m))
        integer_ok = all(t.denominator == 1 for t in target)
        fs_ok = all(Fraction(fs[c]) == target[c] for c in range(m))
        # all units u (for m <= 49), a sample otherwise
        us = units if m <= 49 else [1, 2, units[len(units) // 2], units[-1]]
        det_ok = True
        for u in us:
            dc = det_counts(m, u)
            if any(Fraction(dc[c]) != target[c] for c in range(m)):
                det_ok = False
        mean = sum(L_def(m, c) for c in range(m)) / m
        ok = closed_ok and integer_ok and fs_ok and det_ok and mean == 1
        allok &= ok
        P(f"m = {m:4d}: closed form = definition: {closed_ok}; m^3 L integer: {integer_ok}; "
          f"four squares = m^3 L_m: {fs_ok}; det form (units tested: {len(us)}) = m^3 L_m: {det_ok}; mean L = {mean}")
    P("")
    P("Direct enumeration of all m^4 tuples (numpy), independent of the convolution code:")
    for m in [9, 25, 27, 49]:
        fs_conv = foursq_counts(m)
        fs_dir = direct_foursq(m)
        okf = fs_conv == fs_dir
        okd = True
        for u in [1, 2, m - 1]:
            okd &= (det_counts(m, u) == direct_det(m, u))
        allok &= okf and okd
        P(f"  m = {m}: four squares direct = convolution: {okf}; det form direct = convolution (u = 1, 2, -1): {okd}")
    P("")
    P("Tables (class counts by valuation; count / m^3 = L_m):")
    for m in [9, 25, 27, 49]:
        p = min(factor(m))
        k = factor(m)[p]
        fs = foursq_counts(m)
        dc = det_counts(m, 1)
        P(f"  m = {m} = {p}^{k}")
        for v in range(k + 1):
            cls = [c for c in range(m) if (c % m == 0 and v == k) or (c % m != 0 and v < k and c % p ** v == 0 and c % p ** (v + 1) != 0)]
            vals_fs = sorted(set(fs[c] for c in cls))
            vals_dc = sorted(set(dc[c] for c in cls))
            lim = L_closed(m, cls[0])
            P(f"    v = {v}{'+' if v == k else ' '}: {len(cls):3d} classes; four squares {vals_fs}; det form {vals_dc}; "
              f"m^3 L = {lim * m ** 3}; L = {lim} = {float(lim):.6f}")
    P("")
    # ---- the prime 2 ----
    P("The prime 2: number of solutions modulo 2^k divided by 2^(3k)")
    for k in range(1, 7):
        m = 2 ** k
        fs = foursq_counts(m)
        dc = det_counts(m, 1)
        odd_fs = sorted(set(Fraction(fs[c], m ** 3) for c in range(1, m, 2)))
        odd_dc = sorted(set(Fraction(dc[c], m ** 3) for c in range(1, m, 2)))
        ev_fs = {c: str(Fraction(fs[c], m ** 3)) for c in range(0, m, 2)} if k <= 3 else "..."
        P(f"  2^{k}: four squares, odd classes: {[str(x) for x in odd_fs]}; det form a1a2+b1b2, odd classes: {[str(x) for x in odd_dc]}"
          f"; four squares even classes: {ev_fs}")
    P("  (four squares: density 1 at every odd class; determinant form: 3/4 = 1 - 1/4 at every odd class)")
    P("")
    # ---- Jacobi ----
    P("Jacobi: r_4(n) by direct count against 8 sigma(n) (n odd) and 24 sigma(odd part) (n even), n <= 2000")
    NMAX = 2000
    R = int(NMAX ** 0.5) + 1
    sq = np.zeros(NMAX + 1, dtype=np.int64)
    for a in range(-R, R + 1):
        if a * a <= NMAX:
            sq[a * a] += 1
    r2 = np.convolve(sq, sq)[: NMAX + 1]
    r4 = np.convolve(r2, r2)[: NMAX + 1]

    def sigma(n):
        return sum(d for d in range(1, n + 1) if n % d == 0)

    okj = True
    for n in range(1, NMAX + 1):
        if n % 2:
            okj &= (r4[n] == 8 * sigma(n))
        else:
            o = n
            while o % 2 == 0:
                o //= 2
            okj &= (r4[n] == 24 * sigma(o))
    allok &= bool(okj)
    P(f"  r_4(n) = 8 sigma(n) for all odd n <= {NMAX} and 24 sigma(odd part) for even n: {bool(okj)}")
    P("")
    # ---- real density ----
    P("Real density: sum_{n <= N} r_4(n) against pi^2 N^2 / 2 (volume of the 4-ball of radius sqrt N)")
    cum = np.cumsum(r4)
    for N in [500, 1000, 2000]:
        P(f"  N = {N}: sum = {cum[N]} ; pi^2 N^2/2 = {np.pi ** 2 * N * N / 2:.1f}; ratio = {cum[N] / (np.pi ** 2 * N * N / 2):.5f}")
    P("")
    # ---- singular series identity for four squares ----
    P("Singular series check: r_4(n) = (pi^2 n) * 1 * prod_{p odd} (1 - p^-2) sigma(p^v)/p^v  = 8 sigma(n) for odd n")
    P("  prod_{p odd}(1 - p^-2) = (6/pi^2)/(1 - 1/4) = 8/pi^2 ;  numerically:")
    import math
    P(f"  6/pi^2 = {6 / math.pi ** 2:.10f};  (6/pi^2)/(3/4) = {6 / math.pi ** 2 / 0.75:.10f};  8/pi^2 = {8 / math.pi ** 2:.10f}")
    pr = 1.0
    N = 2000000
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:: i] = False
    primes = np.nonzero(sieve)[0]
    pr = float(np.prod(1 - 1.0 / primes[1:].astype(np.float64) ** 2))
    P(f"  prod_{{3 <= p <= {N}}} (1 - p^-2) = {pr:.10f}")
    P("")
    P("ALL EXACT CHECKS PASS" if allok else "SOME CHECK FAILED")
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
