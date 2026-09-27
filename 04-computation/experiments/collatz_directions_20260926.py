#!/usr/bin/env python3
"""collatz_directions_20260926.py -- computations behind the directions note (session collatz-directions-20260926,
opus, 2026-09-26).

 (1) The no-descent fractal E_inf in Z_2: W_k = number of T-parity words of length k with no coefficient descent
     (3^(o_j) >= 2^j for all prefixes j <= k) = number of length-k cylinders meeting E_inf; the box-dimension
     estimate log_2(W_k)/k for k up to 2000 (exact integer DP), converging to h* = h(log_3 2) = 0.9499555.
 (2) The map rho(r) = 3r + lsb(r) (lsb = lowest set bit as a power of two) on positive integers: r_l = 2^(d_l) m_l
     where m_l is the Syracuse orbit and d_l the cumulative valuation (checked on all odd n <= 2000 for 200
     steps); r reaches a power of two iff the orbit reaches 1 (checked n <= 100000); r_l = 3^l n + S_(l-1).
 (3) Best lower rational approximations A/p of log_2 3 (A/p < log_2 3 with A/p > A'/p' for all p' < p), the
     gaps 3^p - 2^A, and the known integer negative cycles at the first three.
 (4) E_inf cap Z_<0 for |n| <= 10^6: negative odd n whose orbit never decreases in absolute value (for negative
     integers actual and coefficient descent coincide, since |U(m)| < |m| iff v >= 2).
 (5) The identity 2^(d_l) m_l / 3^l = n C_l, C_l = prod_(j<l) (1 + 1/(3 m_j)), and the real partial sums
     S_(l-1)/3^l = n (C_l - 1) of the Bernstein series, on n = 27 (reaching 1: divergent series) and on the
     cycle points -1, -5 (convergent, equal to -n).
Usage: python3 collatz_directions_20260926.py
"""
import math
from fractions import Fraction

LOG23 = math.log2(3)
H = lambda p: -p * math.log2(p) - (1 - p) * math.log2(1 - p)
HSTAR = H(math.log(2) / math.log(3))


def part1():
    print("== (1) cylinders of E_inf: W_k and log_2(W_k)/k ==")
    cur = {0: 1}
    out = {}
    for j in range(1, 2001):
        nxt = {}
        for o, c in cur.items():
            for letter in (0, 1):
                oo = o + letter
                if 3 ** oo >= 2 ** j:  # no coefficient descent: 3^o >= 2^j (equality impossible for j >= 1)
                    nxt[oo] = nxt.get(oo, 0) + c
        cur = nxt
        W = sum(cur.values())
        if j in (10, 20, 41, 60, 100, 200, 500, 1000, 1500, 2000):
            out[j] = W
            print(" k=%4d  W_k has %d digits  log2(W_k)/k = %.6f  (h* = %.7f; log2(W_k)/k + 1.5 log2(k)/k = %.6f)" % (j, len(str(W)), math.log2(W) / j, HSTAR, math.log2(W) / j + 1.5 * math.log2(j) / j))
    return out


def lsb(r):
    return r & -r


def rho(r):
    return 3 * r + lsb(r)


def syracuse(m):
    m = 3 * m + 1
    v = 0
    while m % 2 == 0:
        m //= 2; v += 1
    return m, v


def part2():
    print("== (2) the map rho(r) = 3r + lsb(r) ==")
    for n in range(1, 2001, 2):
        r = n; m = n; d = 0; S = 0; A = 0
        for l in range(200):
            assert r == (1 << d) * m and r == 3 ** l * n + S, (n, l)
            r = rho(r)
            S = 3 * S + (1 << A)
            m, v = syracuse(m); d += v; A += v
    print(" r_l = 2^(d_l) m_l = 3^l n + S_(l-1) for all odd n <= 2000 and l <= 200: OK")
    maxsteps = 0
    for n in range(1, 100001):
        r = n; steps = 0
        while r & (r - 1):
            r = rho(r); steps += 1
        maxsteps = max(maxsteps, steps)
    print(" every n <= 100000 reaches a power of two under rho (max steps %d); once r is a power of two, rho(r) = 4r stays one" % maxsteps)
    r = 27; seq = []
    while r & (r - 1):
        seq.append(r); r = rho(r)
    print(" n = 27: r_l = %s ... reaches %d = 2^%d after %d steps; odd parts %s ..." % (seq[:6], r, r.bit_length() - 1, len(seq), [x // lsb(x) for x in seq[:12]]))


def part3():
    print("== (3) best lower approximations of log_2 3 and the gaps 3^p - 2^A ==")
    best = []
    bestval = 0.0
    for p in range(1, 400):
        A = math.floor(p * LOG23)
        val = A / p
        if val > bestval:
            bestval = val
            gap = 3 ** p - 2 ** A
            best.append((A, p, gap))
    for A, p, gap in best[:12]:
        print(" A/p = %3d/%3d = %.7f   3^p - 2^A = %s" % (A, p, A / p, gap if gap < 10 ** 12 else "%.3e" % gap))
    print(" integer negative cycles occur at 1/1 (gap 1, x = -1), 3/2 (gap 1, x = -5), 11/7 (gap 139, x = -17); 19/12 (gap 7153) has none (necklace search p <= 12, S12 audit)")


def part4():
    print("== (4) E_inf cap Z_<0 for |n| <= 10^6: negative odd n whose orbit never decreases in absolute value ==")
    found = []
    for n in range(-1, -10 ** 6, -2):
        m = n; seen = set(); ok = True
        while True:
            if abs(m) < abs(n):
                ok = False; break
            if m in seen:
                break
            seen.add(m)
            m, v = syracuse(m)
            if len(seen) > 10 ** 4:
                ok = None; break
        if ok:
            found.append(n)
    print(" members:", found)


def part5():
    print("== (5) the two-place identity on 27, -1, -5 ==")
    for n in (27, -1, -5):
        m = n; C = Fraction(1); S = 0; A = 0
        rows = []
        for l in range(1, 25):
            C *= 1 + Fraction(1, 3 * m)
            S = 3 * S + 2 ** A
            m, v = syracuse(m); A += v
            assert Fraction(2 ** A * m, 3 ** l) == n * C, (n, l)
            assert Fraction(S, 3 ** l) == n * (C - 1), (n, l)
            rows.append(float(S) / 3 ** l)
        print(" n = %d: 2^(d_l) m_l / 3^l = n C_l and S_(l-1)/3^l = n (C_l - 1) exact for l <= 24; real partial sums S/3^l: %s ... %s" % (n, ["%.4f" % x for x in rows[:5]], ["%.4f" % x for x in rows[-3:]]))
    print(" (27 reaches 1: C_l -> infinity, the real series diverges; -1 and -5 are cycle points: the real sums converge to -n, the 2-adic sums too)")


def main():
    part1(); part2(); part3(); part4(); part5()


if __name__ == '__main__':
    main()
