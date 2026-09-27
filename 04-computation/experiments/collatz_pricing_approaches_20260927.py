#!/usr/bin/env python3
"""collatz_pricing_approaches_20260927.py -- D12, pricing the approaches (session collatz-posets-zeta5-20260927,
opus, 2026-09-27, fifth note). Syracuse map U(m) = (3m+1)/2^v on odd m.

 (1) The price sheet of a cycle shadow. For a cycle word w (valuations v_1..v_p, A = sum v), rational cycle point
     x_w = S_w/(2^A - 3^p), multiplier mu = 3^p/2^A: for m = x_w + 2^K u (u odd, K >= A+1) the orbit follows w for
     r = floor((K-1)/A) periods and U^(pr)(m) - x_w = mu^r (m - x_w) = 3^(pr) 2^(K-Ar) u exactly. So the deviation loses A bits
     of 2-adic precision and gains p ternary digits per period; the value grows by c_w = (p log_2 3 - A)/A bits per bit of
     precision spent; for a negative integer cycle point the precision K is at most log_2(m + |x_w|), so one shadow multiplies
     the value by at most (m + |x_w|)^(c_w). Checked exactly for -1, -5, -17 and the trivial cycle, K <= 60.
 (2) The three-place multiplier: |3/2^v|_oo |3/2^v|_2 |3/2^v|_3 = 1 per odd step; along a shadow the 2-adic distance to the
     cycle point is multiplied by 2^A, the 3-adic by 3^(-p), the real by 3^p/2^A: real growth = 2-adic expansion x 3-adic
     contraction, exactly. Checked on same-word pairs of integers.
 (3) Regeneration depths along orbits: for n < 2*10^5, the maximal depth v_2(m_j - x) of approaches to x in {-1, -5, -17} along
     the orbit, against log_2 n and against the size bound log_2(m_j + |x|) (never violated), and the empirical law
     P(depth >= k) against the uniform 2^(1-k) (odd m) of regeneration.
 (4) Growth attribution: the share of the positive height increments of an orbit that occur inside shadows of the three negative
     cycle points at depth >= A + 2, for n < 2*10^5 and for the record orbits 27, 703, 6171, 77031, 837799.
 (5) Explicit double-excursion integers (the reset lane's obstruction, priced): the source class of the word (1,2)^r1 + bridge +
     (1,2)^r2 with r2 > r1, its least positive member n, the ledger (K_1, growth_1, dip, K_2, growth_2), the check that both
     excursions stay above n, and the source price: c_(-5)(K_1 + K_2) against log_2 n.
Usage: python3 collatz_pricing_approaches_20260927.py
"""
import math
from fractions import Fraction

LOG23 = math.log2(3)
CYCLES = {
    "-1": (1,),
    "-5": (1, 2),
    "-17": (1, 1, 1, 2, 1, 1, 4),
    "1 (trivial)": (2,),
}


def U(m):
    m = 3 * m + 1; v = 0
    while m % 2 == 0:
        m //= 2; v += 1
    return m, v


def v2(x):
    x = abs(x); v = 0
    while x % 2 == 0 and x:
        x //= 2; v += 1
    return v


def v3(x):
    x = abs(x); v = 0
    while x % 3 == 0 and x:
        x //= 3; v += 1
    return v


def cycle_point(w):
    p = len(w); A = sum(w); S = 0; d = 0
    for t in range(p):
        S = 3 * S + 2 ** d; d += w[t]
    # x with U^p(x) = x along w: 2^A x = 3^p x + S  -> x = S/(2^A - 3^p)
    return Fraction(S, 2 ** A - 3 ** p), p, A


def part1():
    print("== (1) the price sheet of a cycle shadow ==")
    print(" cycle, word, x_w, mu = 3^p/2^A, c_w = (p log_2 3 - A)/A bits of growth per bit of precision")
    for name, w in CYCLES.items():
        x, p, A = cycle_point(w)
        mu = Fraction(3 ** p, 2 ** A)
        c = (p * LOG23 - A) / A
        print("  %-12s %-22s x_w = %-6s mu = %-10s c_w = %+.5f" % (name, w, x, mu, c))
    ok = True; rows = []
    for name, w in CYCLES.items():
        x, p, A = cycle_point(w)
        mu = Fraction(3 ** p, 2 ** A)
        for K in (12, 24, 36, 48, 60):
            for u in (1, 3, 7):
                m = x + 2 ** K * u
                if m <= 0 or m.denominator != 1:
                    continue
                m = int(m); r = (K - 1) // A
                y = m; word = []
                for _ in range(p * r):
                    y, v = U(y); word.append(v)
                ok &= (word == list(w) * r)
                dev = y - x
                ok &= (dev == mu ** r * (m - x))
                ok &= (dev == Fraction(3 ** (p * r) * 2 ** (K - A * r) * u))
                if u == 1 and K in (24, 60):
                    growth = math.log2(y) - math.log2(m)
                    c = (p * LOG23 - A) / A
                    rows.append((name, K, r, round(growth, 3), round(c * K, 3), round(math.log2(m + abs(x)), 3)))
    print(" exactness U^(pr)(m) - x_w = mu^r (m - x_w) = 3^(pr) 2^(K - Ar) u, word = w^r, for K <= 60, u in {1,3,7}: %s" % ok)
    print(" (cycle, K, periods r, growth in bits, c_w K, log_2(m + |x_w|)): the growth is at most c_w K and K <= log_2(m + |x_w|)")
    for row in rows:
        print("  ", row)
    print(" reading: the -1 shadow (all ones) attains the universal maximum log_2 3 - 1 = 0.585 growth bits per bit of precision;")
    print(" the -5 and -17 shadows are cheap approaches (0.057 and 0.0086 bits per bit); the trivial cycle charges -0.2075: descent")


def part2():
    print("== (2) the three-place multiplier identity along a shadow ==")
    ok = True
    for name, w in CYCLES.items():
        if name.startswith("1"):
            continue
        x, p, A = cycle_point(w)
        for K in (20, 30):
            for u, u2 in ((1, 5), (3, 11)):
                m1 = int(x + 2 ** K * u); m2 = int(x + 2 ** K * u2)
                r = (K - 1) // A
                y1, y2 = m1, m2
                for _ in range(p * r):
                    y1, _v = U(y1); y2, _v2 = U(y2)
                d0 = m1 - m2; d1 = y1 - y2
                # real ratio, 2-adic ratio, 3-adic ratio of the pair's distances
                real = Fraction(d1, d0)
                ok &= (real == Fraction(3 ** (p * r), 2 ** (A * r)))
                ok &= (v2(d1) - v2(d0) == -A * r) and (v3(d1) - v3(d0) == p * r)
    print(" for same-word pairs over r periods: real distance x 3^(pr)/2^(Ar), 2-adic distance x 2^(Ar), 3-adic distance x 3^(-pr); product 1: %s" % ok)
    print(" reading: the Syracuse map expands 2-adically by 2^v, contracts 3-adically by 1/3 (within a valuation class) and the real")
    print(" multiplier 3/2^v is the reciprocal of their product; real growth is exactly the excess of 3-adic contraction over 2-adic expansion")


def orbit_depths(n, xs):
    m = n; depths = {x: 0 for x in xs}; viol = 0
    while m != 1:
        for x in xs:
            k = v2(m - x)
            depths[x] = max(depths[x], k)
            if k > math.log2(m + abs(x)) + 1e-9:
                viol += 1
        m, _ = U(m)
    return depths, viol


def part3():
    print("== (3) regeneration depths along orbits (n < 2*10^5) ==")
    xs = (-1, -5, -17)
    N = 2 * 10 ** 5
    maxdepth = {x: (0, 0) for x in xs}; excess = {x: 0.0 for x in xs}; viol = 0
    hist = {x: {} for x in xs}
    for n in range(1, N, 2):
        m = n
        while m != 1:
            for x in xs:
                k = v2(m - x)
                hist[x][k] = hist[x].get(k, 0) + 1
                if k > maxdepth[x][0]:
                    maxdepth[x] = (k, n)
                if k > math.log2(m + abs(x)) + 1e-9:
                    viol += 1
                if k - math.log2(n) > excess[x]:
                    excess[x] = k - math.log2(n)
            m, _ = U(m)
    print(" size bound v_2(m_j - x) <= log_2(m_j + |x|) violated: %d times (a theorem: a positive integer = x mod 2^k with x < 0 is >= 2^k - |x|)" % viol)
    for x in xs:
        tot = sum(hist[x].values())
        tail = [round(sum(c for k, c in hist[x].items() if k >= kk) / tot, 5) for kk in (2, 4, 6, 8, 10, 12)]
        print(" x = %d: max depth %d (at n = %d); max (depth - log_2 n) = %.2f; P(depth >= 2,4,...,12) = %s vs the uniform law 2^(1-k) = %s" % (x, maxdepth[x][0], maxdepth[x][1], excess[x], tail, [round(2.0 ** (1 - kk), 5) for kk in (2, 4, 6, 8, 10, 12)]))
    print(" reading: approaches to a negative cycle point deeper than log_2 n do occur (regenerated precision), but never deeper than log_2 of the current value")


def attribution(n, thresholds):
    m = n; total_up = 0.0; shadow_up = {x: 0.0 for x in thresholds}
    while m != 1:
        m2, v = U(m)
        inc = LOG23 - v
        if inc > 0:
            total_up += inc
            for x, K in thresholds.items():
                if v2(m - x) >= K:
                    shadow_up[x] += inc
                    break
        m = m2
    return total_up, shadow_up


def part4():
    print("== (4) growth attribution: share of positive height increments inside shadows of -1, -5, -17 ==")
    N = 2 * 10 ** 5
    for label, thr in (("shallow (depth >= A + 2: -1 >= 3, -5 >= 5, -17 >= 13)", {-1: 3, -5: 5, -17: 13}),
                       ("deep (-1 >= 8, -5 >= 10, -17 >= 16)", {-1: 8, -5: 10, -17: 16})):
        T = 0.0; S = {x: 0.0 for x in thr}
        for n in range(3, N, 2):
            t, s_ = attribution(n, thr)
            T += t
            for x in thr:
                S[x] += s_[x]
        print(" %s: all odd n < %d, total positive growth %.0f bits; shares -1: %.1f%%, -5: %.1f%%, -17: %.2f%%" % (label, N, T, 100 * S[-1] / T, 100 * S[-5] / T, 100 * S[-17] / T))
        for n in (27, 703, 6171, 77031, 837799):
            t, s_ = attribution(n, thr)
            print("   n = %d: positive growth %.1f bits; shares -1: %.1f%%, -5: %.1f%%, -17: %.1f%%" % (n, t, 100 * s_[-1] / t, 100 * s_[-5] / t, 100 * s_[-17] / t))
    print(" reading: runs of two or more ones (shallow -1 shadows) carry about half of all growth; deep shadows carry a few percent;")
    print(" the rest is generic no-descent words, which no cycle label prices")


def source_of_word(word):
    # least positive n whose exact Syracuse valuation word begins with `word` (THM-4512's exact class mod 2^(A+1))
    L = len(word); A = sum(word); S = 0; d = 0
    for t in range(L):
        S = 3 * S + 2 ** d; d += word[t]
    mod = 2 ** (A + 1)
    n = ((2 ** A - S) * pow(3, -L, mod)) % mod
    return n, mod


def part5():
    print("== (5) explicit double-excursion integers: the reset lane's obstruction, priced ==")
    configs = (("-5", (1, 2), 8, (2, 2), 16), ("-5", (1, 2), 12, (2, 2, 2), 24), ("-5", (1, 2), 16, (2, 2, 2), 40),
               ("-1", (1,), 10, (2, 2), 20), ("-1", (1,), 16, (2, 2, 2), 32))
    for name, per, r1, bridge, r2 in configs:
        x, p, A = cycle_point(per)
        c = (p * LOG23 - A) / A
        word = per * r1 + bridge + per * r2
        n, mod = source_of_word(word)
        m = n; vals = [n]
        for v in word:
            m2, vv = U(m)
            assert vv == v, (n, word)
            m = m2; vals.append(m)
        L1 = p * r1; L2 = L1 + len(bridge)
        K1 = v2(n - x); K2 = v2(vals[L2] - x)
        g1 = math.log2(vals[L1]) - math.log2(n)
        dip = math.log2(vals[L2]) - math.log2(vals[L1])
        g2 = math.log2(vals[-1]) - math.log2(vals[L2])
        above = min(vals[1:]) > n
        print(" cycle %s: word w^%d + %s + w^%d: least source n has %.1f bits (class mod 2^%d)" % (name, r1, bridge, r2, math.log2(n), sum(word) + 1))
        print("   K_1 = %d, growth_1 = %.2f bits; bridge dip = %+.2f bits; re-entry depth K_2 = %d (> K_1: %s), growth_2 = %.2f bits; whole orbit above n: %s" % (K1, g1, dip, K2, K2 > K1, g2, above))
        print("   ledger: dip %.2f < growth_2 %.2f, so log m + (prepaid growth) jumps up at the re-entry; source price c_w (K_1 + K_2) = %.2f bits against log_2 n = %.1f; total growth %.2f = %.3f log_2 n" % (-dip, g2, c * (K1 + K2), math.log2(n), g1 + dip + g2, (g1 + dip + g2) / math.log2(n)))
    print(" reading: the second approach is deeper than the first (regenerated precision), the dip between them is smaller than the")
    print(" second growth, and the whole orbit stays above the source: no bounded-lookahead potential is monotone on these; every such")
    print(" configuration is bought with a source of size 2^(precision consumed), so the combined growth is at most c_w times log_2 n")


def A_of(word):
    return sum(word)


def main():
    part1(); part2(); part3(); part4(); part5()


if __name__ == '__main__':
    main()
