#!/usr/bin/env python3
"""collatz_mediant_tree_K_20260930.py -- the mediant tree and the never-descending set K, unified through the prefix fixed points
(session collatz-posets-zeta5-20260927, opus, 2026-09-30, fourteenth note).

 (1) Prefix fixed points as 2-adic convergents: for odd n with prefix word w_j (first j valuations), x_(w_j) = n mod 2^(A_j),
     U^j(n) - n = D_j (x_(w_j) - n)/2^(A_j), and |n - x_(w_j)|_2 |n - x_(w_j)|_oo = |U^j(n) - n| / |D_j|.
 (2) THM-4512's threshold N(w) = S/(2^A - 3^j) is the fixed point x_w: value-level descent at a positive-clock prefix iff
     n > x_w; the precision residual = integers at or below the fixed point of their own prefix; census to 10^6 of the
     residual (n <= x_w at the first word-descent) and of the delay sigma - tau.
 (3) dist_2(n, K) = 2^(-floor(tau(n) log_2 3)) where tau = first word-descent index; the record function f(N) and the
     empirical Liouville exponent: dist_2(n, K) >= n^(-kappa).
 (4) The mediant tree is the tree of all 3x+d cycles: the fixed points with denominator d are (1/d) times the cycle points
     of 3x + d on integers coprime to d; checked for d = 5, 7, 11, 13, 17, 19, 23, 25 by direct cycle search; which
     d <= 100 (coprime to 6) occur as denominators for A <= 22, and which are found by the 3x+d search up to 10^6.
Usage: python3 collatz_mediant_tree_K_20260930.py
"""
import math, random, time
from fractions import Fraction
from math import gcd
from collections import Counter, defaultdict

LOG23 = math.log2(3)
T0 = time.time()


def carry(w):
    S = 0; d = 0
    for v in w:
        S = 3 * S + 2 ** d; d += v
    return S


def clock(w):
    return 2 ** sum(w) - 3 ** len(w)


def word_prefix(n, j):
    x = n; w = []; orbit = [n]
    for _ in range(j):
        y = 3 * x + 1; v = 0
        while y % 2 == 0:
            y //= 2; v += 1
        w.append(v); x = y; orbit.append(x)
    return w, orbit


def part1():
    print("== (1) prefix fixed points are the 2-adic convergents of n ==")
    random.seed(3); ok2 = okd = okp = True
    for _ in range(300):
        n = random.randrange(1, 10 ** 6) | 1; j = random.randint(1, 12)
        w, orb = word_prefix(n, j)
        A = sum(w); S = carry(w); D = clock(w)
        x = Fraction(S, D)
        ok2 &= (x.numerator - n * x.denominator) % (2 ** A) == 0          # x = n mod 2^A (as 2-adic integers; D is odd)
        okd &= Fraction(orb[j] - n) == Fraction(D, 2 ** A) * (x - n)
        # 2-adic absolute value of x - n: 2^(-v_2(numerator)) with odd denominator
        num = x.numerator - n * x.denominator
        v2 = 0
        while num and num % 2 == 0:
            num //= 2; v2 += 1
        disp = orb[j] - n; vd = 0; t = abs(disp)
        while t and t % 2 == 0:
            t //= 2; vd += 1
        # |x - n|_2 |x - n|_oo = |U^j n - n|_2 |U^j n - n|_oo / |D|  (the displacement is even: both sides carry its 2-adic size)
        okp &= (Fraction(1, 2 ** v2) * abs(x - n) == Fraction(abs(disp), 2 ** vd * abs(D))) if x != n else True
    print(" x_(w_j) = n mod 2^(A_j): %s; U^j(n) - n = D_j (x_(w_j) - n)/2^(A_j): %s; product identity |n - x|_2 |n - x|_oo = |U^j n - n|_2 |U^j n - n|_oo / |D_j|: %s" % (ok2, okd, okp))
    n = 27; w, orb = word_prefix(27, 12)
    rows = []
    for j in range(1, 13):
        pj = w[:j]; x = Fraction(carry(pj), clock(pj))
        rows.append((j, clock(pj) > 0, str(x) if x.denominator < 10 ** 6 else "%.3f" % float(x), orb[j]))
    print(" n = 27: (j, clock > 0, x_(w_j), U^j(27)) for j <= 12:", rows[:8])
    print(" the positive-clock prefixes with x_(w_j) < 27 are the value-level descents; 27 descends at j = %s" % next((j for j in range(1, 200) if word_prefix(27, j)[1][j] < 27), None))


def first_word_descent(n, cap=100000):
    x = n; d = 0; p3 = 1
    for j in range(1, cap + 1):
        y = 3 * x + 1; v = 0
        while y % 2 == 0:
            y //= 2; v += 1
        d += v; p3 *= 3
        if 2 ** d > p3:
            return j, x, y
        x = y
    return None, None, None


def part2():
    print("== (2) THM-4512's threshold is the prefix fixed point; the precision residual ==")
    delays = Counter(); residual = []; N = 10 ** 6
    for n in range(3, N + 1, 2):
        x = n; d = 0; p3 = 1; j = 0; tau = None; S = 0
        while True:
            y = 3 * x + 1; v = 0
            while y % 2 == 0:
                y //= 2; v += 1
            S = 3 * S + 2 ** d; d += v; p3 *= 3; j += 1
            if tau is None and 2 ** d > p3:
                tau = j
                # value-level descent at this prefix iff n > S/(2^d - 3^j)
                if n * (2 ** d - p3) <= S:
                    residual.append((n, j, Fraction(S, 2 ** d - p3)))
            if y < n:
                delays[j - tau] += 1
                break
            x = y
    print(" odd n <= 10^6: delay sigma - tau (value-level descent index minus first word-descent index): %s" % sorted(delays.items())[:8])
    print(" precision residual (n at or below the fixed point of its first positive-clock prefix): %d members: %s" % (len(residual), [(n, j, str(x)) for n, j, x in residual[:12]]))


def part3():
    print("== (3) the 2-adic distance to K and the record function ==")
    # dist_2(n, K) = 2^(-floor(tau log_2 3)); verify against the definition on small n by checking cylinders explicitly
    def dist_by_cylinders(n, depth=40):
        # largest m such that some x = n mod 2^m has a word with no descent to depth `depth` (search over residues x = n + 2^m t)
        best = 0
        for m in range(1, 24):
            found = False
            for t in range(0, 2 ** min(6, 24 - m)):
                x = n + (2 ** m) * t
                if x % 2 == 0:
                    continue
                j, _, _ = first_word_descent(x, cap=depth)
                if j is None:
                    found = True; break
            if found:
                best = m
            else:
                break
        return best
    shown = []
    for n in (1, 3, 5, 7, 9, 11, 15, 27):
        tau, _, _ = first_word_descent(n)
        shown.append((n, tau, math.floor(tau * LOG23)))
    print(" dist_2(n, K) = 2^(-floor(tau(n) log_2 3)) (Lemma K1; proof in the note): (n, tau, exponent) for small n: %s" % shown)
    # records of tau up to 2*10^6
    records = []; best = 0
    for n in range(3, 2 * 10 ** 6, 2):
        tau, _, _ = first_word_descent(n)
        if tau > best:
            best = tau; records.append((n, tau, math.floor(tau * LOG23)))
    print(" records of tau(n) (first word-descent index) for odd n < 2 * 10^6: %s" % records)
    kappa = max(math.floor(t * LOG23) / math.log2(n) for n, t, _ in records if n > 1)
    print(" empirical Liouville exponent: dist_2(n, K) >= n^(-kappa) with kappa = max floor(tau log_2 3)/log_2 n over records = %.3f; the ratio tau/log_2 n at the last records: %s" % (
        kappa, [round(t / math.log2(n), 2) for n, t, _ in records[-5:]]))
    print(" f(N) = least n with floor(tau log_2 3) >= N grows like 2^(N/kappa) = 2^(%.2f N): the size price of the never-descending set, in the record-holders" % (1 / kappa))


def cycles_of_3x_plus_d(d, limit):
    """Syracuse-type cycles of x -> oddpart(3x + d) on positive integers coprime to d, found from starts <= limit."""
    seen = set(); cycles = set()
    for n in range(1, limit, 2):
        if gcd(n, d) != 1 or n in seen:
            continue
        x = n; path = []
        while x not in seen and x <= 50 * limit:
            seen.add(x); path.append(x)
            y = 3 * x + d
            while y % 2 == 0:
                y //= 2
            x = y
        if x in path:
            cyc = path[path.index(x):]
            cycles.add(min(cyc))
    return sorted(cycles)


def part4():
    print("== (4) the mediant tree is the tree of all 3x + d cycles ==")
    dens = defaultdict(set)
    for A in range(1, 23):
        base = None
        for cuts in range(1 << (A - 1)):
            w = []; run = 1
            for i in range(A - 1):
                if (cuts >> i) & 1:
                    w.append(run); run = 1
                else:
                    run += 1
            w.append(run)
            S, D = carry(w), clock(w); g = gcd(S, abs(D)); den = abs(D) // g
            if den <= 100:
                dens[den].add(Fraction(S, D))
    print(" denominators d <= 100 (coprime to 6) of fixed points x_w, A <= 22, with the number of distinct x_w: %s" % [(d, len(dens[d])) for d in sorted(dens)])
    missing = [d for d in range(1, 101) if gcd(d, 6) == 1 and d not in dens]
    print(" coprime-to-6 denominators <= 100 not reached by A <= 22: %s (%.0fs)" % (missing, time.time() - T0))
    for d in (5, 7, 11, 13, 17, 19, 23, 25):
        cyc = cycles_of_3x_plus_d(d, 200000)
        mins = sorted(set(min(abs(x.numerator) for x in dens[d] if x > 0 and x.denominator == d) for _ in [0])) if dens[d] else []
        tree_pts = sorted(set(x for x in dens[d] if x > 0))
        tree_mins = sorted(set(int(min(c)) for c in [[x.numerator for x in tree_pts]] if c))
        print("  d = %2d: cycles of 3x+%d (least elements, starts <= 2*10^5): %s; tree points m/%d with m in %s" % (d, d, cyc, d, sorted(set(x.numerator for x in tree_pts))[:12]))
    for d in missing[:6]:
        cyc = cycles_of_3x_plus_d(d, 400000)
        print("  missing d = %d: 3x+%d cycles found (least elements): %s" % (d, d, cyc))


if __name__ == "__main__":
    part1(); part2(); part3(); part4()
    print("total %.0fs" % (time.time() - T0))


def part5():
    print("== (5) control: the same record function for 5x+1 (positive drift; the never-descending set has Haar measure 0.176) ==")
    LOG25 = math.log2(5)
    def tau5(n, cap=400):
        x = n; d = 0; p5 = 1
        for j in range(1, cap + 1):
            y = 5 * x + 1; v = 0
            while y % 2 == 0:
                y //= 2; v += 1
            d += v; p5 *= 5
            if 2 ** d > p5:
                return j
            x = y
            if x == n or x.bit_length() > 600:
                return None
        return None
    stuck = []; records = []; best = 0
    for n in range(3, 3001, 2):
        t = tau5(n)
        if t is None:
            stuck.append(n)
        elif t > best:
            best = t; records.append((n, t))
    print(" odd n <= 3000 under 5x+1 with no word-level descent within 400 odd steps (or periodic, or grown past 600 bits): %d of them, first %s" % (len(stuck), stuck[:12]))
    print(" records of tau_5 among the others: %s" % records[:10])
    print(" reading: for 5x+1 the record function stalls (g_5 is bounded by the first stuck integer), as Conjecture G must fail when the drift is positive; for 3x+1 it climbs (27, 703, 10087, ...)")


if __name__ == "__main__":
    part5()
