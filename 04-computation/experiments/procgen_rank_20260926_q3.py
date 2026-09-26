#!/usr/bin/env python3
"""
procgen_rank_20260926_q3.py -- Q3 of the rank lane: heights and the product formula (short, typed).

  (a) shadow identity: inside the shadow of a periodic x, y_j = T^j n - x_j = (3^(a_j)/2^j) y_0 exactly, so
      the numerator loses j factors of 2, gains a_j factors of 3, keeps its prime-to-6 part; the local
      contributions at {infinity, 2, 3} of log|y_j| - log|y_0| sum to zero (product formula for 3^(a_j)/2^j)
  (b) the forced charge is a ratio of local expansion rates: chi = log 2 * log|lam|_inf / log|lam|_2
  (c) the backward tree of -1 is binary (2^(j-1) points at depth j >= 1) with log2 H(z) = O(j): the forced
      charges kappa on it have infinite mass and infinite height moment
"""
import math
from fractions import Fraction

from procgen_rank_20260926_lib import (KAPPA, LN2, LN3, check, say, v2, height, T, cycle_points, preimages,
                                       good_integer, chi_of_word)


def vp_int(n, p):
    n = abs(n)
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def prime_to_6(n):
    n = abs(n)
    while n % 2 == 0:
        n //= 2
    while n % 3 == 0:
        n //= 3
    return n


def part_a():
    say("== Q3.a the shadow identity: 2-adic loss + 3-adic gain = archimedean growth ==")
    words = [[1], [1, 1, 0], [1, 1, 1, 1, 0], [1, 1, 1, 1, 1, 1, 0]]
    x = -17
    w17 = []
    for _ in range(11):
        w17.append(x & 1)
        x = T(x)
    words.append(w17)
    nchk = 0
    for w in words:
        pts = cycle_points(w)
        for D in (40, 100, 200):
            n = good_integer(pts[0], D)
            y0 = Fraction(n) - pts[0]
            N0 = y0.numerator
            m = n
            a = 0
            for j in range(D - 1):
                xj = pts[j % len(w)]
                y = Fraction(m) - xj
                if y != Fraction(3 ** a, 2 ** j) * y0:
                    check(False, "shadow identity fails")
                N = y.numerator
                ok = (v2(y) == D - j and vp_int(N, 3) == vp_int(N0, 3) + a and prime_to_6(N) == prime_to_6(N0)
                      and y.denominator == y0.denominator)
                if not ok:
                    check(False, "valuation bookkeeping fails at word %s D %d j %d" % (w, D, j))
                # product formula over {inf, 2, 3} for the ratio y/y0 = 3^a / 2^j
                s = (math.log(abs(y)) - math.log(abs(y0))) - (v2(y) - v2(y0)) * LN2 - (vp_int(N, 3) - vp_int(N0, 3)) * LN3
                if abs(s) > 1e-9:
                    check(False, "product formula fails")
                a += m & 1
                m = T(m)
                nchk += 1
    check(nchk > 1000, "in %d shadow steps (cycles -1, -5, O_5, O_7, -17): y_j = (3^a_j/2^j) y_0, v2 drops by j, "
          "v3 of the numerator rises by a_j, the prime-to-6 part and the denominator are invariant, and the "
          "local log-changes at infinity, 2, 3 sum to 0" % nchk)


def part_b():
    say("== Q3.b the forced charge as a ratio of local expansion rates ==")
    worst = 0.0
    cnt = 0
    for p in range(1, 13):
        for mm in range(1 << p):
            w = [(mm >> i) & 1 for i in range(p)]
            a = sum(w)
            if 3 ** a <= 2 ** p:
                continue
            lam_inf = a * LN3 - p * LN2          # log|lambda|_inf
            lam_2 = p * LN2                      # log|lambda|_2 (|3^a/2^p|_2 = 2^p)
            lam_3 = -a * LN3                     # log|lambda|_3
            worst = max(worst, abs(lam_inf + lam_2 + lam_3), abs(chi_of_word(w) - LN2 * lam_inf / lam_2))
            cnt += 1
    check(worst < 1e-12, "for all %d expanding words of length <= 12: log|lam|_inf + log|lam|_2 + log|lam|_3 = 0 and "
          "chi = log 2 * log|lam|_inf / log|lam|_2" % cnt)


def part_c():
    say("== Q3.c the backward tree of -1: binary, heights linear in the depth ==")
    level = [Fraction(-2)]
    counts = []
    maxrate = 0.0
    total_forced = KAPPA
    for j in range(1, 15):
        counts.append(len(level))
        for z in level:
            if j >= 3:
                maxrate = max(maxrate, math.log2(height(z)) / j)
        total_forced += KAPPA * len(level)
        nxt = []
        for z in level:
            nxt.extend(preimages(z))
        level = nxt
    check(counts == [2 ** (j - 1) for j in range(1, 15)],
          "the backward tree of -1 has exactly 2^(j-1) points at depth j = 1..14 (every point of Z_(2) has the two "
          "preimages 2x and (2x-1)/3; only -1 is its own)")
    check(maxrate < math.log2(6), "log2 H(z) <= %.3f j for every z at depth 3 <= j <= 14 (< log2 6)" % maxrate)
    check(total_forced > 5000, "forced charge kappa on the tree up to depth 14 already totals %.0f: the forced "
          "charges on one backward orbit are not summable" % total_forced)


def run():
    part_a()
    part_b()
    part_c()


if __name__ == "__main__":
    run()
