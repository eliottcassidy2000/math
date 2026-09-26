#!/usr/bin/env python3
"""collatz_thin_20260926_movingbarrier.py -- controls for the o(X^(h*)) form of thin divergence
(session collatz-exponent-atlas-20260926, opus, 2026-09-26).

Walk: letters 0 (step -1) and 1 (step alpha - 1, alpha = log_2 3); S_i = o_i alpha - i; no ties (alpha irrational).
h = h(log_3 2), lambda* = log_2(rho*/(1-rho*))/alpha the tilt at zero drift.

(A) Weighted Spitzer identity (exact, rational letter weights x_0, x_1):
      P_n(x) = sum over words of length n with all prefix sums > 0 (resp. < 0) of x_0^(#0) x_1^(#1),
      B_n(x) = sum over words of length n with total > 0 (resp. < 0) of the same weight,
    then n P_n = sum_(m=1)^n B_m P_(n-m) (P_0 = 1). Checked with Fractions for n <= NMAX_ID, for the weights
    (1,1), (1/2, 2), (3, 1/3), both signs.
(B) Moving-barrier count: M_k(y) = #{ w in {0,1}^k : S_i(w) > -y for all 1 <= i <= k }, exact DP;
    the claim is M_k(y) <= C_eps 2^(hk) k^(-3/2) 2^(lambda* y) e^(eps y), and the observed behaviour is
    M_k(y) ~ (y+1) * const * 2^(hk) k^(-3/2) 2^(lambda* y) for y << sqrt(k). Printed: the normalised counts
    and their quotient by (y+1).
(C) The decomposition identity M_k(y) = W_k + sum_(m=1)^k N_m(y) W_(k-m) with N_m(y) = #{negative words of
    length m with total > -y} (all prefix sums < 0), checked exactly for k <= KDEC and y in {0.5, 1, 2, 4, 8}
    (y non-integer to avoid any question of ties at the barrier: S_i = -y is impossible anyway).
Usage: python3 collatz_thin_20260926_movingbarrier.py [KMAX=800] [NMAX_ID=40] [KDEC=60]
"""
import math, sys
from fractions import Fraction

ALPHA = math.log2(3)
RHO = math.log(2) / math.log(3)
H = -(RHO * math.log2(RHO) + (1 - RHO) * math.log2(1 - RHO))
LAM = math.log2(RHO / (1 - RHO)) / ALPHA


def gt(o, j, y):
    """exact test o*alpha - j > -y (y a Fraction). For o = 0 the value -j is an integer and can equal -y
    (all-zero prefixes at an integer barrier): compare exactly. For o >= 1 the value is irrational, so no tie
    with the rational -y is possible and a float comparison with a safety margin is exact (asserted)."""
    if o == 0:
        return j < y
    v = o * ALPHA - j + float(y)
    assert abs(v) > 1e-9, (o, j, y)
    return v > 0


def words_dp(k, y, sign=+1):
    """sign=+1: all prefix sums > -y; sign=-1: all prefix sums < -y (used with y = 0 for negative words)."""
    cur = {0: 1}
    for j in range(1, k + 1):
        nxt = {}
        for o, c in cur.items():
            for st in (0, 1):
                o2 = o + st
                ok = gt(o2, j, y) if sign > 0 else (not gt(o2, j, y))
                if ok:
                    nxt[o2] = nxt.get(o2, 0) + c
        cur = nxt
    return cur


def weighted_spitzer_check(nmax, x0, x1, sign):
    """exact: n P_n = sum B_m P_(n-m) with rational weights; P counts words with all prefix sums > 0 (sign +) or < 0 (sign -)."""
    P = [Fraction(1)]
    B = [Fraction(0)]
    for n in range(1, nmax + 1):
        # P_n: DP with weights
        cur = {0: Fraction(1)}
        for j in range(1, n + 1):
            nxt = {}
            for o, c in cur.items():
                for st, w in ((0, x0), (1, x1)):
                    o2 = o + st
                    v = o2 * ALPHA - j
                    if (v > 0) if sign > 0 else (v < 0):
                        nxt[o2] = nxt.get(o2, Fraction(0)) + c * w
            cur = nxt
        P.append(sum(cur.values(), Fraction(0)))
        # B_n
        tot = Fraction(0)
        for o in range(n + 1):
            v = o * ALPHA - n
            if (v > 0) if sign > 0 else (v < 0):
                tot += math.comb(n, o) * x0 ** (n - o) * x1 ** o
        B.append(tot)
    ok = all(n * P[n] == sum(B[m] * P[n - m] for m in range(1, n + 1)) for n in range(1, nmax + 1))
    return ok, P[:8]


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 800
    NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    KDEC = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    print("h = %.7f, lambda* = %.6f, lambda*/h - 3/2 = %.6f" % (H, LAM, LAM / H - 1.5))
    print("(A) weighted Spitzer identity n P_n = sum_m B_m P_(n-m), exact rationals, n <= %d:" % NMAX)
    for (x0, x1) in ((Fraction(1), Fraction(1)), (Fraction(1, 2), Fraction(2)), (Fraction(3), Fraction(1, 3))):
        for sign in (+1, -1):
            ok, head = weighted_spitzer_check(NMAX, x0, x1, sign)
            print("    weights (x0, x1) = (%s, %s), %s words: %s; P_1..P_7 = %s" % (x0, x1, "positive" if sign > 0 else "negative", ok, [str(v) for v in head[1:]]))
    print("(B) moving barrier: M_k(y) / (2^(hk) k^(-3/2) 2^(lambda* y)) and the same divided by (y+1)")
    ys = [Fraction(0), Fraction(1), Fraction(2), Fraction(4), Fraction(8), Fraction(16), Fraction(32)]
    print("    y:       " + "  ".join("%7s" % str(y) for y in ys))
    for k in (100, 200, 400, 800, 1200, 1600):
        if k > KMAX:
            break
        row = []
        for y in ys:
            M = sum(words_dp(k, y).values())
            row.append(2 ** (math.log2(M) - H * k + 1.5 * math.log2(k) - LAM * float(y)))
        print("    k=%5d: " % k + "  ".join("%7.3f" % r for r in row) + "  | /(y+1): " + "  ".join("%6.3f" % (r / (float(y) + 1)) for r, y in zip(row, ys)))
    print("(C) decomposition M_k(y) = W_k + sum_m N_m(y) W_(k-m), exact, k <= %d:" % KDEC)
    W = [1] + [sum(words_dp(k, Fraction(0)).values()) for k in range(1, KDEC + 1)]
    allok = True
    for y in (Fraction(1, 2), Fraction(1), Fraction(2), Fraction(4), Fraction(8)):
        N = [0]
        for m in range(1, KDEC + 1):
            neg = words_dp(m, Fraction(0), sign=-1)          # all prefix sums < 0
            N.append(sum(c for o, c in neg.items() if gt(o, m, y)))   # total > -y
        for k in range(1, KDEC + 1):
            M = sum(words_dp(k, y).values())
            rhs = W[k] + sum(N[m] * W[k - m] for m in range(1, k + 1))
            if M != rhs:
                allok = False
                print("    MISMATCH at y=%s k=%d: %d vs %d" % (y, k, M, rhs))
    print("    identity holds for all tested y and k: %s" % allok)


if __name__ == '__main__':
    main()
