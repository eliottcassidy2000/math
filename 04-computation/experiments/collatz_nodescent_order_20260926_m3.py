#!/usr/bin/env python3
"""collatz_nodescent_order_20260926_m3.py -- the Spitzer identity for a three-letter Conway map
(session collatz-exponent-atlas-20260926, opus, 2026-09-26; remark to THM-4495 / Theorem 4 of THM-4487's note).

Map g = x/3, (2x+1)/3, (4x+1)/3 (multipliers 1/3, 2/3, 4/3; letters 0, 1, 2 with log_3-weights
-1, log_3 2 - 1, 2 log_3 2 - 1). A word's partial sum after j letters is s_j log_3 2 - j with s_j = #1 + 2 #2
among the first j letters; no nonempty segment has zero sum (log_3 2 irrational), so THM-4495's proof of
Theorem A applies verbatim:  k W_k = sum_(n=1)^k B_n W_(k-n),  W_k = #{words of length k with all partial
sums > 0} (no k-step descent for large n), B_n = #{words of length n with positive total} = sum over
(#1, #2) with (#1 + 2 #2) log_3 2 > n of the trinomial coefficient. Exact tests: 2^s > 3^j.
Also prints W_k k^(3/2) 3^(-k E) with E = E_g(1) = log_3 min_lambda (1/3)(1 + 2^lambda + 4^lambda)... i.e. the
constrained max entropy at gamma = 1 (Theorem 4), to show the same k^(-3/2) order.
Usage: python3 collatz_nodescent_order_20260926_m3.py [KMAX=400] [KDP=120]
"""
import math, sys
from math import comb

LOG3_2 = math.log(2) / math.log(3)


def positive(s, j):
    return 2 ** s > 3 ** j


def W_dp(K):
    W = [1]
    cur = {0: 1}          # s = #1 + 2#2 -> count
    for j in range(1, K + 1):
        nxt = {}
        for s, c in cur.items():
            for step in (0, 1, 2):
                s2 = s + step
                if positive(s2, j):
                    nxt[s2] = nxt.get(s2, 0) + c
        cur = nxt
        W.append(sum(cur.values()))
    return W


def B_list(K):
    B = [0]
    for n in range(1, K + 1):
        tot = 0
        for n2 in range(n + 1):
            for n1 in range(n - n2 + 1):
                if positive(n1 + 2 * n2, n):
                    tot += comb(n, n2) * comb(n - n2, n1)
        B.append(tot)
    return B


def E_gamma1():
    # E_g(1) = log_3 min_lambda (1/3) sum a_i^lambda + 1 = log_3 min_lambda (1 + 2^lambda + 4^lambda) - 1 + 1 ... :
    # Z(lambda) = sum a_i^lambda with a = (1/3, 2/3, 4/3); E_g(1) = log_3 min Z.
    best = 10
    for k in range(0, 40001):
        lam = 5 * k / 40000
        Z = (1 / 3) ** lam + (2 / 3) ** lam + (4 / 3) ** lam
        best = min(best, Z)
    return math.log(best) / math.log(3)


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    KDP = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    B = B_list(KMAX)
    Wd = W_dp(KDP)
    Wr = [1]
    for k in range(1, KMAX + 1):
        s = sum(B[n] * Wr[k - n] for n in range(1, k + 1))
        assert s % k == 0, k
        Wr.append(s // k)
    ok = all(Wd[k] == Wr[k] for k in range(KDP + 1))
    E = E_gamma1()
    print("m = 3 map x/3, (2x+1)/3, (4x+1)/3: Spitzer identity k W_k = sum B_n W_(k-n): DP == recurrence for k <= %d: %s; integral to k = %d" % (KDP, ok, KMAX))
    print("   W_1..W_15:", Wr[1:16])
    print("   B_1..B_15:", B[1:16])
    print("   E_g(1) = log_3 min_lambda Z(lambda) = %.6f (Theorem 4; thin-divergence exponent 1 - I(g), I(g) = %.4f)" % (E, 1 - E))
    print("   W_k k^1.5 3^(-kE):", "  ".join("k=%d: %.3f" % (k, 3 ** (math.log(Wr[k]) / math.log(3) - E * k + 1.5 * math.log(k) / math.log(3))) for k in (10, 20, 50, 100, 200, 300, 400) if k <= KMAX))


if __name__ == '__main__':
    main()
