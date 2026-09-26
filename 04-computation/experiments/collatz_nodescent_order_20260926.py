#!/usr/bin/env python3
"""collatz_nodescent_order_20260926.py -- the exact order of the no-descent residue count
(session collatz-exponent-atlas-20260926, opus, 2026-09-26).

Objects (Collatz shortcut T, plus sheet; the minus sheet is identical by negation):
  W_k = #{ residue words w in {0,1}^k : 3^(o_j) > 2^j for every 1 <= j <= k }   (o_j = ones among the first j letters)
      = |Bad_k| of THM-4479 = residues mod 2^k with no k-step descent (M_j(r) = 3^(o_j)/2^j > 1 for all j <= k);
  B_n = #{ words of length n with 3^(o_n) > 2^n } = sum_{j : 3^j > 2^n} C(n, j);
  N_k = # binary necklaces of length k with more than k log_3 2 ones (THM-4479's lower bound for delta_k).

Claims checked exactly:
  (S) Spitzer's combinatorial identity  k W_k = sum_{n=1}^k B_n W_{k-n}  (W_0 = 1), i.e.
      sum_k W_k t^k = exp( sum_n B_n t^n / n ): compared with a direct ballot DP for k <= KDP.
  (B) binomial bounds 2^{n h(p)}/sqrt(8 n p (1-p)) <= C(n, pn) <= 2^{n h(p)}/sqrt(2 pi n p (1-p)),
      the geometric tail B_n <= (rho/(2 rho - 1)) C(n, j0), j0 = floor(n log_3 2) + 1, and h(j0/n) <= h(rho).
  (O) the orders: W_k k^{3/2} 2^{-hk}, N_k k^{3/2} 2^{-hk}, B_k k^{1/2} 2^{-hk} stay in bounded windows
      away from 0 and infinity for k up to KMAX (the windows oscillate with the fractional part of k log_3 2),
      and the sandwich N_k <= delta_k <= W_k of THM-4479 therefore pins delta_k to constants.
Usage: python3 collatz_nodescent_order_20260926.py [KMAX=3000] [KDP=300]
"""
import math, sys
from math import comb, gcd

LOG3_2 = math.log(2) / math.log(3)
RHO = LOG3_2
H = -(RHO * math.log2(RHO) + (1 - RHO) * math.log2(1 - RHO))


def positive(o, j):
    return 3 ** o > 2 ** j


def W_dp(K):
    """ballot DP: W[k] for k <= K, exact."""
    W = [1]
    cur = {0: 1}   # ones -> count, all prefixes positive, current length 0
    pow3 = [3 ** o for o in range(K + 2)]
    for j in range(1, K + 1):
        nxt = {}
        p2 = 2 ** j
        for o, c in cur.items():
            for step in (0, 1):
                o2 = o + step
                if pow3[o2] > p2:
                    nxt[o2] = nxt.get(o2, 0) + c
        cur = nxt
        W.append(sum(cur.values()))
    return W


def B_list(K):
    B = [0]
    for n in range(1, K + 1):
        p2 = 2 ** n
        j0 = 0
        while 3 ** j0 <= p2:
            j0 += 1
        B.append(sum(comb(n, j) for j in range(j0, n + 1)))
    return B


def W_rec(K, B):
    W = [1]
    for k in range(1, K + 1):
        s = sum(B[n] * W[k - n] for n in range(1, k + 1))
        assert s % k == 0, ("Spitzer recurrence not integral at k =", k)
        W.append(s // k)
    return W


def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def necklaces_above(k):
    """number of binary necklaces of length k with j ones, summed over j with 3^j > 2^k (Burnside)."""
    p2 = 2 ** k
    total = 0
    for j in range(k + 1):
        if 3 ** j > p2:
            g = gcd(k, j)
            s = 0
            for d in range(1, g + 1):
                if g % d == 0:
                    s += phi(d) * comb(k // d, j // d)
            assert s % k == 0
            total += s // k
    return total


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    KDP = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    print("rho = log_3 2 = %.10f, h = h(rho) = %.10f, 1 - h = %.10f, lambda* = log_2(rho/(1-rho))/log_2 3 = %.6f"
          % (RHO, H, 1 - H, math.log2(RHO / (1 - RHO)) / math.log2(3)))
    B = B_list(KMAX)
    Wd = W_dp(KDP)
    Wr = W_rec(KMAX, B)
    ok = all(Wd[k] == Wr[k] for k in range(KDP + 1))
    print("(S) Spitzer identity k W_k = sum_n B_n W_(k-n): ballot DP == recurrence for all k <= %d: %s; recurrence integral for all k <= %d" % (KDP, ok, KMAX))
    print("    W_k for k = 1..20:", [Wr[k] for k in range(1, 21)])
    print("    B_n for n = 1..20:", [B[n] for n in range(1, 21)])
    # (B) binomial bounds and the geometric tail
    worst_low, worst_up, worst_tail, worst_h = 1e9, 0.0, 0.0, -1e9
    for n in range(3, KMAX + 1):        # n = 1, 2 have j0 = n (p = 1): excluded from the binomial bounds
        p2 = 2 ** n
        j0 = 0
        while 3 ** j0 <= p2:
            j0 += 1
        p = j0 / n
        hp = -(p * math.log2(p) + (1 - p) * math.log2(1 - p))
        c = comb(n, j0)
        lc = math.log2(c)
        low = n * hp - 0.5 * math.log2(8 * n * p * (1 - p))
        up = n * hp - 0.5 * math.log2(2 * math.pi * n * p * (1 - p))
        worst_low = min(worst_low, lc - low)
        worst_up = max(worst_up, lc - up)
        worst_tail = max(worst_tail, B[n] / (c * RHO / (2 * RHO - 1)))
        worst_h = max(worst_h, hp - H)
    print("(B) for 3 <= n <= %d and j0 = floor(n log_3 2) + 1: min(log2 C(n,j0) - lower bound) = %.4f (>= 0), max(log2 C - upper bound) = %.4f (<= 0), "
          "max B_n / ((rho/(2rho-1)) C(n,j0)) = %.6f (<= 1), max h(j0/n) - h = %.2e (<= 0)" % (KMAX, worst_low, worst_up, worst_tail, worst_h))
    # (O) orders
    print("(O) k, W_k k^1.5/2^(hk), N_k k^1.5/2^(hk), B_k k^0.5/2^(hk), W_k/N_k, frac(k log_3 2)")
    Ns = {}
    for k in [5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 700, 1000, 1500, 2000, 2500, 3000]:
        if k > KMAX:
            break
        Nk = necklaces_above(k)
        Ns[k] = Nk
        rw = 2 ** (math.log2(Wr[k]) - H * k + 1.5 * math.log2(k))
        rn = 2 ** (math.log2(Nk) - H * k + 1.5 * math.log2(k))
        rb = 2 ** (math.log2(B[k]) - H * k + 0.5 * math.log2(k))
        print("    %5d  %8.4f  %8.4f  %8.4f  %8.3f  %.4f" % (k, rw, rn, rb, Wr[k] / Nk, (k * RHO) % 1))
    lo, hi = 1e9, 0.0
    lob, hib = 1e9, 0.0
    for k in range(max(2, KMAX // 6), KMAX + 1):
        rw = 2 ** (math.log2(Wr[k]) - H * k + 1.5 * math.log2(k))
        lo, hi = min(lo, rw), max(hi, rw)
        rb = 2 ** (math.log2(B[k]) - H * k + 0.5 * math.log2(k))
        lob, hib = min(lob, rb), max(hib, rb)
    print("    over %d <= k <= %d: W_k k^1.5/2^(hk) in [%.4f, %.4f];  B_k k^0.5/2^(hk) in [%.4f, %.4f]" % (max(2, KMAX // 6), KMAX, lo, hi, lob, hib))
    lo2, hi2 = 1e9, 0.0
    for k in range(max(2, KMAX // 6), KMAX // 6 + 60):
        Nk = necklaces_above(k)
        rn = 2 ** (math.log2(Nk) - H * k + 1.5 * math.log2(k))
        lo2, hi2 = min(lo2, rn), max(hi2, rn)
    print("    over %d <= k < %d: N_k k^1.5/2^(hk) in [%.4f, %.4f]" % (max(2, KMAX // 6), KMAX // 6 + 60, lo2, hi2))
    # sigma and the explicit constant of the elementary upper bound
    b = [0.0] + [2 ** (math.log2(B[n]) - H * n) for n in range(1, KMAX + 1)]
    sigma = sum(b[n] / n for n in range(1, KMAX + 1))
    tail = 2 * 1.99 / math.sqrt(KMAX)      # sum_{n > KMAX} 1.99 n^{-3/2} <= 2 * 1.99 / sqrt(KMAX)
    Cg = max(b[n] * math.sqrt(n) for n in range(1, KMAX + 1))
    print("    b_n sqrt(n) = B_n n^0.5 2^(-hn) <= %.4f for n <= %d (the elementary bound is 2.41/sqrt(2 pi rho(1-rho)) = %.4f); "
          "sigma = sum b_n/n = %.6f (+ tail <= %.4f); explicit upper constant C e^(2 sqrt2 sigma) = %.2f"
          % (Cg, KMAX, 2.4094 / math.sqrt(2 * math.pi * RHO * (1 - RHO)), sigma, tail, 1.99 * math.exp(2 * math.sqrt(2) * (sigma + tail))))
    # cumulative sums for the dip count D(2^T, 1): sum_{t < T} W_t
    print("    sum_(t<T) W_t for T = 8, 12, 16, 20, 24:", [sum(Wr[t] for t in range(1, T)) for T in (8, 12, 16, 20, 24)])


if __name__ == '__main__':
    main()
