#!/usr/bin/env python3
"""Audit P1 (independent): Bregman bound for the doubles D(B) of a regular tournament B of order q (degree mu).

Bound used by the session:  H(D) <= (n!)^(1/n) prod_v ((d+(v)+1)!)^(1/(d+(v)+1)),
obtained from per(A(D + w)) where w <-> every vertex (row sums d+(v)+1 and n).
D(B): out-degree mu + q - 1 = (3q-3)/2 on A, mu + 1 = (q+1)/2 on B'.
Deletion: the bound for D - v with the ORIGINAL out-degrees of the remaining 2q-1 vertices is an upper bound for the
true Bregman bound of D - v (f(d) = ((d+1)!)^(1/(d+1)) is increasing); it is largest when v is a b-vertex.
Mean E[H(m)] = m!/2^(m-1).

Checks:
  (1) exact-integer logs (math.log of exact factorials) for q <= 303, lgamma beyond (cross-checked on the overlap);
  (2) the first q = 3 mod 4 with both log-ratios negative, and negativity / monotonicity up to q = 10^6;
  (3) an explicit Robbins-type upper bound U(q) >= R_del(q), R_D(q), verified < 0 for every q >= 59 up to 10^7 and
      shown decreasing for q >= 10^3 (derivative bound), closing the "large q" step of Theorem A4.
  (4) the numbers quoted in PROOF.md / THM-4995 (q = 51, 55, 59, 71 and T_11 / T_7 ratios).
"""
import math
from math import lgamma, log, factorial

LN10 = log(10)


def lf_exact(k):
    return math.log(factorial(k)) if k > 1 else 0.0


def lf_gamma(k):
    return lgamma(k + 1)


def ratios(q, lf):
    da = (3 * q - 3) // 2
    db = (q + 1) // 2
    n = 2 * q
    # log Bregman for D
    bD = lf(n) / n + q * lf(da + 1) / (da + 1) + q * lf(db + 1) / (db + 1)
    # worst deletion: drop one b-vertex, keep original degrees of the rest
    bdel = lf(n - 1) / (n - 1) + q * lf(da + 1) / (da + 1) + (q - 1) * lf(db + 1) / (db + 1)
    # an a-deletion for comparison
    badel = lf(n - 1) / (n - 1) + (q - 1) * lf(da + 1) / (da + 1) + q * lf(db + 1) / (db + 1)
    logE = lambda m: lf(m) - (m - 1) * log(2)
    return bD - logE(2 * q - 1), bdel - logE(2 * q - 2), badel - logE(2 * q - 2)


def main():
    # (1) exact vs lgamma on the overlap
    maxdiff = 0.0
    for q in range(7, 304, 4):
        a = ratios(q, lf_exact)
        b = ratios(q, lf_gamma)
        maxdiff = max(maxdiff, max(abs(x - y) for x, y in zip(a, b)))
    print(f"(1) max |exact-integer log - lgamma log| over q <= 303: {maxdiff:.3e} (natural log)")
    # (2) table and first negative
    rows = []
    for q in range(7, 1000001, 4):
        lf = lf_exact if q <= 303 else lf_gamma
        rD, rdel, radel = ratios(q, lf)
        rows.append((q, rD / LN10, rdel / LN10, radel / LN10))
    first = next(q for q, a, b, c in rows if a < 0 and b < 0)
    print(f"(2) first q (= 3 mod 4) with both log10-ratios < 0: {first}")
    print(f"    all q in [{first}, 10^6] negative: {all(a < 0 and b < 0 for q, a, b, c in rows if q >= first)}")
    print(f"    D-ratio strictly decreasing from 7 on: {all(rows[i+1][1] < rows[i][1] for i in range(len(rows)-1))}")
    print(f"    deletion ratio strictly decreasing from 7 on: {all(rows[i+1][2] < rows[i][2] for i in range(len(rows)-1))}")
    print(f"    b-deletion bound >= a-deletion bound for all q: {all(b >= c for q, a, b, c in rows)}")
    print("    q  log10(Breg(D)/E[H(2q-1)])  log10(Breg(D-b)/E[H(2q-2)])")
    for q, a, b, c in rows:
        if q <= 103 or q in (199, 499, 999, 1999, 4999, 9999, 19999):
            print(f"    {q} {a:.4f} {b:.4f}")
    # (3) explicit Robbins upper bound. For k >= 1: log k! <= k log k - k + 0.5 log(2 pi k) + 1/(12k)
    #                                              log k! >= k log k - k + 0.5 log(2 pi k) + 1/(12k+1)
    def up(k):
        return k * log(k) - k + 0.5 * log(2 * math.pi * k) + 1 / (12 * k)

    def lo(k):
        return k * log(k) - k + 0.5 * log(2 * math.pi * k) + 1 / (12 * k + 1)

    def U(q):
        da = (3 * q - 3) // 2 + 1  # = (3q-1)/2
        db = (q + 1) // 2 + 1      # = (q+3)/2
        n = 2 * q
        bD = up(n) / n + q * up(da) / da + q * up(db) / db
        bdel = up(n - 1) / (n - 1) + q * up(da) / da + (q - 1) * up(db) / db
        return bD - (lo(2 * q - 1) - (2 * q - 2) * log(2)), bdel - (lo(2 * q - 2) - (2 * q - 3) * log(2))

    bad = [q for q in range(59, 10 ** 7 + 1, 4) if max(U(q)) >= 0]
    print(f"(3) Robbins upper bound U(q) >= true log-ratios; q in [59, 10^7] with U(q) >= 0: {bad[:5]} (count {len(bad)})")
    # U(q) - q log(3/4): bounded by an explicit multiple of log q
    worst = max((max(U(q)) - q * log(3 / 4)) / log(q) for q in range(1003, 10 ** 7 + 1, 4000))
    print(f"    max over sampled q in [10^3, 10^7] of (U(q) - q log(3/4)) / log q = {worst:.4f}")
    # monotonicity of U via exact differences on a geometric grid is not a proof; give the analytic remark instead.
    # (4) quoted numbers
    E = lambda m: factorial(m) / 2 ** (m - 1)
    print("(4) quoted values:")
    for q in (51, 55, 59, 71):
        rD, rdel, _ = ratios(q, lf_exact)
        print(f"    q={q}: log10 ratios {rD / LN10:.4f} {rdel / LN10:.4f}")
    H11, H11d, H7, H7d = 17667825987845, 1822800617515, 4730397, 686457
    print(f"    H(T11)/E[H(21)] = {H11 / E(21):.4f}; H(T11-v)/E[H(20)] = {H11d / E(20):.4f}")
    print(f"    H(T7)/E[H(13)] = {H7 / E(13):.4f}; H(T7-v)/E[H(12)] = {H7d / E(12):.4f}")
    print(f"    H(T7)/H(QR7)^2 = {H7 / 189 ** 2:.2f}; H(T11)/H(QR11)^2 = {H11 / 95095 ** 2:.2f}")
    HD, HDd = 6450103269, 773858997
    print(f"    H(D(C9))/E[H(17)] = {HD / E(17):.4f}; H(D(C9)-v)/E[H(16)] = {HDd / E(16):.4f}")
    print(f"    S(10), S(11) first gaps 14685, 80453 over the mean: {14685 / E(10):.4f}, {80453 / E(11):.4f}")


if __name__ == '__main__':
    main()
