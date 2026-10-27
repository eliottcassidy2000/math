"""Independent census to 2^28 (fable audit of opus S64).  Own code; int16 tau.

Adds to audit_tau_census.py: block stats k = 25..28, the MEDIAN and SKEWNESS of z per block
(the note reports only mean and sd), exact G_k(C) at the note's C grid, exact cumulative
medians, quarter-dyadic floor densities, 4-adic and 3-adic block ratios, two-block exponents,
max tau and max orbit value; and the finite-X evaluation of the independent-3-adic word model
S(X, C) = sum_d C(k_d, d) 3^-d at the same scales, for comparison with the measured exponents.
"""
import math
import sys
import time
from math import comb, log

import numpy as np
from numba import njit

K_MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 28
N = 1 << K_MAX
C_STAR = 3.0 / math.log(4.0 / 3.0)
SIGMA = math.sqrt(2.0 * (1.0 - C_STAR * math.log(2.0)) ** 2 / (2.0 * math.log(2.0) - math.log(3.0)))
ln2, ln3 = math.log(2.0), math.log(3.0)
out = []


def p(*a):
    line = " ".join(str(x) for x in a)
    print(line)
    out.append(line)


@njit(cache=True)
def census16(N):
    tau = np.zeros(N, dtype=np.int16)
    peak = np.int64(0)
    tmax = 0
    nmax = 0
    for n in range(2, N):
        x = np.int64(n)
        steps = 0
        while x >= n:
            if x & 1:
                x = 3 * x + 1
                if x > peak:
                    peak = x
            else:
                x >>= 1
            steps += 1
        t = steps + tau[x]
        tau[n] = t
        if t > tmax:
            tmax = t
            nmax = n
    return tau, peak, tmax, nmax


@njit(cache=True)
def range_stats(tau, lo, hi, cstar, sigma, c_list, floor_c):
    cnt = hi - lo
    s = 0.0
    s2 = 0.0
    s3 = 0.0
    dens = np.zeros(len(c_list), dtype=np.int64)
    fl = 0
    zs = np.empty(cnt, dtype=np.float32)
    rs = np.empty(cnt, dtype=np.float32)
    for i in range(cnt):
        n = lo + i
        t = tau[n]
        L = math.log(n)
        z = (t - cstar * L) / (sigma * math.sqrt(L))
        zs[i] = z
        rs[i] = t / L
        s += z
        s2 += z * z
        s3 += z * z * z
        for j in range(len(c_list)):
            if t <= c_list[j] * L:
                dens[j] += 1
        if t < floor_c * L:
            fl += 1
    mean = s / cnt
    var = s2 / cnt - mean * mean
    sd = math.sqrt(var)
    skew = (s3 / cnt - 3 * mean * var - mean ** 3) / sd ** 3
    zmed = np.median(zs)
    rmed = np.median(rs)
    return cnt, mean, sd, skew, zmed, rmed, dens, fl


def main():
    p("# audit_tau_census_28.py (fable), N = 2^%d" % K_MAX)
    t0 = time.time()
    tau, peak, tmax, nmax = census16(N)
    p("census in %.1f s; max tau = %d at n = %d; max orbit value = %d = 2^%.2f" % (time.time() - t0, tmax, nmax, peak, math.log2(peak)))
    c_list = np.array([3.307, 6.0, 8.0, 9.0, 10.0, C_STAR, 523 / 50, 262 / 25, 11.0, 12.0, 13.0, 14.0, 16.0, 20.0, 30.0])
    p("")
    p("## dyadic blocks [2^(k-1), 2^k): exact G_k(C) (tau <= C ln n), z mean/sd/skew/median, median of tau/ln n, floor count (< 3.307 ln n)")
    p("k  | " + " ".join("%7.3f" % c for c in c_list) + " | zmean  zsd   zskew  zmed  | r_med | floor")
    floors = {}
    for k in range(16, K_MAX + 1):
        cnt, mean, sd, skew, zmed, rmed, dens, fl = range_stats(tau, 1 << (k - 1), 1 << k, C_STAR, SIGMA, c_list, 3.307)
        floors[k] = fl
        p("%2d | " % k + " ".join("%7.4f" % (d / cnt) for d in dens) + " | %+.3f %.3f %+.3f %+.3f | %.3f | %d" % (mean, sd, skew, zmed, rmed, fl))
    p("Gaussian model for the block: zskew = 0, zmed = zmean.  Model median of tau/ln n = C* + zmean*sigma/sqrt(ln N_k).")
    p("")
    p("## floor-count ratios and two-block exponents log_4(count_k / count_(k-2))")
    for k in range(18, K_MAX + 1):
        r1 = floors[k] / floors[k - 1]
        r2 = floors[k] / floors[k - 2]
        p("k = %2d: count %7d  ratio %.3f  two-block exponent %.3f" % (k, floors[k], r1, math.log(r2, 4)))
    ex = [math.log(floors[k] / floors[k - 2], 4) for k in range(21, K_MAX + 1)]
    p("two-block exponents k = 21..%d: mean %.3f, sd %.3f" % (K_MAX, float(np.mean(ex)), float(np.std(ex))))
    p("")
    p("## cumulative F_X(523/50) = #{1 <= n < X : tau <= (523/50) ln n}/(X-1), exact; exact cumulative median of tau/ln n (n >= 2)")
    allr = tau[2:].astype(np.float32) / np.log(np.arange(2, N, dtype=np.float32))
    for k in range(16, K_MAX + 1):
        X = 1 << k
        seg = allr[: X - 2]
        cnt_fast = int(np.count_nonzero(seg <= 523 / 50)) + 1
        p("k = %2d: F = %.4f ; cumulative median tau/ln n = %.4f" % (k, cnt_fast / (X - 1), float(np.median(seg))))
    del allr
    p("")
    p("## quarter-dyadic floor densities x 1e4 (the note's 5.3 numbers)")
    for k in range(20, K_MAX + 1):
        base = 1 << (k - 1)
        row = []
        for i in range(4):
            lo = base + (base * i) // 4
            hi = base + (base * (i + 1)) // 4
            _, _, _, _, _, _, _, fl = range_stats(tau, lo, hi, C_STAR, SIGMA, c_list[:1], 3.307)
            row.append(1e4 * fl / (hi - lo))
        p("k = %2d: %s  (max/min = %.2f)" % (k, ", ".join("%.1f" % v for v in row), max(row) / min(row)))
    p("")
    p("## blocks of ratio 4 and 3: floor counts and ratios")
    prev = None
    j = 8
    while 4 ** (j + 1) <= N:
        _, _, _, _, _, _, _, fl = range_stats(tau, 4 ** j, 4 ** (j + 1), C_STAR, SIGMA, c_list[:1], 3.307)
        p("4-adic j = %2d: floor %7d %s" % (j, fl, "" if prev is None else "ratio %.3f exponent log_4 %.3f" % (fl / prev, math.log(fl / prev, 4))))
        prev = fl
        j += 1
    prev = None
    j = 10
    while 3 ** (j + 1) <= N:
        _, _, _, _, _, _, _, fl = range_stats(tau, 3 ** j, 3 ** (j + 1), C_STAR, SIGMA, c_list[:1], 3.307)
        p("3-adic j = %2d: floor %7d %s" % (j, fl, "" if prev is None else "ratio %.3f exponent log_3 %.3f" % (fl / prev, math.log(fl / prev, 3))))
        prev = fl
        j += 1

    # finite-X word model at the same scales
    p("")
    p("## independent-3-adic word model at finite X: S(X,C) = sum_d C(k_d,d) 3^-d, k_d = round((ln X + d ln3)/ln2), k_d + d <= C ln X; X = 2^k")
    for C in (3.307,):
        prevS = {}
        for k in range(16, 41):
            lnX = k * ln2
            S = 0.0
            for d in range(0, 2 * k):
                kk = round((lnX + d * ln3) / ln2)
                if kk + d <= C * lnX:
                    S += comb(kk, d) * 3.0 ** (-d)
            prevS[k] = S
            ex2 = "" if k - 2 not in prevS else "two-block exponent %.3f" % (math.log(S / prevS[k - 2], 4))
            p("C = %.3f, X = 2^%2d: S = %.4e  ln S/ln X = %.4f  %s" % (C, k, S, math.log(S) / lnX, ex2))

    with open("audit_tau_census_28.out", "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
