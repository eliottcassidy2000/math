"""Independent audit census (fable, 2026-10-27, audit of opus S64).

Ordinary total stopping time tau(n): n -> n/2 (even), n -> 3n+1 (odd), every
application counted, tau(1) = 0.  Computed by a direct loop with memoisation
by first descent for all n < 2^24.  Written BEFORE reading the session's
scripts, so every number here is independent of them.

Outputs (all my own definitions, both <= and < variants where ambiguous):
  F_X(C)      = #{1 <= n < X : tau(n) <= C ln n} / (X - 1)    [and the n>=2, '<' variants]
  block k     = [2^(k-1), 2^k)
  z(n)        = (tau(n) - C* ln n) / (sigma sqrt(ln n)),  C* = 3/ln(4/3),
                sigma = sqrt(2 (1 - C* ln 2)^2 / (2 ln 2 - ln 3))
  G_k(C)      = block density of tau(n) <= C ln n
  floor count = #{n in block k : tau(n) < 3.307 ln n}
"""
import math
import sys
import time

import numpy as np
from numba import njit

K_MAX = 24
N = 1 << K_MAX

C_STAR = 3.0 / math.log(4.0 / 3.0)
SIGMA = math.sqrt(2.0 * (1.0 - C_STAR * math.log(2.0)) ** 2 / (2.0 * math.log(2.0) - math.log(3.0)))
C_MAZUR = 523.0 / 50.0


@njit(cache=True)
def census(N):
    tau = np.zeros(N, dtype=np.int32)
    peak = np.int64(0)
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
        tau[n] = steps + tau[x]
    return tau, peak


@njit(cache=True)
def block_stats(tau, k, cstar, sigma, c_list, floor_c):
    lo = 1 << (k - 1)
    hi = 1 << k
    cnt = hi - lo
    s = 0.0
    s2 = 0.0
    dens_le = np.zeros(len(c_list), dtype=np.int64)
    dens_lt = np.zeros(len(c_list), dtype=np.int64)
    floor_lt = 0
    floor_le = 0
    tmin = 10 ** 9
    tmax = 0
    for n in range(lo, hi):
        t = tau[n]
        L = math.log(n)
        z = (t - cstar * L) / (sigma * math.sqrt(L))
        s += z
        s2 += z * z
        for i in range(len(c_list)):
            if t <= c_list[i] * L:
                dens_le[i] += 1
            if t < c_list[i] * L:
                dens_lt[i] += 1
        if t < floor_c * L:
            floor_lt += 1
        if t <= floor_c * L:
            floor_le += 1
        if t < tmin:
            tmin = t
        if t > tmax:
            tmax = t
    mean = s / cnt
    var_pop = s2 / cnt - mean * mean
    var_samp = var_pop * cnt / (cnt - 1)
    return cnt, mean, math.sqrt(var_pop), math.sqrt(var_samp), dens_le, dens_lt, floor_lt, floor_le, tmin, tmax


@njit(cache=True)
def cumulative(tau, X, C, start):
    le = 0
    lt = 0
    for n in range(start, X):
        t = tau[n]
        L = math.log(n)
        if t <= C * L:
            le += 1
        if t < C * L:
            lt += 1
    return le, lt


def main():
    out = []

    def p(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        out.append(line)

    p("# audit_tau_census.py (fable audit of opus S64), N = 2^%d" % K_MAX)
    p("C* = 3/ln(4/3) = %.6f ; sigma = %.6f ; sigma^2 = %.4f ; 523/50 = %.4f ; 523/50 / C* - 1 = %.5f" % (C_STAR, SIGMA, SIGMA ** 2, C_MAZUR, C_MAZUR / C_STAR - 1))
    t0 = time.time()
    tau, peak = census(N)
    p("census done in %.1f s ; peak value on any trajectory = %d (< 2^63 = %d: no overflow)" % (time.time() - t0, peak, 2 ** 63))
    # sanity checks against known values
    known = {1: 0, 2: 1, 3: 7, 6: 8, 7: 16, 9: 19, 27: 111, 97: 118, 871: 178, 6171: 261, 77031: 350, 837799: 524}
    for n, v in known.items():
        p("  tau(%d) = %d (expected %d) %s" % (n, tau[n], v, "OK" if tau[n] == v else "MISMATCH"))
    # argmax below 10^6 should be 837799 with 524 steps
    am = int(np.argmax(tau[:10 ** 6]))
    p("  argmax tau on n < 10^6: n = %d, tau = %d (expected 837799, 524)" % (am, tau[am]))
    am = int(np.argmax(tau[:N]))
    p("  argmax tau on n < 2^24: n = %d, tau = %d" % (am, tau[am]))

    p("")
    p("## Cumulative density F_X(523/50) = #{n < X : tau(n) <= (523/50) ln n} / count")
    p("X      | n>=1,<= | n>=1,<  | n>=2,<= | n>=2,<  | (denominators X-1 resp. X-2)")
    for k in (16, 20, 24):
        X = 1 << k
        le1, lt1 = cumulative(tau, X, C_MAZUR, 1)
        le2, lt2 = cumulative(tau, X, C_MAZUR, 2)
        p("2^%2d   | %.6f | %.6f | %.6f | %.6f | counts %d %d %d %d" % (k, le1 / (X - 1), lt1 / (X - 1), le2 / (X - 2), lt2 / (X - 2), le1, lt1, le2, lt2))
    # digest figure: 54.22% of n < 10^6 satisfy tau(n) <= 10.46 ln n
    X = 10 ** 6
    le1, lt1 = cumulative(tau, X, 10.46, 1)
    le2, lt2 = cumulative(tau, X, 10.46, 2)
    p("X = 10^6, C = 10.46: n>=1 '<=' %d/%d = %.4f%% ; '<' %d -> %.4f%% ; n>=2 '<=' %d/%d = %.4f%%" % (le1, X - 1, 100 * le1 / (X - 1), lt1, 100 * lt1 / (X - 1), le2, X - 2, 100 * le2 / (X - 2)))
    le1b, _ = cumulative(tau, X + 1, 10.46, 1)
    p("X = 10^6 inclusive (n <= 10^6): '<=' %d/%d = %.4f%%" % (le1b, X, 100 * le1b / X))

    p("")
    p("## Block statistics on [2^(k-1), 2^k), z(n) = (tau - C* ln n)/(sigma sqrt(ln n))")
    p("k  | count    | mean z   | sd(pop)  | sd(samp) | G_k(8) <= | G_k(13) <= | G_k(523/50) <= | G_k(8) < | G_k(13) < | floor<3.307 | floor<=3.307 | min tau/ln n | max tau")
    c_list = np.array([8.0, 13.0, C_MAZUR])
    for k in range(16, K_MAX + 1):
        cnt, mean, sdp, sds, dle, dlt, fl_lt, fl_le, tmin, tmax = block_stats(tau, k, C_STAR, SIGMA, c_list, 3.307)
        p("%2d | %8d | %+.5f | %.5f | %.5f | %.6f | %.6f | %.6f | %.6f | %.6f | %d | %d | %d | %d" % (
            k, cnt, mean, sdp, sds, dle[0] / cnt, dle[1] / cnt, dle[2] / cnt, dlt[0] / cnt, dlt[1] / cnt, fl_lt, fl_le, tmin, tmax))

    p("")
    p("## Floor counts #{n in block k : tau(n) < 3.307 ln n} and consecutive ratios")
    prev = None
    for k in range(17, K_MAX + 1):
        cnt, mean, sdp, sds, dle, dlt, fl_lt, fl_le, tmin, tmax = block_stats(tau, k, C_STAR, SIGMA, c_list, 3.307)
        # also with C_F exact = (1 + (1+log2 3)/2)/ln 2
        cF = (1.0 + (1.0 + math.log2(3.0)) / 2.0) / math.log(2.0)
        _, _, _, _, _, _, fl_lt_F, fl_le_F, _, _ = block_stats(tau, k, C_STAR, SIGMA, c_list, cF)
        ratio = "" if prev is None else "%.4f" % (fl_lt / prev)
        p("k = %2d : floor(<3.307) = %6d  ratio = %s | (<=3.307) = %d | with C_F = %.6f: < %d, <= %d" % (k, fl_lt, ratio, fl_le, cF, fl_lt_F, fl_le_F))
        prev = fl_lt
    p("expected ratio if count ~ X^E with E = 0.738: 2^0.738 = %.4f ; with E = 1: 2" % (2 ** 0.738))

    # minimum of tau(n)/ln n over block for the floor question
    p("")
    p("## Smallest tau(n)/ln n per block (k = 17..24) and the n attaining it")
    for k in range(17, K_MAX + 1):
        lo, hi = 1 << (k - 1), 1 << k
        seg = tau[lo:hi].astype(np.float64) / np.log(np.arange(lo, hi, dtype=np.float64))
        i = int(np.argmin(seg))
        p("k = %2d : min tau/ln n = %.4f at n = %d (tau = %d)" % (k, seg[i], lo + i, tau[lo + i]))

    np.save(sys.argv[1] if len(sys.argv) > 1 else "tau_2p24.npy", tau)
    with open("audit_tau_census.out", "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
