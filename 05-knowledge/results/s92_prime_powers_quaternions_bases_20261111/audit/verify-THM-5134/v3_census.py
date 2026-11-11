#!/usr/bin/env python3
"""Verifier of THM-5134 (after the audit corrections): independent complete censuses by value, own code.

  C1  F1 (base 2), X = 2^22 and 2^24, both sheets: triple counts, the exact statement of (B) (a triple of type d has
      d | N - c_sigma), the remainder of (B) at 9, 27, 81, 25, 125, 49, 45, 75 (table of "Numerical support"),
      the class factors at 9 and 25 by valuation.
  C2  base b, sheet s >= 1, N_b <= 2^24, b = 2, 3, 4, 5, 6, 10: triple counts, remainder of the first-order law at
      the primes 3..13 not dividing b, F p^2, the second digit.
  C3  base b, sheet s = 0, N_b <= 2^22, b = 2, 3, 5, 10: triple counts, violations of 2 N_b + b = 0 (mod p) for a
      p-divisible smaller state, remainder at each prime.
  C4  second digit for b = 3, both sheets, X = 2^16, 2^20, 2^24: max_c |b share - 1|.
Every enumeration uses only sound bounds: psi(h) >= h // 2 gives a lower bound of N_b that is monotone in h2.
Exit status 1 on failure.
"""
import math
import sys
import time

import numpy as np
from numba import njit

FAIL = []


def check(name, ok, extra=""):
    print(("  ok   " if ok else "  FAIL ") + name + ((" : " + extra) if extra else ""), flush=True)
    if not ok:
        FAIL.append(name)


def psi_table(n):
    ps = np.zeros(n, dtype=np.int64)
    for h in range(1, n):
        k = h >> 1
        ps[h] = ps[k] + (k if h % 2 == 0 else (k >> 1) + 1)
    return ps


@njit(cache=False)
def psi_table_fast(n):
    ps = np.zeros(n, dtype=np.int64)
    for h in range(1, n):
        k = h >> 1
        if h % 2 == 0:
            ps[h] = ps[k] + k
        else:
            ps[h] = ps[k] + (k >> 1) + 1
    return ps


@njit(cache=False)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(cache=False)
def census(ps, b, X, sheet, LM, hist, typecount, digit2):
    """All ordered triples (h1, h2, s) of the sheet with N_b <= X (pairs h1 <= h2, weight 2 off the diagonal).
    Value V = N_b on sheet 1, 2 N_b + b on sheet 0 (for b = 2: V = 2 M on sheet 0, M = N + 1).
    hist[V mod LM]; typecount[gcd(g_b(smaller word), LM)] for h1 != h2; digit2[(N_b - 1)/b mod b].
    Returns (number of triples, number of diagonal triples, violations of type | V, 1 if the table was too short)."""
    H = ps.shape[0]
    total = 0
    diag = 0
    viol = 0
    short = 0
    smax = 0 if sheet == 0 else 62
    smin = 0 if sheet == 0 else 1
    for s in range(smin, smax + 1):
        # N_b = x1 x2 + b (2^s y1 y2 - [s = 0]) / 2 ; twoc = b 2^s (so that 2 N_b = 2 x1 x2 + twoc y1 y2 - b [s = 0])
        twoc = b << s
        e0 = b if s == 0 else 0
        if 2 + twoc - e0 > 2 * X:
            break
        for h1 in range(H):
            y1 = 2 * h1 + 1
            lb1 = 1 + b * (h1 >> 1)
            if 2 * lb1 * lb1 + twoc * y1 * y1 - e0 > 2 * X:
                break
            x1 = 1 + b * ps[h1]
            t1 = gcd(gcd(x1, y1), LM)
            closed = False
            for h2 in range(h1, H):
                y2 = 2 * h2 + 1
                if 2 * x1 * (1 + b * (h2 >> 1)) + twoc * y1 * y2 - e0 > 2 * X:
                    closed = True
                    break
                x2 = 1 + b * ps[h2]
                twoN = 2 * x1 * x2 + twoc * y1 * y2 - e0
                if twoN > 2 * X:
                    continue
                N = twoN // 2
                w = 2
                if h2 == h1:
                    w = 1
                    diag += 1
                total += w
                if sheet == 0:
                    V = 2 * N + b
                else:
                    V = N
                hist[V % LM] += w
                digit2[((N - 1) // b) % b] += w
                if h2 != h1:
                    typecount[t1] += w
                    if V % t1 != 0:
                        viol += w
            if not closed:
                short = 1
    return total, diag, viol, short


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def fold(hist, LM, m):
    return hist.reshape(LM // m, m).sum(axis=0)


def remainder_B(hist, typecount, LM, m, total):
    hm = fold(hist, LM, m)
    D = {d: 0 for d in divisors(m)}
    for t in range(1, LM + 1):
        if typecount[t]:
            D[math.gcd(t, m)] += int(typecount[t])
    worst = 0.0
    pred = np.zeros(m)
    for c in range(m):
        g0 = math.gcd(m, c) if c else m
        pred[c] = sum((d / m) * D[d] for d in divisors(g0))
        worst = max(worst, abs(hm[c] - pred[c]))
    return m * worst / total, hm, D


def main():
    t0 = time.time()
    ps = psi_table_fast((1 << 23) + 8)
    print("C1 F1 (base 2): remainder of (B)")
    LM = 81 * 125 * 49
    table = {  # THM-5134, table of Numerical support: (2^22, 2^24), sheet 0 / sheet 1
        9: (("0.035", "0.028"), ("0.021", "0.024")), 27: (("0.046", "0.044"), ("0.022", "0.029")),
        81: (("0.062", "0.055"), ("0.024", "0.042")), 25: (("0.004", "0.008"), ("0.004", "0.005")),
        125: (("0.011", "0.015"), ("0.008", "0.008")), 49: (("0.182", "0.140"), ("0.163", "0.141")),
        45: (("0.083", "0.057"), ("0.040", "0.056")), 75: (("0.097", "0.065"), ("0.050", "0.056")),
    }
    want_T = {22: (10532579, 13218857), 24: (45886662, 57818155)}
    for col, K in enumerate((22, 24)):
        X = 1 << K
        res = []
        for sheet in (0, 1):
            hist = np.zeros(LM, np.int64)
            tc = np.zeros(LM + 1, np.int64)
            d2 = np.zeros(2, np.int64)
            total, diag, viol, short = census(ps, 2, X, sheet, LM, hist, tc, d2)
            res.append((total, diag, viol, short, hist, tc))
        check("X = 2^%d: triples %d / %d (file and cum_%d.out: %d / %d), tables long enough"
              % (K, res[0][0], res[1][0], K, want_T[K][0], want_T[K][1]),
              (res[0][0], res[1][0]) == want_T[K] and res[0][3] == 0 and res[1][3] == 0)
        check("X = 2^%d: every triple of type d has d | N - c_sigma (0 violations on both sheets, modulus %d)" % (K, LM),
              res[0][2] == 0 and res[1][2] == 0)
        okt = True
        cells = []
        for m in (9, 27, 81, 25, 125, 49, 45, 75):
            r = []
            for sheet in (0, 1):
                total, diag, viol, short, hist, tc = res[sheet]
                # sheet 0: V = 2 M, a unit multiple of M: same types, classes permuted with the same gcd
                rem, hm, D = remainder_B(hist, tc, LM, m, total)
                r.append(rem)
            got = ("%.3f" % r[0], "%.3f" % r[1])
            cells.append("%d: %.5f / %.5f" % (m, r[0], r[1]))
            if got != table[m][col]:
                okt = False
        check("X = 2^%d: remainder of (B) equals the table to 3 decimals at the eight moduli" % K, okt, "; ".join(cells))
        if K == 24:
            for m, vals in ((9, None), (25, None)):
                out = []
                for sheet in (0, 1):
                    total, diag, viol, short, hist, tc = res[sheet]
                    hm = fold(hist, LM, m)
                    groups = {}
                    for c in range(m):
                        g0 = math.gcd(m, c) if c else m
                        groups.setdefault(g0, []).append(m * hm[c] / total)
                    out.append(" ".join("g0=%d: %.4f" % (g, sum(v) / len(v)) for g, v in sorted(groups.items())))
                print("       class factors at %d, X = 2^24: sheet 0: %s | sheet 1: %s" % (m, out[0], out[1]))

    print("C2 base b, sheet s >= 1, N_b <= 2^24")
    X = 1 << 24
    primes = (3, 5, 7, 11, 13)
    LMp = 3 * 5 * 7 * 11 * 13
    want = {2: 57818155, 3: 35564915, 4: 24947862, 5: 18838778, 6: 14919260, 10: 7634213}
    want_rem7 = {2: "0.133", 3: "0.067", 4: "0.125", 5: "0.042", 6: "0.105", 10: "0.060"}
    want_d2 = {2: "0.000", 3: "0.040", 4: "0.001", 5: "0.002", 6: "0.052", 10: "0.003"}
    want_F = {(3, 7): "3.37", (10, 7): "3.28", (5, 13): "3.38"}
    okc = okr = okd = okF = okx = True
    others = 0.0
    for b in (2, 3, 4, 5, 6, 10):
        hist = np.zeros(LMp, np.int64)
        tc = np.zeros(LMp + 1, np.int64)
        d2 = np.zeros(b, np.int64)
        total, diag, viol, short = census(ps, b, X, 1, LMp, hist, tc, d2)
        if total != want[b] or short:
            okc = False
        dev2 = max(abs(b * d2[c] / total - 1) for c in range(b))
        if "%.3f" % dev2 != want_d2[b]:
            okd = False
        cells = []
        for p in primes:
            if b % p == 0:
                continue
            hp = fold(hist, LMp, p)
            Dp = sum(int(tc[t]) for t in range(1, LMp + 1) if tc[t] and t % p == 0)
            # exactness at p
            fs = Dp / total
            f = p * hp / total
            rem = max(abs(f[0] - (1 + (p - 1) * fs)), max(abs(f[c] - (1 - fs)) for c in range(1, p)))
            cells.append("p=%d: rem %.3f, F p^2 %.2f" % (p, rem, fs * p * p))
            if p == 7:
                if "%.3f" % rem != want_rem7[b]:
                    okr = False
            else:
                others = max(others, rem)
            if (b, p) in want_F and "%.2f" % (fs * p * p) != want_F[(b, p)]:
                okF = False
        # type | V for every p not dividing b: violations counted with gcd(g_b, LMp), restricted to primes prime to b
        if viol and all(b % p for p in primes):
            okx = False
        print("       b = %2d: %9d triples; second digit flat to %.4f; %s" % (b, total, dev2, "; ".join(cells)), flush=True)
    check("triple counts 57,818,155 / 35,564,915 / 24,947,862 / 18,838,778 / 14,919,260 / 7,634,213", okc)
    check("remainder at 7: 0.133, 0.067, 0.125, 0.042, 0.105, 0.060", okr)
    check("remainder at the other primes at most 0.051", "%.3f" % others == "0.051", "%.4f" % others)
    check("second digit flat to 0.000, 0.040, 0.001, 0.002, 0.052, 0.003", okd)
    check("F p^2 = 3.37 (b = 3, p = 7), 3.28 (b = 10, p = 7), 3.38 (b = 5, p = 13)", okF)

    print("C3 base b, sheet s = 0, N_b <= 2^22")
    X = 1 << 22
    want0 = {2: 10532579, 3: 6030660, 5: 2896626, 10: 1042874}
    want_rem7 = {2: "0.167", 3: "0.101", 5: "0.078", 10: "0.118"}
    okc = okr = okv = True
    others = 0.0
    for b in (2, 3, 5, 10):
        hist = np.zeros(LMp, np.int64)
        tc = np.zeros(LMp + 1, np.int64)
        d2 = np.zeros(b, np.int64)
        total, diag, viol, short = census(ps, b, X, 0, LMp, hist, tc, d2)
        if total != want0[b] or short:
            okc = False
        cells = []
        for p in primes:
            if b % p == 0:
                continue
            hp = fold(hist, LMp, p)
            # violations at p: recount by a direct pass is not needed, V % gcd(g_b, LMp) == 0 was tested for the
            # full gcd; for b with a prime of the list the gcd can contain that prime (g_b is prime to b, see below)
            Dp = sum(int(tc[t]) for t in range(1, LMp + 1) if tc[t] and t % p == 0)
            fs = Dp / total
            f = p * hp / total
            rem = max(abs(f[0] - (1 + (p - 1) * fs)), max(abs(f[c] - (1 - fs)) for c in range(1, p)))
            cells.append("p=%d: rem %.3f" % (p, rem))
            if p == 7:
                if "%.3f" % rem != want_rem7[b]:
                    okr = False
            else:
                others = max(others, rem)
        if viol:
            okv = False
        print("       b = %2d: %9d triples; violations %d; %s" % (b, total, viol, "; ".join(cells)), flush=True)
    check("triple counts 10,532,579 / 6,030,660 / 2,896,626 / 1,042,874", okc)
    check("0 violations: a p-divisible smaller state gives 2 N_b + b = 0 (mod p)  (x_b = 1 mod b, so g_b is prime "
          "to b and the test with gcd(g_b, 15015) covers exactly the primes not dividing b)", okv)
    check("remainder at 7: 0.167, 0.101, 0.078, 0.118", okr)
    check("remainder at the other primes at most 0.057", "%.3f" % others == "0.057", "%.4f" % others)

    print("C4 second digit, b = 3")
    want = {16: ("0.131", "0.127"), 20: ("0.071", "0.073"), 24: ("0.041", "0.040")}
    wantT = {16: (69184, 90768), 20: (1374205, 1837677), 24: (26257636, 35564915)}
    okk = True
    for K in (16, 20, 24):
        out = []
        tots = []
        for sheet in (0, 1):
            hist = np.zeros(LMp, np.int64)
            tc = np.zeros(LMp + 1, np.int64)
            d2 = np.zeros(3, np.int64)
            total, diag, viol, short = census(ps, 3, 1 << K, sheet, LMp, hist, tc, d2)
            out.append(max(abs(3 * d2[c] / total - 1) for c in range(3)))
            tots.append(total)
        print("       X = 2^%d: #T = %d / %d; max_c |b share - 1| = %.5f / %.5f" % (K, tots[0], tots[1], out[0], out[1]))
        if ("%.3f" % out[0], "%.3f" % out[1]) != want[K] or tuple(tots) != wantT[K]:
            okk = False
    check("b = 3: 0.131 / 0.127, 0.071 / 0.073, 0.041 / 0.040 at 2^16, 2^20, 2^24, with the triple counts of the "
          "audit (69184 / 90768, 1374205 / 1837677, 26257636 / 35564915)", okk)
    print("(%.0f s)" % (time.time() - t0))
    print()
    if FAIL:
        print("FAILURES: %d" % len(FAIL))
        for f in FAIL:
            print("   " + f)
        sys.exit(1)
    print("ALL CENSUS CHECKS PASSED")


main()
