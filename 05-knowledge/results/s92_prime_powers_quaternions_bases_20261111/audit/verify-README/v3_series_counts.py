#!/usr/bin/env python3
"""Verification of the corrected S92 note (README.md), section 1, last paragraph: the multiplicity r of the odd N of
[2^24, 2^25) on each sheet of F1 against the series S(N) = (8/pi^2) sigma_1(n)/n, n the odd part of N (sheet 1) or
of N + 1 (sheet 0).  Independent program (own enumeration, own sieve).  Exit status 1 if a quoted number is not
reproduced.
"""
import math
import sys
import time

import numpy as np
from numba import njit

OK = True


def check(cond, msg):
    global OK
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        OK = False


@njit(cache=True)
def build_q(size):
    q = np.zeros(size, np.int32)
    for h in range(1, size):
        k = h >> 1
        q[h] = q[k] + (((k + 1) >> 1) if (h & 1) else 0)
    return q


@njit(cache=True)
def multiplicities(q, A, B, r):
    """r[sheet, N - A] = number of ordered triples (h1, h2, s) with N = x1 x2 + 2^s y1 y2 - [s = 0], A <= N < B."""
    H = q.shape[0]
    short = 0
    for s in range(0, 40):
        d0 = 1 if s == 0 else 0
        sheet = 0 if s == 0 else 1
        if (1 << s) > B:
            break
        for h1 in range(H):
            y1 = 2 * h1 + 1
            lowx = 1 + 2 * (h1 // 2)
            if lowx * lowx + (y1 * y1 << s) - d0 >= B:
                break
            x1 = y1 - 2 * np.int64(q[h1])
            h2 = h1
            while True:
                if h2 >= H:
                    short += 1
                    break
                y2 = 2 * h2 + 1
                if x1 * (1 + 2 * (h2 // 2)) + (y1 * y2 << s) - d0 >= B:
                    break
                N = x1 * (y2 - 2 * np.int64(q[h2])) + (y1 * y2 << s) - d0
                if A <= N < B:
                    r[sheet, N - A] += 1 if h2 == h1 else 2
                h2 += 1
    return short


@njit(cache=True)
def sigma_sieve(n):
    sg = np.zeros(n + 1, np.int64)
    for d in range(1, n + 1):
        for k in range(d, n + 1, d):
            sg[k] += d
    return sg


def main():
    t0 = time.time()
    KLO, KHI = 24, 25
    A, B = 1 << KLO, 1 << KHI
    q = build_q((B >> 1) + 4)
    r = np.zeros((2, B - A), np.int32)
    short = multiplicities(q, A, B, r)
    check(short == 0, "the table of states is long enough")
    sg = sigma_sieve(B)
    N = np.arange(A + 1, B, 2, dtype=np.int64)                 # the odd N of the window
    print(f"odd N in [2^{KLO}, 2^{KHI}): {len(N)} values  ({time.time() - t0:.0f} s)")
    want = {0: dict(tot=49644073, mean=5.9180, miss=8.07, corr=0.40, slope=1.44, r2=16, q=(4.301, 8.550)),
            1: dict(tot=62761167, mean=7.4817, miss=4.86, corr=0.43, slope=1.43, r2=18, q=(5.438, 10.766))}
    for sheet in (0, 1):
        rr = r[sheet, N - A].astype(np.float64)
        check(int(r[sheet].sum()) == int(rr.sum()), f"sheet {sheet}: every value of the sheet in the window is odd")
        arg = N + 1 if sheet == 0 else N.copy()
        n = arg.copy()
        while True:
            ev = (n & 1) == 0
            if not ev.any():
                break
            n[ev] >>= 1
        ab = sg[n] / n
        S = (8 / math.pi ** 2) * ab
        tot = int(rr.sum())
        mean = rr.mean()
        miss = 100 * float((rr == 0).mean())
        corr = float(np.corrcoef(rr, ab)[0, 1])
        slope = float(np.polyfit(S, rr / mean, 1)[0])
        order = np.argsort(ab, kind="stable")
        qs = np.array_split(order, 5)
        qr = [float(rr[ix].mean()) for ix in qs]
        qa = [float(ab[ix].mean()) for ix in qs]
        print(f"sheet {sheet}: {tot} triples, mean r {mean:.4f}, r = 0 for {miss:.3f} %, correlation with sigma_1(n)/n {corr:.4f} "
              f"(r^2 = {100 * corr * corr:.1f} %), slope of r/mean(r) against S {slope:.4f}, mean of S {S.mean():.5f}")
        print(f"    quintiles of sigma_1(n)/n: abundancy {', '.join(f'{v:.4f}' for v in qa)} (ratio {qa[-1] / qa[0]:.3f}); "
              f"mean r {', '.join(f'{v:.3f}' for v in qr)} (ratio {qr[-1] / qr[0]:.3f}); min S {S.min():.4f}")
        w = want[sheet]
        check(tot == w["tot"] and round(mean, 4) == w["mean"], f"sheet {sheet}: {w['tot']} triples, mean multiplicity {w['mean']}")
        check(round(corr, 2) == w["corr"] and round(100 * corr * corr) == w["r2"], f"sheet {sheet}: correlation {w['corr']}, {w['r2']} % of the variance")
        check(round(slope, 2) == w["slope"], f"sheet {sheet}: slope {w['slope']}")
        check(round(miss, 2) == w["miss"], f"sheet {sheet}: {w['miss']} % of the odd N have r = 0")
        check((round(qr[0], 3), round(qr[-1], 3)) == w["q"] and round(qa[-1] / qa[0], 2) == 1.61 and 1.95 < qr[-1] / qr[0] < 2.0,
              f"sheet {sheet}: mean r doubles across the quintiles ({qr[-1] / qr[0]:.3f}) while the series varies by 1.61")
    print(f"({time.time() - t0:.0f} s)")
    print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
    sys.exit(0 if OK else 1)


if __name__ == "__main__":
    main()
