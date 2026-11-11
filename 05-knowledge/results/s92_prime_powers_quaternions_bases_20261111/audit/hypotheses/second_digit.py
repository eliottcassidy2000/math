#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9665 (C). The second base-b digit of N_b, by value, on each sheet.
#   N_b = (1 + b psi1)(1 + b psi2) + b (2^s y1 y2 - [s=0])/2,   alpha_1 = psi1 + psi2 + (2^s y1 y2 - [s=0])/2,
#   second digit = alpha_1 mod b  (N_b = 1 + alpha_1 b + alpha_2 b^2).
# Census over ORDERED triples (h1,h2,s) with N_b <= X = 2^E.
# Also: histogram split by whether the smaller state h_min is 0, and the joint histogram of
#   (q2 mod b, h2 mod b)-free test "fixed smaller state": for h_min in {0,1,2,3} the digit law of the triples
#   with that smaller state (the first-order mechanism: a fixed smaller state should already give a flat digit).
# usage: python3 second_digit.py b E
import sys, time, math
import numpy as np
from numba import njit, prange, get_num_threads, get_thread_id
from f1core import build_psi


@njit(parallel=True, cache=True)
def census(X, b, psi, nthreads, chunk, nsmall):
    nH = psi.shape[0]
    # hist[thread, sheet, digit]; hsmall[thread, sheet, smaller state (< nsmall), digit]
    hist = np.zeros((nthreads, 2, b), dtype=np.int64)
    hsmall = np.zeros((nthreads, 2, nsmall, b), dtype=np.int64)
    nchunks = (nH + chunk - 1) // chunk
    for ci in prange(nchunks):
        tid = get_thread_id()
        lo = ci * chunk
        hi = min(nH, lo + chunk)
        for h1 in range(lo, hi):
            X1 = 1 + b * psi[h1]
            y1 = 2 * h1 + 1
            if X1 > X:
                continue
            s = 0
            while True:
                # value = X1 * X2 + b * (2^s y1 y2 - [s=0]) / 2
                # sheet 0: (y1 y2 - 1)/2 = y1 h2 + h1 ;  sheet s>=1: 2^(s-1) y1 y2
                if s == 0:
                    sheet = 0
                    cy = 0
                else:
                    sheet = 1
                    cy = (np.int64(1) << (s - 1)) * y1
                # lower bound of the value at h2 = 0 : X1 * 1 + b * (s==0 ? h1 : cy)
                if s == 0:
                    v0 = X1 + b * h1
                else:
                    v0 = X1 + b * cy
                if v0 > X:
                    break
                # psi(h2) >= ceil(h2/2) ... use psi2 >= h2/2: X1 (1 + b h2/2) + b * (...)(2 h2 + 1) <= X  -> crude bound
                # sheet 0: X1 (1 + b h2 / 2) + b (y1 h2 + h1) <= X ; sheet 1: X1 (1 + b h2/2) + b cy (2 h2 + 1) <= X
                if s == 0:
                    num = 2 * (X - X1 - b * h1)
                    den = X1 * b + 2 * b * y1
                else:
                    num = 2 * (X - X1 - b * cy)
                    den = X1 * b + 4 * b * cy
                h2max = num // den + 1
                if h2max >= nH:
                    h2max = nH - 1
                for h2 in range(h2max + 1):
                    p2 = psi[h2]
                    if s == 0:
                        half = y1 * h2 + h1
                    else:
                        half = cy * (2 * h2 + 1)
                    Nb = X1 * (1 + b * p2) + b * half
                    if Nb <= X:
                        a1 = psi[h1] + p2 + half
                        d = a1 % b
                        hist[tid, sheet, d] += 1
                        hs = h1 if h1 <= h2 else h2
                        if hs < nsmall:
                            hsmall[tid, sheet, hs, d] += 1
                s += 1
    return hist, hsmall


def main():
    b = int(sys.argv[1])
    E = int(sys.argv[2])
    X = 1 << E
    t0 = time.time()
    # h bound: smallest value for a given h2 is with h1 = 0, s = 0: 1 + b psi2 + b h2 >= 1 + 1.5 b h2
    nH = int(2 * X // (3 * b)) + 4
    psi = build_psi(nH)
    # sanity: psi >= h/2
    hh = np.arange(nH, dtype=np.int64)
    assert np.all(2 * psi >= hh)
    hist, hsmall = census(X, b, psi, get_num_threads(), 64, 6)
    hist = hist.sum(axis=0)
    hsmall = hsmall.sum(axis=0)
    print("# second_digit.py  b = %d  X = 2^%d  (auditor hypotheses; own code)   %.1f s" % (b, E, time.time() - t0))
    for sh in (0, 1):
        T = hist[sh].sum()
        shares = hist[sh] / T
        dev = np.max(np.abs(b * shares - 1))
        noise = math.sqrt(b / T)
        print("sheet %d  #T = %d  b*share by digit: %s   max|b*share - 1| = %.5f   (Poisson scale sqrt(b/#T) = %.5f)" % (
            sh, T, " ".join("%.4f" % (b * v) for v in shares), dev, noise))
        for hs in range(6):
            Ts = hsmall[sh, hs].sum()
            if Ts == 0:
                continue
            ss = hsmall[sh, hs] / Ts
            print("     smaller state h = %d: #T = %9d  max|b*share - 1| = %.5f   b*share: %s" % (
                hs, Ts, np.max(np.abs(b * ss - 1)), " ".join("%.4f" % (b * v) for v in ss)))
    # prime-power components of b
    if b in (6, 10, 12):
        for sh in (0, 1):
            T = hist[sh].sum()
            for pp in [d for d in range(2, b) if b % d == 0]:
                c = np.zeros(pp)
                for d in range(b):
                    c[d % pp] += hist[sh, d]
                print("sheet %d  digit mod %d: max|%d*share - 1| = %.5f" % (sh, pp, pp, np.max(np.abs(pp * c / T - 1))))


if __name__ == "__main__":
    main()
