#!/usr/bin/env python3
# Auditor "hypotheses": consistency check of the base-b censuses quoted in THM-5134 (H) "Numerical support"
# (every_base.out (4): sheet s >= 1, N_b <= 2^24;  bridges_exact.out S10: sheet s = 0, N_b <= 2^22).
# Own COMPLETE census by value (ordered triples), and the same census with the words capped at h < Hcap,
# to see which set the session's programs actually enumerated
# (every_base.py uses q[: 2^21]; bridges_exact.py s10 uses size = 2^12).
# usage: python3 base_census_check.py
import sys, time, math
import numpy as np
from numba import njit, prange, get_num_threads, get_thread_id
from f1core import build_psi

PR = np.array([3, 5, 7, 11, 13], dtype=np.int64)


@njit(parallel=True, cache=True)
def census(X, b, psi, hcap, nthreads, chunk):
    # out[thread, sheet, prime index, t (p-divisible smaller state), class]
    nH = min(psi.shape[0], hcap)
    out = np.zeros((nthreads, 2, 5, 2, 16), dtype=np.int64)
    tot = np.zeros((nthreads, 2), dtype=np.int64)
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
                if s == 0:
                    sheet = 0
                    cy = 0
                    v0 = X1 + b * h1
                else:
                    sheet = 1
                    cy = (np.int64(1) << (s - 1)) * y1
                    v0 = X1 + b * cy
                if v0 > X:
                    break
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
                        tot[tid, sheet] += 1
                        hs = h1 if h1 <= h2 else h2
                        Xs = 1 + b * psi[hs]
                        ys = 2 * hs + 1
                        for j in range(5):
                            p = PR[j]
                            t = 1 if (h1 != h2 and Xs % p == 0 and ys % p == 0) else 0
                            out[tid, sheet, j, t, Nb % p] += 1
                s += 1
    return out, tot


def report(X, b, psi, hcap, label):
    out, tot = census(X, b, psi, hcap, get_num_threads(), 64)
    out = out.sum(axis=0)
    tot = tot.sum(axis=0)
    for sh in (0, 1):
        cells = []
        for j, p in enumerate(PR.tolist()):
            if b % p == 0:
                continue
            c0 = 0 if sh == 1 else (-b * pow(2, -1, p)) % p
            D1 = int(out[sh, j, 0].sum())
            Dp = int(out[sh, j, 1].sum())
            # exactness: divisible smaller state => class c0
            ex = int(out[sh, j, 1].sum() - out[sh, j, 1, c0])
            rem = max(abs(int(out[sh, j, 0, c] + out[sh, j, 1, c]) - (D1 / p + (Dp if c == c0 else 0))) for c in range(p)) * p / tot[sh]
            fac0 = (out[sh, j, 0, c0] + out[sh, j, 1, c0]) * p / tot[sh]
            cells.append("p=%d: class %d factor %.3f, F p^2 = %.2f, rem %.3f%s" % (p, c0, fac0, Dp / tot[sh] * p * p, rem, "" if ex == 0 else " EXACTNESS FAILS"))
        print("  %s  b = %2d  sheet %d: %10d triples; %s" % (label, b, sh, tot[sh], "; ".join(cells)))
    return tot


def main():
    t0 = time.time()
    print("# base_census_check.py (auditor hypotheses; own code)")
    print("# Session, bridges_exact.out S10 (sheet 0, 'N_beta <= 2^22'): beta=2: 2845759 triples, p=3 factor 1.107 F p^2 0.58 rem 0.036; p=7 rem 0.158;")
    print("#   beta=3: 1832488; beta=5: 985666; beta=10: 379472 triples.")
    print("# Session, every_base.out (4) (sheet 1, 'N_b <= 2^24'): b=2: 55892593; b=3: 35516685; b=4: 24945520; b=5: 18836356; b=6: 14916908; b=10: 7632149 triples.")
    print()
    X = 1 << 22
    print("X = 2^22")
    for b in (2, 3, 5, 10):
        nH = int(2 * X // (3 * b)) + 4
        psi = build_psi(nH)
        report(X, b, psi, 1 << 62, "COMPLETE      ")
        report(X, b, psi, 1 << 12, "words < 2^12  ")
    X = 1 << 24
    print()
    print("X = 2^24")
    for b in (2, 3, 4, 5, 6, 10):
        nH = int(2 * X // (3 * b)) + 4
        psi = build_psi(nH)
        report(X, b, psi, 1 << 62, "COMPLETE      ")
        report(X, b, psi, 1 << 21, "words < 2^21  ")
    print("# %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
