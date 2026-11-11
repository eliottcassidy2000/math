#!/usr/bin/env python3
"""Auditor proof-5128-CF, check A10: the series A_p^sigma = sum_h chi(h) w^sigma(h) of THM-5128 (D2), (D4).

usage: python3 a10_ap.py NB      (states below 2^NB; default 26)

Streaming numba implementation (int32 table of psi; own gcd), independent of the session's numpy code.
Printed for each prime:
  - the partial sums over h < 2^n (block summation), n = 16..NB;
  - the block sums b_n = sum over B_n, and b_n 2^(theta_p n) (the proof claims b_n = O(2^(-theta_p n)));
  - the largest excursion of the natural-order partial sums inside each of the last blocks (does the order of
    summation inside a block matter?);
  - sum over h < 2^NB of |chi(h)| w(h) (to show that the convergence is conditional: it grows like ln).
"""
import sys
import time
import numpy as np
from numba import njit

NB = int(sys.argv[1]) if len(sys.argv) > 1 else 26
PR = np.array([3, 5, 7, 11, 13, 17, 31, 43], dtype=np.int64)
RATE = {3: 0.8660, 5: 0.7418, 7: 0.9320, 11: 0.8587, 13: 0.8272, 17: 0.8289, 31: 0.9002, 43: 0.8342}


@njit(cache=True)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(cache=True)
def run(NB, PR):
    H = 1 << NB
    psi = np.zeros(H, dtype=np.int32)
    npr = PR.shape[0]
    part = np.zeros((npr, 2, NB + 1))        # partial sums at h < 2^n
    exc = np.zeros((npr, 2, NB + 1))         # max |partial(h) - partial(start of block)| inside block n
    absum = np.zeros((npr, 2, NB + 1))       # partial sums of |chi| w
    acc = np.zeros((npr, 2))
    comp = np.zeros((npr, 2))                # Kahan compensation
    aacc = np.zeros((npr, 2))
    start = np.zeros((npr, 2))
    n = 0                                    # current block index: h in [2^(n-1), 2^n)
    for h in range(H):
        if h > 0:
            k = h >> 1
            if h & 1:
                psi[h] = psi[k] + (k >> 1) + 1
            else:
                psi[h] = psi[k] + k
        if h == (1 << n):                    # h starts block n + 1; record sums over h < 2^n
            for i in range(npr):
                for s in range(2):
                    part[i, s, n] = acc[i, s]
                    absum[i, s, n] = aacc[i, s]
                    start[i, s] = acc[i, s]
            n += 1
        x = 2 * np.int64(psi[h]) + 1
        y = 2 * np.int64(h) + 1
        g = gcd(y, x)
        xf = float(x)
        yf = float(y)
        w0 = np.log(2 * (xf + yf) / (xf + 2 * yf)) / xf
        w1 = np.log1p(xf / (2 * yf)) / xf
        for i in range(npr):
            p = PR[i]
            chi = (1.0 if g % p == 0 else 0.0) - 1.0 / (p * p)
            for s in range(2):
                term = chi * (w0 if s == 0 else w1)
                yk = term - comp[i, s]
                t = acc[i, s] + yk
                comp[i, s] = (t - acc[i, s]) - yk
                acc[i, s] = t
                aacc[i, s] += abs(term)
                d = abs(acc[i, s] - start[i, s])
                if d > exc[i, s, n]:
                    exc[i, s, n] = d
    for i in range(npr):
        for s in range(2):
            part[i, s, NB] = acc[i, s]
            absum[i, s, NB] = aacc[i, s]
    return part, exc, absum


t0 = time.time()
part, exc, absum = run(NB, PR)
print("states h < 2^%d; time %.1f s" % (NB, time.time() - t0))
for i, p in enumerate(PR.tolist()):
    r = RATE[p]
    th = np.log(1 / r) / np.log(2 / r ** 2)
    print("\np = %d   (theta_p <= %.4f)" % (p, th))
    for s in (0, 1):
        ns = list(range(16, NB + 1, 2))
        print("  A_p^%d partial sums at n = %s: " % (s, ns) + " ".join("%+.5f" % part[i, s, n] for n in ns))
        b = [part[i, s, n] - part[i, s, n - 1] for n in range(1, NB + 1)]
        print("     block sums b_n, n = %d..%d: " % (NB - 9, NB) + " ".join("%+.2e" % b[n - 1] for n in range(NB - 9, NB + 1)))
        print("     b_n 2^(theta n), n = %d..%d: " % (NB - 9, NB)
              + " ".join("%+.4f" % (b[n - 1] * 2 ** (th * n)) for n in range(NB - 9, NB + 1)))
        print("     max natural-order excursion inside block n = %d..%d: " % (NB - 5, NB)
              + " ".join("%.1e" % exc[i, s, n] for n in range(NB - 5, NB + 1)))
        print("     sum |chi| w over h < 2^n, n = 16, 20, 24, %d: " % NB
              + " ".join("%.4f" % absum[i, s, n] for n in (16, 20, 24, NB)) + "   (grows linearly in n: not absolutely convergent)")
