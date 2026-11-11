# Auditor recompute-5128: README section 5 item 2, "the misses sit where the singular series is smallest":
# mean of r(N) (ordered F1 triples) over odd N of a window, restricted to N coprime to 255255 = 3*5*7*11*13*17 and to
# N with N(N+1) coprime to 255255, relative to the mean over all odd N; also with the word-pair weight (s + 1);
# and the miss rates on these sets.
import sys
import math
import numpy as np
from numba import njit
from a0_defs import psi_table32


@njit(cache=True)
def rep_counts_w(psi, B):
    n = np.int64(1) << (B - 1)
    r = np.zeros(n, dtype=np.uint16)       # triples
    rw = np.zeros(n, dtype=np.uint32)      # triples weighted by s + 1
    for s in range(0, B):
        corr = 1 if s == 0 else 0
        Xs = (np.int64(1) << B) - 1 + corr
        h1 = np.int64(0)
        while True:
            y1 = 2 * h1 + 1
            y2crude = (4 * Xs) // ((1 + (np.int64(1) << (s + 2))) * y1)
            if (y2crude - 1) // 2 < h1:
                break
            x1 = 2 * np.int64(psi[h1]) + 1
            b = y1 << s
            y2max = (2 * Xs) // (x1 + 2 * b)
            h2max = (y2max - 1) // 2
            for h2 in range(h1, h2max + 1):
                M = x1 * (2 * np.int64(psi[h2]) + 1) + b * (2 * h2 + 1)
                if M <= Xs:
                    i = (M - corr - 1) >> 1
                    w = 2 if h2 > h1 else 1
                    r[i] += w
                    rw[i] += w * (s + 1)
            h1 += 1
    return r, rw


def main():
    B = 27
    psi = psi_table32(((1 << B) + 1) // 3 + 4)
    r, rw = rep_counts_w(psi, B)
    Q = 255255
    for (a, b) in [(18, 21), (21, 23), (24, 27)]:
        i0, i1 = 1 << (a - 1), 1 << (b - 1)
        N = 2 * np.arange(i0, i1, dtype=np.int64) + 1
        seg = r[i0:i1].astype(np.float64)
        segw = rw[i0:i1].astype(np.float64)
        g1 = np.gcd(N, Q) == 1
        g2 = g1 & (np.gcd(N + 1, Q) == 1)
        g3 = np.gcd(N + 1, Q) == 1
        print(f"[2^{a}, 2^{b}): mean r = {seg.mean():.4f}; word-pair weighted mean = {segw.mean():.4f}; miss rate {100 * np.mean(seg == 0):.3f} %")
        for name, sel in (("N coprime to 255255", g1), ("N + 1 coprime to 255255", g3), ("N (N + 1) coprime to 255255", g2)):
            print(f"    {name:30s}: share of odd N {sel.mean():.4f}; mean r / overall = {seg[sel].mean() / seg.mean():.4f}; weighted (s+1) mean / overall = {segw[sel].mean() / segw.mean():.4f}; miss rate {100 * np.mean(seg[sel] == 0):.3f} % ({100 * np.sum(seg[sel] == 0) / max(1, np.sum(seg == 0)):.1f} % of all misses)")
    print(f"limit values: prod over p in 3..17 of (1 - 1/p^2) = {np.prod([1 - 1 / p ** 2 for p in (3, 5, 7, 11, 13, 17)]):.4f}; 8/pi^2 = {8 / math.pi ** 2:.4f}")


if __name__ == "__main__":
    main()
