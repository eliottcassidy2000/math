# Auditor recompute-5128: F1 parts of README section 5 (and two small items of section 3):
# representation counts r_sigma(N) = #{ordered triples (h1, h2, s) on the sheet sigma with value N} for odd N < 2^B,
# miss rates on the stated windows, mean multiplicity, relative variance, misses by class modulo 3.
import sys
import math
import numpy as np
from numba import njit
from a0_defs import psi_table32


@njit(cache=True)
def rep_counts(psi, B):
    n = np.int64(1) << (B - 1)  # odd N < 2^B  <->  index (N - 1)/2
    r0 = np.zeros(n, dtype=np.uint16)
    r1 = np.zeros(n, dtype=np.uint16)
    sat = 0
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
                    N = M - corr
                    w = 2 if h2 > h1 else 1
                    i = (N - 1) >> 1
                    if s == 0:
                        if r0[i] > 65000:
                            sat += 1
                        else:
                            r0[i] += w
                    else:
                        if r1[i] > 65000:
                            sat += 1
                        else:
                            r1[i] += w
            h1 += 1
    return r0, r1, sat


def main():
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 27
    psi = psi_table32(((1 << B) + 1) // 3 + 4)
    r0, r1, sat = rep_counts(psi, B)
    r = r0.astype(np.int64) + r1.astype(np.int64)
    print(f"# representation counts for odd N < 2^{B}; saturated increments: {sat}; total triples {int(r.sum())} (sheet 0 {int(r0.astype(np.int64).sum())}, sheet 1 {int(r1.astype(np.int64).sum())})")
    print("== small values (README section 3 table: r_F1 = 1, 3, 0, 0, 4, 5, 0 for N = 3, 5, 7, 13, 19, 27, 1039) ==")
    for N in [1, 3, 5, 7, 9, 11, 13, 19, 21, 27, 1039]:
        i = (N - 1) // 2
        print(f"   N = {N:5d}: r = {int(r[i])} (sheet 0: {int(r0[i])}, sheet 1: {int(r1[i])})")
    lo, hi = (1025 - 1) // 2, (2049 - 1) // 2
    print(f"   odd N in [1025, 2049): {hi - lo} numbers, F1 misses: {int(np.sum(r[lo:hi] == 0))}   (README: 65 of 512)")
    print(f"   least odd N missed by F1: {[int(2 * i + 1) for i in np.nonzero(r[:200] == 0)[0][:12]]}")
    print()
    print("== miss rates and moments on windows [2^a, 2^b) ==")
    gbar = math.pi ** 2 / 6 + 2 * float(__import__('mpmath').polylog(2, -0.5))
    for (a, b) in [(18, 21), (21, 23), (23, 24), (24, 27), (24, 25), (25, 26), (26, 27)]:
        if b > B:
            continue
        i0, i1 = 1 << (a - 1), 1 << (b - 1)
        seg = r[i0:i1]
        s0 = r0[i0:i1].astype(np.float64)
        s1 = r1[i0:i1].astype(np.float64)
        miss = int(np.sum(seg == 0))
        n = i1 - i0
        Ns = 2.0 * np.arange(i0, i1) + 1
        meanlog = np.mean(np.log(Ns))
        print(f" [2^{a}, 2^{b}): odd N {n}, misses {miss} = {100 * miss / n:.3f} %; mean r = {seg.mean():.4f} (g_bar * mean ln N = {gbar * meanlog:.4f}, difference {seg.mean() - gbar * meanlog:+.4f}); "
              f"Var/mean^2: all {seg.var() / seg.mean() ** 2:.4f}, sheet 0 {s0.var() / s0.mean() ** 2:.4f}, sheet 1 {s1.var() / s1.mean() ** 2:.4f}; Var/mean: all {seg.var() / seg.mean():.3f}, sheet 0 {s0.var() / s0.mean():.3f}, sheet 1 {s1.var() / s1.mean():.3f}")
        # misses by class modulo 3
        cls = (2 * np.arange(i0, i1) + 1) % 3
        mm = [int(np.sum((seg == 0) & (cls == c))) for c in range(3)]
        nn = [int(np.sum(cls == c)) for c in range(3)]
        rates = [mm[c] / nn[c] for c in range(3)]
        print(f"      misses by N mod 3 (0, 1, 2): {mm}; rates {100 * rates[0]:.3f} %, {100 * rates[1]:.3f} %, {100 * rates[2]:.3f} %; rate(1) / mean rate(0, 2) = {rates[1] / ((mm[0] + mm[2]) / (nn[0] + nn[2])):.2f}")
        mr = [float(seg[cls == c].mean()) for c in range(3)]
        print(f"      mean r by N mod 3 (0, 1, 2): {mr[0]:.4f}, {mr[1]:.4f}, {mr[2]:.4f}")


if __name__ == "__main__":
    main()
