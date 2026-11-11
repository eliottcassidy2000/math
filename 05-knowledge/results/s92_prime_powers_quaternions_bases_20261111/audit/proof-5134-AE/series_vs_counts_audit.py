# series_vs_counts_audit.py -- auditor proof-5134-AE.  Written before reading series_vs_counts.py / its output.
# Multiplicity r_sigma(N) = number of ordered F1 triples of the sheet sigma with value N, for the odd N of a window
# [2^KLO, 2^KHI), against the series S_sigma(N) = (8/pi^2) sigma(n)/n, n the odd part of N + [sigma = 0].
#   - Pearson correlation of r_sigma with sigma(n)/n
#   - mean multiplicity by quintile of sigma(n)/n (lowest .. highest), and the mean of sigma(n)/n in each quintile
#   - the same with the series normalised: mean of r / (mean r * S) by quintile (1 would be "the series is the density")
# usage: python3 series_vs_counts_audit.py KHI KLO
import sys, math, time
import numpy as np
import numba as nb


@nb.njit(cache=True)
def build_psi(H):
    psi = np.zeros(H + 1, np.int32)
    for h in range(1, H + 1):
        k = h >> 1
        if (h & 1) == 0:
            psi[h] = psi[k] + k
        else:
            psi[h] = psi[k] + (k >> 1) + 1
    return psi


@nb.njit(cache=True)
def multiplicities(A, B, psi):
    # r[sheet, (N - A) // 2] for odd N in [A, B)   (A even)
    half = (B - A) // 2
    r = np.zeros((2, half), np.int32)
    smax = 0
    while (np.int64(1) << (smax + 1)) < B:
        smax += 1
    for s in range(0, smax + 1):
        sh = 1
        d0 = 0
        if s == 0:
            sh = 0
            d0 = 1
        p2 = np.int64(1) << np.int64(s)
        h1 = np.int64(0)
        while True:
            y1 = np.int64(2 * h1 + 1)
            x1 = np.int64(2 * psi[h1] + 1)
            if p2 * y1 * y1 - d0 >= B:
                break
            Nd = x1 * x1 + p2 * y1 * y1 - d0
            if Nd >= A and Nd < B:
                r[sh, (Nd - A) // 2] += 1
            den = x1 + 2 * p2 * y1
            y2max = (2 * (B + d0)) // den + 1
            h2max = y2max // 2 + 1
            Ay = p2 * y1
            for h2 in range(h1 + 1, h2max + 1):
                N = x1 * np.int64(2 * psi[h2] + 1) + Ay * np.int64(2 * h2 + 1) - d0
                if N >= A and N < B:
                    r[sh, (N - A) // 2] += 2
            h1 += 1
    return r


@nb.njit(cache=True)
def abundancy_odd(K):
    half = (K + 1) // 2
    ab = np.zeros(half, np.float64)
    for j in range(half):
        d = 2 * j + 1
        inv = 1.0 / d
        n = d
        while n <= K:
            ab[(n - 1) // 2] += inv
            n += 2 * d
    return ab


def main():
    KHI, KLO = int(sys.argv[1]), int(sys.argv[2])
    A, B = 1 << KLO, 1 << KHI
    t0 = time.time()
    H = (2 * (B + 1)) // 3 // 2 + 8
    psi = build_psi(H)
    r = multiplicities(np.int64(A), np.int64(B), psi)
    ab = abundancy_odd(B + 2)
    Ns = np.arange(A + 1, B, 2, dtype=np.int64)
    assert len(Ns) == r.shape[1]
    # sheet 1: n = N ; sheet 0: n = odd part of N + 1
    a1 = ab[(Ns - 1) // 2]
    t = Ns + 1
    while True:
        ev = (t % 2 == 0)
        if not ev.any():
            break
        t = np.where(ev, t // 2, t)
    a0 = ab[(t - 1) // 2]
    c = 8 / math.pi ** 2
    print("window odd N in [2^%d, 2^%d): %d odd N; triples sheet 0: %d, sheet 1: %d  (%.0f s)" % (KLO, KHI, len(Ns), int(r[0].sum()), int(r[1].sum()), time.time() - t0))
    for sh, a in ((0, a0), (1, a1)):
        rr = r[sh].astype(np.float64)
        corr = np.corrcoef(rr, a)[0, 1]
        order = np.argsort(a, kind="stable")
        q = len(order) // 5
        means = [rr[order[i * q:(i + 1) * q]].mean() for i in range(5)]
        amean = [a[order[i * q:(i + 1) * q]].mean() for i in range(5)]
        ratio = [(rr[order[i * q:(i + 1) * q]] / (rr.mean() * c * a[order[i * q:(i + 1) * q]])).mean() for i in range(5)]
        print("sheet %d: mean multiplicity %.4f  variance %.3f  misses (r = 0): %d (%.3f %%)" % (sh, rr.mean(), rr.var(), int((rr == 0).sum()), 100 * (rr == 0).mean()))
        print("   Pearson correlation of r with sigma(n)/n: %.4f" % corr)
        print("   mean multiplicity by quintile of sigma(n)/n: " + "  ".join("%.3f" % v for v in means) + "   highest/lowest = %.3f" % (means[4] / means[0]))
        print("   mean sigma(n)/n by quintile:                 " + "  ".join("%.4f" % v for v in amean) + "   highest/lowest = %.3f" % (amean[4] / amean[0]))
        print("   mean of r / (mean r * S(N)) by quintile:     " + "  ".join("%.4f" % v for v in ratio) + "   (1.000 everywhere would mean: density = series at this size)")
        # regression slope of r/mean on S
        S = c * a
        slope = np.cov(rr / rr.mean(), S)[0, 1] / S.var()
        print("   least-squares slope of r/mean(r) against S(N): %.4f  (1 if the series were the conditional mean)   mean of S on the window: %.5f" % (slope, S.mean()))
        # by exact value of n small: classes n = 1 (sheet 0) etc.
        if sh == 0:
            for nn in (1, 3, 5, 7, 9, 15, 45, 105, 315):
                sel = (t == nn)
                if sel.any():
                    print("      sheet 0, odd part of N + 1 equal to %4d: %3d values, mean r %.3f, r/mean %.3f, S = %.4f" % (nn, int(sel.sum()), rr[sel].mean(), rr[sel].mean() / rr.mean(), c * ab[(nn - 1) // 2]))
        else:
            for dd in (3, 9, 27, 5, 25, 15, 45, 105, 315, 945):
                sel = (Ns % dd == 0)
                print("      sheet 1, %4d | N: mean r/mean %.4f   mean S on these N: %.4f" % (dd, rr[sel].mean() / rr.mean(), (c * a1[sel]).mean()))


if __name__ == "__main__":
    main()
