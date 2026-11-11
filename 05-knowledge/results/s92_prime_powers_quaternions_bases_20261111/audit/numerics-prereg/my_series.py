# Auditor numerics-prereg: independent recomputation of the "series vs counts" numbers
# (HYP-9665 (B), THM-5134 (E)); written before reading series_vs_counts.py.
# usage: python3 -B my_series.py KHI KLO
# r_sigma(N) = number of ordered triples of the sheet sigma with value N, odd N in [2^KLO, 2^KHI);
# n = odd part of N + [sigma = 0];  abundancy sigma(n)/n = sum_{d | n} 1/d.
import sys
import math
import time
import numpy as np
import numba as nb
from aud_common import psi_table


@nb.njit
def multiplicities(psi, A, B):
    # r[sheet, i] for N = A + 1 + 2 i, A even, A < N < B
    n = (B - A) // 2
    r = np.zeros((2, n), np.int32)
    s = 0
    while True:
        e = 1 if s == 0 else 0
        sh = 0 if s == 0 else 1
        if 1 + (1 << s) - e >= B:
            break
        h1 = 0
        while True:
            x1 = 2 * psi[h1] + 1
            y1 = 2 * h1 + 1
            b = y1 << s
            xmin = (y1 + 1) // 2
            # no partner h2 >= h1 for this or any later h1: N >= x1 x2 + b y2 - e >= xmin^2 + b y1 - e
            # (x1 >= xmin = (y1+1)/2 and x2 >= (y2+1)/2 >= xmin).  [A first version of this audit script used
            # y1 * xmin here, which stops too early and lost 0.01 % of the triples; fixed.]
            if xmin * xmin + b * y1 - e >= B:
                break
            # partners h2 >= h1: N >= x1 (y2+1)/2 + b y2 - e
            h2 = h1
            while True:
                y2 = 2 * h2 + 1
                if x1 * ((y2 + 1) // 2) + b * y2 - e >= B:
                    break
                N = x1 * (2 * psi[h2] + 1) + b * y2 - e
                if N > A and N < B:
                    r[sh, (N - A - 1) // 2] += 2 if h2 > h1 else 1
                h2 += 1
            h1 += 1
        s += 1
    return r


@nb.njit
def abundancy_odd(limit):
    # ab[i] = sigma(n)/n for n = 2 i + 1 <= limit
    m = (limit + 1) // 2
    ab = np.zeros(m, np.float64)
    d = 1
    while d <= limit:
        inv = 1.0 / d
        k = d
        while k <= limit:
            ab[(k - 1) // 2] += inv
            k += 2 * d
        d += 2
    return ab


@nb.njit
def sqfree_product_odd(limit):
    # prod over odd primes p | n of (1 + p/(p^2-1)), n = 2 i + 1 <= limit
    m = (limit + 1) // 2
    pr = np.ones(m, np.float64)
    isc = np.zeros(m, np.bool_)
    for i in range(1, m):
        p = 2 * i + 1
        if not isc[i]:
            f = 1.0 + p / (p * p - 1.0)
            k = p
            while k <= limit:
                pr[(k - 1) // 2] *= f
                if k > p:
                    isc[(k - 1) // 2] = True
                k += 2 * p
    return pr


def odd_part(a):
    a = a.copy()
    while True:
        ev = (a % 2 == 0)
        if not ev.any():
            return a
        a[ev] //= 2


def vp(n, p, cap=12):
    v = np.zeros(n.shape, np.int64)
    t = n.copy()
    for _ in range(cap):
        dv = (t % p == 0)
        if not dv.any():
            break
        v[dv] += 1
        t[dv] //= p
    return v


def main():
    KHI = int(sys.argv[1])
    KLO = int(sys.argv[2])
    A = 1 << KLO
    B = 1 << KHI
    t0 = time.time()
    psi = psi_table(B // 3 + 8)
    r = multiplicities(psi, A, B)
    N = A + 1 + 2 * np.arange((B - A) // 2, dtype=np.int64)
    ab = abundancy_odd(B + 1)
    sq = sqfree_product_odd(B + 1)
    print("odd N in [2^%d, 2^%d): %d values; triples sheet 0: %d, sheet 1: %d; mean r_0 = %.4f, mean r_1 = %.4f  (%.0f s)" % (
        KLO, KHI, len(N), int(r[0].sum()), int(r[1].sum()), r[0].mean(), r[1].mean(), time.time() - t0))
    for sh in range(2):
        n = odd_part(N + 1) if sh == 0 else N.copy()
        a = ab[(n - 1) // 2]
        q = sq[(n - 1) // 2]
        rr = r[sh].astype(np.float64)
        c1 = np.corrcoef(rr, a)[0, 1]
        c2 = np.corrcoef(rr, q)[0, 1]
        print("sheet %d: mean of sigma(n)/n over the window %.4f (pi^2/8 = %.4f); correlation of r with sigma(n)/n: %.4f; with the squarefree-level product: %.4f" % (
            sh, a.mean(), math.pi ** 2 / 8, c1, c2))
        for p in (3, 5):
            v = vp(n, p)
            base = rr[v == 0].mean()
            cells = []
            for k in range(5):
                sel = v == k
                lim = (1 - p ** (-(k + 1.0))) / (1 - 1.0 / p)
                cells.append("v=%d: %.3f (limit %.3f, %d N)" % (k, rr[sel].mean() / base, lim, int(sel.sum())))
            print("    mean of r by v_%d(n), relative to v = 0: %s" % (p, "; ".join(cells)))
        order = np.argsort(a, kind="stable")
        qs = np.array_split(order, 5)
        print("    quintiles of sigma(n)/n: mean abundancy -> mean r: " + "; ".join("%.3f -> %.3f" % (a[ix].mean(), rr[ix].mean()) for ix in qs))
        print("    ratio highest / lowest quintile of mean r: %.3f" % (rr[qs[4]].mean() / rr[qs[0]].mean()))
        # the conjectured single-value law of HYP-9665 (B): r ~ 2 kappa_sigma (8/pi^2) sigma(n)/n ln N
        kap = 0.1615389173 if sh == 0 else 0.2125139092
        model = 2 * kap * (8 / math.pi ** 2) * a * np.log(N.astype(np.float64))
        print("    HYP-9665 (B) model 2 kappa (8/pi^2) sigma(n)/n ln N: mean model %.4f against mean r %.4f; mean |r - model| = %.4f; mean |r - mean r| = %.4f" % (
            model.mean(), rr.mean(), np.abs(rr - model).mean(), np.abs(rr - rr.mean()).mean()))


if __name__ == "__main__":
    main()
