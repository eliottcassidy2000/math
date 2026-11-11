#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9665 (B) / THM-5134 (E) numerical paragraph. Own recomputation of
# r_sigma(N) for the odd N of [2^(E-1), 2^E) and its relation to sigma(n)/n, n = odd part of N + [sigma = 0].
# usage: python3 series_counts.py E     (the session's file is E = 25)
import sys, time, math
import numpy as np
from numba import njit, prange
from f1core import build_psi


@njit(cache=True)
def fill(X, A, xs, ys, r0, r1):
    # r_sigma[(N - A) / 2] for odd N in [A, X)
    nH = xs.shape[0]
    tot0 = 0
    tot1 = 0
    for s in range(0, 64):
        p2 = np.int64(1) << s
        sub = 1 if s == 0 else 0
        if 1 + p2 - sub >= X:
            break
        # parallel over h1; increments on a shared array need care -> use per-s sequential over h1 blocks by
        # collecting into thread-safe manner: here we just run h1 sequentially (numba atomic adds are not available)
        for h1 in range(nH):
            x1 = xs[h1]
            y1 = ys[h1]
            num = X - 1 + sub - x1 - p2 * y1
            if num < 0:
                # x1 + p2*y1 grows with h1 only on average; a safe exit needs the lower bound x1 >= h1 + 1
                if (h1 + 1) + p2 * y1 - sub >= X:
                    break
                continue
            h2max = num // (x1 + 2 * p2 * y1)
            if h2max >= nH:
                h2max = nH - 1
            b = p2 * y1
            for h2 in range(h2max + 1):
                N = x1 * xs[h2] + b * ys[h2] - sub
                if N < X:
                    if s == 0:
                        tot0 += 1
                        if N >= A:
                            r0[(N - A) >> 1] += 1
                    else:
                        tot1 += 1
                        if N >= A:
                            r1[(N - A) >> 1] += 1
    return tot0, tot1


@njit(cache=True)
def abundancy_odd(L):
    # ab[i] = sigma(n)/n for odd n = 2i+1 < L
    m = L // 2
    ab = np.zeros(m, dtype=np.float64)
    for i in range(m):
        d = 2 * i + 1
        inv = 1.0 / d
        k = d
        while k < L:
            ab[(k - 1) >> 1] += inv
            k += 2 * d
    return ab


@njit(cache=True)
def sqfree_level(L):
    # prod over odd primes p | n of (1 + p/(p^2-1)) for odd n < L (the squarefree-level product without 8/pi^2)
    m = L // 2
    pr = np.ones(m, dtype=np.float64)
    isp = np.ones(L, dtype=np.uint8)
    for p in range(3, L, 2):
        if isp[p]:
            f = 1.0 + p / (p * p - 1.0)
            k = p
            while k < L:
                pr[(k - 1) >> 1] *= f
                if k > p:
                    isp[k] = 0
                k += 2 * p
    return pr


def corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    return float((a * b).sum() / math.sqrt((a * a).sum() * (b * b).sum()))


def main():
    E = int(sys.argv[1])
    X = 1 << E
    A = 1 << (E - 1)
    t0 = time.time()
    nH = (X - 1) // 3 + 2
    psi = build_psi(nH)
    h = np.arange(nH, dtype=np.int64)
    xs = 2 * psi + 1
    ys = 2 * h + 1
    nodd = (X - A) // 2
    r0 = np.zeros(nodd, dtype=np.int32)
    r1 = np.zeros(nodd, dtype=np.int32)
    tot0, tot1 = fill(X, A, xs, ys, r0, r1)
    print("# series_counts.py E = %d (auditor hypotheses; own code)  %.1f s" % (E, time.time() - t0))
    print("triples with N < 2^%d: sheet 0 %d, sheet 1 %d;  in the window: %d, %d" % (E, tot0, tot1, r0.sum(), r1.sum()))
    N = A + 1 + 2 * np.arange(nodd, dtype=np.int64)
    ab = abundancy_odd(X + 2)
    sq = sqfree_level(X + 2)
    # sheet 1: n = N ; sheet 0: n = odd part of N + 1
    n1 = N
    n0 = N + 1
    while True:
        ev = (n0 & 1) == 0
        if not ev.any():
            break
        n0 = np.where(ev, n0 >> 1, n0)
    for sh, r, n in ((0, r0, n0), (1, r1, n1)):
        a = ab[(n - 1) >> 1]
        s2 = sq[(n - 1) >> 1]
        rf = r.astype(np.float64)
        lnN = np.log(N.astype(np.float64))
        kappa = (0.1615389173, 0.2125139092)[sh]
        print("sheet %d: mean r = %.4f; mean of sigma(n)/n = %.4f (pi^2/8 = %.4f); misses (r = 0): %d (%.4f %%)" % (
            sh, rf.mean(), a.mean(), math.pi ** 2 / 8, int((r == 0).sum()), 100.0 * (r == 0).mean()))
        print("    corr(r, sigma(n)/n) = %.4f ; corr(r, squarefree-level product) = %.4f" % (corr(rf, a), corr(rf, s2)))
        # by v_3, v_5
        for p in (3, 5):
            v = np.zeros(nodd, dtype=np.int64)
            t = n.copy()
            while True:
                dv = (t % p) == 0
                if not dv.any():
                    break
                v += dv
                t = np.where(dv, t // p, t)
            base = rf[v == 0].mean()
            cells = []
            for vv in range(0, 5):
                sel = v == vv
                lim = sum(p ** -i for i in range(vv + 1))
                cells.append("v=%d: %.3f (lim %.3f, %d N)" % (vv, rf[sel].mean() / base, lim, int(sel.sum())))
            sel = v >= 2
            cells.append("v>=2: %.3f" % (rf[sel].mean() / base))
            print("    mean r by v_%d(n) relative to v=0: " % p + "; ".join(cells))
        order = np.argsort(a, kind="stable")
        qs = np.array_split(order, 5)
        print("    quintiles of sigma(n)/n: " + "; ".join("%.3f -> %.3f" % (a[ix].mean(), rf[ix].mean()) for ix in qs))
        # HYP-9665 (B): the L1 distance to 2 kappa S ln N, and to variants
        S = (8 / math.pi ** 2) * a
        pred = 2 * kappa * S * lnN
        l1 = np.abs(rf - pred).mean()
        pred2 = rf.mean() * a / a.mean()
        l1b = np.abs(rf - pred2).mean()
        l1c = np.abs(rf - rf.mean()).mean()
        print("    (B): mean|r - 2 kappa S ln N| = %.4f ; mean|r - rbar S/Sbar| = %.4f ; mean|r - rbar| = %.4f ; mean(2 kappa ln N) = %.4f ; rbar = %.4f" % (
            l1, l1b, l1c, (2 * kappa * lnN).mean(), rf.mean()))
        print("         ratios to rbar: %.4f, %.4f, %.4f ;  Var(r)/mean(r) = %.3f ; Var(r)/mean^2 = %.4f ; Poisson E|r-mean|/mean at this mean = %.4f" % (
            l1 / rf.mean(), l1b / rf.mean(), l1c / rf.mean(), rf.var() / rf.mean(), rf.var() / rf.mean() ** 2,
            poisson_mad(rf.mean()) / rf.mean()))


def poisson_mad(lam):
    # E|X - lam| for Poisson(lam) = 2 lam^(floor(lam)+1) e^-lam / floor(lam)!
    k = int(math.floor(lam))
    return 2 * math.exp(-lam + (k + 1) * math.log(lam) - math.lgamma(k + 1))


if __name__ == "__main__":
    main()
