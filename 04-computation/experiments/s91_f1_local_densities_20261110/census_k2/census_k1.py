#!/usr/bin/env python3
"""S91 helper (census_k2), step 2: the K = 1 census of THM-4916 by value, in THM-4916's normalisation (M_1 = words with
at most one J = F1 triples weighted s + 1), against
  (a) THM-4916's letter-count limits   L(0) = 1 + (3p-4)/(4p^2), L(-1) = 1 - 1/p^2 + 1/(4p), others 1 - 1/p^2;
  (b) the by-value limits of S91        L(0) = 1 - 1/p^2 + w_1/p, L(-1) = 1 - 1/p^2 + w_0/p, others 1 - 1/p^2,
                                        (w_0, w_1) = (0.19641, 0.80359) (pair normalisation);
  (c) the first-order finite-size prediction of window_split.py, adapted to the weights s + 1:
        W(h1, s) = (1/x1) int_{f0}^{f1} w(f) df  (formula of window_split.py, reused verbatim),
        T_0 = sum_h1 W(h1, 0),  T_1 = sum_h1 sum_{s >= 1} (s + 1) W(h1, s),  f_sheet = same sums over the p-divisible h1,
        D_sheet = 2 f_sheet - 1/p^2,  class law on a sheet: p P(c) = 1 - D + p D [c = c_sheet],
        L(c) = sum over sheets of share * p P(c).   No window data enters (c).
Usage: python3 census_k1.py [LOG_LO LOG_HI]      (default 20 24)
Checks (exit 1 on failure): the word census equals the weighted triple count for EVERY odd N of the window; the
confinement lemma for every confined word; #S-words = #T-words on each sheet (h1 <-> h2 symmetry); (c) within the
stated tolerances; the dilogarithm identities.
"""
import math
import sys
import time

import numpy as np
from numba import njit, prange

import rc_census as rc

PRIMES = [3, 5, 7, 11, 13, 17]


@njit(cache=True)
def triples_weighted(Q, A, B):
    """sum over F1 triples (h1, h2, s) with value N in [A, B) of the weight s + 1, per odd N; also plain counts by sheet."""
    size = (B - A) // 2
    r = np.zeros(size, np.int32)
    tot = np.zeros(2, np.int64)
    H = Q.shape[0]
    s = 0
    while (1 << s) < B:
        d0 = 1 if s == 0 else 0
        sheet = 0 if s == 0 else 1
        for h1 in range(H):
            y1 = 2 * h1 + 1
            Y = y1 << s
            if Y >= B:
                break
            x1 = y1 - Q[h1]
            y2max = (2 * (B - 1 + d0)) // (2 * Y + x1)
            if y2max < 1:
                continue
            h2max = (y2max - 1) // 2 + 1
            if h2max >= H:
                h2max = H - 1
            h2min = ((A // (Y + x1)) - 1) // 2 - 1
            if h2min < 0:
                h2min = 0
            for h2 in range(h2min, h2max + 1):
                y2 = 2 * h2 + 1
                x2 = y2 - Q[h2]
                N = x1 * x2 + Y * y2 - d0
                if N >= A and N < B:
                    r[(N - A) >> 1] += s + 1
                    tot[sheet] += 1
        s += 1
    return r, tot


@njit(cache=True)
def w_closed(x1, Y, A, B):
    """window_split.py: number of second states with A <= N < B for the first state x1 and Y = 2^s y1 (rho uniform)."""
    f0 = 0.5 * x1 + Y
    f1 = 1.0 * x1 + Y
    w = 0.0
    if f0 < A:
        top = f1 if f1 < A else A
        w += (B - A) * math.log(top / f0)
    lo = f0 if f0 > A else A
    hi = f1 if f1 < B else B
    if lo < hi:
        w += B * math.log(hi / lo)
    return w / x1


@njit(parallel=True, cache=True)
def predicted(Q, A, B, primes, nchunks):
    """T[sheet, 0] = triple mass, T[sheet, 1] = pair mass (weights s + 1); F the same over the p-divisible h1."""
    size = Q.shape[0]
    npr = primes.shape[0]
    T = np.zeros((nchunks, 2, 2))
    F = np.zeros((nchunks, npr, 2, 2))
    step = size // nchunks + 1
    for c in prange(nchunks):
        for h1 in range(c * step, min(size, (c + 1) * step)):
            y1 = 2 * h1 + 1
            x1 = y1 - Q[h1]
            mk = 0
            for j in range(npr):
                p = primes[j]
                if y1 % p == 0 and x1 % p == 0:
                    mk |= (1 << j)
            s = 0
            while True:
                Y = y1 << s
                if Y >= B:
                    break
                w = w_closed(x1, Y, A, B)
                sheet = 0 if s == 0 else 1
                T[c, sheet, 0] += w
                T[c, sheet, 1] += (s + 1) * w
                if mk != 0:
                    for j in range(npr):
                        if (mk >> j) & 1:
                            F[c, j, sheet, 0] += w
                            F[c, j, sheet, 1] += (s + 1) * w
                s += 1
    return T.sum(axis=0), F.sum(axis=0)


def dilog_constants():
    import mpmath as mp
    mp.mp.dps = 30
    L2 = mp.polylog(2, mp.mpf(-1) / 2)
    L4 = mp.polylog(2, mp.mpf(-1) / 4)
    gbar = mp.pi ** 2 / 6 + 2 * L2
    g0 = 2 * (mp.pi ** 2 / 12 + 2 * L2 - L4)
    g1 = 2 * (L4 - L2)
    pair = mp.pi ** 2 / 6
    # independent quadrature of the defining integrals (rho uniform on [1/2, 1], density 2)
    q0 = mp.quad(lambda r: 2 * mp.log((1 + r) / (r / 2 + 1)) / r, [mp.mpf(1) / 2, 1])
    qall = mp.quad(lambda r: 2 * mp.log(1 + r) / r, [mp.mpf(1) / 2, 1])
    qpair = mp.quad(lambda r: 2 * mp.fsum(mp.log(1 + r * mp.mpf(2) ** (-s)) for s in range(0, 110)) / r, [mp.mpf(1) / 2, 1])
    ok = abs(q0 - g0) < mp.mpf(10) ** (-20) and abs(qall - gbar) < mp.mpf(10) ** (-20) and abs(qpair - pair) < mp.mpf(10) ** (-12)
    ok = ok and abs(g0 + g1 - gbar) < mp.mpf(10) ** (-25)
    return dict(L2=L2, L4=L4, gbar=gbar, g0=g0, g1=g1, pair=pair, w0_tri=g0 / gbar, w1_tri=g1 / gbar,
                w0_pair=g0 / pair, w1_pair=1 - g0 / pair, ok=bool(ok))


def main():
    lo = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    hi = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    A, B = 1 << lo, 1 << hi
    ok = True
    c = dilog_constants()
    print("dilogarithm constants (mpmath, 30 digits; closed forms checked against quadrature of the defining integrals):")
    print(f"  Li2(-1/2) = {float(c['L2']):.12f}   Li2(-1/4) = {float(c['L4']):.12f}")
    print(f"  triples: g_bar = pi^2/6 + 2 Li2(-1/2) = {float(c['gbar']):.10f};  g_0 = pi^2/6 + 4 Li2(-1/2) - 2 Li2(-1/4) = "
          f"{float(c['g0']):.10f};  g_1 = 2 (Li2(-1/4) - Li2(-1/2)) = {float(c['g1']):.10f}")
    print(f"  triple shares  w_0 = {float(c['w0_tri']):.10f}  w_1 = {float(c['w1_tri']):.10f}")
    print(f"  pairs (weight s + 1): total pi^2/6 = {float(c['pair']):.10f};  w_0 = 1 + (12/pi^2)(2 Li2(-1/2) - Li2(-1/4)) = "
          f"{float(c['w0_pair']):.10f}   w_1 = {float(c['w1_pair']):.10f}      (letter-count: 1/4, 3/4)")
    ok &= c["ok"]
    w0, w1 = float(c["w0_pair"]), float(c["w1_pair"])

    t0 = time.time()
    res = rc.census(1, A, B, PRIMES, nthreads=3)
    t1 = time.time()
    tot = res["tot"]
    nodd = (B - A) // 2
    print(f"\nK = 1 word census on odd N in [2^{lo}, 2^{hi}): {res['nwords']} counted words, mean M_1 = {res['nwords'] / nodd:.4f}, "
          f"max M_1 = {res['maxM']}, {res['ntasks']} tasks, {t1 - t0:.1f} s")
    ok &= res["overflow_ok"] and res["lemma_violations"] == 0
    print(f"  confinement lemma (value = -[c = 0] E mod p for every word with S or T): violations {res['lemma_violations']}")

    Q = rc.build_Q(B // 2 + 4)
    r_tri, tri_tot = triples_weighted(Q, A, B)
    same = bool(np.array_equal(r_tri.astype(np.int64), tot))
    print(f"  dictionary on the whole window: M_1(N) = sum over triples of (s + 1) for every odd N: {same}   "
          f"({int(tri_tot.sum())} triples; sheet s = 0: {int(tri_tot[0])}, s >= 1: {int(tri_tot[1])})")
    ok &= same

    allcnt = res["allcnt"]          # [nJ, parity, z]
    words = int(allcnt.sum())
    sheet0 = int(allcnt[:2, 0, 0].sum())       # even tail and E = 1 (the J-free words carry the virtual E = 1)
    sheet1 = words - sheet0
    free = int(allcnt[0].sum())
    primes = np.array(PRIMES, dtype=np.int64)
    Tp, Fp = predicted(Q, A, B, primes, 48)
    tp0, tp1 = Tp[0, 1], Tp[1, 1]
    print(f"  sheets (pair normalisation): measured s = 0: {sheet0} ({sheet0 / words:.5f}), s >= 1: {sheet1} ({sheet1 / words:.5f});  "
          f"J-free words {free} ({free / words:.5f}; (B - A) ln 2 = {(B - A) * math.log(2):.0f})")
    print(f"     predicted finite-size (W):  s = 0: {tp0:.0f} ({tp0 / (tp0 + tp1):.5f}), s >= 1: {tp1:.0f} ({tp1 / (tp0 + tp1):.5f});  "
          f"total {tp0 + tp1:.0f} against {words} (ratio {(tp0 + tp1) / words:.5f});  limit shares {w0:.5f}, {w1:.5f}")
    print(f"     triples (unweighted): measured share s = 0: {tri_tot[0] / tri_tot.sum():.5f}, predicted {Tp[0, 0] / Tp[:, 0].sum():.5f}, "
          f"limit {float(c['w0_tri']):.5f}")
    # window_split.py's own tolerance is 2 % per sheet; in pair weights the two sheet errors (+0.8 %, -0.9 %) no longer
    # cancel in the total.  W is exact for small h1 and biased for large h1 (few, small second states: rho not uniform).
    ok &= abs((tp0 + tp1) / words - 1) < 0.01 and abs(tp0 / (tp0 + tp1) - sheet0 / words) < 0.004
    ok &= abs(tp0 / sheet0 - 1) < 0.02 and abs(tp1 / sheet1 - 1) < 0.02
    ok &= abs(Tp[0, 0] / Tp[:, 0].sum() - tri_tot[0] / tri_tot.sum()) < 0.004

    cnt = res["cnt"]                # [prime, nJ, S, T, C1, st]
    print("\n  p-divisible shares by sheet (pair weights): measured first-state share = measured second-state share exactly "
          "(h1 <-> h2 symmetry), against f predicted from W")
    print(f"  {'p':>3} {'sheet':>6} {'f first':>10} {'f second':>10} {'f pred':>10} {'1/p^2':>9} {'both':>10} {'f^2':>10}")
    fmeas = {}
    for j, p in enumerate(PRIMES):
        for sheet in (0, 1):
            sts = (2,) if sheet == 0 else (0, 1)
            nS = sum(int(cnt[j, :2, 1, :, :, st].sum()) for st in sts)
            nT = sum(int(cnt[j, :2, :, 1, :, st].sum()) for st in sts)
            nST = sum(int(cnt[j, :2, 1, 1, :, st].sum()) for st in sts)
            den = sheet0 if sheet == 0 else sheet1
            fpred = Fp[j, sheet, 1] / Tp[sheet, 1]
            fmeas[(p, sheet)] = nS / den
            ok &= nS == nT
            print(f"  {p:3d} {('s=0' if sheet == 0 else 's>=1'):>6} {nS / den:10.6f} {nT / den:10.6f} {fpred:10.6f} {1 / p ** 2:9.6f} "
                  f"{nST / den:10.6f} {(nS / den) ** 2:10.6f}")
            ok &= abs(fpred - nS / den) < 0.02 * (nS / den) + 0.0005

    print("\n  class factors L(c) = (mean of M_1 over N = c mod p) / (overall mean);  'other' = mean over the classes != 0, -1")
    print(f"  {'p':>3} {'class':>6} {'measured':>9} {'(a) letter':>10} {'(b) value':>10} {'(c) finite':>10} {'(c) - meas':>10}")
    worst = 0.0
    for j, p in enumerate(PRIMES):
        means, n = rc.class_means(tot, A, B, p)
        overall = tot.sum() / nodd
        L = means / overall
        others = [r for r in range(p) if r not in (0, p - 1)]
        meas = {"0": L[0], "-1": L[p - 1], "other": float(np.mean(L[others]))}
        a = {"0": 1 + (3 * p - 4) / (4 * p * p), "-1": 1 - 1 / p ** 2 + 1 / (4 * p), "other": 1 - 1 / p ** 2}
        b = {"0": 1 - 1 / p ** 2 + w1 / p, "-1": 1 - 1 / p ** 2 + w0 / p, "other": 1 - 1 / p ** 2}
        f0 = Fp[j, 0, 1] / Tp[0, 1]
        f1 = Fp[j, 1, 1] / Tp[1, 1]
        D0, D1 = 2 * f0 - 1 / p ** 2, 2 * f1 - 1 / p ** 2
        s0, s1 = tp0 / (tp0 + tp1), tp1 / (tp0 + tp1)
        cpred = {"0": s0 * (1 - D0) + s1 * (1 + (p - 1) * D1), "-1": s0 * (1 + (p - 1) * D0) + s1 * (1 - D1),
                 "other": s0 * (1 - D0) + s1 * (1 - D1)}
        for key in ("0", "-1", "other"):
            print(f"  {p:3d} {key:>6} {meas[key]:9.4f} {a[key]:10.4f} {b[key]:10.4f} {cpred[key]:10.4f} {cpred[key] - meas[key]:+10.4f}")
            if p not in (3, 7):
                worst = max(worst, abs(cpred[key] - meas[key]))
        spread = float(np.max(np.abs(L[others] - meas["other"]))) if len(others) > 1 else 0.0
        print(f"      ripple of the other classes (max |L - mean|): {spread:.4f};  measured excess of 0 over other "
              f"{meas['0'] / meas['other'] - 1:+.4f}, of -1 over other {meas['-1'] / meas['other'] - 1:+.4f}")
    print(f"  worst |(c) - measured| over p = 5, 11, 13, 17 (fast primes): {worst:.4f}   (3 and 7 are the slow primes of "
          f"twisted_spectrum.py: rates 0.866, 0.932)")
    ok &= worst < 0.012
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
