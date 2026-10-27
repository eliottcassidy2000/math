#!/usr/bin/env python3
"""
fast_convergence_density_20261027.py  (opus-2026-10-27-S64)

Where Mazur's positive-density set sits. Mazur (ProofAtlas, v2.1, 2026-09-06; Lean-checked, review pending) proves
that F_C = { n : tau(n) <= C ln n } has positive lower natural density for C = 523/50 = 10.46, where tau is the
ordinary total stopping time (n -> n/2 for even n, n -> 3n+1 for odd n; every application counted). The constant is
the typical-word constant

    C* = 3 / ln(4/3) = 10.4281...

(one Syracuse step of mean valuation 2 costs three ordinary steps and contracts by 3/4; the digest
05-knowledge/results/mazur_positive_density_20260928.md, section 1). This script measures the empirical density
F_X(C) = #{1 <= n < X : tau(n) <= C ln n} / (X - 1) for X = 2^k and compares the fluctuations with the
Lagarias-Weiss renewal model, which has NO free parameter:

    tau(n) ~ C* ln n + sigma sqrt(ln n) Z,   Z standard normal,
    sigma^2 = Var(c - C* delta) / E[delta],  c = 1 + a (ordinary steps per Syracuse step),
    delta = a ln 2 - ln 3 (the drift per Syracuse step), a ~ Geometric(1/2) on {1, 2, ...}:
    Var(a) = 2, E[delta] = 2 ln 2 - ln 3, c - C* delta = const + a (1 - C* ln 2), hence
    sigma^2 = 2 (1 - C* ln 2)^2 / (2 ln 2 - ln 3) = 269.7...,  sigma = 16.42...

So the model predicts, for n in a dyadic block around N, the block density G(C) = Phi((C - C*) sqrt(ln N) / sigma):
a sigmoid centred at C* that sharpens only like sqrt(ln N). Mazur's set at C = 10.46 is the lower half of the
population, and its density should tend to 1/2 (not to 1 and not to 0) under the model.

PRE-REGISTERED PREDICTIONS (fixed before any tau was computed in this session; opus-2026-10-27-S64, 2026-10-27):
  P1. z(n) := (tau(n) - C* ln n) / (sigma sqrt(ln n)) over dyadic blocks [2^(k-1), 2^k): the block standard
      deviation lies in [0.85, 1.15] for every k in 20..28, and the block mean satisfies |mean| <= 0.35 for k >= 20.
  P2. The cumulative density F_X(523/50) lies in [0.45, 0.60] for every X = 2^k, 16 <= k <= 28, and the median of
      tau(n)/ln n over n < 2^k lies within 1.0 of C* for every k >= 20.
  P3. Sigmoid sharpening: the block density G_k(8) is non-increasing and G_k(13) non-decreasing in k for k = 20..28,
      with G_28(8) <= 0.33 and G_28(13) >= 0.70.
  P4. The elementary almost-all floor: the count of n in [2^(k-1), 2^k) with tau(n) < 3.307 ln n has consecutive-block
      ratio < 2 for every k in 20..28 (Proposition F of the note: the count below X is O(X^(1-c))).
  P5. (3x-1 seeds, the sign-blind analogue of Mazur's Lemma 4.1) r_j = (2^(2j+1) + 1)/3 is an odd integer with
      3 r_j - 1 = 2^(2j+1), v_3(2^(2j+1) + 1) = 1 + v_3(2j+1), and r_0, ..., r_(3^q - 1) permute Z/3^q for q = 1..6.
Outcomes are printed as PASS/FAIL against these thresholds.

Universe: all 1 <= n < 2^KMAX (default KMAX = 28); exact 64-bit integer arithmetic (the maximal orbit value met is
printed; it must stay far below 2^63). tau(1) = 0; n = 1 counts as fast for every C >= 0.

Usage: python fast_convergence_density_20261027.py [KMAX]     (KMAX = 28 takes a few minutes and ~0.6 GB)

Resolution note (audit, 2026-10-27): the density columns are read off a histogram of r = tau/ln n with bins of width
0.01 (bin centres <= C), exact for the integer, 10.46 and 10.48 columns; the 3.307 column therefore counts r < 3.31
(4718 against the exact 4699 at k = 20), the floor counts (P4) are exact, and medians have resolution 0.01.
"""
import math
import sys
import time

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover
    def njit(*a, **k):
        def wrap(f):
            return f
        return wrap if not a or not callable(a[0]) else a[0]

LN2 = math.log(2.0)
LN3 = math.log(3.0)
CSTAR = 3.0 / math.log(4.0 / 3.0)
SIGMA = math.sqrt(2.0 * (1.0 - CSTAR * LN2) ** 2 / (2.0 * LN2 - LN3))
FLOOR_C = 3.307          # Proposition F: elementary almost-all lower bound tau(n) >= (3.307 - eps) ln n
C_GRID = [FLOOR_C, 6.0, 8.0, 9.0, 10.0, CSTAR, 523 / 50, 262 / 25, 11.0, 12.0, 13.0, 14.0, 16.0, 20.0, 30.0]
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 28
CHUNK = 1 << 22

RESULTS = []


def verdict(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("  [" + str(detail) + "]") if detail != "" else ""))


@njit(cache=True)
def tau_all(X):
    """Ordinary total stopping times for all n < X, by first descent below n (int16; max tau below 2^28 is < 2^15)."""
    tau = np.zeros(X, dtype=np.int16)
    maxm = 0
    for n in range(2, X):
        m = n
        s = 0
        while m >= n:
            if m & 1:
                m = 3 * m + 1
                if m > maxm:
                    maxm = m
            else:
                m >>= 1
            s += 1
        tau[n] = s + tau[m]
    return tau, maxm


def phi(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def main():
    print("C* = 3/ln(4/3) = %.6f   sigma = %.4f   (model: tau ~ C* ln n + sigma sqrt(ln n) Z)" % (CSTAR, SIGMA))
    X = 1 << KMAX
    t0 = time.time()
    tau, maxm = tau_all(X)
    print("tau(n) for all n < 2^%d in %.1f s; max tau = %d at n = %d; max orbit value %d = 2^%.2f"
          % (KMAX, time.time() - t0, int(tau.max()), int(tau.argmax()), maxm, math.log2(maxm)))
    assert maxm < (1 << 62)

    # histograms of r = tau / ln n per dyadic block (bins of width 0.01 on [0, 400)), z-statistics, floor counts
    nb = 40000
    edges = np.linspace(0.0, 400.0, nb + 1)
    blocks = {}
    for k in range(2, KMAX + 1):
        lo, hi = 1 << (k - 1), 1 << k
        hist = np.zeros(nb, dtype=np.int64)
        zs = 0.0
        zss = 0.0
        floor_cnt = 0
        cnt = 0
        for a in range(lo, hi, CHUNK):
            b = min(hi, a + CHUNK)
            n = np.arange(a, b, dtype=np.float64)
            ln = np.log(n)
            t = tau[a:b].astype(np.float64)
            r = t / ln
            h, _ = np.histogram(r, bins=edges)
            hist += h
            z = (t - CSTAR * ln) / (SIGMA * np.sqrt(ln))
            zs += float(z.sum())
            zss += float((z * z).sum())
            floor_cnt += int(np.count_nonzero(t < FLOOR_C * ln))
            cnt += b - a
        blocks[k] = dict(hist=hist, zmean=zs / cnt, zsd=math.sqrt(max(zss / cnt - (zs / cnt) ** 2, 0.0)),
                         floor=floor_cnt, cnt=cnt)
    centers = 0.5 * (edges[1:] + edges[:-1])

    def frac_le(hist, C):
        # fraction of the histogram mass in bins whose centre is <= C (bin width 0.01)
        return float(hist[centers <= C].sum()) / float(hist.sum())

    print()
    print("Block densities G_k(C) = #{n in [2^(k-1), 2^k) : tau(n) <= C ln n} / 2^(k-1), with the model "
          "Phi((C - C*) sqrt(ln N_k)/sigma), N_k = 2^(k - 1/2):")
    hdr = "  k  " + " ".join("%7.3f" % C for C in C_GRID) + "   median   zmean   zsd   floor#  ratio"
    print(hdr)
    prev_floor = None
    for k in range(10, KMAX + 1):
        bl = blocks[k]
        hist = bl["hist"]
        cum = np.cumsum(hist)
        med_idx = int(np.searchsorted(cum, cum[-1] / 2.0))
        median = centers[med_idx]
        row = "%3d  " % k + " ".join("%7.4f" % frac_le(hist, C) for C in C_GRID)
        ratio = (bl["floor"] / prev_floor) if prev_floor else float("nan")
        row += "   %6.3f  %6.3f  %5.3f  %7d  %5.3f" % (median, bl["zmean"], bl["zsd"], bl["floor"], ratio)
        print(row)
        lnN = (k - 0.5) * LN2
        print("mod  " + " ".join("%7.4f" % phi((C - CSTAR) * math.sqrt(lnN) / SIGMA) for C in C_GRID))
        prev_floor = bl["floor"]

    print()
    print("Cumulative densities F_X(C) = #{1 <= n < X : tau(n) <= C ln n}/(X - 1), X = 2^k (n = 1 counted as fast):")
    print("  k  " + " ".join("%7.3f" % C for C in C_GRID) + "   median r (n < X)")
    cum_hist = np.zeros(nb, dtype=np.int64)
    FX = {}
    for k in range(2, KMAX + 1):
        cum_hist += blocks[k]["hist"]
        if k >= 10:
            total = X_k = (1 << k) - 1
            vals = []
            for C in C_GRID:
                vals.append((float(cum_hist[centers <= C].sum()) + 1.0) / total)
            FX[k] = dict(zip(C_GRID, vals))
            cc = np.cumsum(cum_hist)
            med = centers[int(np.searchsorted(cc, cc[-1] / 2.0))]
            FX[k]["median"] = med
            print("%3d  " % k + " ".join("%7.4f" % v for v in vals) + "   %6.3f" % med)

    print()
    # P1
    ok1 = all(0.85 <= blocks[k]["zsd"] <= 1.15 and abs(blocks[k]["zmean"]) <= 0.35 for k in range(20, KMAX + 1))
    verdict("P1 block sd of z in [0.85, 1.15] and |mean| <= 0.35 for k = 20..%d" % KMAX, ok1,
            "sd: " + ", ".join("%.3f" % blocks[k]["zsd"] for k in range(20, KMAX + 1)) +
            "; mean: " + ", ".join("%.3f" % blocks[k]["zmean"] for k in range(20, KMAX + 1)))
    # P2
    ok2 = all(0.45 <= FX[k][523 / 50] <= 0.60 for k in range(16, KMAX + 1)) and \
        all(abs(FX[k]["median"] - CSTAR) <= 1.0 for k in range(20, KMAX + 1))
    verdict("P2 F_X(523/50) in [0.45, 0.60] for k = 16..%d and |median - C*| <= 1 for k >= 20" % KMAX, ok2,
            "F: " + ", ".join("%.4f" % FX[k][523 / 50] for k in range(16, KMAX + 1)) +
            "; median: " + ", ".join("%.3f" % FX[k]["median"] for k in range(20, KMAX + 1)))
    # P3
    g8 = [frac_le(blocks[k]["hist"], 8.0) for k in range(20, KMAX + 1)]
    g13 = [frac_le(blocks[k]["hist"], 13.0) for k in range(20, KMAX + 1)]
    ok3 = all(g8[i + 1] <= g8[i] + 1e-12 for i in range(len(g8) - 1)) and \
        all(g13[i + 1] >= g13[i] - 1e-12 for i in range(len(g13) - 1)) and g8[-1] <= 0.33 and g13[-1] >= 0.70
    verdict("P3 G_k(8) non-increasing, G_k(13) non-decreasing for k = 20..%d; G_%d(8) <= 0.33, G_%d(13) >= 0.70"
            % (KMAX, KMAX, KMAX), ok3, "G(8): " + ", ".join("%.4f" % v for v in g8) +
            "; G(13): " + ", ".join("%.4f" % v for v in g13))
    # P4
    ratios = [blocks[k]["floor"] / max(blocks[k - 1]["floor"], 1) for k in range(20, KMAX + 1)]
    ok4 = all(r < 2.0 for r in ratios)
    verdict("P4 floor-count block ratio < 2 for k = 20..%d (tau < 3.307 ln n)" % KMAX, ok4,
            "counts: " + ", ".join(str(blocks[k]["floor"]) for k in range(19, KMAX + 1)) +
            "; ratios: " + ", ".join("%.3f" % r for r in ratios))
    # P5 (pure Python, exact)
    ok5 = True
    detail5 = []
    for q in range(1, 7):
        mod = 3 ** q
        seen = set()
        for j in range(mod):
            num = (1 << (2 * j + 1)) + 1
            assert num % 3 == 0
            r = num // 3
            assert r % 2 == 1 and 3 * r - 1 == (1 << (2 * j + 1))
            v = 0
            t = num
            while t % 3 == 0:
                t //= 3
                v += 1
            w = 0
            u = 2 * j + 1
            while u % 3 == 0:
                u //= 3
                w += 1
            if v != 1 + w:
                ok5 = False
            seen.add(r % mod)
        detail5.append("q=%d: %d classes" % (q, len(seen)))
        if len(seen) != mod:
            ok5 = False
    verdict("P5 3x-1 seeds r_j = (2^(2j+1)+1)/3 permute Z/3^q (q = 1..6), v_3(2^(2j+1)+1) = 1 + v_3(2j+1)", ok5,
            "; ".join(detail5))
    print()
    nfail = sum(1 for _, ok in RESULTS if not ok)
    print("%d pre-registered predictions, %d failed" % (len(RESULTS), nfail))
    print("ALL PASSED" if nfail == 0 else "SOME FAILED")


if __name__ == "__main__":
    main()
