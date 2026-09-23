#!/usr/bin/env python3
"""collatz_procgen_20260922_ladder_ballot.py -- no-choice exceptional sets of q x + 1.

Map (q odd): T_q(x) = x/2 (x even), (q x + 1)/2 (x odd), on Z_2.
Bad(q) = {x : q^{a_s(x)} > 2^s for every s >= 1}, a_s = number of odd terms among
x, T x, ..., T^{s-1} x (multiplicative non-descent = infinite coefficient stopping time).

By the parity-vector isometry (Lemma 1 of the sibling-dimension-ladder note) the number
N_m(q) of classes mod 2^m that do not descend within m steps is EXACTLY the number of
0/1 words of length m all of whose prefixes (length s = 1..m) contain more than
s*log_q(2) ones.  This script

  A. computes N_m(q) exactly (Python big integers) and checks it against the q=3 counts
     of lane one (collatz_procgen_20260922_choice_ladder.out, section 1) and a direct
     parity-vector bijection check mod 2^k;
  B. q=3: float DP (renormalised) up to m = M3 and a fit of
     log2 N_m = alpha*m + beta*log2(m) + gamma  (prediction alpha = h(log_3 2), beta = -3/2),
     with the residual collapse onto a function of theta_m = frac(m*log_3 2);
  C. q>=5: exact f_m = N_m/2^m and a rigorous enclosure of mu_q = Haar(Bad(q)) = lim f_m,
     from 0 <= f_m - mu_q <= rho^{m+1}/(1-rho), rho = 2^{t-1}(1+q^{-t}) at a fixed t>0;
  D. an independent evaluation of mu_q by the Sparre Andersen / Spitzer identity
     mu_q = exp(-sum_{n>=1} P(S_n<0)/n),  P(S_n<0) = 2^{-n} sum_{k<n log_q 2} C(n,k).

Usage: python3 collatz_procgen_20260922_ladder_ballot.py [M3]   (default M3 = 100000)
"""
import math
import os
import re
import sys
import time

import mpmath
import numpy as np

mpmath.mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
LANE1_OUT = os.path.join(ROOT, "05-knowledge", "results",
                         "collatz_procgen_20260922_choice_ladder.out")


def binary_entropy(p):
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def amin_table(q, M):
    """amin[s] = least a with q^a > 2^s (exact, via bit lengths), s = 1..M; amin[0] = 0."""
    res = [0] * (M + 1)
    a, qa = 0, 1
    for s in range(1, M + 1):
        while qa.bit_length() <= s:  # q^a < 2^s  (q^a = 2^s is impossible for s >= 1)
            qa *= q
            a += 1
        res[s] = a
    return res


def exact_counts(q, M):
    """N[m] for m = 0..M (N[0] = 1), exact big integers."""
    am = amin_table(q, M)
    lo, vec = 0, [1]  # vec[i] = number of surviving words with a = lo + i ones
    N = [1]
    for s in range(1, M + 1):
        new = vec + [0]
        for i in range(len(vec)):
            new[i + 1] += vec[i]
        cut = am[s] - lo
        if cut > 0:
            new = new[cut:]
            lo = am[s]
        vec = new
        N.append(sum(vec))
    return N


def float_log2_counts(q, M):
    """log2 N[m], m = 0..M, by a renormalised float64 DP (all additions positive)."""
    am = amin_table(q, M)
    lo, vec, logscale = 0, np.array([1.0]), 0.0
    out = np.zeros(M + 1)
    for s in range(1, M + 1):
        new = np.empty(len(vec) + 1)
        new[:-1] = vec
        new[-1] = 0.0
        new[1:] += vec
        cut = am[s] - lo
        if cut > 0:
            new = new[cut:]
            lo = am[s]
        vec = new
        if s % 32 == 0:
            mx = vec.max()
            vec = vec / mx
            logscale += math.log2(mx)
        out[s] = logscale + math.log2(vec.sum())
    return out


def parity_bijection_check(q, K):
    """Lemma 1 check: for k <= K, x mod 2^k -> first k parities of T_q is a bijection."""
    for k in range(1, K + 1):
        seen = bytearray(1 << k)
        for x in range(1 << k):
            y, w = x, 0
            for i in range(k):
                b = y & 1
                w |= b << i
                y = (q * y + 1) // 2 if b else y // 2
            if seen[w]:
                return False, k
            seen[w] = 1
    return True, K


def lane1_q3_counts():
    """Parse the Collatz-T counts (section 1) from lane one's output file."""
    vals = {}
    if not os.path.exists(LANE1_OUT):
        return vals
    with open(LANE1_OUT) as fh:
        txt = fh.read()
    block = txt.split("### 1.")[1].split("Greedy G")[0]
    for m, bad in re.findall(r"m=\s*(\d+) bad=(\d+)", block):
        vals[int(m)] = int(bad)
    return vals


def rho_of(q, t):
    t = mpmath.mpf(t)
    return mpmath.power(2, t - 1) * (1 + mpmath.power(q, -t))


def best_t(q):
    f = lambda t: float(rho_of(q, t))
    lo, hi = 1e-6, 5.0
    for _ in range(200):  # golden-section search (rho is log-convex in t)
        a = lo + (hi - lo) * 0.381966
        b = lo + (hi - lo) * 0.618034
        if f(a) < f(b):
            hi = b
        else:
            lo = a
    return round((lo + hi) / 2, 4)  # a fixed rational t; any t>0 gives a valid bound


def spitzer_mu(q, Nmax):
    """exp(-sum_{n<=Nmax} P(S_n<0)/n) with exact binomial partial sums; returns (value, partial sum)."""
    am = amin_table(q, Nmax)
    # F = sum_{j<=k} C(n,j) with k = k_n = am[n]-1, maintained exactly; Cnk = C(n,k_n)
    # start n = 1: k_1 = am[1]-1
    n = 1
    k = am[1] - 1
    Cnk = math.comb(1, k) if k >= 0 else 0
    F = sum(math.comb(1, j) for j in range(0, k + 1))
    total = mpmath.mpf(0)
    total += mpmath.mpf(F) / mpmath.mpf(2) ** n / n
    for n in range(1, Nmax):
        # move to n+1 keeping k: F(n+1,k) = 2F(n,k) - C(n,k); C(n+1,k) = C(n,k)(n+1)/(n+1-k)
        F = 2 * F - Cnk
        Cnk = Cnk * (n + 1) // (n + 1 - k)
        knew = am[n + 1] - 1
        while k < knew:  # C(n+1,k+1) = C(n+1,k)(n+1-k)/(k+1)
            Cnk = Cnk * (n + 1 - k) // (k + 1)
            k += 1
            F += Cnk
        assert k == knew
        total += mpmath.mpf(F) / mpmath.mpf(2) ** (n + 1) / (n + 1)
    return mpmath.e ** (-total), total


def main():
    M3 = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    t0 = time.time()
    print("=" * 78)
    print("LADDER-A. exact no-choice counts N_m(q) (ballot DP) and Lemma 1 checks")
    print("=" * 78)
    for q in (3, 5, 7):
        ok, K = parity_bijection_check(q, 14)
        print(f"Lemma 1 (parity-vector bijection mod 2^k, k<=14) for q={q}: {'PASS' if ok else 'FAIL at k=%d' % K}")
    N3 = exact_counts(3, 2000)
    lane = lane1_q3_counts()
    if lane:
        agree = all(N3[m] == lane[m] for m in lane)
        print(f"q=3: ballot DP N_m vs lane-one Collatz-T counts, m=1..{max(lane)}: "
              f"{'IDENTICAL' if agree else 'MISMATCH'} ({len(lane)} levels)")
        if not agree:
            sys.exit("mismatch with lane one")
    else:
        print("q=3: lane-one output not found; skipping comparison")
    for q in (3, 5, 7, 9):
        N = exact_counts(q, 26)
        row = " ".join(f"{m}:{N[m]}" for m in (16, 20, 24, 26))
        print(f"q={q}: N_m for m=16,20,24,26 -> {row}")
    print("(these are compared with the lane-one C program in collatz_procgen_20260922_ladder_choice.py)")

    print()
    print("=" * 78)
    print(f"LADDER-B. q=3 dimension: float DP to m={M3}; model log2 N_m = alpha m + beta log2 m + gamma")
    print("=" * 78)
    p3 = math.log(2) / math.log(3)
    h3 = binary_entropy(p3)
    print(f"p = log_3 2 = {p3:.12f};  h(p) = {h3:.12f} bits (predicted Hausdorff dimension)")
    L = float_log2_counts(3, M3)
    # float vs exact cross-check
    worst = max(abs(L[m] - math.log2(N3[m])) for m in range(1, 2001))
    print(f"float DP vs exact big-int DP, m=1..2000: max |log2 difference| = {worst:.2e}")
    print(f"{'m':>7} {'log2 N_m':>14} {'log2N/m':>9} {'(log2N+1.5log2m)/m':>19} {'R_m=log2N-hm+1.5log2m':>22} {'theta_m':>8}")
    for m in (10, 20, 26, 50, 100, 300, 1000, 3000, 10000, 30000, 100000):
        if m <= M3:
            R = L[m] - h3 * m + 1.5 * math.log2(m)
            th = (m * p3) % 1.0
            print(f"{m:>7} {L[m]:>14.4f} {L[m]/m:>9.5f} {(L[m]+1.5*math.log2(m))/m:>19.5f} {R:>22.5f} {th:>8.4f}")

    def fit(ms, fix_beta=None):
        ms = np.array(ms, dtype=float)
        y = np.array([L[int(m)] for m in ms])
        if fix_beta is None:
            A = np.vstack([ms, np.log2(ms), np.ones_like(ms)]).T
            coef, *_ = np.linalg.lstsq(A, y, rcond=None)
            return coef
        A = np.vstack([ms, np.ones_like(ms)]).T
        coef, *_ = np.linalg.lstsq(A, y - fix_beta * np.log2(ms), rcond=None)
        return np.array([coef[0], fix_beta, coef[1]])

    print("least-squares fits (alpha should be h = %.5f, beta should be -1.5):" % h3)
    for lo_, hi_, fb in ((5, 26, None), (5, 26, -1.5), (10, 26, 0.0), (100, 1000, None),
                         (1000, 10000, None), (10000, M3, None), (10000, M3, -1.5)):
        if hi_ <= M3:
            c = fit(range(lo_, hi_ + 1), fb)
            tag = "free beta" if fb is None else f"beta fixed {fb:+.1f}"
            print(f"  m in [{lo_},{hi_}] ({tag}): alpha={c[0]:.5f} beta={c[1]:+.4f} gamma={c[2]:+.4f}")
    # residual collapse onto a function C(theta) of theta_m = frac(m log_3 2)
    windows = [(1000, 5000), (5000, 20000)] + ([(M3 // 2, M3)] if M3 >= 40000 else [])
    for (lo_, hi_) in windows:
        ms = np.arange(lo_, hi_ + 1)
        R = L[ms] - h3 * ms + 1.5 * np.log2(ms)
        th = (ms * p3) % 1.0
        for nb in (20, 400):
            bins = np.minimum((th * nb).astype(int), nb - 1)
            spreads = [R[bins == b].max() - R[bins == b].min() for b in range(nb) if (bins == b).any()]
            print(f"  residual R_m, m in [{lo_},{hi_}], {nb:>3} theta-bins: overall range {R.max()-R.min():.4f}; "
                  f"within-bin range median {float(np.median(spreads)):.4f}, max {max(spreads):.4f}")
    # the large within-bin ranges sit at jump discontinuities of C(theta) at theta = frac(j log_3 2)
    lo_ = M3 // 2 if M3 >= 40000 else 5000
    ms = np.arange(lo_, M3 + 1)
    R = L[ms] - h3 * ms + 1.5 * np.log2(ms)
    th = (ms * p3) % 1.0
    print(f"  jumps of C(theta) (m in [{lo_},{M3}], means over theta-windows of width 0.002 on each side):")
    for j in range(1, 7):
        c = (j * p3) % 1.0
        below = R[(th > c - 0.002) & (th < c)]
        above = R[(th > c) & (th < c + 0.002)]
        print(f"    theta = frac({j} log_3 2) = {c:.4f}: C- = {below.mean():.4f}, C+ = {above.mean():.4f}, "
              f"jump {above.mean() - below.mean():+.4f}")
    nb = 400
    bins = np.minimum((th * nb).astype(int), nb - 1)
    prof = np.array([R[bins == b].mean() for b in range(nb)])
    gam = float(prof.mean())
    print(f"  theta-average of C(theta) (log2) = {gam:.4f}; range of the 400-bin profile "
          f"[{prof.min():.4f}, {prof.max():.4f}]")
    print("  finite-size part F_m = R_m - C(theta_m) (C from the 400-bin profile above):")
    for m in (26, 50, 100, 300, 1000, 3000, 10000):
        b = min(int(((m * p3) % 1.0) * nb), nb - 1)
        Rm = L[m] - h3 * m + 1.5 * math.log2(m)
        print(f"    m={m:>6}: R_m={Rm:.4f}  C(theta_m)={prof[b]:.4f}  F_m={Rm - prof[b]:+.4f}")
    R26 = L[26] - h3 * 26 + 1.5 * math.log2(26)
    print(f"  decomposition at m=26: log2(N)/m = {L[26]/26:.4f} = h {h3:.4f} - 1.5 log2(26)/26 "
          f"{-1.5*math.log2(26)/26:.4f} + R_26/26 {R26/26:+.4f}")
    for m in (26, 100, 1000, 10000, M3):
        print(f"    corrected exponent (log2 N_m + 1.5 log2 m - {gam:.3f})/m at m={m}: "
              f"{(L[m] + 1.5*math.log2(m) - gam)/m:.6f}   (naive log2 N_m/m = {L[m]/m:.6f})")
    # exponent along fixed-theta subsequences
    ms = np.arange(2000, M3 + 1)
    th = (ms * p3) % 1.0
    for c0 in (0.1, 0.5, 0.9):
        sel = ms[np.abs(th - c0) < 0.01]
        A = np.vstack([sel, np.log2(sel), np.ones_like(sel, dtype=float)]).T
        coef, *_ = np.linalg.lstsq(A, L[sel], rcond=None)
        A4 = np.vstack([sel, np.log2(sel), np.ones_like(sel, dtype=float), 1.0 / sel]).T
        coef4, *_ = np.linalg.lstsq(A4, L[sel], rcond=None)
        print(f"  fixed theta ~ {c0} (|theta-{c0}|<0.01, m in [2000,{M3}], {len(sel)} pts): "
              f"alpha={coef[0]:.8f} beta={coef[1]:+.4f};  with a delta/m term: alpha={coef4[0]:.8f} "
              f"beta={coef4[1]:+.4f} delta={coef4[3]:+.1f}")
    print(f"  exceptional Haar mass at level m: N_m/2^m = 2^(log2N_m - m): "
          f"m=26 {2**(L[26]-26):.4g}, m=1000 {2**(L[1000]-1000):.3e}, m={M3} 2^{L[M3]-M3:.1f}")

    print()
    print("=" * 78)
    print("LADDER-C. q>=5: mu_q = Haar(Bad(q)) = lim N_m/2^m, rigorous enclosure")
    print("=" * 78)
    print("tail bound: 0 <= f_m - mu_q <= sum_{n>m} P(S_n<0) <= rho^(m+1)/(1-rho),")
    print("            rho = 2^(t-1)(1+q^(-t)) at the stated fixed t (Chernoff; any t>0 is valid)")
    for q in (3, 5, 7):
        print(f"drift E[S_1] = log2({q})/2 - 1 = {math.log2(q)/2 - 1:+.6f}"
              + ("  (< 0: mu_3 = 0 by the strong law)" if q == 3 else "  (> 0: mu_q > 0)"))
    results = {}
    for q in (5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31):
        t = best_t(q)
        rho = rho_of(q, t)
        target = mpmath.mpf("1e-22")
        m = int(mpmath.ceil(mpmath.log(target * (1 - rho)) / mpmath.log(rho)))
        N = exact_counts(q, m)
        f_m = mpmath.mpf(N[m]) / mpmath.mpf(2) ** m
        tail = rho ** (m + 1) / (1 - rho)
        f20 = N[20] / 2 ** 20
        f24 = N[24] / 2 ** 24
        results[q] = (f_m, tail)
        print(f"q={q:>2}: t={t} rho={mpmath.nstr(rho, 10)} m={m} tail<={mpmath.nstr(tail, 3)}  "
              f"mu_q in [{mpmath.nstr(f_m - tail, 24)}, {mpmath.nstr(f_m, 24)}]   [f_20={f20:.6f}, f_24={f24:.6f}]")
    print()
    print("=" * 78)
    print("LADDER-D. independent check: Spitzer / Sparre Andersen  mu_q = exp(-sum P(S_n<0)/n)")
    print("=" * 78)
    for q in (5, 7, 9, 11, 13):
        t = best_t(q)
        rho = rho_of(q, t)
        Nmax = int(mpmath.ceil(mpmath.log(mpmath.mpf("1e-24") * (1 - rho)) / mpmath.log(rho)))
        val, tot = spitzer_mu(q, Nmax)
        tail = rho ** (Nmax + 1) / ((Nmax + 1) * (1 - rho))
        lo_s = mpmath.e ** (-(tot + tail))
        f_m, tl = results[q]
        lo_d, hi_d = f_m - tl, f_m
        overlap = max(lo_s, lo_d) <= min(val, hi_d)
        print(f"q={q:>2}: Spitzer (n<={Nmax}) mu in [{mpmath.nstr(lo_s, 20)}, {mpmath.nstr(val, 20)}]; "
              f"DP enclosure [{mpmath.nstr(lo_d, 20)}, {mpmath.nstr(hi_d, 20)}]; "
              f"{'CONSISTENT' if overlap else 'INCONSISTENT'}; |diff|={mpmath.nstr(abs(val - f_m), 3)}")
    print(f"(ballot script time {time.time() - t0:.1f}s)")


if __name__ == "__main__":
    main()
