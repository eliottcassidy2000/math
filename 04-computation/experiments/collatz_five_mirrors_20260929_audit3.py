#!/usr/bin/env python3
"""Independent audit (3) of section 2c of `collatz_five_mirrors_20260929.md` (S21): the closed recursion on the
powers of two to level 300, the rate/prefactor analysis, the split of the resonant coefficient by the excursion of
the prefix sums above the critical line, and the random units at h = 30.  Own code throughout; nothing imported
from the session's scripts (their outputs are parsed only for comparison).

Sections (tags refer to the output):
  [A] the recursion to h = 300 with AMAX = 40 (the session's setting) and with AMAX = 60 (certified: the a priori
      truncation bound 300 2^-40 = 2.7e-10 exceeds M(300) = 7.4e-11, so AMAX = 40 is not self-certifying); the
      rigorous rounding bound 50 eps sum_k maxwin|m_k| (all summed terms lie inside the computed window); a
      reversed-summation run and a noise-injection run for the empirical rounding sensitivity; comparison with the
      session's level-300 output at every level; doubling exponents, geometric-mean ratios, s - h log2 3.
  [B] P_h to 300 (own DP); free and pinned fits C h^-beta r^h over 100..300 for M and P_h; the (beta, r) trade-off
      curve (r fitted for fixed beta, with residuals) that the note calls the prefactor degeneracy; M/P_h.
  [C] own DP for the excursion split: validation at c = infinity against the closed recursion (h = 20, 40), mass
      = P_h at c = 0, completeness of the A-window; the numbers at h = 40, 80, 120 for c = 0, 4, 6, 16 (24 at 120).
  [D] random units at h = 30 (100 own random exponents): median, rms, against 0.84 3^-15; the exact distribution of
      |mu_hat_h| over all units at h = 12 and 14 (own FFT law): median/rms, to judge the twelve-unit sample.
Run: python 04-computation/experiments/collatz_five_mirrors_20260929_audit3.py   (about 5 minutes, < 3 GB)
"""
from __future__ import annotations

import cmath
import math
import random
import re
import time

import numpy as np

T0 = time.time()
L23 = math.log2(3)
I23 = L23 * math.log(2) + (L23 - 1) * math.log(L23 - 1) - L23 * math.log(L23)
R0 = math.exp(-I23)  # = 3^(h* - 1) = 0.94650


def stamp() -> str:
    return f"[{time.time() - T0:.0f}s]"


# ----------------------------------------------------------------------------------------------------------
# the closed recursion on the powers of two (numpy over an exponent window)
# ----------------------------------------------------------------------------------------------------------
def powers_recursion(hmax: int, amax: int, jmax: int, reverse: bool = False, noise: float = 0.0, seed: int = 0):
    """m_n(j) on j in [-amax (hmax - n), jmax], n = 1..hmax.  Returns per level (lo, max|m|, argmax j, sum over the
    window of nothing else); keeps only the last two levels in memory.  reverse: sum a downward (different rounding);
    noise: multiply every computed value by (1 + noise * gaussian) at every level (rounding-sensitivity probe)."""
    rng = np.random.default_rng(seed)
    w = np.array([2.0 ** -a for a in range(1, amax + 1)])
    jlo0 = -amax * hmax
    prev = np.ones(jmax - jlo0 + 1, dtype=np.complex128)
    prev_lo = jlo0
    res = {}
    order = range(amax, 0, -1) if reverse else range(1, amax + 1)
    for n in range(1, hmax + 1):
        mod = 3 ** n
        lo = -amax * (hmax - n)
        kmin = lo - amax
        cnt = jmax - kmin + 1
        pw = [0] * cnt
        p = pow(2, kmin, mod)
        for i in range(cnt):
            pw[i] = p
            p = (p * 2) % mod
        phase = np.exp(2j * np.pi * np.fromiter((r / mod for r in pw), dtype=np.float64, count=cnt))
        cur = np.zeros(jmax - lo + 1, dtype=np.complex128)
        for a in order:
            cur += w[a - 1] * phase[lo - a - kmin: jmax - a - kmin + 1] * prev[lo - a - prev_lo: jmax - a - prev_lo + 1]
        if noise:
            cur *= 1 + noise * (rng.standard_normal(cur.shape) + 1j * rng.standard_normal(cur.shape))
        a_ = np.abs(cur)
        k = int(np.argmax(a_))
        res[n] = (float(a_[k]), k + lo)
        prev, prev_lo = cur, lo
    return res


def cone_value(h: int, j: int, amax: int = 40) -> complex:
    """m_h(j) for a single exponent through its dependency cone (numpy, float64)."""
    lo0, hi0 = j - h * amax, j - h
    cur = np.ones(hi0 - lo0 + 1, dtype=np.complex128)
    cur_lo = lo0
    w = np.array([2.0 ** -a for a in range(1, amax + 1)])
    for n in range(1, h + 1):
        mod = 3 ** n
        lo, hi = j - (h - n) * amax, j - (h - n)
        kmin = lo - amax
        cnt = hi - kmin + 1
        pw = [0] * cnt
        p = pow(2, kmin, mod)
        for i in range(cnt):
            pw[i] = p
            p = (p * 2) % mod
        phase = np.exp(2j * np.pi * np.fromiter((r / mod for r in pw), dtype=np.float64, count=cnt))
        new = np.zeros(hi - lo + 1, dtype=np.complex128)
        for a in range(1, amax + 1):
            new += w[a - 1] * phase[lo - a - kmin: hi - a - kmin + 1] * cur[lo - a - cur_lo: hi - a - cur_lo + 1]
        cur, cur_lo = new, lo
    return complex(cur[j - cur_lo])


HMAX = 300
JMAX = 2 * HMAX + 40
r40 = powers_recursion(HMAX, 40, JMAX)
print(f"[A] recursion to h = {HMAX}, AMAX = 40 {stamp()}")
r60 = powers_recursion(HMAX, 60, JMAX)
print(f"[A] recursion to h = {HMAX}, AMAX = 60 {stamp()}")
rrev = powers_recursion(HMAX, 40, JMAX, reverse=True)
rnoise = powers_recursion(HMAX, 40, JMAX, noise=1e-12, seed=1)
print(f"[A] reversed-order and noise-injected (1e-12 per level) runs {stamp()}")
M = {n: r60[n][0] for n in r60}
S = {n: r60[n][1] for n in r60}
# the session's level-300 output
sess = {}
try:
    for line in open("05-knowledge/results/collatz_five_mirrors_powers_of_two300_20260929.out", encoding="utf-8"):
        mm = re.match(r"\s+h=\s*(\d+): ([0-9.e+-]+)\s+argmax j=(-?\d+)", line)
        if mm:
            sess[int(mm.group(1))] = (float(mm.group(2)), int(mm.group(3)))
except OSError:
    pass
print("[A1] h, M(h) [AMAX 60], AMAX 40, session's value, s - h, s - h log2 3, level ratio:")
for n in (20, 50, 100, 150, 180, 200, 240, 250, 280, 300):
    print(f"      h={n:3d}: {M[n]:.6e}  {r40[n][0]:.6e}  {sess.get(n, (float('nan'),))[0]:.6e}  s-h={S[n] - n:+d}  s-h log2 3={S[n] - n * L23:+.2f}  ratio={M[n] / M[n - 1]:.4f}")
print(f"[A2] AMAX 40 vs 60: max relative |difference| over h <= 300: {max(abs(r40[n][0] / M[n] - 1) for n in M):.1e}; levels where the argmax differs: {[n for n in M if r40[n][1] != S[n]]} (full-period ties for h <= 8: the window holds several periods 2 3^(h-1) and the first maximum is reported)")
print(f"[A3] against the session's output (AMAX 40, dict/cmath implementation): max relative |difference| over the {len(sess)} levels: {max(abs(sess[n][0] / r40[n][0] - 1) for n in sess):.1e}; levels where the argmax differs: {[n for n in sess if sess[n][1] != r40[n][1]]}")
print(f"[A4] reversed summation order: max relative |difference| {max(abs(rrev[n][0] / r40[n][0] - 1) for n in M):.1e}; noise 1e-12 injected at every level: relative |difference| at h = 100, 200, 300: "
      + ", ".join(f"{abs(rnoise[n][0] / r40[n][0] - 1):.1e}" for n in (100, 200, 300)) + "  (linear accumulation would give 3e-10 at 300; errors are averaged down, not accumulated)")
eps = np.finfo(float).eps
summax = sum(M[n] for n in range(1, HMAX))
print("[A5] error bounds. Truncation: the dropped terms sum_(a>AMAX) 2^-a m_(n-1)(j-a) involve exponents OUTSIDE the computed window, so only |m| <= 1 is")
print(f"     available a priori: <= n 2^-AMAX = {HMAX * 2.0 ** -40:.1e} (AMAX 40, > M(300) = {M[300]:.1e}: not self-certifying) or {HMAX * 2.0 ** -60:.1e} (AMAX 60).")
print(f"     Rounding: every summed term lies inside the window, so per level <= 50 eps maxwin|m_(n-1)|; total <= 50 eps sum_k maxwin|m_k| = 50 eps * {summax:.3f} = {50 * eps * summax:.1e}")
print(f"     Certified (AMAX 60): |error at h = 300| <= {HMAX * 2.0 ** -60 + 50 * eps * summax:.1e}, i.e. {(HMAX * 2.0 ** -60 + 50 * eps * summax) / M[300]:.1e} relative. The AMAX-40 values agree with it to {max(abs(r40[n][0] / M[n] - 1) for n in M):.0e}.")
print("[A6] local doubling exponents log2(M(h)/M(2h)): " + ", ".join(f"{h}->{2*h}: {math.log2(M[h] / M[2*h]):.3f}" for h in (25, 50, 75, 100, 125, 150)))
dh = np.array([25, 50, 75, 100, 125, 150])
de = np.array([math.log2(M[h] / M[2 * h]) for h in dh])
sl, ic = np.polyfit(dh, de, 1)
print(f"      linear fit: slope {sl:.4f} per level (endpoints: {(de[-1] - de[0]) / 125:.4f}), intercept {ic:.2f}; |log2 r0| = {-math.log2(R0):.4f}; the slope gives r = {2 ** -sl:.4f}")
print("[A7] geometric-mean level ratios of M: " + ", ".join(f"{a}-{b}: {(M[b] / M[a]) ** (1 / (b - a)):.5f}" for a, b in ((100, 200), (150, 300), (200, 300), (250, 300))))
print(f"[A8] s - h log2 3 over 20..300: min {min(S[n] - n * L23 for n in range(20, 301)):+.2f}, max {max(S[n] - n * L23 for n in range(20, 301)):+.2f}; at 300: {S[300] - 300 * L23:+.2f}")

# ----------------------------------------------------------------------------------------------------------
# [B] P_h and the fits
# ----------------------------------------------------------------------------------------------------------
def no_descent(hmax: int, amax: int = 80):
    dist = {0: 1.0}
    out = {}
    for j in range(1, hmax + 1):
        bound = j * L23
        nd = {}
        for s_, p in dist.items():
            for a in range(1, amax + 1):
                q = s_ + a
                if q < bound:
                    nd[q] = nd.get(q, 0.0) + p * 2.0 ** (-a)
                else:
                    break
        dist = nd
        out[j] = sum(dist.values())
    return out


P = no_descent(HMAX)
print(f"[B1] P_h at 100, 200, 300: {P[100]:.4e}, {P[200]:.4e}, {P[300]:.4e}; geometric-mean ratios: " + ", ".join(f"{a}-{b}: {(P[b] / P[a]) ** (1 / (b - a)):.5f}" for a, b in ((100, 200), (150, 300), (200, 300))))
print("[B2] M/P_h: " + ", ".join(f"{h}: {M[h] / P[h]:.3f}" for h in (20, 60, 100, 120, 140, 180, 200, 240, 280, 300)) + f"; range over 20..180: {min(M[h] / P[h] for h in range(20, 181)):.3f}-{max(M[h] / P[h] for h in range(20, 181)):.3f}; over 20..300: {min(M[h] / P[h] for h in range(20, 301)):.3f}-{max(M[h] / P[h] for h in range(20, 301)):.3f}")
hs = np.arange(100, HMAX + 1)


def fit_free(series):
    ys = np.log(np.array([series[h] for h in hs]))
    X = np.column_stack([np.ones(len(hs)), -np.log(hs), hs])
    c, *_ = np.linalg.lstsq(X, ys, rcond=None)
    return c[1], math.exp(c[2]), float(np.abs(X @ c - ys).max())


def fit_pinned(series, r):
    ys = np.log(np.array([series[h] for h in hs])) - hs * math.log(r)
    X = np.column_stack([np.ones(len(hs)), -np.log(hs)])
    c, *_ = np.linalg.lstsq(X, ys, rcond=None)
    return c[1], float(np.abs(X @ c - ys).max())


def fit_beta_fixed(series, beta):
    ys = np.log(np.array([series[h] for h in hs])) + beta * np.log(hs)
    b, a = np.polyfit(hs, ys, 1)
    return math.exp(b), float(np.abs(a + b * hs - ys).max())


for name, ser in (("M", M), ("P_h", P)):
    beta, r, res = fit_free(ser)
    bp, rp = fit_pinned(ser, R0)
    print(f"[B3] {name} ~ C h^-beta r^h over 100..300: free beta = {beta:.3f}, r = {r:.5f} (max log-residual {res:.4f}); r pinned to {R0:.5f}: beta = {bp:.3f} (residual {rp:.4f})")
    print(f"      trade-off (beta fixed -> fitted r, residual): " + ", ".join(f"beta={b:.2f}: r={fit_beta_fixed(ser, b)[0]:.4f} ({fit_beta_fixed(ser, b)[1]:.3f})" for b in (1.0, 1.25, 1.5, 1.75, 2.0, 2.5)))
b40, a40 = np.polyfit(hs, np.log(np.array([M[h] for h in hs])), 1)
print(f"[B4] pure geometric fit of M over 100..300: ratio {math.exp(b40):.5f}, max log-residual {float(np.abs(a40 + b40 * hs - np.log(np.array([M[h] for h in hs]))).max()):.4f}")
print(f"[B5] is the rate of P_h 'measured' better than that of M? free fit of P_h gives r = {fit_free(P)[1]:.5f} against the theorem {R0:.5f} (error {fit_free(P)[1] / R0 - 1:+.1e}); the same free fit applied to M gives {fit_free(M)[1]:.5f}")

# ----------------------------------------------------------------------------------------------------------
# [C] the excursion split (own DP)
# ----------------------------------------------------------------------------------------------------------
def excursion_sum(h: int, s: int, c, amax: int = 40, awin: int = 45, wide: bool = False):
    """Sum of 2^-A e(2^s Y_h(w)/3^h) over words with P_j < j log2 3 + c for all j (c = None: no constraint), by a DP over
    (A, j, P) with phases e((2^(s - A + P_(j-1)) mod 3^j)/3^j).  Returns (complex sum, mass, number of A values)."""
    if c is None:
        lim = None
        A_lo, A_hi = h, int(2 * h + 8 * math.sqrt(2 * h)) + 40
    else:
        lim = [0] + [math.floor(j * L23 + c) for j in range(1, h + 1)]  # largest integer < j log2 3 + c (irrational)
        if wide:
            A_lo, A_hi = h, lim[h]
        else:
            A_lo, A_hi = max(h, s - awin - c), min(lim[h], s + awin + c)
    tot = 0j
    mass = 0.0
    nA = 0
    for A in range(A_lo, A_hi + 1):
        if lim is None:
            lm = [min(A - (h - j), A) for j in range(h + 1)]
            lm[0] = 0
        else:
            lm = [min(lim[j], A - (h - j)) for j in range(h + 1)]  # room to reach A with a_i >= 1
        if any(lm[j] < j for j in range(h + 1)):
            continue
        nA += 1
        dp = np.zeros(lm[h] + 1, dtype=np.complex128)
        dq = np.zeros(lm[h] + 1)
        dp[0] = 1.0
        dq[0] = 1.0
        for j in range(1, h + 1):
            mod = 3 ** j
            pmax = lm[j - 1]
            # 2^(s - A + P) mod 3^j for P = 0..pmax by repeated doubling from the inverse power
            base = pow(2, s - A, mod)
            ph = np.empty(pmax + 1)
            for Pp in range(pmax + 1):
                ph[Pp] = base / mod
                base = (base * 2) % mod
            src = dp[: pmax + 1] * np.exp(2j * np.pi * ph)
            srq = dq[: pmax + 1]
            new = np.zeros(lm[j] + 1, dtype=np.complex128)
            nq = np.zeros(lm[j] + 1)
            for a in range(1, amax + 1):
                hi = min(pmax, lm[j] - a)
                if hi < 0:
                    break
                new[a: a + hi + 1] += (2.0 ** -a) * src[: hi + 1]
                nq[a: a + hi + 1] += (2.0 ** -a) * srq[: hi + 1]
            dp, dq = new, nq
        tot += dp[A]
        mass += dq[A]
    return tot, mass, nA


print("[C1] validation of the DP: c = infinity (all words, A up to 2h + 8 sqrt(2h) + 40) against the closed recursion:")
for h, s_ in ((20, 26), (40, 57)):
    full = cone_value(h, s_)
    tot, mass, nA = excursion_sum(h, s_, None)
    print(f"      h={h}, s={s_}: DP {abs(tot):.6e} (mass covered {mass:.6f}, {nA} values of A) vs recursion {abs(full):.6e}; |difference| = {abs(tot - full):.1e}")
tot0, mass0, _ = excursion_sum(40, 57, 0)
print(f"[C2] c = 0 at h = 40: DP mass {mass0:.6e} = P_40 = {P[40]:.6e} (difference {abs(mass0 - P[40]):.1e}); |coh_0| = {abs(tot0):.4e}")
tot6, mass6, nA6 = excursion_sum(40, 57, 6)
tot6w, mass6w, nA6w = excursion_sum(40, 57, 6, wide=True)
print(f"[C3] completeness of the A-window at h = 40, c = 6: window |A - s| <= 51 ({nA6} values) vs all A in [h, lim_h] ({nA6w} values): |sum| {abs(tot6):.6e} vs {abs(tot6w):.6e}, mass {mass6:.6e} vs {mass6w:.6e}")
print("[C4] the split at h = 40, 80, 120 (s = 57, 120, 184): c, family mass, mass/P_h, |coh_c|/|full|, |full - coh_c|/|full|, coherent fraction |coh_c|/mass")
note = {40: {0: (0.5440, 0.4810, 0.2520), 4: (0.9235, 0.1625, 0.0357), 6: (1.0082, 0.0121, 0.0198), 16: (0.9996, 0.0010, 0.0031)},
        80: {0: (0.4081, 0.6019, 0.1887), 4: (0.9284, 0.0961, 0.0311), 6: (1.0067, 0.1081, 0.0153), 16: (1.0128, 0.0129, 0.0011)},
        120: {0: (0.3457, 0.6802, 0.1530), 4: (0.8883, 0.1288, 0.0270), 6: (0.9969, 0.2735, 0.0131), 16: (1.0333, 0.0645, 0.0007), 24: (1.0002, 0.0024, 0.0001)}}
worst = 0.0
for h, s_ in ((40, 57), (80, 120), (120, 184)):
    full = cone_value(h, s_)
    for c in (0, 4, 6, 16, 24):
        if c == 24 and h != 120:
            continue
        tot, mass, _ = excursion_sum(h, s_, c)
        vals = (abs(tot) / abs(full), abs(full - tot) / abs(full), abs(tot) / mass)
        ref = note[h].get(c)
        if ref:
            worst = max(worst, max(abs(v - r_) for v, r_ in zip(vals, ref)))
        print(f"      h={h:3d} c={c:2d}: mass {mass:.4e} (= {mass / P[h]:6.2f} P_h)  |coh|/|full| = {vals[0]:.4f}  |full-coh|/|full| = {vals[1]:.4f}  coherent {vals[2]:.4f}   |full| = {abs(full):.4e} {stamp()}")
print(f"[C5] max |difference| against the note's quoted values (three ratios, all listed (h, c)): {worst:.1e}")
# is the c = 6 magnitude match stable?  the magnitude ratio and the vector remainder as functions of c at the three levels
print("[C6] magnitude ratio |coh_c|/|full| and vector remainder |full - coh_c|/|full| for c = 0..24:")
for h, s_ in ((40, 57), (80, 120), (120, 184)):
    full = cone_value(h, s_)
    rows = []
    for c in range(0, 25):
        t, m_, _ = excursion_sum(h, s_, c)
        rows.append((c, abs(t) / abs(full), abs(full - t) / abs(full)))
    stable10 = next((c for c in range(25) if all(r[2] < 0.10 for r in rows[c:])), None)
    stable05 = next((c for c in range(25) if all(r[2] < 0.05 for r in rows[c:])), None)
    print(f"      h={h:3d}: " + " ".join(f"c{c}:{mr:.2f}/{rm:.2f}" for c, mr, rm in rows) + f"  -> remainder stays below 10% from c = {stable10}, below 5% from c = {stable05} {stamp()}")

# ----------------------------------------------------------------------------------------------------------
# [D] random units at h = 30 and the exact coefficient distribution at h = 12, 14
# ----------------------------------------------------------------------------------------------------------
h = 30
Lh = 2 * 3 ** (h - 1)
rng = random.Random(2024)
vals = np.array([abs(cone_value(h, rng.randrange(Lh))) for _ in range(100)])
scale = math.sqrt(0.70) * 3 ** (-h / 2)
print(f"[D1] 100 random units at h = 30: median {np.median(vals):.2e}, rms {math.sqrt((vals ** 2).mean()):.2e}, mean {vals.mean():.2e}, min {vals.min():.2e}, max {vals.max():.2e}; square-root scale sqrt(0.70) 3^-15 = {scale:.2e}; "
      f"share below 0.5 x scale: {(vals < 0.5 * scale).mean():.2f}; mean |mu_hat|^2 3^h = {(vals ** 2).mean() * 3 ** h:.3f} (level-h average over all units: 0.70) {stamp()}")


def forward_law(N, AMAX=64):
    mu = np.array([1.0])
    for n in range(1, N + 1):
        mod = 3 ** n
        y = np.arange(3 ** (n - 1), dtype=np.int64)
        base = (3 * y + 1) % mod
        new = np.zeros(mod)
        for a in range(1, AMAX + 1):
            new[(pow(2, -a, mod) * base) % mod] += (2.0 ** -a) * mu
        mu = new
    return mu


for n in (12, 14):
    mu = forward_law(n)
    amp = np.abs(np.fft.fft(mu))
    units = np.arange(3 ** n) % 3 != 0
    au = amp[units]
    rms = math.sqrt((au ** 2).mean())
    print(f"[D2] exact distribution over all units at h = {n}: rms {rms:.3e} (= {rms * 3 ** (n / 2):.3f} 3^-h/2), median {np.median(au):.3e} = {np.median(au) / rms:.3f} rms, mean = {au.mean() / rms:.3f} rms, "
          f"share below 0.5 rms {(au < 0.5 * rms).mean():.3f}, share below 0.3 rms {(au < 0.3 * rms).mean():.3f}; Rayleigh would give median 0.83 rms, mean 0.89 rms, share below 0.5 rms 0.22")
print(f"DONE {stamp()}")
