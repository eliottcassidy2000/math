#!/usr/bin/env python3
"""Independent audit (2) of section 2b of `collatz_five_mirrors_20260929.md`: the Fourier coefficients of the
3-adic Syracuse law at the powers of two followed to level 120 by the closed frequency recursion, and the
no-descent probability P_h.  Own implementation; nothing imported from the session's scripts.

Recursion (from Y_n = 2^-a (3 Y_(n-1) + 1) mod 3^n):  mu_hat_n(t) = sum_(a>=1) 2^-a e((t 2^-a mod 3^n)/3^n)
mu_hat_(n-1)(t 2^-a mod 3^(n-1)).  On t = 2^j (j in Z) it closes: m_n(j) = sum_a 2^-a e((2^(j-a) mod 3^n)/3^n)
m_(n-1)(j-a), m_0 = 1.  Implemented as a numpy recursion over an exponent window, the phases from the exact
integers 2^k mod 3^n (Python ints, correctly rounded to float64), valuations truncated at a <= AMAX.

Sections: [A] closure and exactness (the recursion against a brute-force word sum at small n, exact rationals of
the phases, window bookkeeping); [B] the FFT values M(h), h <= 18, and the argmax offsets; [C] the values,
argmax offsets, ratios, local exponents and fits to h = 120; float64 error analysis: the recursion is an
l^inf contraction, so |error| <= n (2^-AMAX + 41 eps); AMAX = 40 against AMAX = 50; a 30-digit mpmath
recomputation of the argmax coefficient at h = 60 (and h = 120 if time allows) through its full dependency cone;
[D] the no-descent probability P_h by an own DP with the strict and the non-strict inequality, its ratios, the
identity e^(-I(log_2 3)) = 3^(h* - 1), the polynomial prefactor of P_h, and the fit C h^-beta r^h of the
coefficient (the ballot form) against the pure geometric and the shifted power law.
Run: python 04-computation/experiments/collatz_five_mirrors_20260929_audit2.py [HMAX]   (default 120; ~4 minutes)
"""
from __future__ import annotations

import cmath
import math
import sys
import time
from fractions import Fraction

import numpy as np

T0 = time.time()


def stamp() -> str:
    return f"[{time.time() - T0:.0f}s]"


HMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 120
L23 = math.log2(3)


# ----------------------------------------------------------------------------------------------------------
# the closed recursion on the powers of two: numpy over an exponent window
# ----------------------------------------------------------------------------------------------------------
def powers_recursion(hmax: int, amax: int, jmax: int, dtype=np.complex128):
    """m_n(j) for n = 1..hmax on the window j in [-amax (hmax - n), jmax]; returns {n: (jlo, array)}.
    Every value returned at level n depends only on values inside the level-(n-1) window (bookkeeping checked)."""
    w = np.array([2.0 ** -a for a in range(1, amax + 1)])
    jlo0 = -amax * hmax
    prev = np.ones(jmax - jlo0 + 1, dtype=dtype)  # level 0: mu_hat_0 = 1
    prev_lo = jlo0
    out = {}
    inv2 = None
    for n in range(1, hmax + 1):
        mod = 3 ** n
        lo = -amax * (hmax - n)
        # exact 2^k mod 3^n for k in [lo - amax, jmax] by repeated multiplication from 2^(lo - amax)
        kmin = lo - amax
        pw = np.empty(jmax - kmin + 1, dtype=object)
        p = pow(2, kmin, mod)
        for i in range(jmax - kmin + 1):
            pw[i] = p
            p = (p * 2) % mod
        phase = np.exp(2j * np.pi * np.array([r / mod for r in pw], dtype=np.float64))  # correctly rounded r/mod
        cur = np.zeros(jmax - lo + 1, dtype=dtype)
        for a in range(1, amax + 1):
            # cur[j] += w_a phase(j - a) prev[j - a] for j in [lo, jmax]; k = j - a in [lo - a, jmax - a]
            ks = slice(lo - a - kmin, jmax - a - kmin + 1)
            ps = slice(lo - a - prev_lo, jmax - a - prev_lo + 1)
            assert lo - a - prev_lo >= 0
            cur += w[a - 1] * phase[ks] * prev[ps]
        out[n] = (lo, cur)
        prev, prev_lo = cur, lo
    return out


def brute_m(n: int, j: int, amax: int) -> complex:
    """mu_hat_n(2^j mod 3^n) = sum over words of length n of 2^-A e(2^j Y_n(w)/3^n), valuations <= amax (tiny n)."""
    mod = 3 ** n
    t = pow(2, j, mod)
    total = 0j
    # iterate over words by odometer
    word = [1] * n
    while True:
        y = 0
        for a in word:
            y = (pow(2, -a, mod) * (3 * y + 1)) % mod
        total += 2.0 ** (-sum(word)) * cmath.exp(2j * math.pi * ((t * y) % mod) / mod)
        i = n - 1
        while i >= 0 and word[i] == amax:
            word[i] = 1
            i -= 1
        if i < 0:
            break
        word[i] += 1
    return total


# ----------------------------------------------------------------------------------------------------------
# [A] closure and exactness
# ----------------------------------------------------------------------------------------------------------
print("[A1] closure: t = 2^j gives t 2^-a = 2^(j-a); the family {2^j mod 3^n : j in Z} is the whole unit group (2 generates it, order 2 3^(n-1)),")
print("     so max over a full period of j IS M(n); over a window shorter than the period it is a lower bound for M(n).")
# the recursion against the brute-force word sum at n = 1..4 with a small amax (both truncated identically)
AM = 6
rec = powers_recursion(4, AM, 12)
worst = 0.0
for n in range(1, 5):
    lo, arr = rec[n]
    for j in (0, 1, 2, 5, 8, 12):  # inside the level-4 window [0, 12]; negative j are exact only at the lower levels
        worst = max(worst, abs(arr[j - lo] - brute_m(n, j, AM)))
lo3, arr3 = rec[3]
worst = max(worst, abs(arr3[-3 - lo3] - brute_m(3, -3, AM)))
print(f"[A2] recursion vs brute-force word sum over all {AM}^n words (n <= 4, valuations <= {AM}, six exponents each, plus j = -3 at n = 3): max deviation {worst:.1e}")
# level 1 exactly: m_1(j) = (1/3) e(1/3) + (2/3) e(2/3) for j even (|.| = 1/sqrt 3), conjugate for j odd
lo, arr = powers_recursion(1, 40, 6)[1]
v_even = arr[0 - lo]
print(f"[A3] m_1(0) = {v_even:.12f}; (1/3)e(1/3)+(2/3)e(2/3) = {(cmath.exp(2j*math.pi/3)/3 + 2*cmath.exp(4j*math.pi/3)/3):.12f}; |m_1| = {abs(v_even):.12f} = 1/sqrt3 = {1/math.sqrt(3):.12f}")
# the phase integers: 2^(j-a) mod 3^n for negative exponents is the inverse power; check pow(2, -k, mod) * 2^k = 1
ok = all((pow(2, -k, 3 ** n) * pow(2, k, 3 ** n)) % 3 ** n == 1 for n in (5, 30, 120) for k in (1, 7, 100, 4800))
print(f"[A4] inverse powers 2^-k mod 3^n well defined (2^-k 2^k = 1 mod 3^n at n = 5, 30, 120): {ok}")
print("[A5] bookkeeping: level n computes j in [-AMAX (HMAX - n), JMAX] from level n-1 values at j - a, a <= AMAX, i.e. from")
print("     [-AMAX (HMAX - n + 1), JMAX - 1], inside the level-(n-1) window (asserted in the code); truncation a <= AMAX drops mass 2^-AMAX per level.")

# ----------------------------------------------------------------------------------------------------------
# [B], [C] the main run
# ----------------------------------------------------------------------------------------------------------
AMAX = 40
JMAX = 2 * HMAX + 40
rec = powers_recursion(HMAX, AMAX, JMAX)
print(f"[B] main run: HMAX = {HMAX}, AMAX = {AMAX}, JMAX = {JMAX} {stamp()}")
known = {1: 0.577350, 2: 0.377924, 3: 0.252237, 4: 0.176999, 5: 0.129274, 6: 0.096106, 7: 0.075870, 8: 0.060891,
         9: 0.048026, 10: 0.038278, 11: 0.031944, 12: 0.026458, 13: 0.022052, 14: 0.019128, 15: 0.016284,
         16: 0.014409, 17: 0.012511, 18: 0.011187}
mx = {}
arg = {}
for n in range(1, HMAX + 1):
    lo, arr = rec[n]
    a = np.abs(arr)
    k = int(np.argmax(a))
    mx[n] = float(a[k])
    arg[n] = k + lo
period_ok = {n: (JMAX - rec[n][0] + 1) >= 2 * 3 ** (n - 1) for n in range(1, 12)}
print("[B1] window covers a full period of the unit group (max = M(h) by definition) for h <= " + str(max(n for n, v in period_ok.items() if v)))
print("[B2] h, max_j |m_h(j)|, FFT M(h), difference, argmax offset s - h (h = 1..18):")
for n in range(1, 19):
    print(f"      h={n:2d}: {mx[n]:.6e}  M={known[n]:.6f}  diff={mx[n] - known[n]:+.1e}  s-h={arg[n] - n:+d}")
print("[B3] max |max_j |m_h(j)| - M(h)| over h <= 18: " + f"{max(abs(mx[n] - known[n]) for n in range(1, 19)):.1e}; offsets 9..18: " + ", ".join(f"{arg[n] - n:+d}" for n in range(9, 19)))
print("[C1] h = 20, 30, ..., 120: value, argmax j, s - h, ratio to previous level:")
for n in range(20, HMAX + 1, 10):
    print(f"      h={n:3d}: {mx[n]:.6e}  j={arg[n]}  s-h={arg[n] - n:+d}  ratio={mx[n] / mx[n - 1]:.4f}")
note_vals = {20: 8.88e-3, 30: 3.29e-3, 40: 1.38e-3, 50: 6.16e-4, 60: 2.83e-4, 70: 1.33e-4, 80: 6.45e-5, 90: 3.19e-5, 100: 1.59e-5, 110: 8.06e-6, 120: 4.12e-6}
note_off = {20: 6, 30: 12, 40: 17, 50: 23, 60: 29, 70: 35, 80: 40, 90: 46, 100: 52, 110: 58, 120: 64}
if HMAX >= 120:
    print("[C2] against the note's table: max relative deviation of the values " + f"{max(abs(mx[n] / v - 1) for n, v in note_vals.items()):.1e}; offsets agree: {all(arg[n] - n == o for n, o in note_off.items())}")
    print("[C3] argmax exponent s against h log_2 3 (s - h log2 3): " + ", ".join(f"h={n}: {arg[n] - n * L23:+.1f}" for n in range(20, 121, 10)))
    print("[C3b] the resonant frequency as a fraction of the modulus, 2^s / 3^h: " + ", ".join(f"h={n}: {2.0 ** (arg[n] - n * L23):.4f}" for n in range(20, 121, 20)) + "  (h <= 18: 0.12 -> 0.022; it stops decreasing near 2^-6)")
    print("[C4] local exponents log2(m(h)/m(2h)): " + ", ".join(f"{h}->{2*h}: {math.log2(mx[h] / mx[2*h]):.3f}" for h in (10, 15, 20, 25, 30, 40, 50, 60)))
    dh = np.array([10, 15, 20, 25, 30, 40, 50, 60])
    de = np.array([math.log2(mx[h] / mx[2 * h]) for h in dh])
    sl, ic = np.polyfit(dh, de, 1)
    print(f"[C4b] for C h^-beta r^h the doubling exponent is h |log2 r| + beta: linear fit slope {sl:.4f} per level (|log2 0.94650| = {-math.log2(0.94650):.4f}; |log2 0.930| = {-math.log2(0.930):.4f}), intercept beta = {ic:.2f}; "
          f"a shifted power law C (h + c)^-alpha would give exponents bounded by alpha and saturating")
    hs = np.arange(40, HMAX + 1)
    ys = np.log(np.array([mx[h] for h in hs]))
    b, a = np.polyfit(hs, ys, 1)
    res_geo = np.abs(ys - (a + b * hs)).max()
    best = None
    for c in np.arange(0, 40.01, 0.1):
        bb, aa = np.polyfit(np.log(hs + c), ys, 1)
        r = np.abs(ys - (aa + bb * np.log(hs + c))).max()
        if best is None or r < best[0]:
            best = (r, c, -bb)
    print(f"[C5] fits over h = 40..{HMAX}: geometric ratio {math.exp(b):.5f} (max log-residual {res_geo:.4f}); best shifted power law C (h + {best[1]:.1f})^(-{best[2]:.3f}) (max log-residual {best[0]:.4f})")
    # the ballot form C h^-beta r^h (see [D]) and a pure geometric fit over 80..120
    X = np.column_stack([hs, np.log(hs), np.ones(len(hs))])
    coef, *_ = np.linalg.lstsq(X, ys, rcond=None)
    res_b = np.abs(ys - X @ coef).max()
    b2, a2 = np.polyfit(hs[hs >= 80], ys[hs >= 80], 1)
    print(f"[C6] C h^-beta r^h over 40..{HMAX}: r = {math.exp(coef[0]):.5f}, beta = {-coef[1]:.3f} (max log-residual {res_b:.4f}); pure geometric over 80..{HMAX}: ratio {math.exp(b2):.5f}; asymptotic no-descent rate 3^(h*-1) = 0.94650")
    print("[C7] max_j |m_h(j)| e^(I h) with I = I(log2 3) (the prefactor of the coefficient, if the rate is e^-I): " + ", ".join(f"h={n}: {mx[n] * math.exp((L23 * math.log(2) + (L23 - 1) * math.log(L23 - 1) - L23 * math.log(L23)) * n):.4f}" for n in range(20, 121, 20)))

# ----------------------------------------------------------------------------------------------------------
# float64 error analysis
# ----------------------------------------------------------------------------------------------------------
eps = np.finfo(float).eps
print("[C8] error analysis: the truncated recursion is a weighted average with weights 2^-a (sum 1 - 2^-AMAX) times unimodular phases,")
print("     hence an l^inf contraction; per level the computed value differs from the exact one by at most 2^-AMAX (dropped mass, |m| <= 1)")
print("     plus the rounding of one 40-term complex sum with correctly rounded phases (<= ~41 eps); errors add, never amplify:")
print(f"     |error at level n| <= n (2^-{AMAX} + 41 eps) = {HMAX} * {2.0 ** -AMAX + 41 * eps:.2e} = {HMAX * (2.0 ** -AMAX + 41 * eps):.2e} at n = {HMAX}, against the value {mx[HMAX]:.2e} (relative {HMAX * (2.0 ** -AMAX + 41 * eps) / mx[HMAX]:.1e}).")
rec50 = powers_recursion(HMAX, 50, JMAX)
d = {}
for n in (20, 60, 100, HMAX):
    lo40, a40 = rec[n]
    lo50, a50 = rec50[n]
    # compare on the common window [lo40, JMAX]
    d[n] = float(np.abs(a40[:] - a50[lo40 - lo50: lo40 - lo50 + len(a40)]).max())
print("[C9] AMAX = 40 against AMAX = 50 (float64), max |difference| on the window: " + ", ".join(f"h={n}: {v:.1e}" for n, v in d.items()) + f"  (predicted truncation gap <= n 2^-40 = {HMAX * 2.0 ** -40:.1e} at {HMAX}) {stamp()}")
try:
    import mpmath as mp

    mp.mp.dps = 30

    def cone_value(h: int, j: int, amax: int):
        """m_h(j) in 30-digit arithmetic through its dependency cone: level n needs j - (h-n) amax .. j - (h-n)."""
        cur = {k: mp.mpc(1) for k in range(j - h * amax, j - h + 1)}  # level 0 on the cone base
        for n in range(1, h + 1):
            mod = 3 ** n
            lo, hi = j - (h - n) * amax, j - (h - n)
            new = {}
            for jj in range(lo, hi + 1):
                s = mp.mpc(0)
                for a in range(1, amax + 1):
                    r = pow(2, jj - a, mod)
                    s += mp.mpf(2) ** (-a) * mp.expjpi(2 * mp.mpf(r) / mod) * cur[jj - a]
                new[jj] = s
            cur = new
        return cur[j]

    for h in (30, 60):
        t1 = time.time()
        v = cone_value(h, arg[h], AMAX)
        lo, arr = rec[h]
        print(f"[C10] mpmath 30 digits, h={h}, j={arg[h]}: |m| = {mp.nstr(abs(v), 15)}; float64: {mx[h]:.15e}; |difference| = {float(abs(abs(v) - mx[h])):.1e} [{time.time() - t1:.0f}s]")
except ImportError:
    print("[C10] mpmath not available: high-precision recomputation skipped")

# ----------------------------------------------------------------------------------------------------------
# [D] the no-descent probability
# ----------------------------------------------------------------------------------------------------------
def no_descent(hmax: int, strict: bool, amax: int = 64):
    """P_h = P(S_j < j log2 3 for all j <= h) (or <=), i.i.d. geometric(1/2) valuations, by DP over the prefix sum."""
    dist = {0: 1.0}
    out = []
    for j in range(1, hmax + 1):
        bound = j * L23
        nd = {}
        for s, p in dist.items():
            for a in range(1, amax + 1):
                q = s + a
                if (q < bound) if strict else (q <= bound):
                    nd[q] = nd.get(q, 0.0) + p * 2.0 ** (-a)
                else:
                    break  # q only grows with a
        dist = nd
        out.append(sum(dist.values()))
    return out


P = no_descent(HMAX, True)
P2 = no_descent(HMAX, False)
print(f"[D1] strict vs non-strict inequality in P_h: max |difference| over h <= {HMAX} = {max(abs(x - y) for x, y in zip(P, P2)):.1e} (j log2 3 is irrational for j >= 1, so S_j < j log2 3 iff S_j <= floor(j log2 3))")
I23 = L23 * math.log(2) + (L23 - 1) * math.log(L23 - 1) - L23 * math.log(L23)
p3 = math.log(2) / math.log(3)
hstar = -p3 * math.log2(p3) - (1 - p3) * math.log2(1 - p3)
print(f"[D2] I(log2 3) = {I23:.6f}, e^-I = {math.exp(-I23):.5f}; h* = binary entropy of log_3 2 in bits = {hstar:.5f}, 3^(h*-1) = {3 ** (hstar - 1):.5f}; identity: (h*-1) ln 3 = -ln p - ((1-p)/p) ln(1-p) - ln 3 = -I with p = log_3 2 = 1/log_2 3 (algebra in the report)")
print("[D3] P_h at h = 20, 30, ..., 120 with the ratio P_h/P_(h-1), the prefactor P_h e^(I h), and max_j |m_h(j)| / P_h:")
note_P = {20: 1.93e-2, 30: 7.07e-3, 40: 2.99e-3, 50: 1.33e-3, 60: 6.15e-4, 70: 2.83e-4, 80: 1.40e-4, 90: 6.94e-5, 100: 3.45e-5, 110: 1.79e-5, 120: 9.31e-6}
for n in range(20, HMAX + 1, 10):
    print(f"      h={n:3d}: P={P[n-1]:.4e}  ratio={P[n-1]/P[n-2]:.4f}  P e^(Ih)={P[n-1] * math.exp(I23 * n):.4f}  m/P={mx[n] / P[n-1]:.4f}" + (f"  note P={note_P[n]:.2e}" if n in note_P else ""))
if HMAX >= 120:
    ratios_mp = {n: mx[n] / P[n - 1] for n in range(20, 121)}
    viol = [(n, round(r, 3)) for n, r in ratios_mp.items() if abs(r - 0.46) > 0.02]
    print(f"[D4] against the note: max relative deviation of P_h {max(abs(P[n-1] / v - 1) for n, v in note_P.items()):.1e}; m/P over 20..120: min {min(ratios_mp.values()):.3f}, max {max(ratios_mp.values()):.3f}, at h = 120: {ratios_mp[120]:.3f}; "
          f"levels with |m/P - 0.46| > 0.02: {viol}; window means of m/P: 20-40 {np.mean([ratios_mp[n] for n in range(20, 41)]):.4f}, 41-80 {np.mean([ratios_mp[n] for n in range(41, 81)]):.4f}, 81-120 {np.mean([ratios_mp[n] for n in range(81, 121)]):.4f}")
    print(f"[D5] P_h^(1/h) at h = 30, 60, 120: {P[29] ** (1/30):.4f}, {P[59] ** (1/60):.4f}, {P[119] ** (1/120):.4f} (limit e^-I = {math.exp(-I23):.4f}); the ratio P_h/P_(h-1) averaged over 61..120: {(P[119] / P[59]) ** (1/60):.4f}; over 101..120: {(P[119] / P[99]) ** (1/20):.4f}")
    # the prefactor exponent: log(P e^{Ih}) against log h
    hh = np.arange(40, 121)
    pref = np.array([P[n - 1] * math.exp(I23 * n) for n in hh])
    slope = np.polyfit(np.log(hh), np.log(pref), 1)[0]
    print(f"[D6] prefactor P_h e^(I h) ~ h^(-beta): fitted beta over 40..120 = {-slope:.3f} (local: 30->60: {math.log2((P[29] * math.exp(30 * I23)) / (P[59] * math.exp(60 * I23))):.2f}, 60->120: {math.log2((P[59] * math.exp(60 * I23)) / (P[119] * math.exp(120 * I23))):.2f}); the same for the coefficient: 30->60: {math.log2((mx[30] * math.exp(30 * I23)) / (mx[60] * math.exp(60 * I23))):.2f}, 60->120: {math.log2((mx[60] * math.exp(60 * I23)) / (mx[120] * math.exp(120 * I23))):.2f}")
    print("[D7] Chernoff upper bound P_h <= e^(-I h) holds at every h: " + str(all(P[n - 1] <= math.exp(-I23 * n) * (1 + 1e-12) for n in range(1, HMAX + 1))) + f"; P_h e^(Ih) at h = 1, 2, 5, 10: " + ", ".join(f"{P[n-1] * math.exp(I23 * n):.3f}" for n in (1, 2, 5, 10)))
# ----------------------------------------------------------------------------------------------------------
# [E] the multiplier families u 2^j (u odd, 3 !| u, u <= 49): is the maximum on the pure powers of two?
# ----------------------------------------------------------------------------------------------------------
def family_recursion(u: int, hmax: int, amax: int, jmax: int):
    """m_n(u, j) = mu_hat_n(u 2^j mod 3^n): same recursion with phases from u 2^(j-a) mod 3^n."""
    w = np.array([2.0 ** -a for a in range(1, amax + 1)])
    jlo0 = -amax * hmax
    prev = np.ones(jmax - jlo0 + 1, dtype=np.complex128)
    prev_lo = jlo0
    best = {}
    for n in range(1, hmax + 1):
        mod = 3 ** n
        lo = -amax * (hmax - n)
        kmin = lo - amax
        pw = []
        p = (u * pow(2, kmin, mod)) % mod
        for _ in range(jmax - kmin + 1):
            pw.append(p)
            p = (p * 2) % mod
        phase = np.exp(2j * np.pi * np.array([r / mod for r in pw], dtype=np.float64))
        cur = np.zeros(jmax - lo + 1, dtype=np.complex128)
        for a in range(1, amax + 1):
            cur += w[a - 1] * phase[lo - a - kmin: jmax - a - kmin + 1] * prev[lo - a - prev_lo: jmax - a - prev_lo + 1]
        best[n] = float(np.abs(cur).max())
        prev, prev_lo = cur, lo
    return best


HE = min(HMAX, 60)
fam = {}
for u in range(1, 50, 2):
    if u % 3 == 0:
        continue
    fam[u] = family_recursion(u, HE, 36, 2 * HE + 40)
print(f"[E1] seventeen families u 2^j, u <= 49 odd, 3 !| u, to h = {HE} (AMAX = 36, error <= h 2^-36 = {HE * 2.0 ** -36:.0e}): ratio max_(u > 1) / (u = 1) at h = 10, 20, ..., {HE}: "
      + ", ".join(f"h={n}: {max(fam[u][n] for u in fam if u > 1) / fam[1][n]:.4f} (u={max((u for u in fam if u > 1), key=lambda u: fam[u][n])})" for n in range(10, HE + 1, 10)) + f" {stamp()}")
print(f"[E2] is the pure family the largest at EVERY level 1..{HE}? " + str(all(fam[1][n] >= max(fam[u][n] for u in fam if u > 1) for n in range(1, HE + 1)))
      + "; levels where another family ties within 1e-9: " + str([n for n in range(1, HE + 1) if any(abs(fam[u][n] - fam[1][n]) < 1e-9 for u in fam if u > 1)][:20]))
print("[E3] reading: 2 generates the units mod 3^n, so every u is 2^k for a profinite exponent k in Z/2 x Z_3 (its discrete log is not an integer");
print("     a family u 2^j is the pure family shifted by that 3-adic exponent, and for h beyond the window period (h >= 9 here) the seventeen")
print("     families are seventeen windows of ~ (36 (60-h) + 160) residues each out of 2 3^(h-1) units: a sample of the unit group, not a cover.")
print(f"DONE {stamp()}")
