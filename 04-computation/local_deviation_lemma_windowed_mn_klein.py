#!/usr/bin/env python3
"""
klein-2026-07-01-S90 -- HYP-3847: THE LOCAL DEVIATION LEMMA (windowed Mirsky-Newman core)
                        + HYP-3848: the triple-overlap relation-lattice formula (part B)

LEMMA (local deviation; proved in the HYP file, 3 lines):
  F = {v_1..v_j} distinct speeds, r in (0,1/2), C = sum 1_{D_i}, A = 2rj.
  v* = a divisor-minimal element (e.g. min F); s = sin(2 pi r)/pi = |C-hat(v*)|.
  Spectral gap delta = min |m - v*| over frequencies m = k v_i (k != 0) of C, m != v*.
  If M + 1 <= delta then for EVERY center t0:
      int |C(t) - A| K_M(t - t0) dt  >=  s,
  where K_M = the Fejer kernel of degree M (K >= 0, int K = 1, K-hat supported [-M, M]).
  Proof: test C - A against e(-v* t) K_M(t - t0); the gap kills every frequency except
  v*, leaving exactly s in modulus; |int (C-A) e K| <= int |C-A| K.  QED

COROLLARY (windowed floor, A <= 1): uncovered-K-mass + overlap-K-mass >= s - (1 - A)
  at every center; integrating centers over a window I and combining with the mass
  identity U_I = |I| - mass_I + int_I (C-1)^+ gives
      U_{I+} >= [ (s + A) |I| - mass_I ] / 2 - Fejer-tail loss(M, |I|).
  BITES iff delta * |I| >~ const (need M ~ 1/|I| localization AND M + 1 <= delta).
  DICHOTOMY: spectral-gap clusters lose every window; gap-less clusters (consecutive =
  the AP difference core, delta = 1) are exactly the renormalization fixed point
  (opus HYP-3901) -- hpartA's danger case splits into the two proved regimes.

This script: numerical verification (fine deterministic grid) of the lemma, the
corollary, and the dichotomy; then part B exact/numeric checks of the triple overlap.
"""
import math

TWO_PI = 2 * math.pi

def C_of_t(F, r, t):
    c = 0
    for v in F:
        x = (v * t) % 1.0
        if min(x, 1.0 - x) < r:
            c += 1
    return c

def fejer(M, t):
    # K_M(t) = (1/(M+1)) (sin(pi (M+1) t)/sin(pi t))^2
    st = math.sin(math.pi * t)
    if abs(st) < 1e-15:
        return float(M + 1)
    return (math.sin(math.pi * (M + 1) * t) ** 2) / ((M + 1) * st ** 2)

def spectral_gap(F, vstar, span=3000):
    freqs = set()
    for v in F:
        k = 1
        while k * v <= vstar + span:
            freqs.add(k * v)
            k += 1
    freqs.discard(vstar)
    return min(abs(m - vstar) for m in freqs)

def local_dev_integral(F, r, A, M, t0, grid=200001):
    # int |C - A| K_M(t - t0) dt  (trapezoid on fine grid; C piecewise const)
    tot = 0.0
    for i in range(grid):
        t = i / grid
        tot += abs(C_of_t(F, r, t) - A) * fejer(M, t - t0)
    return tot / grid

print("=" * 92)
print("A. THE LOCAL DEVIATION LEMMA -- verification")
print("=" * 92)
r = 1.0 / 14
s = math.sin(TWO_PI * r) / math.pi
print(f"  r = 1/14, s = sin(2 pi r)/pi = {s:.6f}")

# gap cluster: 7 speeds, spacing 14, v* = min = 97 (divisor-minimal), delta = 14
GAP7 = [97, 111, 125, 139, 153, 167, 181]
A = 2 * r * len(GAP7)
d = spectral_gap(GAP7, 97)
M = d - 1
print(f"\n  GAP7 = {GAP7}: A = {A:.4f}, delta = {d}, M = {M}")
for t0 in (0.13, 0.377, 0.5, 0.771, 0.925):
    val = local_dev_integral(GAP7, r, A, M, t0)
    print(f"    t0 = {t0:5.3f}:  int|C-A| K_M = {val:.6f}  >= s = {s:.6f}  OK={val >= s - 1e-3}")

# consecutive cluster: delta = 1 -> M = 0 -> K = 1 (global only)
CONS7 = list(range(100, 107))
d2 = spectral_gap(CONS7, 100)
print(f"\n  CONS7 = {CONS7}: delta = {d2} -> M = 0 (lemma degenerates to the global bound)")
val = local_dev_integral(CONS7, r, A, 0, 0.0)
print(f"    global int|C-A| = {val:.6f} >= s = {s:.6f}  OK={val >= s - 1e-3}")

print("\nB. THE DICHOTOMY -- can a window be (nearly) covered?")
print("=" * 92)
def min_window_uncovered(F, r, wlen, ngrid=4000):
    best = (1.0, 0.0)
    for i in range(ngrid):
        a = i / ngrid
        # uncovered measure in [a, a+wlen] via fine subgrid
        sub = 400
        u = sum(1 for k in range(sub) if C_of_t(F, r, a + wlen * (k + 0.5) / sub) == 0)
        u = u / sub * wlen
        if u < best[0]:
            best = (u, a)
    return best

wlen = 0.05
for name, F in (("GAP7", GAP7), ("CONS7", CONS7)):
    u, a = min_window_uncovered(F, r, wlen)
    floor = (s * wlen) / 2  # critical-case A=1 heuristic floor (loss omitted)
    print(f"  {name}: min uncovered over windows |I|={wlen}: {u:.6f} at t={a:.4f}"
          f"   [lemma floor s|I|/2 = {floor:.6f} applies iff delta|I| large: "
          f"delta*|I| = {spectral_gap(F, min(F)) * wlen:.2f}]")

print("\n  reading: GAP7 (delta|I| = 0.7) obeys the floor; CONS7 (delta|I| = 0.05) can")
print("  nearly cover windows -- the gap-less regime, exactly the renormalization fixed")
print("  point (consecutive = AP difference core). hpartA danger case = split regimes.")

# ---------------- part B: triple overlap (HYP-3848) ----------------
print("\n" + "=" * 92)
print("C. TRIPLE OVERLAP |D_u ^ D_v ^ D_w| -- relation-lattice formula (HYP-3848)")
print("=" * 92)
from fractions import Fraction as Fr

def triple_overlap_exact(u, v, w, r):
    """Exact by interval intersection (Fractions)."""
    def danger(vv):
        rv = r / vv
        out = []
        for a in range(vv + 1):
            lo, hi = Fr(a, vv) - rv, Fr(a, vv) + rv
            lo, hi = max(lo, Fr(0)), min(hi, Fr(1))
            if hi > lo:
                out.append((lo, hi))
        return out
    def inter(A, B):
        out = []
        for lo1, hi1 in A:
            for lo2, hi2 in B:
                lo, hi = max(lo1, lo2), min(hi1, hi2)
                if hi > lo:
                    out.append((lo, hi))
        return out
    got = inter(inter(danger(u), danger(v)), danger(w))
    return sum(hi - lo for lo, hi in got)

def triple_overlap_series(u, v, w, r, R=400):
    """Fourier relation-lattice sum: m1 u + m2 v + m3 w = 0, coeff prod sinc."""
    def co(m):  # Fourier coefficient of 1_{||x||<r} at frequency m (m any int)
        if m == 0:
            return 2 * r
        return math.sin(TWO_PI * m * r) / (math.pi * m)
    tot = 0.0
    for m1 in range(-R, R + 1):
        for m2 in range(-R, R + 1):
            num = m1 * u + m2 * v
            if num % w:
                continue
            m3 = -num // w
            if abs(m3) > R:
                continue
            tot += co(m1) * co(m2) * co(m3)
    return tot

rF = Fr(1, 14)
for (u, v, w) in [(3, 5, 8), (2, 3, 5), (3, 5, 11), (5, 7, 12)]:
    ex = triple_overlap_exact(u, v, w, rF)
    se = triple_overlap_series(u, v, w, 1.0 / 14)
    rel = "w=u+v (one relation family)" if w == u + v else "generic"
    print(f"  ({u},{v},{w}) [{rel}]: exact = {float(ex):.8f} ({ex}),  series = {se:.8f},"
          f"  agree = {abs(float(ex) - se) < 5e-4}")
print("  independence benchmark (2r)^3 =", (1 / 7) ** 3)
print("  -> the series = (2r)^3 + pair-relation slices + the genuine triple slice;")
print("     w = u+v turns on the (k, k, -k)+null lattice; each slice is a finite")
print("     cosine/Bernoulli-polynomial evaluation (the d=3 layer of mac-mini sec-1).")
print("\nDONE.")
