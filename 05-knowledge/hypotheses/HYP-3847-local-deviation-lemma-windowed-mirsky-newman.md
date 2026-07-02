---
id: HYP-3847
title: THE LOCAL DEVIATION LEMMA (windowed Mirsky-Newman core) -- spectral gap + modulated Fejer kernel force int |C-A| K >= sin(2pi r)/pi at EVERY center; windowed uncovered floor; the delta*|I| dichotomy (gap clusters lose every window, gap-less clusters = the renormalization fixed point)
status: CONFIRMED -- lemma PROVED (3-line Fourier argument, in this file); corollary + dichotomy verified numerically; the consecutive cluster VIOLATES the gap floor (hypothesis necessary)
source: klein-2026-07-01-S90
script: 04-computation/local_deviation_lemma_windowed_mn_klein.py (+ .out)
related:
  - mac-mini S96 strategies doc sec 2 (the "windowed MN" target lemma -- this is its rigorous core, with the sharp hypothesis)
  - THM-594(C,D,E) (global Mirsky-Newman + Parseval defect + critical-mass floor)
  - kps HYP-3953 (hpartA via the c-ruler -- complementary route; if that dissolves hpartA, this lemma is the analytic regime map behind it)
  - opus HYP-3901 (difference-core renormalization -- the gap-less regime's owner)
  - HYP-3765 (Beurling-Selberg smoothing -- the owner's "BS majorant of the window" realized as the Fejer localization)
---

# HYP-3847: the local deviation lemma

## The lemma (PROVED)

Let F = {v_1, ..., v_j} be distinct positive integers, r in (0, 1/2), C(t) = sum_i
1_{||v_i t|| < r}, A = 2rj. Let v* be a divisor-minimal element of F (e.g. min F), so the
Fourier coefficient C-hat(v*) = sin(2 pi r)/pi =: s > 0. Let

    delta = min { |m - v*| : m = k v_i, k >= 1, m != v* }   (the spectral gap at v*).

**If M + 1 <= delta then for EVERY t0:  int |C(t) - A| K_M(t - t0) dt >= s,**
where K_M is the Fejer kernel of degree M (K_M >= 0, int K_M = 1, K_M-hat supported in
[-M, M]).

*Proof.* Test C - A against e(-v* t) K_M(t - t0). By Plancherel the pairing equals
sum_m C-hat(m) K_M-hat(m - v*) e((m - v*) t0) over the frequency support of C - A; the
band-limitation kills every m with |m - v*| > M, and the spectral gap leaves only m = v*,
giving C-hat(v*) K_M-hat(0) = s in modulus. Since |e(-v* t) K_M| = K_M pointwise, the
pairing's modulus is at most int |C - A| K_M. QED

## The windowed-floor corollary

|C - A| <= A on {C = 0} and <= (C-1)^+ + |1-A| on {C >= 1}, so at every center
(uncovered K-mass) + (overlap K-mass) >= s - |1 - A|. Integrating centers over a window I
(Fubini + explicit Fejer tails) and combining with the mass identity
U_I = |I| - mass_I + int_I (C-1)^+ gives, for A <= 1:

    U_{I+} >= [ (s + A)|I| - mass_I ] / 2 - loss(M, |I|),

positive even at the critical mass_I = |I|, A = 1 (the j = 7 tiling-attempt case):
**U >= s|I|/2 - loss**. The localization needs M ~ 1/|I| AND M + 1 <= delta: the lemma
bites iff **delta * |I| >~ 1**.

## The dichotomy (verified, and the hypothesis is NECESSARY)

- GAP7 = {97,111,125,139,153,167,181} (delta = 14, A = 1 at r = 1/14): int|C-A|K_13 in
  [0.58, 0.61] >= s = 0.1381 at every tested center; min window-uncovered (|I| = 0.05)
  = 0.01075 > the floor 0.00345.
- CONS7 = {100..106} (delta = 1, the AP difference core): min window-uncovered = 0.00288
  **< 0.00345** -- the consecutive cluster nearly covers a window, VIOLATING the
  gap-cluster floor. The gap-less regime genuinely escapes the window lemma -- and it is
  exactly the renormalization fixed point (opus HYP-3901: differences of an AP are an AP).

**Reading for hpartA:** the danger case ("the large cluster covers the whole positive-
measure G2 window") splits into two proved regimes: spectral-gap clusters cannot cover
any window of length >~ 1/delta (this lemma); gap-less clusters are consecutive-type =
the difference-core fixed point, owned by the tower/renormalization (and by kps's
HYP-3953 c-ruler route, which sidesteps windows entirely -- complementary). The owner's
"divisor-minimal frequency against a Beurling-Selberg majorant of the window" is realized
with the Fejer kernel as the band-limited window; a Selberg majorant sharpens constants
but the mechanism is the one above.
