---
id: THM-711
title: The hit-empty product form of the k=9 base — the deg-2 requirement is EXACTLY E[N(7−N)] = 6m1 − m2 = Σ_{s≠s′} P(s empty ∧ s′ hit) ≥ 432/91 (exact identity E[N(7−N)] = 7m1 − (m2+m1); threshold 12(1−cap₁₀)); adversarial hunt over 57 families in every enemy class (consec/shifted, doubling, near-AP, mod-7-aligned, random, far-mixed) gives global min 4465/882 = 5.0624 AT SHIFTED-CONSEC {1..9}, margin +0.3151 — evidence the base holds UNCONDITIONALLY over all 9-element integer sets, killing the core-bound bookkeeping
status: IDENTITY PROVED (two lines); threshold exact; the unconditional inf ≥ 432/91 with minimizer {1..9} is the SHARP CONJECTURE (0 violations, 57 adversarial families; dilate-invariant; far elements RAISE J via the THM-710 eigen-transfer, so wide directions are covered — the conjecture's content is compact). Partial floors: p₇ = 0 always (some sector is hit); J ≥ 6(1−p₀) pointwise (insufficient alone: needs the N ≥ 2 middle mass, which consec retains).
source: mac-mini-2026-07-09-S65 (cont.35, 2026-07-11)
depends_on:
  - THM-705 (the linear form), THM-710 (the eigen-transfer that propagates it)
related:
  - kps 3D box (the k=8/m3 analog), opus-S220 ladder. COLLISION FLAG: opus-S227 reused THM-709 (mine = the doubling singleton, cont.33, prior — please renumber theirs).
---
# THM-711 — the hit-empty product form
h(N) = N/2 − N(N−1)/12 = N(7−N)/12, so the THM-705 requirement m1/2 − m2/12 ≥ 1 − cap₁₀ is
E[N(7−N)] ≥ 432/91: the mass of (empty, hit) sector pairs. Decorrelated value 42[(6/7)⁹−(5/7)⁹]
= 8.18 (73% headroom); the correlations at consec eat it down to 5.06, still +0.315 clear.
The k=9 wide-direction base = ONE compact inequality: **inf over 9-sets of E[N(7−N)] = 4465/882,
attained at {1..9}** (conjectured exact; verified across every adversarial class).
Files: 04-computation/lrc14_hit_empty_product_macmini_S65cont35.py (+ .out).
