# LEM-017 — The interior-augmentation cascade: μ is pointwise-monotone under adding runners, so the k=11 and k=12 hfloor bands reduce to the k=13 band per diameter

**Status:** PROVED (one paragraph, below). Companion to HYP-5775's exhaustive band sweep.
**Source:** mac-mini-2026-07-09-S65 (cont. 6).
**Depends on:** THM-657 (W/μ dictionary). **Serves:** THM-661/LEM-005's diam ∈ [19,35] band —
the last unproved step under hfloor (opus-S186).

**Setup.** For a finite shape `E ⊂ Z≥0` with `0 ∈ E`, let
`μ(E) = meas{x ∈ [0,1) : maxCircGap{frac(e·x) : e ∈ E} > 1/7}` (THM-661's density-floor object).

> **Lemma.** (i) If `E ⊆ E'` then `μ(E') ≤ μ(E)`.
> (ii) Consequently, for every diameter `d ≥ 13`,
> `min_{|E|=11, diam E = d} μ(E) ≥ min_{|E|=12, diam E = d} μ(E) ≥ min_{|E|=13, diam E = d} μ(E)`,
> where the minima are over shapes with `0, d ∈ E` (augment by interior points, which exist since
> `|[1, d−1] \ E| ≥ d − 1 − (|E| − 2) > 0`).

*Proof.* (i) Adding a phase point to a circular configuration either lands on an existing point
or splits one gap into two shorter ones; every gap length is nonincreasing, hence the maximal
gap is nonincreasing POINTWISE in `x`, hence the region `{maxgap > 1/7}` shrinks. (ii) Iterate
(i) with any one or two interior augmentation points, which change neither `min = 0` nor
`max = d`. ∎

## Consequence — the hfloor band architecture

The k=11 leg needs `μ ≥ bar₁₁ = 0.3312` on diam ∈ [19,35]; the k=12 leg needs `≥ bar₁₂ =
0.1993` there. By the cascade both follow from the SINGLE statement
> `min_{|E|=13, diam E = d} μ(E) ≥ 0.3312` for each `d ∈ [25, 35]`
(k=11 checked directly for d ≤ 24 by HYP-5775's exhaustive sweep — min μ = 0.71+, margin > 2×
the bar), leaving ONE finite computation (the k=13 band, C-scale) in place of three. If the
k=13 band minimum dips below 0.3312 at some d, the direct k=11/12 sweeps at that d remain as
the fallback — the reduction is conservative, never lossy.

**Also proved by (i):** the reflection symmetry `μ(E) = μ(d − E)` (the phase configuration of
`d − E` at `x` is a rotation-plus-reflection of that of `E`, preserving all gaps) — the factor-2
reduction used in the sweeps.

**Related:** THM-661, LEM-005, THM-660, HYP-5775, opus-S186 (hfloor terminal surface),
kps-S88 (k=11 compact), mac-mini-S58 (k=12,13 legs).
