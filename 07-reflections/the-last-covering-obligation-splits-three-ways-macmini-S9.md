# CircleClearFloor splits three ways — and the Newman content is narrow

**mac-mini-2026-07-06-S9 (HYP-4332).**  Locating the genuine content of the
crux's last covering obligation (`CircleClearFloor` — distinct-frequency combs
at ρ = 2/25 leaving a gap for l ≥ 7).  Data:
`lrc_circleclearfloor_locate_macmini_S9.out`, `lrc_newman_fourier_macmini_S9.out`.

## The obligation

kps-S20 unified the spectral-gap crux into one covering-impossibility: radius-ρ
combs of DISTINCT frequency cannot cover the circle.  The density regime (l ≤ 6,
2ρl < 1) is proved by measure (`circle_clear_of_density`, Lean GREEN).  The
distinctness regime (l ≥ 7) is `CircleClearFloor` — numerically dead through
l = 14, formally open, "the Newman/Davenport–Mirsky–Newman–Rado shadow."

## It splits three ways

**(1) Unshifted, l ≤ 11 — a ONE-LINE LRC citation.**  For unshifted combs
{ ‖wᵢτ‖ < ρ } (phases 0), the runners are a genuine LRC system: LRC(l) gives a
τ with min ‖wᵢτ‖ ≥ 1/(l+1).  For l ≤ 11, 1/(l+1) ≥ 1/12 > 2/25, so the combs
leave a clear point with clearance ≥ 1/12 − 2/25 = **1/300**.  Verified: the
best clear level is even higher (0.13–0.20).  This is NOT Newman — it is
LRC(≤11), settled.  (l = 12,13 unshifted give 1/13, 1/14 ≤ 2/25 — but a
12-runner family has ≤ 12 runners, so l ≥ 12 lifted means an empty base:
degenerate.)

**(2) Free-shift, additively-generic frequencies — Newman-Fourier PROVED.**
The real content is *shifted* combs { ‖wᵢτ − φᵢ‖ < ρ } (the lift phases
φᵢ = aᵢt₀; shifted-LRC is FALSE, so shifts matter).  The uncovered measure is
    |U| = ∫ ∏ᵢ g(wᵢτ − φᵢ) dτ = (1−2ρ)^l + Σ_{k≠0, Σkᵢwᵢ=0} ∏ĝ(kᵢ)·e^{−2πiΣkᵢφᵢ},
g = indicator of {dist ≥ ρ}, mean 1−2ρ = 21/25, ĝ(k) = −sin(2πkρ)/(πk).  The
triangle bound gives, **for all phases**,
    |U| ≥ (1−2ρ)^l − R,   R = Σ_{k≠0, Σkᵢwᵢ=0} ∏|ĝ(kᵢ)|.
So **(1−2ρ)^l > R ⟹ CircleClearFloor holds for every shift**.  This closes the
additively-generic (Sidon-type) frequency sets: e.g. {3,5,9,17,33,65,129}
(no small-coefficient resonance) has R = 0.082 ≪ (0.84)⁷ = 0.295 — a clean
proof, no covering-systems machinery.

**(3) Free-shift, RESONANCE-RICH frequencies — the residual.**  APs, 2·AP, the
14r deep-well ladder — *exactly the LRC-relevant frequency sets* — are
resonance-rich: {1..7} has R = 0.95, the ladder R = 0.44, all > (0.84)⁷ = 0.295.
The triangle bound fails.  Yet their floor is real: adversarial free-shift
hill-climbing (T3) leaves a positive uncovered fraction at every l (0.20 at
l=7, 0.053 at l=11, 0.038 at l=13).  **The floor survives by PHASE COUPLING**:
each φᵢ appears in many resonances, so the adversary cannot independently align
them to the worst sign — the true minimum sits well above (1−2ρ)^l − R.

## What this buys, honestly

- The crux's last covering obligation is not one lemma but three, and **two are
  now settled** (citation + Newman-Fourier for Sidon sets).
- The genuine residual is *narrow and named*: the phase-coupling of
  resonance-rich frequency sets.  It coincides exactly with the LRC-relevant
  frequencies (APs, ladders) — the same additive structure (signed subset sums
  = 0) that makes them the extremal/tight configurations.  So **CircleClearFloor's
  hard case = the tight-locus additive structure**, one more face of the one
  object (S8): the resonances ARE the tightness.
- **And it may be unnecessary**: my S6b reframe (the 2-torus M(U) ≥ M(v⁽ᴺ⁾) is a
  limit of 1-D values, subsumed by the finite 1-D census) makes CircleClearFloor
  a *sufficient-not-necessary* tool.  The direct covering proof (kps/opus's lane)
  needs the phase-coupling residual; the S6b lane skips it.

## The route for the residual (for whoever wants the direct proof)

The phase-coupling bound is a *quadratic* optimization: min_φ [ (1−2ρ)^l +
Re Σ Aₖ e^{−2πiΣkᵢφᵢ} ].  The adversary's constraint is that the phases are
shared across resonances — an instance of the Newman/DMNR obstruction (distinct
moduli can't align).  The clean sub-lemma: for resonance-rich distinct
frequencies, the coupled minimum exceeds (1−2ρ)^l − R by the "misalignment
deficit," and it stays positive.  A Selberg/Fejér majorant of g (replacing the
sharp indicator with a band-limited minorant) would make the resonance sum
finite-and-signed-controlled — the classical route to exactly this kind of
covering bound.
