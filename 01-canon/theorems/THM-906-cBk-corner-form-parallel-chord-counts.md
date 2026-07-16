---
id: THM-906
title: THE c_B(k) CORNER FORM + THE PARALLEL-CHORD COUNTING LAW — (I) the residue-6 kernel's per-relation plateau in closed form; c_B(k) = Σ_{j≠0}∏Êᵢ(jkᵢ) = −(1/(24·∏kᵢ))·Σ_{sections,corners}(−1)^{#upper}B₄({Σkᵢ(cᵢ+εᵢ)/7}) — a Bernoulli-B₄ alternating corner sum (the quartic analog of THM-880's B₂ law); validated R/c_B = 0.998 at the plateau-dominant regime (wide boxes, exact-ℚ strata), constants of the form N/7⁴ (e.g. 45/2401). (II) Additive quadruples ARE the Guy/Zarankiewicz species: circular parallel-chord count Qcirc(Z_n) = n·⌊(n−2)²/4⌋/2 EXACTLY (= n(n−1)(n−3)/8 odd n; proof: pairs in a sum class are automatically disjoint ⟹ Qcirc = Σ_s C(r(s),2)); linear count = partial sums of quarter-squares (Zarankiewicz atoms A002620); the near-AP quadruple ladder 125 > 95 > 83 ≫ 7 (AP/deep-well/GW/random) quantifies the kernel's near-AP concentration (the E₄ face of THM-730's AP-extremality)
status: (I) DERIVED (three-line Fourier: interval transforms → 16-corner expansion → Σe(jx)/j⁴ = −(2π)⁴B₄({x})/24) + REFEREED at the plateau-dominant regime (0.998, 1.079 at wide boxes; singleton-box ladders honest-noisy: plateau 1/14406 ≈ 7e-5 sits under O(1/v) at v ≤ 215 — larger ladders or wider boxes separate; boxeph's measured ~5e-3 was a different box/family phase, consistent with the form's corner-phase dependence). (II) PROVED (class-count identity, verified n = 5..15 exact; Guy Z(12) = Qcirc(Z₁₂) = 150 is a coincidence of the two floor-products)
source: death-star-2026-07-16-S25 (owner: near-AP quadruples vs Guy/A000241 + the residue-6 kernel); completes boxeph-S33 THM-898's evidence-log item "evaluate c_B(k) exactly" and supplies codex-HYP-7062's "universal relation-stratified proofs" machinery
depends_on: [boxeph THM-898 (S33), codex THM-891/904/905, THM-880 (B₂ analog), THM-730]
verification: 04-computation/cBk_closed_form_and_guy_quadruples_deathstar_S25.py -> 05-knowledge/results/cBk_closed_form_and_guy_quadruples_deathstar_S25.out
---

# THM-906 — the corner form and the parallel-chord law

## (I) c_B(k) in closed form

For runners v with unique primitive relation k ⊥ v and section boxes Bᵢ ⊆ Z₇, the orbit
equidistributes on the subtorus {Σkᵢxᵢ ∈ Z}, whose Haar measure of the box product exceeds
the independent product by exactly the nonzero dual modes:

> c_B(k) = Σ_{j≠0} ∏ᵢ Êᵢ(jkᵢ)
>       = −(1/(24·∏ᵢkᵢ)) · Σ_{cᵢ∈Bᵢ} Σ_{ε∈{0,1}⁴} (−1)^{Σε} B₄({Σᵢ kᵢ(cᵢ+εᵢ)/7}),

B₄(x) = x⁴−2x³+x²−1/30. Three-line proof: 1̂_{[a,b)}(m) = (e(−ma)−e(−mb))/(2πim); expand the
product over the 16 corner choices; (2πij)⁴ = (2π)⁴j⁴ and Σ_{j≠0}e(jx)/j⁴ = −(2π)⁴B₄({x})/24.
The 1/∏kᵢ factor + corner-argument spread PROVE boxeph's "taller relations ⟹ smaller
plateau" observation; all constants are rationals with denominator dividing 30·7⁴-type
(e.g. 45/2401 for k = (1,1,−1,−1) at [0..2]⁴-sectors, refereed at R/c_B = 0.998).
This is the named tool for codex-HYP-7062's sector-box caps: their candidate constants are
now derivable corner sums; the per-relation uniform law |R| ≤ Σ c_B(k) + O(1/v) has its
constants.

## (II) The parallel-chord counting law (the Guy/A000241 face)

Kernel quadruples v₁+v₂ = v₃+v₄ are PARALLEL CHORDS (equal midpoint) on the residue
circle — the pairing-complement of Guy's crossing (interleaved) quadruples. Exactly:

- **Circular:** on Z_n, two distinct pairs in one sum class never share an element, so
  Qcirc(Z_n) = Σ_s C(r(s),2) = **n·⌊(n−2)²/4⌋/2** (= n(n−1)(n−3)/8 for odd n) — a
  floor-product of exactly Guy's species Z(n) = (1/4)⌊n/2⌋⌊(n−1)/2⌋⌊(n−2)/2⌋⌊(n−3)/2⌋;
  the two agree at n = 12 (both 150), a floor-product coincidence.
- **Linear:** Qlin({1..n}) = Σ_k ⌊k²/4⌋-partial-sums — the QUARTER-SQUARE atoms of
  Zarankiewicz's cr(K_{m,n}) conjecture.
- **Near-AP ladder:** quadruple counts 125 (AP {1..13}) > 95 (deep well) > 83 (GW-ish)
  ≫ 7 (random 13-set): the kernel's near-AP concentration is the E₄ (additive-energy)
  face of THM-730's AP-extremality, with counting functions literally of the
  crossing-number species — the T-crossing-numbers tangent's "same species" intuition is
  now a pair of exact formulas.
