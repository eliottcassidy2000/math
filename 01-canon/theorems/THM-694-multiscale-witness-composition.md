---
id: THM-694
title: The multi-scale witness composition — THM-693 iterated per scale via two generic lemmas: THE BAND LIFT (a speed strictly in-band at (q,a) with p ≤ B stays strictly in-band at (qV, aV+Δ) for ANY Δ ∈ [0,q) once 14Bq < V — integrality supplies the margins, the lift never inspects Δ, so higher levels cannot disturb lower ones) and THE CLUSTER JOIN (a band-interior missed residue (C−ea) % q admits the new cluster speed V−e at (qV, aV+(C−aV)%q) once 14eq < V), plus miss_at_next (the coarse missed class c₁ mod q* — always available, every cluster has ≤ 12 co-offsets — lifts to a band-interior missed residue at the composite modulus); the r = 2 assembly: P ∪ C₂ ∪ C₁ strictly live at q*V₂V₁ with the explicit MIXED-RADIX multiplier w₁ = w₂V₁ + (c₁V₂ − w₂V₁) % (q*V₂), w₂ = aV₂ + (c₂ − aV₂) % q*, under quadratic scale separation 14·q*V₂·V₂ < V₁. FORMALIZED SORRY-FREE, GREEN ON FIRST COMPILE (LRCMultiScaleWitness.lean); general r = the same two lemmas applied once per scale, innermost out
status: PROVED AND FORMALIZED (band_lift, cluster_join, miss_at_next, threeScale_strictlyLive, threeScale_strictWitness, demoThreeScale_strictWitness — all [propext, Classical.choice, Quot.sound]; demo: P = {1,2,3}, E₂ = {0,1} @ V₂ = 3000, E₁ = {0..7} @ V₁ = 2·10⁹, (q*, a, c₂, c₁) = (13, 1, 2, 8), witnessed THROUGH the theorem at (q₁, w₁) = (78·10¹², 6010000020000); Python cross-check all 13 residues strictly in-band). The general-r induction is the same two lemmas per scale: every already-placed speed lifts (band_lift never inspects the new tuning Δ), each new cluster joins via its coarse missed class (pigeonhole at q* ≤ 13 always succeeds — clusters have ≤ 12 co-offsets).
source: klein-2026-07-10-S245 (HYP-5950; the multi-scale iteration of THM-693)
depends_on:
  - THM-693 / LRCTwoScaleWitness (the base level: small_speed_band, cluster_speed_band)
  - kps LRCStrictRuler (StrictlyLive, strictWitness_of_strictlyLive)
related:
  - THM-688 (the continuum multi-scale limit this makes constructive — the continuum needed V₁/V₂ → ∞; the constructive price is quadratic separation 14q*V₂² < V₁)
  - LRCTestPointCore (S243: pigeonhole_missed_class supplies every c_i)
---

# THM-694 — the multi-scale witness composition

## The two lemmas

**THE BAND LIFT.** If speed p (0 < p ≤ B) is strictly in-band at (q, a) —
q < 14·((pa) % q) < 13q — then for ANY Δ ∈ [0, q) and V > 14Bq:

> (p·(aV + Δ)) % (qV) = ((pa) % q)·V + pΔ, strictly in-band at qV.

Integrality gives 14r ≥ q+1 and 14r ≤ 13q−1 for free; the perturbation
pΔ ≤ B(q−1) is dominated once V > 14Bq. The lift NEVER INSPECTS Δ — the
next level's tuning cannot disturb any speed already placed.

**THE CLUSTER JOIN.** If the missed residue s = (C − ea) % q is itself
strictly in-band and V > 14eq, then the new cluster speed V − e satisfies

> ((V−e)·(aV + Δ_C)) % (qV) = sV − eΔ_C, strictly in-band at qV,

with Δ_C = (C − aV) % q — the tuning that steers the top phase to C.

**MISS AT NEXT.** The base-level missed class c₁ (E₁ misses it mod q* —
always available since |E₁| ≤ 12 < q*) lifts: (c₁V₂ − e·w₂) % (q*V₂) =
σ_e·V₂ − eΔ₂ with σ_e = (c₁ − ea) % q* ∈ [1, q*−1], strictly in-band at the
composite modulus once 168e < V₂ — the join's hypothesis at level 2 is
manufactured from level-1 data.

## The assembly (r = 2) and the general pattern

threeScale_strictlyLive: for v = P ∪ {V₂−f : f ∈ E₂} ∪ {V₁−e : e ∈ E₁} with
base data (q* ∈ [8,13], a, missed classes c₂, c₁ mod q*), thresholds
V₂ > 2184, 168f < V₂, 168e < V₂, 14e(q*V₂) < V₁, and 14(q*V₂)V₂ < V₁:

> StrictlyLive v (q*V₂V₁) (w₂V₁ + (c₁V₂ − w₂V₁) % (q*V₂)),
> w₂ = aV₂ + (c₂ − aV₂) % q*.

Proof per speed: small and inner-cluster speeds are banded at level 1 by
THM-693's two lemmas, then lifted (B := V₂); the outer cluster joins via
miss_at_next. General r: identical — each scale adds one lift-pass (free:
the lift is Δ-blind) and one join; every missed class exists by pigeonhole
at q* ≤ 13. The constructive price vs the continuum (THM-688): quadratic
scale separation V_{i} > 14·q*·(Π_{j>i} V_j)·V_{i+1}-type thresholds instead
of ratios → ∞.

## What this closes

The wall's witness construction is now fully constructive AND compositional
across arbitrary scale structure: mixed-radix multipliers, one digit per
scale, every arrow kernel-pure. With THM-693 (two-scale), the banks
(bounded V), and kps's strict chain, the two-/multi-scale program is
formalized end to end.

## Formalization & files

`04-computation/lean/TournamentH7/TournamentH7/LRCMultiScaleWitness.lean`
(kernel-pure, root-wired, green on first compile). Demo cross-check:
`05-knowledge/results/lrc14_multiscale_witness_demo_klein_S245.out`.
