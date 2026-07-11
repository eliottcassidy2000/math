---
id: THM-688
title: The multi-scale limit and the rational-ray cover — middle speeds completed. (A) SEPARATED SCALES: S = P ∪ C₁ ∪ … ∪ C_r (clusters C_i = {V_i − e : e ∈ E_i}, all ratios V_i/V_{i+1} → ∞): μ(S) → ∫_{G_P} Π_i m_{E_i}(α) dα with rate C·(1/V_r + Σ V_{i+1}/V_i) — independent fibers, iterated two-scale slicing, exact product evaluator (piecewise-polynomial cells); (B) RATIONAL RAYS: a middle cluster at V_mid ≈ (a/b)·V_top is COUPLED through the b-fold cover — the product limit FAILS by a persistent correction (+0.0124 measured at b = 2) and the exact b-cover joint limit is the truth (49379/470400 for the demo class, |err|·V ≤ 0.29); (C) floors: independent clusters contribute (1−k_i/7) factors, ray-coupled clusters merge (k_i + k_j ≤ 6 ⟹ unconditional); (D) with THM-685 + kps's strict chain: THE WALL on every multi-scale slice with positive limit — and the limit is exactly computable in every case of the taxonomy
status: PROVED ((A) the iterated slice argument — each level is THM-687(B) verbatim with the slower structure frozen; rate constants explicit per level; (B) the b = 2 cover derived exactly (fiber γ = 2β on the circumference-2 circle: top-cluster arcs appear twice at width 1/7, mid-cluster arcs once at width 2/7), evaluator exact, convergence verified at V₂ = 300/600/1200; general (a,b) stated with the (β, j mod b)-fiber, implementation mechanical; (C) arc-count bounds per fiber circle). VERIFIED: all convergence tables exact (companion output); the three-scale sanity case included.
source: klein-2026-07-10-S237 (HYP-5910; the named multi-scale follow-up of THM-687)
depends_on:
  - THM-687   # the two-scale base case and the wall composition
  - THM-685   # the modulus transfer
related:
  - kps LRCStrictRuler (the strict chain), opus LRCPrimitivePeel
  - THM-686 (the bounded-V complement)
---

# THM-688 — the multi-scale limit and the rational-ray cover

## (A) Separated scales (proved)

S = P ∪ C₁ ∪ … ∪ C_r with C_i = {V_i − e : e ∈ E_i}, |E_i| = k_i,
V₁ > V₂ > … > V_r ≫ 13. Each cluster carries its own fast phase
β_i = frac(V_i α); when all ratios diverge the phases decouple:

> **μ(S) = ∫_{G_P} Π_{i=1}^r m_{E_i}(α) dα + O(1/V_r + Σ_i V_{i+1}/V_i)**

(m_E as in THM-687). *Proof:* induct on r; slice α at the slowest cluster
scale V_r and apply THM-687(B)'s freeze-and-count with the inner (faster)
structure replaced by its own limit — the inner application contributes the
V_{i+1}/V_i terms (the faster phase sweeps V_i/V_{i+1} times per slice), the
outer the 1/V_r. Each level is the identical three-error bookkeeping
(boundary slices, center Lipschitz, piecewise Riemann). ∎

The limit is EXACTLY computable: on each cell of the common breakpoint
refinement every m_{E_i} is linear (fit from interior thirds); the product
is a degree-r polynomial integrated exactly (Fraction coefficients).
Demo class P = {1..5}, E₂ = {0,1}, E₁ = {0..5}: μ∞ = **191/2058**; the
convergence table (V₂ = 40..100, V₁ = 240..10000) tracks the predicted rate.

## (B) Rational rays — the true middle-speed case (proved at b = 2; general form stated)

A middle cluster at bounded ratio, V_mid = (a·V_top − c)/b (gcd(a,b) = 1,
c bounded), does NOT decouple: frac(V_mid·α) is determined by
β = frac(V_top·α) together with ⌊V_top·α⌋ mod b — the fiber is the b-fold
cover [0,1) × Z_b, and the joint fiber measure replaces the product.

**The b = 2 case, exact** (V_top = 2V₂, the classic "W ≈ V/2" middle speed):
on the circumference-2 circle γ = 2β₂ the mid-cluster (V₂) arcs appear once
with width 2/7 at centers 2eα, and the top-cluster arcs appear TWICE with
width 1/7 at centers e′α and e′α + 1; m_joint = (2 − |union|)/2, and
μ∞ = ∫_{G_P} m_joint — exact evaluator in the companion script. Measured on
the demo class:

> product limit 0.092809 (WRONG: persistent error +0.0124 at every V) vs
> **2-cover limit 49379/470400 ≈ 0.104972** — the exact sweeps converge to
> it at the C/V rate (|err|·V₂ = 0.28, 0.04, 0.29 at V₂ = 300, 600, 1200).

Notably the coupling correction is POSITIVE here — the ray overlap wastes
less of the fiber circle; no sign claim in general (the evaluator decides
per class). As a → ∞ (integer-multiple ratios with growing multiplier) the
a-cover arcs spread and the correction decays — integer-ratio sequences with
ratio → ∞ correctly fall under (A).

**The taxonomy is complete for two clusters** (and pairwise for more):
bounded ratio → rational ray (b-cover limit, exact); diverging ratio →
separated (product limit, exact). Every middle-speed structure has an
exactly computable limit measure.

## (C) Floors

- Independent (separated) clusters: m_{E_i} ≥ 1 − k_i/7 each, so
  **μ∞ ≥ Π_i (1 − k_i/7) · m_P** — unconditional when every k_i ≤ 6
  (m_P ≥ 1/(91·maxP) as in THM-687(C)).
- Ray-coupled pairs MERGE for the floor: on the b-cover circle the total arc
  mass is (k_mid + k_top)·(2/7)·(1/2)-normalized — the pair behaves as one
  cluster of size k_mid + k_top: unconditional when the SUM is ≤ 6; finite
  criterion otherwise (the evaluator).
- At most one merged cluster can exceed 6 (thirteen runners), so every
  multi-scale class is [unconditional] or [one exact positivity check].

## (D) The wall on multi-scale slices

Composition identical to THM-687(E): limit measure > 0 and scales beyond the
explicit thresholds ⟹ μ(S) > 0 ⟹ THM-685 at 14-coprime q ⟹ StrictlyLive
(strictness free) ⟹ kps's strict chain ⟹ strict witness. With (A)+(B) the
limit is computable for EVERY scale structure, so the wall on the unbounded
residual reduces everywhere to per-class exact positivity checks + the
bounded-V bank/census territory (THM-686).

## Honest scope

General (a,b) ray evaluator: derived, not implemented (mechanical — the
(β, s)-fiber sum); three-scale ray mixtures compose pairwise. The universal
positivity of the limit measures (the dead-zone lemma) remains the named
open structural question — per-class decidable today, extremal shapes
verified positive (THM-687(D)).

## Verification & files

`04-computation/lrc14_multiscale_middle_speeds_klein_S237.py` (+ `.out`):
the product evaluator, the ratio study, the b = 2 cover evaluator and its
convergence, the three-scale case, floor checks.

## The general-(a,b) evaluator (S246)

`lrc14_general_ray_evaluator_klein_S246.py`: the exact limit measure for ANY
ray a/b (gcd = 1, shift c) — the joint fiber (β, s′) ∈ [0,1) × Z_b with the
mid-phase φ = (s′ + aβ)/b − (c/b)α; per (α, s′) the mid-cluster bad set is
aβ ∈ (c + bf)α − s′ ± b/14 (mod b), clipped to [0, a); m_joint = (1/b)Σ_s′
meas(safe β); μ∞ = ∫_{G_P} m_joint by piecewise-linear trapezoid (endpoint
rates e and (c+bf)/a). VERIFIED: exact regression on the b = 2 special case
(49379/470400, equality); new rays (2,3,0), (1,3,0), (3,4,1) against exact
finite-V sweeps, |err|·V ≤ 0.94 (the C/V rate); the merged floor; a 36-class
ray positivity census (all μ∞ > 0, including top-block P's). The taxonomy's
"stated, implementation mechanical" item is now implemented. Named follow-up:
the CONSTRUCTIVE ray witness (the mixed-radix machinery of THM-694 with one
more digit at modulus q*·b·V — same two lemmas).
