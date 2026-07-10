---
id: THM-687
title: The two-scale cluster limit — for S_V = P ∪ {V−e : e ∈ E} (small speeds P ⊆ [1,13], fixed co-offsets E, V → ∞): (A) μ(S_V) → μ∞(P,E) = ∫_{G_P} m_E(α) dα with the cluster fiber measure m_E = 1 − meas(∪_e arc(eα, 1/14)), EXACTLY computable per class; (B) THE TWO-SCALE TRANSFER |μ(S_V) − μ∞| ≤ C(P,E)/V with C explicit (measured C ≤ 0.97, worst at the deep well's resonant V = 182); (C) THE k ≤ 6 UNCONDITIONAL FLOOR μ∞ ≥ (1−k/7)·m_P ≥ (1−k/7)/(91·maxP) > 0 free from LRC(≤13); (D) THE APEX POSITIVITY: the consecutive co-offset shape E = {0..k−1} — the measured MINIMIZER — has μ∞ > 0 at every k = 7..13; composed with THM-685 (transfer at 14∤q) and kps's LRCStrictRuler chain, THE WALL (strict live modulus) holds on every two-scale slice with μ∞ > 0 beyond an explicit V₀ — the unbounded direction of the residual reduces to per-class finite positivity checks plus bounded-V bank territory
status: PROVED ((A) dominated convergence along the slice decomposition — the proof of (B) below gives it with rate; (B) complete elementary proof, slice + crossing counting, crude explicit constant, measured C ≤ 0.97 across five families and V up to 2000; (C) two-line composition of the pointwise arc bound with THM-667 Lemma A / LRC(≤13)-strict; (D) exact rational evaluations — per-class decidable criterion; the SHAPE-minimality of consecutive E is measured (41-candidate battery), not proved). VERIFIED: k=1 exact identity μ∞ = (6/7)·m_P (6617/226380 for the deep-well family — the deep well {1..12,182} is the V = 182 sample, its μ = 4637/194040 sits 0.0053 below the limit, the largest resonant deviation observed); all convergence tables in the companion outputs.
source: klein-2026-07-10-S236 (HYP-5905; the measure-floor program's unbounded half, answering kps-S127's WALL on the bounded-co-offset slices)
depends_on:
  - THM-685   # the modulus transfer the floors feed
  - THM-667 (Lemma A)   # m_P ≥ 1/(91·maxP) from LRC(≤13)-strict
related:
  - kps LRCStrictRuler / LRCStrictWitnessFloor (the strict integer chain consuming these floors)
  - opus LRCPrimitivePeel (primitivity; the slices here are primitive for gcd(P)=1)
  - THM-686 (bounded-window censuses = the V ≤ V₀ complement)
  - THM-667 / the mid-band frame (m_E is the drift/fiber object, now exactly integrable)
---

# THM-687 — the two-scale cluster limit and the wall on bounded-co-offset slices

## Setting

S_V = P ∪ {V − e : e ∈ E}, with P ⊆ [1,13] the small speeds (|P| = 13 − k),
E = {0 = e₁ < … < e_k} the fixed co-offsets, V → ∞. Substituting
β = frac(Vα): frac((V−e)α) = frac(β − eα), so the cluster runner V−e is safe
iff β avoids the **bad arc** (eα − 1/14, eα + 1/14) — k arcs of length 1/7
centered at the points eα. Define

> **m_E(α) = 1 − meas(∪_{e∈E} (eα − 1/14, eα + 1/14))** (the fiber measure),
> **μ∞(P,E) = ∫_{G_P} m_E(α) dα**, G_P = {α : frac(pα) ∈ B̃ ∀p ∈ P}.

m_E is piecewise linear with breakpoints where frac((e−e′)α) ∈ {0, 1/7, 6/7};
μ∞ is EXACTLY computable (cell decomposition + trapezoid on interior thirds,
all rational — the evaluator in the companion script).

## (B) The two-scale transfer (proved)

> **|μ(S_V) − μ∞(P,E)| ≤ C(P,E)/V**, C ≤ 2K_P + 2Σ_E e + 2N_cells
> (K_P = components of G_P ≤ Σp; N_cells = the m_E cell count
> ≤ 2Σp + 3Σ_{pairs}|e−e′| + 2).

*Proof.* Partition [0,1) into V slices I_j = [j/V, (j+1)/V); on I_j write
α = (j+β)/V. (i) Freezing the P-conditions at α_j = j/V changes the slice
integral only on slices meeting ∂G_P: ≤ 2K_P slices, error ≤ 1/V each.
(ii) Freezing the arc centers at eα_j moves each center by ≤ e/V across the
slice; the avoided-set measure is Lipschitz in each center with constant 2,
so the per-slice error is ≤ 2Σe/V, and these average (weight 1/V per slice)
to ≤ 2Σe/V total. (iii) The frozen sum (1/V)Σ_j 1_{G_P}(α_j)m_E(α_j) is a
Riemann sum of a piecewise-linear function with N_cells pieces and values in
[0,1]: error ≤ 2N_cells/V. Sum the three. ∎

Measured (five families, V ≤ 2000): **C ≤ 0.97** — the crude constant is
~100–300× conservative; the maximum occurs at the deep-well resonance
(family ({1..12},{0}) at V = 182 = 13·14, error −0.0053·1). The k = 1 limit
is exact in closed form: μ∞(P,{0}) = (6/7)·m_P (verified: 6617/226380).

## (C) The k ≤ 6 unconditional floor (proved)

k arcs of length 1/7 cannot cover the circle, so m_E ≥ 1 − k/7 pointwise,
hence μ∞ ≥ (1 − k/7)·m_P. For |P| = 13 − k ≤ 12, LRC(≤13)-strict gives a
time with all P-clearances ≥ 1/13, and the Lipschitz ball of radius
1/(182·maxP) stays ≥ 1/14 (THM-667 Lemma A): **m_P ≥ 1/(91·maxP)**. So

> **μ∞(P,E) ≥ (1 − k/7)/(91·maxP) > 0 for every k ≤ 6 — no structure
> hypothesis on E whatsoever.**

## (D) k ≥ 7: the finite criterion and the apex positivity (exact; shape-minimality measured)

For k ≥ 7 the arcs can cover (m_E can vanish on a **dead zone**); μ∞ > 0 is
a per-class FINITE exact computation. Evaluations:

- The consecutive apex E = {0,…,6}, P = {1..6}: **μ∞ = 2171/29400 ≈ 0.0738**
  — positive; and it is the MINIMUM over a 41-shape adversarial battery
  (APs d = 1..5, detuned APs, geometric, paired, 30 random) — every other
  shape is larger (spread E = {0,2,5,9,14,20,27}: 0.1559). The 7-structured
  obstruction is confirmed as the extremal SHAPE, and the extremum is
  POSITIVE: the dead zone never exhausts G_P.
- The consecutive apex at every depth: k = 8: 0.0489, k = 9: 0.0447,
  k = 10: 0.0448, k = 11: 0.0586, k = 12: 0.0750, k = 13 (P = {1}): 0.0654 —
  all positive, minimum near k = 9–10.

(The universal statement "μ∞ > 0 for all (P,E)" is NOT claimed — the
criterion is decidable per class and the extremal-shape evidence is strong;
a proof would need the dead-zone/G_P disjointness structurally.)

## (E) The wall on the two-scale slices (composition)

For any class with μ∞(P,E) > 0 and any V > V₀ := 2C(P,E)/μ∞:
μ(S_V) ≥ μ∞/2 > 0; then THM-685 at any modulus q > Σv/(μ∞/2) with 14 ∤ q
gives a live multiplier, which is automatically STRICT (14 ∤ q makes the
closed rounded band strict — kps's int_band_bound), i.e. StrictlyLive; kps's
strictWitness_of_strictlyLive completes:

> **every two-scale family S_V with μ∞(P,E) > 0 and V > V₀(P,E) has a strict
> witness — THE WALL holds on the slice beyond V₀; V ≤ V₀ is finite
> bank/census territory (THM-686 machinery).**

With (C), the k ≤ 6 slices need no per-class check at all: V₀ ≤
2C·91·maxP/(1 − k/7) explicit (measured-C version: V₀ ≈ 70–250; crude-C
version: V₀ ≈ 10⁴–10⁵, still bank-reachable).

## Honest scope

- Middle speeds (between 14 and V − e_max) are outside the two-scale form;
  the natural extension is multi-scale (one β per scale, the same slice
  argument) — named follow-up, not claimed.
- The (P,E) quantification: (C) closes k ≤ 6 wholesale; k ≥ 7 needs the
  per-class finite check (or the structural dead-zone lemma) — the checks
  are exact rationals, one evaluator call each.
- E is fixed as V → ∞; E growing with V re-enters the mid-band/drift regime
  (THM-667's territory).

## Verification & files

`04-computation/lrc14_two_scale_cluster_limit_klein_S236.py` (+ `.out`):
the evaluator, five convergence tables, the k ≤ 6 floor checks, the k = 1
closed form. `04-computation/lrc14_apex_adversarial_klein_S236.py`
(+ `.out`): the 41-shape battery and the consecutive apex at k = 8..13.
