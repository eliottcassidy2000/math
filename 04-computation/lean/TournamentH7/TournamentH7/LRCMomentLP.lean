/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S191)
-/
import Mathlib

/-!
# The one-sided moment-LP core (THM-661, the analytic link for the per-k D3 bars)

The per-k density-floor bars `rhoGlobFloorRat(k) ≤ μ(GOOD)` (`= P(W>0)`, `W = Σ(gapᵢ − 1/7)₊`) are proved
by the moment-LP route (THM-661): pick a polynomial `p` with `p(w) ≤ 1_{w>0}` on `[0, 6/7]`; then
`Σᵢ cᵢ E[Wⁱ] = ∫ p(W) dμ ≤ μ(W>0)`. `LRCD3FloorCert` (kps-S89) certifies the *rational* moment bound
`Σ cᵢ mᵢ ≥ bar` per shape; `GoodSetBridge` (death-star, AP certs) proves `μ(GOOD) ≥ m_P` for bounded
diameter. The MISSING analytic link between them is the moment-LP inequality `∫ p(W) ≤ μ(W>0)` — the
Markov–Krein / Paley–Zygmund core, formalized here abstractly and reusably.

`integral_le_measure_pos` : on a probability space, for measurable `W` and any `p` with
`p(w) ≤ 1_{w>0}` pointwise, `∫ p∘W dμ ≤ (μ{W>0}).toReal`. Instantiated with a feasible degree-`d`
polynomial `p`, and `∫ p(W) = Σ cᵢ E[Wⁱ]`, this is EXACTLY `D_d(E) ≤ μ(GOOD)` — so `bar ≤ D_d(E) ≤
μ(GOOD)`, reducing the per-k bar to (i) the rational moment bound `bar ≤ Σ cᵢ mᵢ` (LRCD3FloorCert,
native_decide) and (ii) the moment identity `E[Wⁱ] = ∫ Wⁱ` (Farey-cell integration). The remaining
analytic content is the min over ALL clusters (compact check + the decorrelation tail, LEM-005).

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LRCMomentLP

open MeasureTheory Set

variable {Ω : Type*} [MeasurableSpace Ω]

/-- **The one-sided moment-LP core (THM-661).** On a probability space, for a measurable `W` and any
function `p : ℝ → ℝ` that is pointwise below the positivity indicator (`p(w) ≤ 1` if `w > 0`, `p(w) ≤ 0`
if `w ≤ 0`), the integral of `p∘W` lower-bounds the good-set measure `μ{W>0}`:

  `∫ p(W x) dμ(x) ≤ (μ {x | 0 < W x}).toReal`.

With `p` a feasible one-sided moment polynomial this is `Σᵢ cᵢ E[Wⁱ] ≤ μ(GOOD)` — the rigorous reduction
of the density floor to a moment bound. -/
theorem integral_le_measure_pos (μ : Measure Ω) [IsProbabilityMeasure μ]
    (W : Ω → ℝ) (hW : Measurable W)
    (p : ℝ → ℝ) (hint : Integrable (fun x => p (W x)) μ)
    (hbound : ∀ w : ℝ, p w ≤ ({y : ℝ | 0 < y}).indicator (fun _ => (1 : ℝ)) w) :
    ∫ x, p (W x) ∂μ ≤ (μ {x | 0 < W x}).toReal := by
  have hmeas : MeasurableSet {x : Ω | 0 < W x} := by
    have : {x : Ω | 0 < W x} = W ⁻¹' (Ioi 0) := by ext x; simp
    rw [this]; exact hW measurableSet_Ioi
  -- the positivity indicator of `W x`, transported to the domain
  have hindint : Integrable (fun x => ({x : Ω | 0 < W x}).indicator (fun _ => (1 : ℝ)) x) μ :=
    (integrable_const (1 : ℝ)).indicator hmeas
  have key : ∀ x, p (W x) ≤ ({x : Ω | 0 < W x}).indicator (fun _ => (1 : ℝ)) x := by
    intro x
    have hb := hbound (W x)
    simp only [Set.indicator_apply, Set.mem_setOf_eq] at hb ⊢
    exact hb
  calc ∫ x, p (W x) ∂μ
      ≤ ∫ x, ({x : Ω | 0 < W x}).indicator (fun _ => (1 : ℝ)) x ∂μ := integral_mono hint hindint key
    _ = (μ {x | 0 < W x}).toReal := integral_indicator_one hmeas

/-- Packaged for the density floor: if a feasible moment functional `L = ∫ p(W)` clears a rational bar
`b`, then so does the good-set measure `μ{W>0} ≥ b`. (Chains `b ≤ L` with the moment-LP core.) -/
theorem measure_pos_ge_of_moment_ge (μ : Measure Ω) [IsProbabilityMeasure μ]
    (W : Ω → ℝ) (hW : Measurable W)
    (p : ℝ → ℝ) (hint : Integrable (fun x => p (W x)) μ)
    (hbound : ∀ w : ℝ, p w ≤ ({y : ℝ | 0 < y}).indicator (fun _ => (1 : ℝ)) w)
    (b : ℝ) (hb : b ≤ ∫ x, p (W x) ∂μ) :
    b ≤ (μ {x | 0 < W x}).toReal :=
  le_trans hb (integral_le_measure_pos μ W hW p hint hbound)

/-- **The moment-LP bound from polynomial coefficients (THM-661, packaged).** For a degree-`d` polynomial
`p(w) = Σ_{i≤d} cᵢ wⁱ` that is feasible (`p(w) ≤ 1_{w>0}` pointwise), the moment functional
`Σ_{i≤d} cᵢ E[Wⁱ]` lower-bounds the good-set measure `μ{W>0}`. This is the general/abstract content of the
D3/B_d density-floor route: it reduces the per-k bar `bar ≤ μ(GOOD)` to (a) the RATIONAL moment bound
`bar ≤ Σ cᵢ mᵢ` (LRCD3FloorCert, native_decide) and (b) the moment identity `mᵢ = E[Wⁱ] = ∫ Wⁱ` (the
Farey-cell integration) — nothing else. -/
theorem momentLP_from_coeffs (μ : Measure Ω) [IsProbabilityMeasure μ]
    (W : Ω → ℝ) (hW : Measurable W) (d : ℕ) (c : ℕ → ℝ)
    (hint : ∀ i ∈ Finset.range (d + 1), Integrable (fun x => c i * W x ^ i) μ)
    (hfeasible : ∀ w : ℝ, (∑ i ∈ Finset.range (d + 1), c i * w ^ i)
        ≤ ({y : ℝ | 0 < y}).indicator (fun _ => (1 : ℝ)) w) :
    (∑ i ∈ Finset.range (d + 1), c i * ∫ x, W x ^ i ∂μ) ≤ (μ {x | 0 < W x}).toReal := by
  have hpint : Integrable (fun x => ∑ i ∈ Finset.range (d + 1), c i * W x ^ i) μ :=
    integrable_finsetSum _ hint
  have hlin : (∑ i ∈ Finset.range (d + 1), c i * ∫ x, W x ^ i ∂μ)
      = ∫ x, (∑ i ∈ Finset.range (d + 1), c i * W x ^ i) ∂μ := by
    rw [integral_finsetSum _ hint]
    exact Finset.sum_congr rfl (fun i _ => (integral_const_mul (c i) _).symm)
  rw [hlin]
  exact integral_le_measure_pos μ W hW (fun w => ∑ i ∈ Finset.range (d + 1), c i * w ^ i)
    hpint hfeasible

end LRCMomentLP
