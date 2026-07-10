/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S193)
-/
import Mathlib

/-!
# The single-cell Farey moment identity (THM-661 / LRCD3FloorCert)

The moment route to the density floor computes `m_i = ∫_0^1 W(x)^i dx` (`W = Σ(gapⱼ − 1/7)₊`, the
uncovered measure) by partitioning `[0,1]` into **Farey cells** on each of which `W` is affine,
`W(x) = A + B·x`, and summing the exact per-cell integrals. This file proves the **single-cell** step —
the atomic exact contribution of one cell — which, summed over cells, is the Farey moment identity
`m_i = Σ_cells ∫_cell W^i` (`opus-S192 momentLP_from_coeffs` is the only other piece of the moment route,
plus the LRCD3FloorCert rational bar).

`cell_moment` : for `B ≠ 0`, `∫_a^b (A + B·x)^i dx = ((A + B·b)^{i+1} − (A + B·a)^{i+1}) / (B·(i+1))`.
`cell_moment_const` : the degenerate `B = 0` (constant `W = A`) cell, `∫_a^b A^i = A^i·(b − a)`.

Proof: FTC with antiderivative `F(x) = (A + B·x)^{i+1} / (B·(i+1))` (`F' = (A + B·x)^i`).

Kernel-pure: no `sorry`, no `native_decide`. Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LRCFareyCell

open intervalIntegral

/-- **The single-cell Farey moment identity** (`B ≠ 0`, the generic affine cell). On a Farey cell `[a,b]`
where `W(x) = A + B·x`, the degree-`i` moment contribution is the exact closed form. -/
theorem cell_moment (A B : ℝ) (hB : B ≠ 0) (i : ℕ) (a b : ℝ) :
    (∫ x in a..b, (A + B * x) ^ i)
      = ((A + B * b) ^ (i + 1) - (A + B * a) ^ (i + 1)) / (B * ((i : ℝ) + 1)) := by
  have hi1 : ((i : ℝ) + 1) ≠ 0 := by positivity
  -- `F(x) = (A + B x)^(i+1) / (B (i+1))` has derivative `(A + B x)^i`
  have hderiv : ∀ x : ℝ,
      HasDerivAt (fun x => (A + B * x) ^ (i + 1) / (B * ((i : ℝ) + 1))) ((A + B * x) ^ i) x := by
    intro x
    have h1 : HasDerivAt (fun x => A + B * x) B x := by
      simpa using ((hasDerivAt_id x).const_mul B).const_add A
    have h2 : HasDerivAt (fun x => (A + B * x) ^ (i + 1))
        ((↑(i + 1)) * (A + B * x) ^ ((i + 1) - 1) * B) x := h1.pow (i + 1)
    have h3 := h2.div_const (B * ((i : ℝ) + 1))
    convert h3 using 1
    rw [Nat.add_sub_cancel]
    push_cast
    field_simp
  have hcont : Continuous fun x : ℝ => (A + B * x) ^ i := by fun_prop
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ => hderiv x)
      (hcont.intervalIntegrable a b), div_sub_div_same]

/-- The degenerate constant cell (`B = 0`, `W ≡ A`): `∫_a^b A^i = (b − a)·A^i`. -/
theorem cell_moment_const (A : ℝ) (i : ℕ) (a b : ℝ) :
    (∫ _x in a..b, (A) ^ i) = (b - a) * A ^ i := by
  rw [intervalIntegral.integral_const, smul_eq_mul]

open MeasureTheory in
/-- **The Farey cell decomposition = the Farey moment identity.** Given a partition
`0 = t₀ < t₁ < ⋯ < t_N = 1` of `[0,1]` on each cell of which the uncovered measure `W` is affine
(`W(x) = Aⱼ + Bⱼ·x` on `[tⱼ, tⱼ₊₁]`, `Bⱼ ≠ 0`), the degree-`i` moment `∫₀¹ Wⁱ` is the sum of the exact
per-cell contributions (`cell_moment`). Combined with `momentLP_from_coeffs` (opus-S192), this closes the
NUMERIC content of the moment route; the ONLY input still tied to the concrete `W` is the per-cell
affineness `haffine` (which the Farey-breakpoint / three-distance structure supplies). -/
theorem farey_moment_decomp (W : ℝ → ℝ) (i N : ℕ) (t : ℕ → ℝ) (A B : ℕ → ℝ)
    (ht0 : t 0 = 0) (htN : t N = 1) (hB : ∀ j, j < N → B j ≠ 0)
    (haffine : ∀ j, j < N → ∀ x ∈ Set.uIcc (t j) (t (j + 1)), W x = A j + B j * x)
    (hint : ∀ j, j < N → IntervalIntegrable (fun x => W x ^ i) volume (t j) (t (j + 1))) :
    (∫ x in (0 : ℝ)..1, W x ^ i)
      = ∑ j ∈ Finset.range N,
          ((A j + B j * t (j + 1)) ^ (i + 1) - (A j + B j * t j) ^ (i + 1)) / (B j * ((i : ℝ) + 1)) := by
  have hbound : (∫ x in (0 : ℝ)..1, W x ^ i) = ∫ x in (t 0)..(t N), W x ^ i := by rw [ht0, htN]
  rw [hbound, ← intervalIntegral.sum_integral_adjacent_intervals hint]
  refine Finset.sum_congr rfl (fun j hj => ?_)
  rw [Finset.mem_range] at hj
  rw [intervalIntegral.integral_congr (g := fun x => (A j + B j * x) ^ i)
      (fun x hx => by rw [haffine j hj x hx])]
  exact cell_moment (A j) (B j) (hB j hj) i (t j) (t (j + 1))

end LRCFareyCell
