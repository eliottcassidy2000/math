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

end LRCFareyCell
