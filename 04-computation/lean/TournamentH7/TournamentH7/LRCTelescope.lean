/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S92)
-/
import Mathlib

/-!
# The telescope bound: the abstract core of the tiling-degradation theorem (HYP-4206)

The fixed-order telescope (S89 Lemma A, the heart of the S88–S91 desert-bound chain) in
its abstract form: `n` signed gaps `g i ≤ 0` (coverage), drifting linearly at rates
`δ i` with the coverage maintained for `N` periods, total deficit `Σ g = −σ` (the mass
identity), and total positive drift `≥ d` (the circle walk: piece speeds rise `0 → D/W`
and return, so `d = D/W` in the folded application). Conclusion: `N·d ≤ σ` — the
universal desert bound `N ≤ σW/D` (per covering sub-cluster).

Pure `Finset.sum` arithmetic: no measure theory, no combs — the folded geometry enters
only through the instantiation. Kernel-pure.
-/

namespace LonelyRunner
namespace Telescope

open Finset

/-- **The telescope bound.** Gaps `g i ≤ 0` at time `0` and at time `N` under linear
drift (`g i + N·δ i ≤ 0`), total deficit `Σ g = −σ`, total positive drift `≥ d`:
then `N·d ≤ σ`. -/
theorem telescope_bound {n : ℕ} (g δ : Fin n → ℝ) (σ d N : ℝ)
    (hN : 0 ≤ N)
    (hcover0 : ∀ i, g i ≤ 0)
    (hcoverN : ∀ i, g i + N * δ i ≤ 0)
    (hsum : ∑ i, g i = -σ)
    (hdrift : d ≤ ∑ i, max (δ i) 0) :
    N * d ≤ σ := by
  have key : ∀ i, N * max (δ i) 0 ≤ -(g i) := by
    intro i
    rcases le_total (δ i) 0 with h | h
    · rw [max_eq_right h]
      have h0 := hcover0 i
      nlinarith
    · rw [max_eq_left h]
      have hcN := hcoverN i
      linarith
  calc N * d ≤ N * ∑ i, max (δ i) 0 := by
        exact mul_le_mul_of_nonneg_left hdrift hN
    _ = ∑ i, N * max (δ i) 0 := Finset.mul_sum _ _ _
    _ ≤ ∑ i, -(g i) := Finset.sum_le_sum (fun i _ => key i)
    _ = -(∑ i, g i) := by rw [← Finset.sum_neg_distrib]
    _ = σ := by rw [hsum, neg_neg]

/-- The folded instantiation shape (documentation form): `d = D/W`, `σ = 2ρc − 1`,
gaps = the circular tiling gaps of a covered strip — `N ≤ σ·W/D`. -/
theorem desert_period_bound {n : ℕ} (g δ : Fin n → ℝ) (σ D W N : ℝ)
    (hN : 0 ≤ N) (hW : 0 < W) (hD : 0 < D)
    (hcover0 : ∀ i, g i ≤ 0)
    (hcoverN : ∀ i, g i + N * δ i ≤ 0)
    (hsum : ∑ i, g i = -σ)
    (hdrift : D / W ≤ ∑ i, max (δ i) 0) :
    N ≤ σ * W / D := by
  have h := telescope_bound g δ σ (D / W) N hN hcover0 hcoverN hsum hdrift
  have hDW : 0 < D / W := div_pos hD hW
  have heq : σ * W / D = σ / (D / W) := by field_simp
  rw [heq, le_div_iff₀ hDW]
  exact h

#print axioms telescope_bound

end Telescope
end LonelyRunner
