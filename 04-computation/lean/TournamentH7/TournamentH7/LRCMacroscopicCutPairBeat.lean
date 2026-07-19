/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Codex (LRC multi-agent project, 2026-07-19)
-/
import Mathlib.Tactic

/-!
# Arithmetic consumers for THM-1238 pair-sum beat dichotomy

The paper proof supplies the consecutive-beat run argument.  This module
checks its ratio/length regimes, the exact seam divisibility, the elementary
third-blocker logic, and the infinite parity-locked residual ledger.
-/

namespace LRC14
namespace MacroscopicCutPairBeat

/-- In the moderate-ratio branch, the pair-sum clock step lies in
`[1/7,1/2)`, so two consecutive residues cannot both fit in the strict danger
arc of length `1/7`. -/
theorem moderate_ratio_step
    {x y : ℝ} (hx : 0 < x) (hxy : x < y) (hy6x : y ≤ 6 * x) :
    1 / 7 ≤ x / (x + y) ∧ x / (x + y) < 1 / 2 := by
  have hq : 0 < x + y := by linarith
  constructor
  · apply (le_div_iff₀ hq).2
    nlinarith
  · apply (div_lt_iff₀ hq).2
    nlinarith

/-- In the far-ratio branch the step is strictly shorter than the danger arc. -/
theorem far_ratio_step
    {x y : ℝ} (hx : 0 < x) (hy6x : 6 * x < y) :
    0 < x / (x + y) ∧ x / (x + y) < 1 / 7 := by
  have hq : 0 < x + y := by linarith
  constructor
  · positivity
  · apply (div_lt_iff₀ hq).2
    nlinarith

/-- Pair sum at least `7c/3` makes the scaled slow gap at least two beat
spacings long. -/
theorem two_beat_length
    {c q : ℝ} (hc : 0 < c) (hq : 7 * c / 3 ≤ q) :
    2 ≤ 6 * q / (7 * c) := by
  apply (le_div_iff₀ (by positivity : 0 < 7 * c)).2
  nlinarith

/-- Every fast-fast pair sum exceeds `2c`; in the residual horn its scaled
gap length lies strictly between `12/7` and `2`. -/
theorem residual_length_horn
    {c x y : ℝ} (hcx : c < x) (hcy : c < y)
    (hq : x + y < 7 * c / 3) (hc : 0 < c) :
    12 / 7 < 6 * (x + y) / (7 * c) ∧
      6 * (x + y) / (7 * c) < 2 := by
  have hden : 0 < 7 * c := by positivity
  constructor
  · apply (div_lt_div_iff₀ (by norm_num : (0 : ℝ) < 7) hden).2
    nlinarith
  · apply (div_lt_iff₀ hden).2
    nlinarith

/-- Zero active-pair curvature has the exact THM-1156 seam divisibility. -/
theorem zero_curvature_divisibility {D q : ℕ} (hzero : 14 * D = q) :
    14 ∣ q := by
  exact ⟨D, hzero.symm⟩

/-- If the defining pair is safe at a point covered by six combs, one of the
other four labels is a third blocker. -/
theorem third_blocker_logic
    {A₁ A₂ A₃ A₄ A₅ A₆ : Prop}
    (hcover : A₁ ∨ A₂ ∨ A₃ ∨ A₄ ∨ A₅ ∨ A₆)
    (h₁ : ¬A₁) (h₂ : ¬A₂) :
    A₃ ∨ A₄ ∨ A₅ ∨ A₆ := by
  tauto

/-- For even defining speeds, the half-beat numerator is the common zero
residue for both pair masks. -/
theorem even_pair_halfbeat (a b : ℕ) :
    ∃ q p : ℕ, q = 2 * a + 2 * b ∧ 2 * p = q ∧
      q ∣ (2 * a) * p ∧ q ∣ (2 * b) * p := by
  refine ⟨2 * (a + b), a + b, by ring, by ring, ?_, ?_⟩
  · exact ⟨a, by ring⟩
  · exact ⟨b, by ring⟩

/-- Exact arithmetic ledger for the infinite parity-locked family. -/
theorem parity_locked_packet_ledger (m : ℕ) (hm : 1 ≤ m) :
    let c := 420 * m + 1
    let d₁ := 428 * m + 2
    let d₂ := 440 * m + 2
    let d₃ := 452 * m + 2
    let d₄ := 464 * m + 2
    let d₅ := 476 * m + 2
    let d₆ := 500 * m + 2
    d₁ < d₂ ∧ d₂ < d₃ ∧ d₃ < d₄ ∧ d₄ < d₅ ∧ d₅ < d₆ ∧
      6 * d₅ < 7 * c ∧ 7 * c < 6 * d₆ ∧
      14 * ((d₆ - d₁) + (d₆ - d₂) + (d₆ - d₃) +
        (d₆ - d₄) + (d₆ - d₅)) > d₆ ∧
      d₆ - d₅ = 24 * m ∧ 180 * (d₆ - d₅) > c := by
  dsimp
  omega

#print axioms moderate_ratio_step
#print axioms far_ratio_step
#print axioms two_beat_length
#print axioms residual_length_horn
#print axioms zero_curvature_divisibility
#print axioms third_blocker_logic
#print axioms even_pair_halfbeat
#print axioms parity_locked_packet_ledger

end MacroscopicCutPairBeat
end LRC14
