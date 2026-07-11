/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-11-S127)
-/
import Mathlib
import TournamentH7.LRCHyperbolaBox
import TournamentH7.LRCFourierCompletionB

set_option maxHeartbeats 800000

/-!
# LEM-022 Fourier completion, Stage B.3 — the off-diagonal aggregation (kps-S127)

The middle brick of opus-S213's Stage B.3, connecting his coefficient bound to death-star's harmonic sum.

The completion identity (the remaining algebraic piece) writes the interval pair-correlation as
`C_w = b²/q + (1/q)·Σ_{h≠0} B̂(h)·conj(B̂(w·h))`.  Bounding `|C_w − b²/q|` therefore reduces to bounding the
**off-diagonal Fourier sum** `Σ_{h≠0} B̂(h)·conj(B̂(w·h))`.  This file supplies that bound, kernel-pure and
abstract in the coefficient function `bc`:

  `‖Σ_{h≠0} bc(h)·conj(bc(w·h))‖ ≤ (q²/4)·Σ_{h≠0} 1/(cdist h · cdist(w·h))`

given only the per-coefficient bound `‖bc(h)‖ ≤ q/(2·cdist h)` — which is exactly opus-S213's
`FourierCompletion.norm_bandCoeff_le` (B.2).  Mechanism: the triangle inequality `‖Σ‖ ≤ Σ‖·‖`, the
multiplicativity `‖bc(h)·conj(bc(w·h))‖ = ‖bc(h)‖·‖bc(w·h)‖`, and the termwise product of the two B.2
bounds; `w` a unit keeps `w·h ≠ 0` so both bounds apply.

**Composition (the constants close).**  The RHS sum `S` is exactly the object death-star's
`HyperbolaBox.harmonic_ratio_sum_mul_le` bounds: `S·P ≤ 20·(log₂q+1)²`.  So
`‖offdiag‖ ≤ (q²/4)·S ≤ (q²/4)·20(log₂q+1)²/P = 5q²(log₂q+1)²/P`, and dividing the completion identity by
`q` gives LEM-022's target `|C_w − b²/q| ≤ 5q(log₂q+1)²/P`.  This brick is the analytic bridge between
opus's per-cell coefficient bound (B.2) and death-star's dyadic harmonic aggregation (Step 3).

Kernel-pure: no `sorry`, no `native_decide`.  Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace FourierCompletion

open Complex Finset
open LonelyRunner.HyperbolaBox

/-- **B.3 off-diagonal aggregation.**  For any coefficient function `bc : ZMod q → ℂ` obeying the per-cell
bound `‖bc h‖ ≤ q/(2·cdist h)` (opus's `norm_bandCoeff_le`), and any unit `w`, the off-diagonal Fourier
correlation is bounded by the harmonic ratio sum:

  `‖Σ_{h≠0} bc(h)·conj(bc(w·h))‖ ≤ (q²/4)·Σ_{h≠0} 1/(cdist h · cdist(w·h))`.

The RHS is `death-star`'s `harmonic_ratio_sum_mul_le` object; composing gives LEM-022's error bound. -/
theorem offDiag_bandSum_le {q : ℕ} [NeZero q] (w : ZMod q) (hw : IsUnit w) (bc : ZMod q → ℂ)
    (hbc : ∀ h : ZMod q, h ≠ 0 → ‖bc h‖ ≤ (q : ℝ) / (2 * (cdist h : ℝ))) :
    ‖∑ h ∈ (univ : Finset (ZMod q)).filter (fun h => h ≠ 0),
        bc h * (starRingEnd ℂ) (bc (w * h))‖
      ≤ (q : ℝ) ^ 2 / 4 *
        ∑ h ∈ (univ : Finset (ZMod q)).filter (fun h => h ≠ 0),
          1 / ((cdist h : ℝ) * (cdist (w * h) : ℝ)) := by
  classical
  -- `w · h ≠ 0` when `h ≠ 0` (multiplication by a unit is injective on `0`)
  have hwh : ∀ h : ZMod q, h ≠ 0 → w * h ≠ 0 := by
    obtain ⟨u, rfl⟩ := hw
    intro h hh hz
    rw [Units.mul_right_eq_zero] at hz
    exact hh hz
  refine le_trans (norm_sum_le _ _) ?_
  rw [Finset.mul_sum]
  apply Finset.sum_le_sum
  intro h hh
  have hh0 : h ≠ 0 := (Finset.mem_filter.mp hh).2
  have hwh0 : w * h ≠ 0 := hwh h hh0
  have hc1 : (1 : ℝ) ≤ (cdist h : ℝ) := by exact_mod_cast one_le_cdist hh0
  have hc2 : (1 : ℝ) ≤ (cdist (w * h) : ℝ) := by exact_mod_cast one_le_cdist hwh0
  have hne1 : (cdist h : ℝ) ≠ 0 := by linarith
  have hne2 : (cdist (w * h) : ℝ) ≠ 0 := by linarith
  have hbnd1 : (0 : ℝ) ≤ (q : ℝ) / (2 * (cdist h : ℝ)) := by positivity
  rw [norm_mul, RCLike.norm_conj]
  calc ‖bc h‖ * ‖bc (w * h)‖
      ≤ ((q : ℝ) / (2 * (cdist h : ℝ))) * ((q : ℝ) / (2 * (cdist (w * h) : ℝ))) :=
        mul_le_mul (hbc h hh0) (hbc (w * h) hwh0) (norm_nonneg _) hbnd1
    _ = (q : ℝ) ^ 2 / 4 * (1 / ((cdist h : ℝ) * (cdist (w * h) : ℝ))) := by
        field_simp
        ring

/-- **B.3 off-diagonal bound in closed form.**  Composing the aggregation with death-star's dyadic
harmonic sum `harmonic_ratio_sum_mul_le` (`S·P ≤ 20(log₂q+1)²`) eliminates the sum: for any unit `w`, a
coefficient bound `‖bc h‖ ≤ q/(2·cdist h)`, and the ratio-lattice floor `P ≤ cdist h · cdist(w·h)`,

  `‖Σ_{h≠0} bc(h)·conj(bc(w·h))‖ ≤ 5·q²·(log₂q+1)² / P`.

Dividing by `q` (the completion identity `C_w − b²/q = (1/q)·offdiag`) gives LEM-022's target
`|C_w − b²/q| ≤ 5q(log₂q+1)²/P`.  This reduces the entire Fourier-completion node (opus-S213 Stage B.3) to
its one remaining algebraic piece — the Parseval completion identity itself. -/
theorem offDiag_bandSum_le_closed {q : ℕ} [NeZero q] (w : ZMod q) (hw : IsUnit w) (bc : ZMod q → ℂ)
    (P : ℕ) (hP : 0 < P)
    (hPmin : ∀ h : ZMod q, h ≠ 0 → P ≤ cdist h * cdist (w * h))
    (hbc : ∀ h : ZMod q, h ≠ 0 → ‖bc h‖ ≤ (q : ℝ) / (2 * (cdist h : ℝ))) :
    ‖∑ h ∈ (univ : Finset (ZMod q)).filter (fun h => h ≠ 0),
        bc h * (starRingEnd ℂ) (bc (w * h))‖
      ≤ 5 * (q : ℝ) ^ 2 * ((Nat.log 2 q : ℝ) + 1) ^ 2 / (P : ℝ) := by
  have hPR : (0 : ℝ) < (P : ℝ) := by exact_mod_cast hP
  have hagg := offDiag_bandSum_le w hw bc hbc
  have hds := LonelyRunner.HyperbolaBox.harmonic_ratio_sum_mul_le w P hPmin
  have hdsR : (∑ h ∈ (univ : Finset (ZMod q)).filter (fun h => h ≠ 0),
        (1 : ℝ) / ((cdist h : ℝ) * (cdist (w * h) : ℝ))) * (P : ℝ)
      ≤ 20 * ((Nat.log 2 q : ℝ) + 1) ^ 2 := by
    calc (∑ h ∈ (univ : Finset (ZMod q)).filter (fun h => h ≠ 0),
            (1 : ℝ) / ((cdist h : ℝ) * (cdist (w * h) : ℝ))) * (P : ℝ)
        = (((∑ h ∈ (univ : Finset (ZMod q)).filter (fun h => h ≠ 0),
            (1 : ℚ) / (cdist h * cdist (w * h))) * (P : ℚ) : ℚ) : ℝ) := by push_cast; ring
      _ ≤ ((20 * ((Nat.log 2 q : ℚ) + 1) ^ 2 : ℚ) : ℝ) := by exact_mod_cast hds
      _ = 20 * ((Nat.log 2 q : ℝ) + 1) ^ 2 := by push_cast; ring
  set S : ℝ := ∑ h ∈ (univ : Finset (ZMod q)).filter (fun h => h ≠ 0),
      (1 : ℝ) / ((cdist h : ℝ) * (cdist (w * h) : ℝ)) with hSdef
  have hSle : S ≤ 20 * ((Nat.log 2 q : ℝ) + 1) ^ 2 / (P : ℝ) := by
    rw [le_div_iff₀ hPR]; linarith [hdsR]
  calc ‖∑ h ∈ (univ : Finset (ZMod q)).filter (fun h => h ≠ 0),
          bc h * (starRingEnd ℂ) (bc (w * h))‖
      ≤ (q : ℝ) ^ 2 / 4 * S := hagg
    _ ≤ (q : ℝ) ^ 2 / 4 * (20 * ((Nat.log 2 q : ℝ) + 1) ^ 2 / (P : ℝ)) :=
        mul_le_mul_of_nonneg_left hSle (by positivity)
    _ = 5 * (q : ℝ) ^ 2 * ((Nat.log 2 q : ℝ) + 1) ^ 2 / (P : ℝ) := by field_simp; ring

end FourierCompletion
end LonelyRunner
