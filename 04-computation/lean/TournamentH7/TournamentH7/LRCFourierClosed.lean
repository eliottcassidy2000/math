/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-11-S127)
-/
import Mathlib
import TournamentH7.LRCFourierCompletionC
import TournamentH7.LRCHyperbolaBox

set_option maxHeartbeats 1200000

/-!
# LEM-022 Fourier completion — the CLOSED-FORM capstone (kps-S127)

opus-S214 (`LRCFourierCompletionC`) proved the completion identity and the difference bound
`‖C_w − b²/q‖ ≤ (1/q)·Σ_{k≠0} ‖B̂(k)‖·‖B̂(w·k)‖`.  This file closes the node: applying the per-coefficient
bound (opus-S213's `norm_bandCoeff_le`, B.2) termwise and composing death-star's dyadic harmonic aggregation
(`HyperbolaBox.harmonic_ratio_sum_mul_le`, `S·P ≤ 20(log₂q+1)²`) collapses the sum to LEM-022's target

  `‖C_w − b²/q‖ ≤ 5·q·(log₂q+1)² / P`,   `P = P(w) = min_{k≠0} cdist(k)·cdist(w·k)`.

The coefficient bound enters as a hypothesis `hcoeff` (discharged verbatim by `norm_bandCoeff_le` when the
band `B` is an interval — the LRC safe arc).  Two bridges: `P > 0` with the ratio floor forces `w·z ≠ 0` for
`z ≠ 0` (so B.2 applies to both factors, no unit hypothesis needed); and the diff-bound sum over `range q`
transports to death-star's sum over `ZMod q` by the `natCast`/`val` bijection.

With this, the OffLine/LEM-022 Fourier route is closed modulo its (proved) inputs: identity (opus C),
coefficient bound (opus B.2), harmonic aggregation (death-star), and this capstone.

Kernel-pure: no `sorry`, no `native_decide`.  Axioms: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace FourierCompletion

open Complex Finset
open LonelyRunner.HyperbolaBox

/-- **LEM-022 closed-form completion bound.**  From opus's difference bound, opus's per-coefficient bound
`hcoeff`, and the ratio-lattice floor `P ≤ cdist z · cdist(w·z)`, the interval pair-correlation error is
`‖C_w − b²/q‖ ≤ 5·q·(log₂q+1)² / P`. -/
theorem completion_closed_of_coeffBound (q : ℕ) (hq : 0 < q) (B : Finset ℕ)
    (hB : B ⊆ Finset.range q) (w : ℕ) (P : ℕ) (hP : 0 < P)
    (hPmin : ∀ z : ZMod q, z ≠ 0 → P ≤ cdist z * cdist ((w : ZMod q) * z))
    (hcoeff : ∀ j : ℤ, (j : ZMod q) ≠ 0 →
        ‖bandDFT q B j‖ ≤ (q : ℝ) / (2 * (cdist (j : ZMod q) : ℝ))) :
    ‖(corrCount q B w : ℂ) - (B.card : ℂ) ^ 2 / q‖
      ≤ 5 * (q : ℝ) * ((Nat.log 2 q : ℝ) + 1) ^ 2 / (P : ℝ) := by
  haveI : NeZero q := ⟨hq.ne'⟩
  have hqR : (0 : ℝ) < q := by exact_mod_cast hq
  have hPR : (0 : ℝ) < (P : ℝ) := by exact_mod_cast hP
  classical
  -- filters
  set sR := (Finset.range q).filter (fun k => k ≠ 0) with hsR
  set sZ := (Finset.univ : Finset (ZMod q)).filter (fun z => z ≠ 0) with hsZ
  -- `w·z ≠ 0` for `z ≠ 0` (else `P ≤ cdist z · 0 = 0`)
  have hwk_ne : ∀ z : ZMod q, z ≠ 0 → (w : ZMod q) * z ≠ 0 := by
    intro z hz hc
    have hh := hPmin z hz
    rw [hc, cdist_zero, Nat.mul_zero] at hh
    omega
  -- `(k : ZMod q) ≠ 0` for `k ∈ range q, k ≠ 0`
  have hkz_ne : ∀ k : ℕ, k < q → k ≠ 0 → ((k : ℕ) : ZMod q) ≠ 0 := by
    intro k hklt hkne hc
    have hv : ((k : ZMod q)).val = k := ZMod.val_natCast_of_lt hklt
    rw [hc, ZMod.val_zero] at hv
    exact hkne hv.symm
  -- cast bridges for the ℤ arguments of `bandDFT`
  have hkzℤ : ∀ k : ℕ, (((k : ℤ)) : ZMod q) = ((k : ℕ) : ZMod q) := fun k => by push_cast; ring
  have hwkzℤ : ∀ k : ℕ, (((w : ℤ) * (k : ℤ)) : ZMod q) = (w : ZMod q) * ((k : ℕ) : ZMod q) := by
    intro k; push_cast; ring
  -- opus's difference bound
  have hdb := completion_diff_bound q hq B hB w
  -- termwise: `‖B̂(k)‖·‖B̂(wk)‖ ≤ (q²/4)·1/(cdist k · cdist(w·k))`
  have hterm : ∀ k ∈ sR,
      ‖bandDFT q B (k : ℤ)‖ * ‖bandDFT q B ((w : ℤ) * (k : ℤ))‖
        ≤ (q : ℝ) ^ 2 / 4 *
          (1 / ((cdist ((k : ℕ) : ZMod q) : ℝ) * (cdist ((w : ZMod q) * ((k : ℕ) : ZMod q)) : ℝ))) := by
    intro k hk
    rw [hsR, Finset.mem_filter, Finset.mem_range] at hk
    obtain ⟨hklt, hkne⟩ := hk
    have hz1 : ((k : ℕ) : ZMod q) ≠ 0 := hkz_ne k hklt hkne
    have hz2 : (w : ZMod q) * ((k : ℕ) : ZMod q) ≠ 0 := hwk_ne _ hz1
    have hb1 := hcoeff (k : ℤ) (by exact_mod_cast hz1)
    have hb2 := hcoeff ((w : ℤ) * (k : ℤ)) (by exact_mod_cast hz2)
    simp only [Int.cast_mul, Int.cast_natCast] at hb1 hb2
    have hc1 : (1 : ℝ) ≤ (cdist ((k : ℕ) : ZMod q) : ℝ) := by exact_mod_cast one_le_cdist hz1
    have hc2 : (1 : ℝ) ≤ (cdist ((w : ZMod q) * ((k : ℕ) : ZMod q)) : ℝ) := by
      exact_mod_cast one_le_cdist hz2
    have hne1 : (cdist ((k : ℕ) : ZMod q) : ℝ) ≠ 0 := by linarith
    have hne2 : (cdist ((w : ZMod q) * ((k : ℕ) : ZMod q)) : ℝ) ≠ 0 := by linarith
    have hbnd1 : (0 : ℝ) ≤ (q : ℝ) / (2 * (cdist ((k : ℕ) : ZMod q) : ℝ)) := by positivity
    calc ‖bandDFT q B (k : ℤ)‖ * ‖bandDFT q B ((w : ℤ) * (k : ℤ))‖
        ≤ ((q : ℝ) / (2 * (cdist ((k : ℕ) : ZMod q) : ℝ)))
            * ((q : ℝ) / (2 * (cdist ((w : ZMod q) * ((k : ℕ) : ZMod q)) : ℝ))) :=
          mul_le_mul hb1 hb2 (norm_nonneg _) hbnd1
      _ = (q : ℝ) ^ 2 / 4 *
            (1 / ((cdist ((k : ℕ) : ZMod q) : ℝ) * (cdist ((w : ZMod q) * ((k : ℕ) : ZMod q)) : ℝ))) := by
          field_simp; ring
  -- the range-`q` harmonic sum
  set Hr : ℝ := ∑ k ∈ sR,
      1 / ((cdist ((k : ℕ) : ZMod q) : ℝ) * (cdist ((w : ZMod q) * ((k : ℕ) : ZMod q)) : ℝ)) with hHr
  -- `(1/q)·Σ ‖·‖‖·‖ ≤ (q/4)·Hr`
  have hsum_le : (1 / (q : ℝ)) * ∑ k ∈ sR,
        ‖bandDFT q B (k : ℤ)‖ * ‖bandDFT q B ((w : ℤ) * (k : ℤ))‖
      ≤ (q : ℝ) / 4 * Hr := by
    have hstep : ∑ k ∈ sR, ‖bandDFT q B (k : ℤ)‖ * ‖bandDFT q B ((w : ℤ) * (k : ℤ))‖
        ≤ ∑ k ∈ sR, (q : ℝ) ^ 2 / 4 *
            (1 / ((cdist ((k : ℕ) : ZMod q) : ℝ) * (cdist ((w : ZMod q) * ((k : ℕ) : ZMod q)) : ℝ))) :=
      Finset.sum_le_sum hterm
    rw [← Finset.mul_sum, ← hHr] at hstep
    calc (1 / (q : ℝ)) * ∑ k ∈ sR, ‖bandDFT q B (k : ℤ)‖ * ‖bandDFT q B ((w : ℤ) * (k : ℤ))‖
        ≤ (1 / (q : ℝ)) * ((q : ℝ) ^ 2 / 4 * Hr) :=
          mul_le_mul_of_nonneg_left hstep (by positivity)
      _ = (q : ℝ) / 4 * Hr := by rw [← mul_assoc]; congr 1; field_simp
  -- bridge: `Hr = S` (death-star's `ZMod q` harmonic sum)
  set S : ℝ := ∑ z ∈ sZ, 1 / ((cdist z : ℝ) * (cdist ((w : ZMod q) * z) : ℝ)) with hS
  have hbridge : Hr = S := by
    rw [hHr, hS]
    refine Finset.sum_bij' (fun k _ => ((k : ℕ) : ZMod q)) (fun z _ => z.val) ?_ ?_ ?_ ?_ ?_
    · intro k hk
      rw [hsR, Finset.mem_filter, Finset.mem_range] at hk
      rw [hsZ, Finset.mem_filter]
      exact ⟨Finset.mem_univ _, hkz_ne k hk.1 hk.2⟩
    · intro z hz
      rw [hsZ, Finset.mem_filter] at hz
      rw [hsR, Finset.mem_filter, Finset.mem_range]
      have hvp := ZMod.val_pos.mpr hz.2
      exact ⟨ZMod.val_lt z, hvp.ne'⟩
    · intro k hk
      rw [hsR, Finset.mem_filter, Finset.mem_range] at hk
      exact ZMod.val_natCast_of_lt hk.1
    · intro z _
      simp
    · intro k _
      rfl
  -- death-star's harmonic aggregation (ℚ → ℝ)
  have hds := LonelyRunner.HyperbolaBox.harmonic_ratio_sum_mul_le (w : ZMod q) P hPmin
  have hdsR : S * (P : ℝ) ≤ 20 * ((Nat.log 2 q : ℝ) + 1) ^ 2 := by
    rw [hS]
    calc (∑ z ∈ sZ, 1 / ((cdist z : ℝ) * (cdist ((w : ZMod q) * z) : ℝ))) * (P : ℝ)
        = (((∑ z ∈ sZ, (1 : ℚ) / (cdist z * cdist ((w : ZMod q) * z))) * (P : ℚ) : ℚ) : ℝ) := by
          push_cast; ring
      _ ≤ ((20 * ((Nat.log 2 q : ℚ) + 1) ^ 2 : ℚ) : ℝ) := by exact_mod_cast hds
      _ = 20 * ((Nat.log 2 q : ℝ) + 1) ^ 2 := by push_cast; ring
  have hSle : S ≤ 20 * ((Nat.log 2 q : ℝ) + 1) ^ 2 / (P : ℝ) := by
    rw [le_div_iff₀ hPR]; linarith [hdsR]
  -- assemble
  calc ‖(corrCount q B w : ℂ) - (B.card : ℂ) ^ 2 / q‖
      ≤ (1 / (q : ℝ)) * ∑ k ∈ sR,
          ‖bandDFT q B (k : ℤ)‖ * ‖bandDFT q B ((w : ℤ) * (k : ℤ))‖ := hdb
    _ ≤ (q : ℝ) / 4 * Hr := hsum_le
    _ = (q : ℝ) / 4 * S := by rw [hbridge]
    _ ≤ (q : ℝ) / 4 * (20 * ((Nat.log 2 q : ℝ) + 1) ^ 2 / (P : ℝ)) :=
        mul_le_mul_of_nonneg_left hSle (by positivity)
    _ = 5 * (q : ℝ) * ((Nat.log 2 q : ℝ) + 1) ^ 2 / (P : ℝ) := by field_simp; ring

end FourierCompletion
end LonelyRunner
