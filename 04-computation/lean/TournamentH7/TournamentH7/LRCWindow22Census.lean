/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S202)
-/
import Mathlib
import TournamentH7.LRCKernelGate
import TournamentH7.LRC14CertRoute
import TournamentH7.LRCSpread13
import TournamentH7.LRC14WindowWiring

/-!
# The ≤22 covering census as a 6-witness pigeonhole (LEM-024) — removing the `winData22` native_decide

Every census-branch tuple (13 distinct integers in `[1,22]`, covering, containing `1`) is lonely at one of
the SIX witnesses `{12/25, 9/26, 7/27, 11/28, 4/23, 11/26}`, whose danger sets over `[1,22]` are
`{2},{3},{4},{5},{6,17},{7,19}`. Failing all six forces `{1,2,3,4,5}` ∪ (one of {6,17}) ∪ (one of {7,19})
∪ (covering: {12,13,14} ∪ one of each of {8,16},{9,18},{10,20},{11,22}) = **14 distinct elements** in a
13-set — impossible. So one witness works; `KernelGate.lonely_of_kernelWitness` finishes.

This replaces the `C(22,13) = 497 420` native_decide census (kernel `decide` infeasible, MISTAKE-135) with
a fixed witness list + a pigeonhole — all facts small and kernel-pure. See LEM-024.

Kernel-pure: no `sorry`, no `native_decide`. Axioms target: `[propext, Classical.choice, Quot.sound]`.
-/

namespace LonelyRunner
namespace Window22Census

open KernelGate

/-! ## The six far-set facts: over `[1,22]`, a speed that is NOT far at witness `w` lies in `D(w)`. -/

theorem danger_12_25 (s : ℤ) (h1 : 1 ≤ s) (h2 : s ≤ 22) (h : ¬ speedOK s 12 25) : s = 2 := by
  interval_cases s <;> revert h <;> decide

theorem danger_9_26 (s : ℤ) (h1 : 1 ≤ s) (h2 : s ≤ 22) (h : ¬ speedOK s 9 26) : s = 3 := by
  interval_cases s <;> revert h <;> decide

theorem danger_7_27 (s : ℤ) (h1 : 1 ≤ s) (h2 : s ≤ 22) (h : ¬ speedOK s 7 27) : s = 4 := by
  interval_cases s <;> revert h <;> decide

theorem danger_11_28 (s : ℤ) (h1 : 1 ≤ s) (h2 : s ≤ 22) (h : ¬ speedOK s 11 28) : s = 5 := by
  interval_cases s <;> revert h <;> decide

theorem danger_4_23 (s : ℤ) (h1 : 1 ≤ s) (h2 : s ≤ 22) (h : ¬ speedOK s 4 23) : s = 6 ∨ s = 17 := by
  interval_cases s <;> revert h <;> decide

theorem danger_11_26 (s : ℤ) (h1 : 1 ≤ s) (h2 : s ≤ 22) (h : ¬ speedOK s 11 26) : s = 7 ∨ s = 19 := by
  interval_cases s <;> revert h <;> decide

/-! ## Covering forces the multiples of `q` in `[1,22]`. -/

theorem cov_12 (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) (h22 : ∀ i, v i ≤ 22)
    (hcov : LRC14.CoveringFamily v) : ∃ i, v i = 12 := by
  obtain ⟨i, hi⟩ := hcov 12 (by norm_num) (by norm_num)
  have hp := hpos i; have hb := h22 i; have hd : (12 : ℤ) ∣ v i := by exact_mod_cast hi
  exact ⟨i, by omega⟩

theorem cov_13 (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) (h22 : ∀ i, v i ≤ 22)
    (hcov : LRC14.CoveringFamily v) : ∃ i, v i = 13 := by
  obtain ⟨i, hi⟩ := hcov 13 (by norm_num) (by norm_num)
  have hp := hpos i; have hb := h22 i; have hd : (13 : ℤ) ∣ v i := by exact_mod_cast hi
  exact ⟨i, by omega⟩

theorem cov_14 (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) (h22 : ∀ i, v i ≤ 22)
    (hcov : LRC14.CoveringFamily v) : ∃ i, v i = 14 := by
  obtain ⟨i, hi⟩ := hcov 14 (by norm_num) (by norm_num)
  have hp := hpos i; have hb := h22 i; have hd : (14 : ℤ) ∣ v i := by exact_mod_cast hi
  exact ⟨i, by omega⟩

theorem cov_8 (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) (h22 : ∀ i, v i ≤ 22)
    (hcov : LRC14.CoveringFamily v) : ∃ i, v i = 8 ∨ v i = 16 := by
  obtain ⟨i, hi⟩ := hcov 8 (by norm_num) (by norm_num)
  have hp := hpos i; have hb := h22 i; have hd : (8 : ℤ) ∣ v i := by exact_mod_cast hi
  exact ⟨i, by omega⟩

theorem cov_9 (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) (h22 : ∀ i, v i ≤ 22)
    (hcov : LRC14.CoveringFamily v) : ∃ i, v i = 9 ∨ v i = 18 := by
  obtain ⟨i, hi⟩ := hcov 9 (by norm_num) (by norm_num)
  have hp := hpos i; have hb := h22 i; have hd : (9 : ℤ) ∣ v i := by exact_mod_cast hi
  exact ⟨i, by omega⟩

theorem cov_10 (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) (h22 : ∀ i, v i ≤ 22)
    (hcov : LRC14.CoveringFamily v) : ∃ i, v i = 10 ∨ v i = 20 := by
  obtain ⟨i, hi⟩ := hcov 10 (by norm_num) (by norm_num)
  have hp := hpos i; have hb := h22 i; have hd : (10 : ℤ) ∣ v i := by exact_mod_cast hi
  exact ⟨i, by omega⟩

theorem cov_11 (v : Fin 13 → ℤ) (hpos : ∀ i, 0 < v i) (h22 : ∀ i, v i ≤ 22)
    (hcov : LRC14.CoveringFamily v) : ∃ i, v i = 11 ∨ v i = 22 := by
  obtain ⟨i, hi⟩ := hcov 11 (by norm_num) (by norm_num)
  have hp := hpos i; have hb := h22 i; have hd : (11 : ℤ) ∣ v i := by exact_mod_cast hi
  exact ⟨i, by omega⟩

/-! ## The pigeonhole: every census-branch tuple is far at one of the six witnesses. -/

/-- **LEM-024 (min=1 case), kernel-pure.** A covering family of 13 distinct speeds in `[1,22]` containing
`1` is lonely — it is far at one of the six witnesses `{12/25, 9/26, 7/27, 11/28, 4/23, 11/26}` (the
14-forced-elements pigeonhole), and `lonely_of_kernelWitness` finishes. No enumeration, no `native_decide`. -/
theorem window22_min1_lonely (v : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < v i) (h22 : ∀ i, v i ≤ 22) (hinj : Function.Injective v)
    (hcov : LRC14.CoveringFamily v) (hmin1 : ∃ i, v i = 1) :
    ∃ t : ℝ, Lonely 14 v t := by
  suffices hdisj : (∀ i, speedOK (v i) 12 25) ∨ (∀ i, speedOK (v i) 9 26) ∨
      (∀ i, speedOK (v i) 7 27) ∨ (∀ i, speedOK (v i) 11 28) ∨
      (∀ i, speedOK (v i) 4 23) ∨ (∀ i, speedOK (v i) 11 26) by
    rcases hdisj with h | h | h | h | h | h
    · exact ⟨_, lonely_of_kernelWitness (by norm_num) h⟩
    · exact ⟨_, lonely_of_kernelWitness (by norm_num) h⟩
    · exact ⟨_, lonely_of_kernelWitness (by norm_num) h⟩
    · exact ⟨_, lonely_of_kernelWitness (by norm_num) h⟩
    · exact ⟨_, lonely_of_kernelWitness (by norm_num) h⟩
    · exact ⟨_, lonely_of_kernelWitness (by norm_num) h⟩
  by_contra hcon
  push_neg at hcon
  obtain ⟨⟨j2, hj2⟩, ⟨j3, hj3⟩, ⟨j4, hj4⟩, ⟨j5, hj5⟩, ⟨j6, hj6⟩, ⟨j7, hj7⟩⟩ := hcon
  have e2 : v j2 = 2 := danger_12_25 _ (hpos j2) (h22 j2) hj2
  have e3 : v j3 = 3 := danger_9_26 _ (hpos j3) (h22 j3) hj3
  have e4 : v j4 = 4 := danger_7_27 _ (hpos j4) (h22 j4) hj4
  have e5 : v j5 = 5 := danger_11_28 _ (hpos j5) (h22 j5) hj5
  have e6 : v j6 = 6 ∨ v j6 = 17 := danger_4_23 _ (hpos j6) (h22 j6) hj6
  have e7 : v j7 = 7 ∨ v j7 = 19 := danger_11_26 _ (hpos j7) (h22 j7) hj7
  obtain ⟨j12, e12⟩ := cov_12 v hpos h22 hcov
  obtain ⟨j13, e13⟩ := cov_13 v hpos h22 hcov
  obtain ⟨j14, e14⟩ := cov_14 v hpos h22 hcov
  obtain ⟨j8, e8⟩ := cov_8 v hpos h22 hcov
  obtain ⟨j9, e9⟩ := cov_9 v hpos h22 hcov
  obtain ⟨j10, e10⟩ := cov_10 v hpos h22 hcov
  obtain ⟨j11, e11⟩ := cov_11 v hpos h22 hcov
  obtain ⟨j1, e1⟩ := hmin1
  have hmem : ∀ x : ℤ, (∃ i, v i = x) → x ∈ Finset.image v Finset.univ :=
    fun x ⟨i, hi⟩ => Finset.mem_image.mpr ⟨i, Finset.mem_univ i, hi⟩
  have hR13 : (Finset.image v Finset.univ).card ≤ 13 :=
    le_trans Finset.card_image_le (by simp)
  set F : Finset ℤ := {1, 2, 3, 4, 5, 12, 13, 14, v j6, v j7, v j8, v j9, v j10, v j11} with hFdef
  have hFsub : F ⊆ Finset.image v Finset.univ := by
    rw [hFdef]
    intro x hx
    simp only [Finset.mem_insert, Finset.mem_singleton] at hx
    rcases hx with rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl|rfl
    · exact hmem 1 ⟨j1, e1⟩
    · exact hmem 2 ⟨j2, e2⟩
    · exact hmem 3 ⟨j3, e3⟩
    · exact hmem 4 ⟨j4, e4⟩
    · exact hmem 5 ⟨j5, e5⟩
    · exact hmem 12 ⟨j12, e12⟩
    · exact hmem 13 ⟨j13, e13⟩
    · exact hmem 14 ⟨j14, e14⟩
    · exact hmem _ ⟨j6, rfl⟩
    · exact hmem _ ⟨j7, rfl⟩
    · exact hmem _ ⟨j8, rfl⟩
    · exact hmem _ ⟨j9, rfl⟩
    · exact hmem _ ⟨j10, rfl⟩
    · exact hmem _ ⟨j11, rfl⟩
  have hFcard : F.card = 14 := by
    rw [hFdef]
    rcases e6 with h6 | h6 <;> rcases e7 with h7 | h7 <;> rcases e8 with h8 | h8 <;>
      rcases e9 with h9 | h9 <;> rcases e10 with h10 | h10 <;> rcases e11 with h11 | h11 <;>
      rw [h6, h7, h8, h9, h10, h11] <;> decide
  have := Finset.card_le_card hFsub
  omega

/-- **LEM-024, full census replacement (kernel-pure).** EVERY covering family of 13 distinct positive
speeds bounded by 22 is lonely — `min = 1` by the six-witness pigeonhole (`window22_min1_lonely`), `min ≥ 2`
by `spread13_lonely` (ratio `≤ 22/2 = 11 ≤ 13`). This is a foundational-axioms-only replacement for the
`winData22` native_decide census (it even drops the `gcd = 1` hypothesis `hdistinct22_from_data` carried).
-/
theorem window22_lonely (v : Fin 13 → ℤ)
    (hpos : ∀ i, 0 < v i) (h22 : ∀ i, v i ≤ 22) (hinj : Function.Injective v)
    (hcov : LRC14.CoveringFamily v) :
    ∃ t : ℝ, Lonely 14 v t := by
  by_cases h1 : ∃ i, v i = 1
  · exact window22_min1_lonely v hpos h22 hinj hcov h1
  · push_neg at h1
    refine ⟨1 / ((2 : ℤ) + 22), LRC14.spread13_lonely v 2 22 (by norm_num) (fun i => ?_) (fun i => ?_)
      (by norm_num)⟩
    · have h0 := hpos i; have hn := h1 i; rw [abs_of_pos h0]; omega
    · rw [abs_of_pos (hpos i)]; exact h22 i

/-- **Kernel-pure drop-in for `WindowData.hdistinct22_from_data`.** Identical signature (it even ignores
the `gcd = 1` hypothesis), proved by `window22_lonely`. Consequently `WindowData.hwindowW_closed 22 cite
hdistinct22_kernel` is a foundational-axioms-only `hwindow22_closed`, and replacing the `hwindow22_closed
cite` call in `lrc14_grand_assembly` (LRC14GrandAssembly.lean, the `bounded window: all |v i| ≤ 22`
branch) with it removes the two `winData22` native_decide axioms from the LRC(14) top theorem — the last
non-foundational axioms on the with-census route (opus-S201's "option (c)"). -/
theorem hdistinct22_kernel : ∀ u : Fin 13 → ℤ, (∀ i, 0 < u i) → StrictMono u →
    (∀ i, u i ≤ 22) → LRC14.CoveringFamily u → LRC14.tupleGcd u = 1 → ∃ t : ℝ, Lonely 14 u t :=
  fun u hpos hsm h22 hcov _ => window22_lonely u hpos h22 hsm.injective hcov

end Window22Census
end LonelyRunner
