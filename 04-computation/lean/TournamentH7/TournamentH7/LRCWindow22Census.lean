/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-09-S202)
-/
import Mathlib
import TournamentH7.LRCKernelGate
import TournamentH7.LRC14CertRoute

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

end Window22Census
end LonelyRunner
