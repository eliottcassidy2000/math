/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S82)
-/
import TournamentH7.LRCKernelGate13

/-!
# Class-level kernel rows for the l ≥ 7 grid stratum (HYP-4119)

The bottom-cluster residual of the corrected l ≥ 7 architecture (HYP-4116) is served by
strict witnesses `t = a/169`.  The key semantics: `speedOK13 s a 169` depends only on
`s mod 169`, so ONE kernel row certifies the ENTIRE congruence class of the family mod 169
— a row per PATTERN `(r, k mod 13)`, not per family.  `speedOK13_congr` is the transport;
the four rows below are the systematic bottom-cluster patterns probed in S81
(`six_top_ceiling_gap_descent_opus_S81.out`, Part 3), including the pure `k = 1` bottom
block `{14..20} ∪ {8..12}`.

The complete row set for the stratum is determined by the shadow-dichotomy search
(`shadow_dichotomy_169_opus_S82.py`): patterns whose kill sets do NOT cover all 156 cells
`a ∈ [1,168] \ 13ℤ` have a strict witness; the witness-free patterns are expected to be
exactly the 156 dilation shadows `λ·{1..12} mod 169` — non-primitive upstream, and any
primitive family on a shadow deviates by ≥ 169 at some position, hence carries a runner
`≥ 170 > 134` = descent-eligible (LRCGapDescent).
-/

namespace LonelyRunner
namespace LiftRowsL7

open KernelGate13

/-- `speedOK13` sees only the speed's residue class mod `den`: the kernel-row transport.
One checked row certifies every family congruent to it. -/
theorem speedOK13_congr {s s' num : ℤ} {den : ℕ}
    (h : s % (den : ℤ) = s' % (den : ℤ)) (hOK : speedOK13 s' num den) :
    speedOK13 s num den := by
  unfold speedOK13 at hOK ⊢
  have hm : (s * num) % (den : ℤ) = (s' * num) % (den : ℤ) :=
    Int.ModEq.mul_right num h
  rw [hm]
  exact hOK

/-- Class-level row schema: a kernel-checked representative certifies its whole
congruence class mod 169 as strictly `1/13`-lonely. -/
theorem strictLonely13_of_row (v v₀ : Fin 12 → ℤ) (a : ℤ)
    (hrow : ∀ i, speedOK13 (v₀ i) a 169)
    (hcong : ∀ i, v i % 169 = v₀ i % 169) :
    ∃ t : ℝ, StrictLonely13 v t :=
  strictLonely13_of_kernelWitness (by norm_num)
    (fun i => speedOK13_congr (hcong i) (hrow i))

/-! ## The four systematic bottom-cluster rows (S81 probe, kernel-checked) -/

/-- Pattern A: the pure `k = 1` bottom block — residues 1..7 lifted once. -/
def rowA : Fin 12 → ℤ := ![14, 15, 16, 17, 18, 19, 20, 8, 9, 10, 11, 12]

theorem rowA_check : ∀ i, speedOK13 (rowA i) 6 169 := by decide

/-- Every family ≡ (14,…,20,8,…,12) mod 169 is strictly loose — witness `t = 6/169`. -/
theorem rowA_class (v : Fin 12 → ℤ) (hcong : ∀ i, v i % 169 = rowA i % 169) :
    ∃ t : ℝ, StrictLonely13 v t :=
  strictLonely13_of_row v rowA 6 rowA_check hcong

/-- Pattern B: residues 6..12 lifted once (`k = 1` on the unique-multiple block). -/
def rowB : Fin 12 → ℤ := ![1, 2, 3, 4, 5, 19, 20, 21, 22, 23, 24, 25]

theorem rowB_check : ∀ i, speedOK13 (rowB i) 19 169 := by decide

theorem rowB_class (v : Fin 12 → ℤ) (hcong : ∀ i, v i % 169 = rowB i % 169) :
    ∃ t : ℝ, StrictLonely13 v t :=
  strictLonely13_of_row v rowB 19 rowB_check hcong

/-- Pattern C: residues 1..6 and 12 lifted once. -/
def rowC : Fin 12 → ℤ := ![14, 15, 16, 17, 18, 19, 7, 8, 9, 10, 11, 25]

theorem rowC_check : ∀ i, speedOK13 (rowC i) 5 169 := by decide

theorem rowC_class (v : Fin 12 → ℤ) (hcong : ∀ i, v i % 169 = rowC i % 169) :
    ∃ t : ℝ, StrictLonely13 v t :=
  strictLonely13_of_row v rowC 5 rowC_check hcong

/-- Pattern D: alternating lift (residues 2,4,6,8,10,11,12 lifted once). -/
def rowD : Fin 12 → ℤ := ![1, 15, 3, 17, 5, 19, 7, 21, 9, 23, 24, 25]

theorem rowD_check : ∀ i, speedOK13 (rowD i) 83 169 := by decide

theorem rowD_class (v : Fin 12 → ℤ) (hcong : ∀ i, v i % 169 = rowD i % 169) :
    ∃ t : ℝ, StrictLonely13 v t :=
  strictLonely13_of_row v rowD 83 rowD_check hcong

#print axioms rowA_class
#print axioms strictLonely13_of_row

end LiftRowsL7
end LonelyRunner
