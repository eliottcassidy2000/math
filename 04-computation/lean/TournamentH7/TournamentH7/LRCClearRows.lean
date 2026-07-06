/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S94)
-/
import TournamentH7.LRCClearCert

/-!
# The end-to-end instantiation demo: the l ≥ 7 chain FIRES (HYP-4226)

One concrete run of the complete formal chain on generated data
(`clear_cert_gen_opus_S94.out`): the |B| = 5 worst base `(3,5,6,8,11)`, its zone-clear
component rationalized to `(3/106, 4/106)`, kernel `toothMiss` certificates for all five
base runners (each an `interval_cases + decide` over a two-element `m`-window), and the
pinned spread chain `[222, 535, 1265, 2997, 7094, 16780, 39662]` (residues
1,2,4,7,9,10,12 mod 13; consecutive ratios ≥ 26/11; entry `2 ≤ 222/106`).  Conclusion:
a single time strictly `1/13`-lonely for all twelve runners.  This is the shape every
production row of the l ≥ 7 table will take.
-/

namespace LonelyRunner
namespace ClearRows

open ClearCert DescentSurface

/-- The five base-runner certificates, kernel-checked over their finite tooth windows. -/
theorem cert_base (w : ℤ) (hw : w = 3 ∨ w = 5 ∨ w = 6 ∨ w = 8 ∨ w = 11) :
    ∀ m : ℤ, w * 3 - 106 ≤ 106 * m → 106 * m ≤ w * (3 + 1) + 106 →
      toothMiss w 3 1 106 m := by
  rcases hw with rfl | rfl | rfl | rfl | rfl <;>
  · intro m h1 h2
    have hlo : 0 ≤ m := by omega
    have hhi : m ≤ 1 := by omega
    interval_cases m <;> decide

/-- The base is strictly clear on `(3/106, 4/106)`, pointwise. -/
theorem demo_base_clear :
    ∀ t : ℝ, (3 : ℝ) / 106 < t → t < (3 : ℝ) / 106 + 1 / 106 →
      ∀ w ∈ ([3, 5, 6, 8, 11] : List ℝ), ∀ m : ℤ, (1 : ℝ) / 13 < |w * t - m| := by
  intro t ht1 ht2 w hw m
  have hglue : ∀ (wz : ℤ), wz = 3 ∨ wz = 5 ∨ wz = 6 ∨ wz = 8 ∨ wz = 11 →
      (1 : ℝ) / 13 < |(wz : ℝ) * t - m| := by
    intro wz hwz
    apply clear_of_cert wz 3 1 106 (by rcases hwz with rfl|rfl|rfl|rfl|rfl <;> norm_num)
      (by norm_num) (by norm_num) (cert_base wz hwz) t
    · push_cast; linarith
    · push_cast; linarith
  fin_cases hw
  · exact_mod_cast hglue 3 (by norm_num)
  · exact_mod_cast hglue 5 (by norm_num)
  · exact_mod_cast hglue 6 (by norm_num)
  · exact_mod_cast hglue 8 (by norm_num)
  · exact_mod_cast hglue 11 (by norm_num)

/-- **The demo**: the full chain fires — twelve concrete runners, one strictly
`1/13`-lonely time. -/
theorem demo_strict_lonely :
    ∃ t : ℝ,
      (∀ w ∈ ([3, 5, 6, 8, 11] : List ℝ), ∀ m : ℤ, (1 : ℝ) / 13 < |w * t - m|) ∧
      (∀ w ∈ ([222, 535, 1265, 2997, 7094, 16780, 39662] : List ℝ), ∀ m : ℤ,
        (1 : ℝ) / 13 < |w * t - m|) := by
  apply strict_lonely_of_clear_component (1/13) (by norm_num) (by norm_num)
    ([3, 5, 6, 8, 11] : List ℝ) ([222, 535, 1265, 2997, 7094, 16780, 39662] : List ℝ)
    ((3 : ℝ) / 106) (1 / 106) (by norm_num) demo_base_clear
  · intro w hw
    fin_cases hw <;> norm_num
  · repeat rw [List.isChain_cons]
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, List.isChain_nil⟩ <;>
      (intro y hy; cases hy <;> norm_num)
  · intro w hw
    cases hw <;> norm_num

#print axioms demo_strict_lonely

end ClearRows
end LonelyRunner
