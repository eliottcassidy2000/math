/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-05-S84)
-/
import TournamentH7.LRCLiftRowsL7

/-!
# The tower lift: kernel rows migrate up the 13-adic tower for free (HYP-4136)

The tower's one-way traffic (HYP-4126, verified numerically S83) at the GATE level:
`speedOK13 s num den → speedOK13 s (13·num) (13·den)`.  A strict witness `num/den` is the
same rational as `13·num/(13·den)`, and the strict-gate inequality survives the rescaling
exactly — so every kernel row at denominator `13^ℓ` is automatically a row at `13^{ℓ+1}`:
covers project DOWN the tower, witnesses lift UP.  Consequently the S77/S82 row corpus at
169 feeds any future level-3 (2197) table without recomputation, and the tower-limit
dichotomy (HYP-4126: only 13-adic dilations cover at every level) is the right asymptotic
form of the shadow question.

`strictLonely13_of_kernelWitness` consumes lifted rows unchanged (it is
denominator-generic), so no new consumer is needed.
-/

namespace LonelyRunner
namespace TowerLift

open KernelGate13

/-- **The gate-level tower lift**: the strict `1/13`-gate at witness `num/den` transports
verbatim to `13·num/(13·den)` — rows migrate up the 13-adic tower for free. -/
theorem speedOK13_lift {s num : ℤ} {den : ℕ} (h : speedOK13 s num den) :
    speedOK13 s (13 * num) (13 * den) := by
  unfold speedOK13 at h ⊢
  have hmod : (s * (13 * num)) % ((13 * den : ℕ) : ℤ) = 13 * ((s * num) % den) := by
    push_cast
    rw [show s * (13 * num) = 13 * (s * num) by ring]
    exact Int.mul_emod_mul_of_pos _ _ (by norm_num)
  rw [hmod]
  push_cast
  omega

/-- Demo: the S82 pattern-A row at 169 is a row at 2197 with no new kernel check. -/
theorem rowA_check_2197 : ∀ i, speedOK13 (LiftRowsL7.rowA i) (13 * 6) (13 * 169) :=
  fun i => speedOK13_lift (LiftRowsL7.rowA_check i)

/-- **The strict gap lemma** (assembly item (iii): hdich needs `> 1/13`, and interior
points of a full inter-tooth gap deliver it -- strictness is free in the nested tower,
since a nested subinterval's interior is strictly inside every enclosing gap). -/
theorem gap_interior_strict (ρ w a L : ℝ) (hρ0 : 0 < ρ) (hρ : ρ < 1/2)
    (hw : 0 < w) (hL : 2 ≤ w * L) :
    ∃ c d : ℝ, a ≤ c ∧ d ≤ a + L ∧ c < d ∧ d - c = (1 - 2*ρ)/w ∧
      ∀ t, c < t → t < d → ∀ m : ℤ, ρ < |w * t - m| := by
  have hk1 : w * a + ρ ≤ ((⌈w * a + ρ⌉ : ℤ) : ℝ) := Int.le_ceil _
  have hk2 : ((⌈w * a + ρ⌉ : ℤ) : ℝ) < w * a + ρ + 1 := Int.ceil_lt_add_one _
  refine ⟨(((⌈w * a + ρ⌉ : ℤ) : ℝ) + ρ) / w, (((⌈w * a + ρ⌉ : ℤ) : ℝ) + 1 - ρ) / w,
    ?_, ?_, ?_, ?_, ?_⟩
  · rw [le_div_iff₀ hw]
    have hmul : a * w = w * a := mul_comm a w
    linarith
  · rw [div_le_iff₀ hw]
    have hmul : (a + L) * w = w * a + w * L := by ring
    linarith
  · rw [← sub_pos, div_sub_div_same]
    apply div_pos ?_ hw
    linarith
  · rw [div_sub_div_same]
    congr 1
    ring
  · intro t htc htd m
    rw [div_lt_iff₀ hw] at htc
    rw [lt_div_iff₀ hw] at htd
    have hcomm : w * t = t * w := mul_comm w t
    by_cases hm : m ≤ ⌈w * a + ρ⌉
    · have hmR : (m : ℝ) ≤ ((⌈w * a + ρ⌉ : ℤ) : ℝ) := by exact_mod_cast hm
      exact lt_abs.mpr (Or.inl (by linarith))
    · have hm1 : (⌈w * a + ρ⌉ : ℤ) + 1 ≤ m := by omega
      have hmR : ((⌈w * a + ρ⌉ : ℤ) : ℝ) + 1 ≤ (m : ℝ) := by exact_mod_cast hm1
      exact lt_abs.mpr (Or.inr (by linarith))

#print axioms speedOK13_lift
#print axioms gap_interior_strict

end TowerLift
end LonelyRunner
