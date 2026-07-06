/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-06-S101)
-/
import TournamentH7.LRC13Citation

/-!
# The unshifted circle-clear floor at l ≤ 11 is a citation corollary (HYP-4336)

kps's `CircleClearFloor (2/25) l` (LRCCircleCover) asks that `l` distinct-frequency combs
of radius `2/25` with ARBITRARY shifts leave a clear point; it is proved for `l ≤ 6` by
density (`2ρl < 1`).  mac-mini-S9 (HYP-4332) located the genuinely open ("Newman") content
in the SHIFTED case: the UNSHIFTED version at `l ≤ 11` is a direct corollary of the
`LRC(≤13)` citation, because `l ≤ 11` distinct nonzero speeds are lonely at margin
`1/(l+1) ≥ 1/12 > 2/25`.

This file discharges that clean part formally, extending the density lane (`l ≤ 6`) to
`l ≤ 11` for the unshifted (`s = 0`) case, so the remaining obligation is precisely the
shifted `7 ≤ l ≤ 11` Mirsky–Newman lemma.
-/

namespace LonelyRunner
namespace CircleUnshifted

/-- **The unshifted floor, l ≤ 11** (citation corollary): any `l ≤ 11` positive
frequencies `r i` (i ∈ S) admit a point `θ` with every `r i · θ` at circle-distance
`≥ 2/25` from `ℤ` — no shifts.  Margin `1/(l+1) ≥ 1/12 > 2/25` from `LRC(≤13)`. -/
theorem circleClear_unshifted_of_le11 (cite : LRCUpTo13)
    (S : Finset (Fin 12)) (r : Fin 12 → ℤ) (hl : S.card ≤ 11)
    (hr : ∀ i ∈ S, 0 < r i) :
    ∃ θ : ℝ, ∀ i ∈ S, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(r i : ℝ) * θ - m| := by
  classical
  -- enumerate S and cite it as an ≤ 12 nonzero family
  let e : Fin S.card ≃ S := S.equivFin.symm
  set w : Fin S.card → ℤ := fun j => r (e j : Fin 12) with hw
  have hwne : ∀ j, w j ≠ 0 := by
    intro j
    have hmem : (e j : Fin 12) ∈ S := (e j).2
    have := hr _ hmem
    rw [hw]; simp only
    exact ne_of_gt this
  obtain ⟨t, hL⟩ := cite S.card (by omega) w hwne
  refine ⟨t, fun i hi m => ?_⟩
  -- transfer the margin at the enumerated index back to i
  set j0 : Fin S.card := e.symm ⟨i, hi⟩ with hj0
  have hwj0 : w j0 = r i := by
    rw [hw, hj0]; simp only [Equiv.apply_symm_apply]
  have hstep : (1 : ℝ) / ((S.card : ℝ) + 1) ≤ |(r i : ℝ) * t - m| := by
    have := hL j0 m
    rw [hwj0] at this
    have hcast : ((S.card + 1 : ℕ) : ℝ) = (S.card : ℝ) + 1 := by push_cast; ring
    rwa [hcast] at this
  -- 1/(l+1) ≥ 1/12 > 2/25 since l ≤ 11
  have hbig : (2 : ℝ) / 25 ≤ 1 / ((S.card : ℝ) + 1) := by
    have hpos : (0 : ℝ) < (S.card : ℝ) + 1 := by positivity
    have hle12 : ((S.card : ℝ) + 1) ≤ 12 := by
      have : (S.card : ℝ) ≤ 11 := by exact_mod_cast hl
      linarith
    calc (2 : ℝ) / 25 ≤ 1 / 12 := by norm_num
      _ ≤ 1 / ((S.card : ℝ) + 1) := one_div_le_one_div_of_le hpos hle12
  linarith

#print axioms circleClear_unshifted_of_le11

end CircleUnshifted
end LonelyRunner
