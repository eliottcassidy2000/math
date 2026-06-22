/-
  TournamentH7.LRCMaxGapPigeonhole  (mac-mini-2026-06-22-S26)

  A self-contained building block for the LRC(14) witness-floor `hnu1` node
  ("k ≤ 7 ⟹ nu = 1", i.e. for ≤ 7 phases the max gap exceeds 1/7 a.e.).

  CORE: the max-gap pigeonhole.  If `k` cyclic gaps are nonnegative and sum to
  `1`, then the largest is `≥ 1/k`.  Hence for `k ≤ 6`, the max gap is
  `≥ 1/6 > 1/7` ALWAYS — the clean (non-a.e.) part of `hnu1`.  We state it as a
  pure averaging fact over `Fin k`, ready to feed a concrete max-gap definition.
-/
import Mathlib

namespace TournamentH7.MaxGapPigeonhole

open Finset

/-- **Averaging / pigeonhole.**  `k` nonnegative reals summing to `1` have one
that is `≥ 1/k`. -/
theorem exists_one_div_card_le {k : ℕ} (hk : 0 < k) (g : Fin k → ℝ)
    (hsum : ∑ i, g i = 1) : ∃ i, (1 : ℝ) / k ≤ g i := by
  by_contra h
  push_neg at h
  have hlt : ∑ i, g i < ∑ _i : Fin k, (1 : ℝ) / k := by
    apply Finset.sum_lt_sum_of_nonempty
    · exact Finset.univ_nonempty_iff.mpr (Fin.pos_iff_nonempty.mp hk)
    · intro i _; exact h i
  rw [hsum, Finset.sum_const, Finset.card_univ, Fintype.card_fin] at hlt
  rw [nsmul_eq_mul, mul_one_div, div_self (by exact_mod_cast hk.ne')] at hlt
  exact lt_irrefl 1 hlt

/-- **Max-gap pigeonhole (the usable form).**  If the cyclic gaps `g : Fin k → ℝ`
are nonnegative, sum to `1`, and `k ≤ 6`, then some gap exceeds `1/7`.  This is
the clean (everywhere, not just a.e.) part of the witness-floor `hnu1` node:
`≤ 6` phases always leave a `> 1/7` runner-free arc. -/
theorem exists_gap_gt_one_seventh {k : ℕ} (hk : 0 < k) (hk6 : k ≤ 6)
    (g : Fin k → ℝ) (hsum : ∑ i, g i = 1) : ∃ i, (1 : ℝ) / 7 < g i := by
  obtain ⟨i, hi⟩ := exists_one_div_card_le hk g hsum
  refine ⟨i, lt_of_lt_of_le ?_ hi⟩
  have hkpos : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
  have hk7 : (k : ℝ) < 7 := by exact_mod_cast Nat.lt_succ_of_le hk6
  exact one_div_lt_one_div_of_lt hkpos hk7

end TournamentH7.MaxGapPigeonhole
