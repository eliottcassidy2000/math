/-
  TournamentH7.LRCMaxGapPigeonhole  (mac-mini-2026-06-22-S26)

  A self-contained building block for the LRC(14) witness-floor `hnu1` node
  ("k ≤ 7 ⟹ nu = 1", i.e. for ≤ 7 phases the max gap exceeds 1/7 a.e.).

  CORE: the max-gap pigeonhole.  If `k` cyclic gaps are nonnegative and sum to
  `1`, then the largest is `≥ 1/k`.  Hence for `k ≤ 6`, the max gap is
  `≥ 1/6 > 1/7` ALWAYS — the clean (non-a.e.) part of `hnu1`.  We state it as a
  pure averaging fact over `Fin k`, ready to feed a concrete max-gap definition.
-/
import Mathlib.Tactic

namespace TournamentH7.MaxGapPigeonhole

open Finset

/-- **Averaging / pigeonhole.**  `k` nonnegative reals summing to `1` have one
that is `≥ 1/k`. -/
theorem exists_one_div_card_le {k : ℕ} (hk : 0 < k) (g : Fin k → ℝ)
    (hsum : ∑ i, g i = 1) : ∃ i, (1 : ℝ) / k ≤ g i := by
  by_contra h
  push Not at h
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

/-- **Seven-gap equality boundary.**  If seven cyclic gaps sum to `1` and every
gap is at most `1/7`, then every gap is exactly `1/7`.  This isolates the
remaining `k = 7` `hnu1` obstruction to the equal-spacing event; the strict
`> 1/7` event follows everywhere off this boundary. -/
theorem all_eq_one_seventh_of_le (g : Fin 7 → ℝ)
    (hsum : ∑ i, g i = 1)
    (hle : ∀ i, g i ≤ (1 : ℝ) / 7) :
    ∀ i, g i = (1 : ℝ) / 7 := by
  intro i
  by_contra hne
  have hlt_i : g i < (1 : ℝ) / 7 := lt_of_le_of_ne (hle i) hne
  have hsum_lt : ∑ j, g j < ∑ _j : Fin 7, (1 : ℝ) / 7 := by
    exact Finset.sum_lt_sum (s := Finset.univ) (fun j _ => hle j)
      ⟨i, Finset.mem_univ i, hlt_i⟩
  rw [hsum, Finset.sum_const, Finset.card_univ, Fintype.card_fin] at hsum_lt
  norm_num at hsum_lt

/-- **Seven-gap strict-or-boundary split.**  Seven cyclic gaps summing to `1`
either contain a gap strictly larger than `1/7`, or all seven gaps are exactly
`1/7`.  This is the finite algebraic part of the a.e. `k = 7` `hnu1` node. -/
theorem exists_gap_gt_or_all_eq_one_seventh (g : Fin 7 → ℝ)
    (hsum : ∑ i, g i = 1) :
    (∃ i, (1 : ℝ) / 7 < g i) ∨ ∀ i, g i = (1 : ℝ) / 7 := by
  by_cases hle : ∀ i, g i ≤ (1 : ℝ) / 7
  · exact Or.inr (all_eq_one_seventh_of_le g hsum hle)
  · push Not at hle
    obtain ⟨i, hi⟩ := hle
    exact Or.inl ⟨i, hi⟩

/-! ## Axiom audit -/

#print axioms exists_one_div_card_le
#print axioms exists_gap_gt_one_seventh
#print axioms all_eq_one_seventh_of_le
#print axioms exists_gap_gt_or_all_eq_one_seventh

end TournamentH7.MaxGapPigeonhole
