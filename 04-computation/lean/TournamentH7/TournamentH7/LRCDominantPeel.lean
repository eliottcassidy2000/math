/-
  TournamentH7.LRCDominantPeel — THE SHARP DOMINANT-RUNNER PEEL
  (kind-pasteur-2026-07-04-S5, HYP-4086).

  Closes the DOMINANT case of the covering dispatch (`covering_lonely_of_dominant_or_compressed`,
  opus LRCHlargeRoute) at the SHARP LINEAR threshold `v_i > 13·(max of the others)` — replacing the
  corpus far-peel's QUADRATIC threshold (`far_peel_lonely_of_cite`: `w > ~(ΣB)²·267`, the V² artifact).
  So EVERY dominant covering family is lonely, not just the quadratically-far tail.

  Argument (LRC(13) + Lipschitz + a covering step, no measure theory):
   1. The 12 runners `≠ i` are 13-lonely at some `tstar` (`M(base) ≥ 1/13`, LRC(13)).
   2. Reverse-triangle: since each `|v_j| ≤ B`, the base stays `≥ 1/14` for `|t − tstar| ≤ 1/(182B)`
      (because `1/13 − B·(1/(182B)) = 1/14`).
   3. COVERING (`far_safe_point`): one of `t = a` or `t = a + 1/(7 v_i)` has `‖v_i t‖ ≥ 1/14`; both lie
      in the base-good interval exactly when `v_i > 13B` (danger gap `1/(7v_i) < 1/(91B)`).
   4. At that `t` all 13 runners are `≥ 1/14`.  Lonely.

  Kernel-pure (propext, Classical.choice, Quot.sound); no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRC13Citation

namespace LonelyRunner
namespace LRC14

/-- **THE COVERING STEP**: for any real `y`, one of `y` or `y + 1/7` is `≥ 1/14` from every integer.
(The danger gap of `‖·‖` has width `1/7`; shifting by `1/7` escapes it.) -/
lemma far_safe_point (y : ℝ) :
    ∃ s : ℝ, (s = y ∨ s = y + 1 / 7) ∧ ∀ m : ℤ, (1 : ℝ) / 14 ≤ |s - m| := by
  by_cases h : ∀ m : ℤ, (1 : ℝ) / 14 ≤ |y - m|
  · exact ⟨y, Or.inl rfl, h⟩
  · push_neg at h
    obtain ⟨m0, hm0⟩ := h
    rw [abs_lt] at hm0
    refine ⟨y + 1 / 7, Or.inr rfl, fun m => ?_⟩
    rcases le_or_gt m m0 with hle | hlt
    · have hmm : (m : ℝ) ≤ (m0 : ℝ) := by exact_mod_cast hle
      rw [abs_of_nonneg (by linarith [hm0.1, hmm])]
      linarith [hm0.1, hmm]
    · have hmm : (m0 : ℝ) + 1 ≤ (m : ℝ) := by exact_mod_cast (by omega : m0 + 1 ≤ m)
      rw [abs_of_nonpos (by linarith [hm0.2, hmm])]
      linarith [hm0.2, hmm]

/-- **THE SHARP DOMINANT PEEL**: a positive 13-family with a runner `i` exceeding `13×` every other,
whose 12 other runners are 13-lonely at `tstar` (LRC(13)), is `Lonely 14`.  Closes `hdom` at the sharp
linear threshold. -/
theorem dominant_lonely (V : Fin 13 → ℤ) (i : Fin 13) (B : ℤ) (tstar : ℝ)
    (hBpos : 0 < B) (hV : ∀ j, 0 < V j) (hB : ∀ j, j ≠ i → V j ≤ B)
    (hdom : 13 * B < V i)
    (hbase : ∀ j, j ≠ i → ∀ m : ℤ, (1 : ℝ) / 13 ≤ |(V j : ℝ) * tstar - m|) :
    ∃ t : ℝ, Lonely 14 V t := by
  have hBR : (0 : ℝ) < B := by exact_mod_cast hBpos
  have hB0 : (B : ℝ) ≠ 0 := hBR.ne'
  have hwpos : (0 : ℝ) < (V i : ℝ) := by exact_mod_cast hV i
  have hwne : (V i : ℝ) ≠ 0 := ne_of_gt hwpos
  have hwR : (13 : ℝ) * B < (V i : ℝ) := by exact_mod_cast hdom
  have h7w : (0 : ℝ) < 7 * (V i : ℝ) := by linarith
  have hp182 : (0 : ℝ) < 1 / (182 * B) := div_pos one_pos (by linarith)
  obtain ⟨s, hs_or, hs_safe⟩ := far_safe_point ((V i : ℝ) * (tstar - 1 / (182 * B)))
  refine ⟨s / (V i : ℝ), ?_⟩
  have hwt : (V i : ℝ) * (s / (V i : ℝ)) = s := by field_simp
  have ht_cases : s / (V i : ℝ) = tstar - 1 / (182 * B) ∨
      s / (V i : ℝ) = (tstar - 1 / (182 * B)) + 1 / (7 * (V i : ℝ)) := by
    rcases hs_or with h | h
    · left; rw [h]; field_simp
    · right; rw [h]; field_simp; try ring
  have hdist : |s / (V i : ℝ) - tstar| ≤ 1 / (182 * B) := by
    rcases ht_cases with h | h
    · rw [h, abs_le]; exact ⟨by linarith [hp182], by linarith [hp182]⟩
    · rw [h]
      have e : tstar - 1 / (182 * B) + 1 / (7 * (V i : ℝ)) - tstar
          = 1 / (7 * (V i : ℝ)) - 1 / (182 * B) := by ring
      rw [e, abs_le]
      have hpos7w : (0 : ℝ) < 1 / (7 * (V i : ℝ)) := div_pos one_pos h7w
      have hle : 1 / (7 * (V i : ℝ)) ≤ 1 / (91 * B) :=
        one_div_le_one_div_of_le (by linarith) (by linarith [hwR])
      have h91 : (1 : ℝ) / (91 * B) = 2 * (1 / (182 * B)) := by
        rw [mul_one_div, div_eq_div_iff (by positivity) (by positivity)]; ring
      exact ⟨by linarith [hpos7w], by linarith [hle, h91]⟩
  have h14 : (1 : ℝ) / (14 : ℕ) = 1 / 14 := by norm_num
  intro idx m
  by_cases hidx : idx = i
  · rw [h14, hidx]
    have key : (V i : ℝ) * (s / (V i : ℝ)) - (m : ℝ) = s - m := by rw [hwt]
    rw [key]; exact hs_safe m
  · have hbi := hbase idx hidx m
    have hBi : (V idx : ℝ) ≤ B := by exact_mod_cast hB idx hidx
    have hVpos : (0 : ℝ) < (V idx : ℝ) := by exact_mod_cast hV idx
    have hstep : |(V idx : ℝ) * tstar - m| - |(V idx : ℝ) * (s / (V i : ℝ)) - m|
        ≤ |(V idx : ℝ)| * |s / (V i : ℝ) - tstar| := by
      calc |(V idx : ℝ) * tstar - m| - |(V idx : ℝ) * (s / (V i : ℝ)) - m|
          ≤ |((V idx : ℝ) * tstar - m) - ((V idx : ℝ) * (s / (V i : ℝ)) - m)| :=
            abs_sub_abs_le_abs_sub _ _
        _ = |(V idx : ℝ)| * |s / (V i : ℝ) - tstar| := by
            rw [show ((V idx : ℝ) * tstar - m) - ((V idx : ℝ) * (s / (V i : ℝ)) - m)
                  = (V idx : ℝ) * (tstar - s / (V i : ℝ)) from by ring,
              abs_mul, abs_sub_comm tstar (s / (V i : ℝ))]
    have hVabs : |(V idx : ℝ)| = (V idx : ℝ) := abs_of_pos hVpos
    have hprod : |(V idx : ℝ)| * |s / (V i : ℝ) - tstar| ≤ 1 / 182 := by
      have hle2 : |(V idx : ℝ)| * |s / (V i : ℝ) - tstar| ≤ (B : ℝ) * (1 / (182 * B)) := by
        rw [hVabs]; exact mul_le_mul hBi hdist (abs_nonneg _) (le_of_lt hBR)
      have he : (B : ℝ) * (1 / (182 * B)) = 1 / 182 := by field_simp
      linarith [hle2, he.le, he.ge]
    rw [h14]
    linarith [hbi, hstep, hprod]

/-- **hdom CLOSED from the LRC(13) citation.**  Any positive covering family with a runner `i`
that exceeds `13×` every other runner is `Lonely 14`.  Self-contained: the 12 other runners are
`13`-lonely by `LRCUpTo13` (`k = 12`), and `dominant_lonely` peels `i` at the sharp threshold.
This discharges the `hdom` obligation of `covering_lonely_of_dominant_or_compressed` — the WHOLE
dominant branch, at the linear threshold, no quadratic blow-up. -/
theorem hdom_closed (cite : LonelyRunner.LRCUpTo13) (V : Fin 13 → ℤ) (hV : ∀ j, 0 < V j)
    (i : Fin 13) (hdom : ∀ j, j ≠ i → 13 * V j < V i) :
    ∃ t : ℝ, Lonely 14 V t := by
  obtain ⟨tstar, htstar⟩ :=
    cite 12 (by norm_num) (fun idx => V (i.succAbove idx)) (fun idx => (hV _).ne')
  have hbase : ∀ j, j ≠ i → ∀ m : ℤ, (1 : ℝ) / 13 ≤ |(V j : ℝ) * tstar - m| := by
    intro j hj m
    obtain ⟨idx, hidx⟩ := Fin.exists_succAbove_eq hj
    simpa only [hidx] using htstar idx m
  obtain ⟨j0, hj0mem, hj0max⟩ := Finset.exists_max_image
    (Finset.univ.filter (fun j => j ≠ i)) V
    (by obtain ⟨j, hj⟩ := exists_ne i;
        exact ⟨j, Finset.mem_filter.mpr ⟨Finset.mem_univ j, hj⟩⟩)
  have hj0ne : j0 ≠ i := (Finset.mem_filter.mp hj0mem).2
  refine dominant_lonely V i (V j0) tstar (hV j0) hV ?_ (hdom j0 hj0ne) hbase
  intro j hj
  exact hj0max j (Finset.mem_filter.mpr ⟨Finset.mem_univ j, hj⟩)

/-- Loneliness is invariant under negating a speed (`‖−y‖ = ‖y‖` on the circle). -/
lemma lonely_neg_arg {c y : ℝ} (h : ∀ m : ℤ, c ≤ |y - m|) (m : ℤ) : c ≤ |(-y) - m| := by
  have hh := h (-m)
  rwa [show (-y) - (m : ℝ) = -(y - ((-m : ℤ) : ℝ)) by push_cast; ring, abs_neg]

/-- **hdom CLOSED, `|·|`-form.**  The exact `hdom` obligation of
`covering_lonely_of_dominant_or_compressed`: any nonzero 13-family with a runner exceeding `13×`
every other in absolute value is `Lonely 14` (from `LRCUpTo13`, sharp linear threshold). -/
theorem hdom_closed_abs (cite : LonelyRunner.LRCUpTo13) (V : Fin 13 → ℤ) (hV : ∀ j, V j ≠ 0)
    (i : Fin 13) (hdom : ∀ j, j ≠ i → 13 * |V j| < |V i|) :
    ∃ t : ℝ, Lonely 14 V t := by
  obtain ⟨t, ht⟩ := hdom_closed cite (fun j => |V j|) (fun j => abs_pos.mpr (hV j)) i hdom
  refine ⟨t, fun idx m => ?_⟩
  have hc : ((|V idx| : ℤ) : ℝ) = |(V idx : ℝ)| := by rw [Int.cast_abs]
  have h := ht idx
  rw [hc] at h
  rcases abs_choice (V idx : ℝ) with he | he
  · rw [he] at h; exact h m
  · rw [he, neg_mul] at h
    have := lonely_neg_arg h m
    rwa [neg_neg] at this

#print axioms far_safe_point
#print axioms dominant_lonely
#print axioms hdom_closed
#print axioms hdom_closed_abs

end LRC14
end LonelyRunner
