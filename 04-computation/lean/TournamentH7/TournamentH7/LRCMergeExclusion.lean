/-
  TournamentH7.LRCMergeExclusion — THE MERGE EXCLUSION, FORMAL
  (kind-pasteur-2026-07-05-S5, HYP-4110).

  The S3 gap-value analysis as a kernel theorem.  By grid attainment (HYP-4108)
  the margin maximum sits at `t* = m/s` with `s = |vᵢ| + |vⱼ|`; at a rational
  point every `distZ` is an integer multiple of `1/s` (`distZ_rat_den`), so the
  maximum VALUE is `k/s` with `k ∈ ℤ`.  If that value lies in the spectral-gap
  window `(1/13, 2/25)` then `13k > s` and `25k < 2s`, and integer arithmetic
  kills the small numerators:

    k = 1:  `2s ∈ (25, 26)` — empty;
    k = 2:  `2s = 51` would be odd — parity kill;
    k ≥ 3:  `s > 37.5`, i.e. `s ≥ 38`.

  Hence (`gap_forces_big_pair`): a family whose margin exceeds the `1/13` floor
  somewhere but never reaches `2/25` has TWO runners with `|vᵢ| + |vⱼ| ≥ 38` —
  any spectral-gap violator carries a large binding pair (`w_max ≥ 19`).  This
  completes the formalization of the S3 merge-denominator exclusion: every claim
  in the HYP-4105 gap analysis is now machine-checked.

  Kernel-pure; no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCGridAttainment

namespace LonelyRunner
namespace MergeExclusion

open TournamentH7.LRCWitness GridAttainment

/-- At a rational point the distance to ℤ is an integer multiple of `1/s`. -/
theorem distZ_rat_den (a s : ℤ) (hs : 0 < s) :
    ∃ k : ℤ, 0 ≤ k ∧ distZ ((a : ℝ) / (s : ℝ)) = (k : ℝ) / (s : ℝ) := by
  have hsR : (0 : ℝ) < (s : ℝ) := by exact_mod_cast hs
  refine ⟨|a - s * round ((a : ℝ) / (s : ℝ))|, abs_nonneg _, ?_⟩
  rw [distZ_eq_round]
  rw [show (a : ℝ) / s - (round ((a : ℝ) / s) : ℝ)
      = ((a - s * round ((a : ℝ) / s) : ℤ) : ℝ) / s by push_cast; field_simp]
  rw [abs_div, abs_of_pos hsR, Int.cast_abs]

variable {k : ℕ} [Nonempty (Fin k)]

/-- The margin at a rational point `a/s` is `κ/s` for some integer `κ ≥ 0`. -/
theorem margin_rat_den (v : Fin k → ℤ) (a s : ℤ) (hs : 0 < s) :
    ∃ κ : ℤ, 0 ≤ κ ∧ margin v ((a : ℝ) / (s : ℝ)) = (κ : ℝ) / (s : ℝ) := by
  -- the margin is attained by some runner; that runner's distZ has the form
  obtain ⟨i0, -, hi0⟩ := Finset.exists_mem_eq_inf' (Finset.univ_nonempty (α := Fin k))
    (fun i => distZ ((v i : ℝ) * ((a : ℝ) / (s : ℝ))))
  have hform : ∃ κ : ℤ, 0 ≤ κ ∧
      distZ ((v i0 : ℝ) * ((a : ℝ) / (s : ℝ))) = (κ : ℝ) / (s : ℝ) := by
    have heq : (v i0 : ℝ) * ((a : ℝ) / (s : ℝ)) = ((v i0 * a : ℤ) : ℝ) / (s : ℝ) := by
      push_cast; ring
    rw [heq]
    exact distZ_rat_den (v i0 * a) s hs
  obtain ⟨κ, hκ0, hκ⟩ := hform
  exact ⟨κ, hκ0, by rw [margin, hi0, hκ]⟩

/-- **THE MERGE EXCLUSION (gap ⟹ big pair).**  If the margin exceeds the LRC(13)
floor `1/13` somewhere but no point reaches `2/25`, then two runners have
`|vᵢ| + |vⱼ| ≥ 38`.  Any spectral-gap violator carries a binding pair summing to
at least 38 (in particular `max |vᵢ| ≥ 19`). -/
theorem gap_forces_big_pair (v : Fin k → ℤ) (hv : ∀ i, v i ≠ 0)
    (habove : ∃ t : ℝ, ∀ i, ∀ m : ℤ, 1 / 13 < |(v i : ℝ) * t - m|)
    (hnl : ∀ t : ℝ, ∃ i, ∃ m : ℤ, |(v i : ℝ) * t - m| < 2 / 25) :
    ∃ i j : Fin k, 38 ≤ |v i| + |v j| := by
  -- global maximizer
  obtain ⟨tstar, hmax⟩ := exists_global_max_margin v
  obtain ⟨t1, ht1⟩ := habove
  -- M > 1/13: the strict floor transports through the max
  have hM13 : 1 / 13 < margin v tstar := by
    have h1 : margin v t1 ≤ margin v tstar := hmax t1
    -- margin v t1 > 1/13: the min over runners of strict bounds
    have h2 : 1 / 13 < margin v t1 := by
      rw [margin, Finset.lt_inf'_iff]
      intro i _
      have h3 := (le_distZ_iff (distZ ((v i : ℝ) * t1)) ((v i : ℝ) * t1)).1 le_rfl
      -- distZ is the inf; strictness: distZ = |y − round y| > 1/13 via ht1
      rw [distZ_eq_round]
      exact ht1 i (round ((v i : ℝ) * t1))
    linarith
  have hM0 : 0 < margin v tstar := lt_trans (by norm_num) hM13
  -- M < 2/25: the no-2/25-point hypothesis at t*
  have hM25 : margin v tstar < 2 / 25 := by
    obtain ⟨i, m, him⟩ := hnl tstar
    have h1 : margin v tstar ≤ distZ ((v i : ℝ) * tstar) := margin_le_distZ v tstar i
    have h2 : distZ ((v i : ℝ) * tstar) ≤ |(v i : ℝ) * tstar - m| :=
      distZ_le_abs_sub _ m
    linarith
  -- the grid form
  obtain ⟨i, j, m, hgrid⟩ := maximizer_on_grid v hv tstar hmax hM0
  set s : ℤ := |v i| + |v j| with hsdef
  have hs0 : 0 < s := by
    have h1 : (0 : ℤ) < |v i| := abs_pos.mpr (hv i)
    have h2 : (0 : ℤ) < |v j| := abs_pos.mpr (hv j)
    omega
  have hsR : (0 : ℝ) < (s : ℝ) := by exact_mod_cast hs0
  -- t* = m/s
  have ht : tstar = (m : ℝ) / (s : ℝ) := by
    rw [eq_div_iff (ne_of_gt hsR)]
    rw [mul_comm]
    exact_mod_cast hgrid
  -- the margin VALUE is κ/s
  obtain ⟨κ, hκ0, hκ⟩ := margin_rat_den v m s hs0
  rw [← ht] at hκ
  -- 1/13 < κ/s < 2/25 with κ, s integers ⟹ s ≥ 38
  rw [hκ] at hM13 hM25
  have h13 : (13 : ℝ) * κ > s := by
    rw [div_lt_div_iff₀ (by norm_num) hsR] at hM13
    linarith
  have h25 : (25 : ℝ) * κ < 2 * s := by
    rw [div_lt_div_iff₀ hsR (by norm_num)] at hM25
    linarith
  have h13' : 13 * κ > s := by exact_mod_cast h13
  have h25' : 25 * κ < 2 * s := by exact_mod_cast h25
  exact ⟨i, j, by omega⟩

#print axioms distZ_rat_den
#print axioms margin_rat_den
#print axioms gap_forces_big_pair

end MergeExclusion
end LonelyRunner
