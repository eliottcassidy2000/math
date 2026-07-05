/-
  TournamentH7.LRCMergeGridAttainment — THE MERGE-GRID ATTAINMENT THEOREM
  (klein-2026-07-05-S137, HYP-4114; THM-592's formalization, kps-S3's suggested-next).

  The max-min profile  M(V) = max_t min_i dist(vᵢt, ℤ)  of a finite positive family is
  ATTAINED, and attained at a MERGE-GRID point  t* = m/(vᵢ + vⱼ)  (i = j giving the
  half-integer peaks).  This is the theorem beneath every exact-M sweep in the project:
  with it, scanning the finite grid `{m/(vᵢ+vⱼ)}` is a KERNEL-CHECKABLE evaluation of M.

  Stage A (this section): the distance-to-ℤ toolkit (`distZ`, nearest-point property,
  1-Lipschitz continuity, integer periodicity), the profile as a finite inf, its
  continuity/periodicity/positivity, and ATTAINMENT via compactness + period transport.
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner
namespace MergeGrid

/-- Distance to the nearest integer, through `round`. -/
noncomputable def distZ (x : ℝ) : ℝ := |x - round x|

theorem distZ_nonneg (x : ℝ) : 0 ≤ distZ x := abs_nonneg _

theorem distZ_le_half (x : ℝ) : distZ x ≤ 1/2 := by
  simpa [distZ] using abs_sub_round x

/-- `round` is the nearest integer: `distZ x ≤ |x − m|` for every integer `m`. -/
theorem distZ_le (x : ℝ) (m : ℤ) : distZ x ≤ |x - m| := by
  rcases eq_or_ne m (round x) with rfl | hne
  · simp [distZ]
  · have h1 : (1 : ℝ) ≤ |(m : ℝ) - round x| := by
      have hne' : m - round x ≠ 0 := sub_ne_zero.mpr hne
      have hz : (1 : ℤ) ≤ |m - round x| := Int.one_le_abs hne'
      have hz' : ((1 : ℤ) : ℝ) ≤ ((|m - round x| : ℤ) : ℝ) := by exact_mod_cast hz
      rw [Int.cast_abs] at hz'
      push_cast at hz'
      linarith
    have h2 : |x - round x| ≤ 1/2 := abs_sub_round x
    have htri : |(m : ℝ) - round x| ≤ |(m : ℝ) - x| + |x - round x| := abs_sub_le _ _ _
    have hsw : |(m : ℝ) - x| = |x - (m : ℝ)| := abs_sub_comm _ _
    unfold distZ
    linarith
  
/-- `distZ` is 1-Lipschitz. -/
theorem distZ_lipschitz (x y : ℝ) : distZ x ≤ |x - y| + distZ y := by
  have h1 : distZ x ≤ |x - round y| := distZ_le x (round y)
  have h2 : |x - (round y : ℝ)| ≤ |x - y| + |y - round y| := abs_sub_le _ _ _
  unfold distZ at *
  linarith

theorem distZ_continuous : Continuous distZ := by
  have hlip : LipschitzWith 1 distZ := by
    apply LipschitzWith.of_dist_le_mul
    intro x y
    rw [Real.dist_eq, Real.dist_eq, NNReal.coe_one, one_mul]
    rw [abs_le]
    constructor
    · have := distZ_lipschitz y x
      have habs : |y - x| = |x - y| := abs_sub_comm _ _
      linarith [habs ▸ this]
    · linarith [distZ_lipschitz x y]
  exact hlip.continuous

/-- Integer shift invariance. -/
theorem distZ_add_int (x : ℝ) (k : ℤ) : distZ (x + k) = distZ x := by
  unfold distZ
  rw [round_add_intCast]
  push_cast
  ring_nf

/-- The min-distance profile of a positive family (nonempty index). -/
noncomputable def profile {n : ℕ} (v : Fin (n+1) → ℤ) (t : ℝ) : ℝ :=
  Finset.univ.inf' Finset.univ_nonempty fun i => distZ ((v i : ℝ) * t)

theorem profile_continuous {n : ℕ} (v : Fin (n+1) → ℤ) : Continuous (profile v) := by
  apply Continuous.finset_inf'_apply
  intro i _
  exact distZ_continuous.comp (continuous_const.mul continuous_id)

theorem profile_periodic {n : ℕ} (v : Fin (n+1) → ℤ) :
    Function.Periodic (profile v) 1 := by
  intro t
  unfold profile
  congr 1
  funext i
  have hexp : (v i : ℝ) * (t + 1) = (v i : ℝ) * t + ((v i : ℤ) : ℝ) := by
    push_cast
    ring
  rw [hexp, distZ_add_int]

/-- The profile is somewhere positive: at `t₀ = 1/(2W+1)` every runner sits in
`(0, 1/2)`. -/
theorem profile_pos {n : ℕ} (v : Fin (n+1) → ℤ) (hv : ∀ i, 0 < v i) (W : ℤ)
    (hW : ∀ i, v i ≤ W) :
    (1 : ℝ) / (2 * (W : ℝ) + 1) ≤ profile v (1 / (2 * (W : ℝ) + 1)) := by
  have hWpos : (0 : ℝ) < (W : ℝ) := by
    have := hv 0
    have := hW 0
    exact_mod_cast lt_of_lt_of_le (hv 0) (hW 0)
  have hden : (0 : ℝ) < 2 * (W : ℝ) + 1 := by linarith
  unfold profile
  rw [Finset.le_inf'_iff]
  intro i _
  have hviR : (0 : ℝ) < (v i : ℝ) := by exact_mod_cast hv i
  have hviW : (v i : ℝ) ≤ (W : ℝ) := by exact_mod_cast hW i
  set x : ℝ := (v i : ℝ) * (1 / (2 * (W : ℝ) + 1)) with hx
  have hxpos : 0 < x := by
    rw [hx]
    positivity
  have hxlt : x < 1/2 := by
    rw [hx]
    rw [mul_one_div, div_lt_iff₀ hden]
    linarith
  have hround : round x = 0 := by
    apply round_eq_zero_iff.mpr
    exact Set.mem_Ico.mpr ⟨by linarith, by linarith⟩
  unfold distZ
  rw [hround]
  push_cast
  rw [sub_zero, abs_of_pos hxpos, hx]
  have h1 : (1 : ℝ) ≤ (v i : ℝ) := by exact_mod_cast hv i
  calc (1 : ℝ) / (2 * (W : ℝ) + 1) = 1 * (1 / (2 * (W : ℝ) + 1)) := by ring
    _ ≤ (v i : ℝ) * (1 / (2 * (W : ℝ) + 1)) :=
        mul_le_mul_of_nonneg_right h1 (le_of_lt (div_pos one_pos hden))

/-- **ATTAINMENT**: the profile attains its global maximum. -/
theorem profile_attains_max {n : ℕ} (v : Fin (n+1) → ℤ) :
    ∃ tstar : ℝ, ∀ t : ℝ, profile v t ≤ profile v tstar := by
  obtain ⟨tstar, -, hmax⟩ :=
    (isCompact_Icc (a := (0:ℝ)) (b := 1)).exists_isMaxOn
      (Set.nonempty_Icc.mpr (by norm_num))
      (profile_continuous v).continuousOn
  refine ⟨tstar, fun t => ?_⟩
  have hper : profile v (t - (⌊t⌋ : ℝ) * 1) = profile v t :=
    (profile_periodic v).sub_int_mul_eq ⌊t⌋
  rw [mul_one] at hper
  rw [← hper]
  have hf : t - (⌊t⌋ : ℝ) = Int.fract t := Int.self_sub_floor t
  refine hmax (Set.mem_Icc.mpr ⟨?_, ?_⟩)
  · rw [hf]
    exact Int.fract_nonneg t
  · rw [hf]
    exact le_of_lt (Int.fract_lt_one t)

end MergeGrid
end LonelyRunner
