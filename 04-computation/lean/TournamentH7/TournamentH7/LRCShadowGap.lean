/-
  TournamentH7.LRCShadowGap — the SHADOW-GAP middle-witness (THM-744, klein-2026-07-13-S297/S298).

  For a covering cluster `C` with smallest even element `e` and every element `≤ m`, if `m < 6e`
  then the family `{1} ∪ C` is `1/14`-lonely at a time `t = 1/2 + δ` in the MIDDLE of the circle,
  for any `δ` in the (non-empty) window `(1/(14e), 3/(7m))`.

  This is the first UNCONDITIONAL middle-reach for tight covering clusters — exactly the residual
  that the isolated-far disc certificate and the simultaneous multi-peel do not cover.

  The proof is elementary: at `t = 1/2` the parities separate.  An even speed `c` lands at
  `‖ct‖ = cδ` (good once `cδ ≥ 1/14`); an odd speed lands at `‖ct‖ = 1/2 − cδ` (good while
  `cδ ≤ 3/7`).  Each reduces to a single fact: a real in `[1/14, 13/14]` is `≥ 1/14` from every
  integer.  Fully sorry-free, all-ℚ, no measure theory.
-/
import TournamentH7.LRC14CertRoute
import Mathlib.Tactic

namespace LonelyRunner
namespace LRC14

/-- **Core distance bound.**  A real `y ∈ [1/14, 13/14]` is at distance `≥ 1/14` from every
integer (the nearest integers are `0` and `1`, at distances `y ≥ 1/14` and `1 − y ≥ 1/14`). -/
theorem dist_ge_of_mem {y : ℝ} (h1 : 1 / 14 ≤ y) (h2 : y ≤ 13 / 14) (j : ℤ) :
    (1 : ℝ) / 14 ≤ |y - (j : ℝ)| := by
  by_cases hj : (j : ℝ) ≤ 0
  · rw [abs_of_nonneg (by linarith)]
    linarith
  · push_neg at hj
    have hjZ : (0 : ℤ) < j := by exact_mod_cast hj
    have hj1 : (1 : ℝ) ≤ (j : ℝ) := by exact_mod_cast hjZ
    rw [abs_of_nonpos (by linarith)]
    linarith

/-- **Reduction to the fractional value.**  If `c·t = N + y` with `N` an integer and
`y ∈ [1/14, 13/14]`, then speed `c` is `1/14`-lonely at `t`. -/
theorem lonely_of_reduce (c N : ℤ) (y : ℝ) {t : ℝ}
    (heq : (c : ℝ) * t = (N : ℝ) + y) (h1 : 1 / 14 ≤ y) (h2 : y ≤ 13 / 14) :
    LonelySpeed 14 c t := by
  intro m
  have hn : (1 : ℝ) / (14 : ℕ) = 1 / 14 := by norm_num
  rw [hn]
  have hrw : (c : ℝ) * t - (m : ℝ) = y - (((m - N : ℤ)) : ℝ) := by
    rw [heq]; push_cast; ring
  rw [hrw]
  exact dist_ge_of_mem h1 h2 (m - N)

/-- **Even speed at `t = 1/2 + δ`.**  For an even `c` with `e ≤ c ≤ m`, given `1 < 14 e δ`
(so `cδ ≥ 1/14`) and `7 m δ < 3` (so `cδ ≤ 3/7 < 13/14`), speed `c` is lonely: `‖ct‖ = cδ`. -/
theorem shadowGap_even (e m : ℤ) (he : 1 ≤ e) {δ : ℝ}
    (hlo : 1 < 14 * (e : ℝ) * δ) (hhi : 7 * (m : ℝ) * δ < 3)
    (c : ℤ) (hev : Even c) (hce : e ≤ c) (hcm : c ≤ m) :
    LonelySpeed 14 c (1 / 2 + δ) := by
  have heR : (1 : ℝ) ≤ (e : ℝ) := by exact_mod_cast he
  have hδ0 : 0 < δ := by nlinarith [heR, hlo]
  have hceR : (e : ℝ) ≤ (c : ℝ) := by exact_mod_cast hce
  have hcmR : (c : ℝ) ≤ (m : ℝ) := by exact_mod_cast hcm
  obtain ⟨N, hN⟩ := hev
  have hcR : (c : ℝ) = 2 * (N : ℝ) := by rw [hN]; push_cast; ring
  refine lonely_of_reduce c N ((c : ℝ) * δ) ?_ ?_ ?_
  · rw [hcR]; ring
  · nlinarith [mul_le_mul_of_nonneg_right hceR hδ0.le, hlo]
  · nlinarith [mul_le_mul_of_nonneg_right hcmR hδ0.le, hhi]

/-- **Odd speed at `t = 1/2 + δ`.**  For an odd `c` with `1 ≤ c ≤ m` and `δ ≥ 0`, given
`7 m δ < 3` (so `cδ < 3/7`), speed `c` is lonely: it lands at `‖ct‖ = 1/2 − cδ`. -/
theorem shadowGap_odd (m : ℤ) {δ : ℝ} (hδ0 : 0 ≤ δ) (hhi : 7 * (m : ℝ) * δ < 3)
    (c : ℤ) (hodd : Odd c) (hc1 : 1 ≤ c) (hcm : c ≤ m) :
    LonelySpeed 14 c (1 / 2 + δ) := by
  have hc1R : (1 : ℝ) ≤ (c : ℝ) := by exact_mod_cast hc1
  have hcmR : (c : ℝ) ≤ (m : ℝ) := by exact_mod_cast hcm
  obtain ⟨N, hN⟩ := hodd
  have hcR : (c : ℝ) = 2 * (N : ℝ) + 1 := by rw [hN]; push_cast; ring
  refine lonely_of_reduce c N (1 / 2 + (c : ℝ) * δ) ?_ ?_ ?_
  · rw [hcR]; ring
  · have : (0 : ℝ) ≤ (c : ℝ) * δ := mul_nonneg (by linarith) hδ0
    linarith
  · nlinarith [mul_le_mul_of_nonneg_right hcmR hδ0, hhi]

/-- **THM-744 (shadow-gap middle-witness).**  If every speed of `v` is either an even integer in
`[e, m]` or an odd integer in `[1, m]` (`1 ≤ e`), then for any `δ` with `1 < 14 e δ` and
`7 m δ < 3`, the family is `1/14`-lonely at the middle time `t = 1/2 + δ`.  The observer speed `1`
(odd, `≤ m`) and every cluster speed are handled uniformly. -/
theorem shadowGap_lonely {ι : Type*} (v : ι → ℤ) (e m : ℤ) (he : 1 ≤ e) {δ : ℝ}
    (hlo : 1 < 14 * (e : ℝ) * δ) (hhi : 7 * (m : ℝ) * δ < 3)
    (hcase : ∀ i, (Even (v i) ∧ e ≤ v i ∧ v i ≤ m) ∨ (Odd (v i) ∧ 1 ≤ v i ∧ v i ≤ m)) :
    Lonely 14 v (1 / 2 + δ) := by
  have heR : (1 : ℝ) ≤ (e : ℝ) := by exact_mod_cast he
  have hδ0 : (0 : ℝ) ≤ δ := by nlinarith [heR, hlo]
  intro i
  rcases hcase i with ⟨hev, hce, hcm⟩ | ⟨hod, hc1, hcm⟩
  · exact shadowGap_even e m he hlo hhi (v i) hev hce hcm
  · exact shadowGap_odd m hδ0 hhi (v i) hod hc1 hcm

/-- The `δ`-window `(1/(14e), 3/(7m))` is non-empty exactly when `m < 6e`, so a lonely time
in the middle exists: `{1} ∪ C` (a tight covering cluster) has a `1/14`-lonely time. -/
theorem shadowGap_exists {ι : Type*} (v : ι → ℤ) (e m : ℤ) (he : 1 ≤ e) (hm : 1 ≤ m)
    (hme : m < 6 * e)
    (hcase : ∀ i, (Even (v i) ∧ e ≤ v i ∧ v i ≤ m) ∨ (Odd (v i) ∧ 1 ≤ v i ∧ v i ≤ m)) :
    ∃ t : ℝ, Lonely 14 v t := by
  have heR : (1 : ℝ) ≤ (e : ℝ) := by exact_mod_cast he
  have hmR : (1 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm
  have hmeR : (m : ℝ) < 6 * (e : ℝ) := by exact_mod_cast hme
  have hepos : (0 : ℝ) < (e : ℝ) := by linarith
  have hmpos : (0 : ℝ) < (m : ℝ) := by linarith
  have hwin : 1 / (14 * (e : ℝ)) < 3 / (7 * (m : ℝ)) := by
    have hnum : (0 : ℝ) < 42 * (e : ℝ) - 7 * (m : ℝ) := by nlinarith [hmeR]
    have hden : (0 : ℝ) < 7 * (m : ℝ) * (14 * (e : ℝ)) := by positivity
    have hsub : 3 / (7 * (m : ℝ)) - 1 / (14 * (e : ℝ))
        = (42 * (e : ℝ) - 7 * (m : ℝ)) / (7 * (m : ℝ) * (14 * (e : ℝ))) := by
      field_simp; ring
    have hpos : 0 < 3 / (7 * (m : ℝ)) - 1 / (14 * (e : ℝ)) := by
      rw [hsub]; exact div_pos hnum hden
    linarith
  obtain ⟨δ, hδ1, hδ2⟩ := exists_between hwin
  have hc1 : (14 * (e : ℝ)) * (1 / (14 * (e : ℝ))) = 1 := by field_simp
  have hc2 : (7 * (m : ℝ)) * (3 / (7 * (m : ℝ))) = 3 := by field_simp
  have hlo : 1 < 14 * (e : ℝ) * δ := by
    have := mul_lt_mul_of_pos_left hδ1 (show (0 : ℝ) < 14 * (e : ℝ) by positivity)
    rw [hc1] at this; linarith
  have hhi : 7 * (m : ℝ) * δ < 3 := by
    have := mul_lt_mul_of_pos_left hδ2 (show (0 : ℝ) < 7 * (m : ℝ) by positivity)
    rw [hc2] at this; linarith
  exact ⟨1 / 2 + δ, shadowGap_lonely v e m he hlo hhi hcase⟩

end LRC14
end LonelyRunner
