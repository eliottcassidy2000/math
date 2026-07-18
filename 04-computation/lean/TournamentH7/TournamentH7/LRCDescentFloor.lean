/-
  TournamentH7.LRCDescentFloor — THM-1008 + THM-1010, the LRC-descent floor
  (boxeph-2026-07-18, S83/S84).

  Remove the largest speed `v⋆`; a `μ`-lonely time `t₀` for the other twelve
  (from LRC(13) with `μ = M(V∖{v⋆}) ≥ 1/13`) lifts to a `1/14`-lonely time for the
  whole family, PROVIDED each kept `v_i` obeys the sweep budget `v_i ≤ 14 v⋆ (μ−1/14)`.
  The minimal kick `s = ±1/(14 v⋆)` lifts `v⋆` into the `1/14` band while every kept
  runner loses at most `μ − 1/14`, staying `≥ 1/14`.

  `descent_general` (THM-1010, S84) is the SHARP recursion: with `μ = M(V∖{v⋆})` it
  gives `M(V) ≥ ρ·M(V∖{v⋆})/(ρ+1)`, strictly beating the factor-2 floor `1/26` once
  the kept set is non-tight (`μ > 1/13`).  `descent_dominant` (THM-1008, S83) is the
  `μ = 1/13` corollary: `ρ = v_max/v_2nd ≥ 13 ⟹ M ≥ 1/14`.

  Conditional on the kept-set lonely time (`hkept`), so kernel-pure and citation-free.
  Elementary: `round` + reverse triangle on ℝ/ℤ.  Companions: `sieve_frac`,
  `fill1_perturbation`.
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner

variable {ι : Type*}

/-- If `1/14 ≤ |w|` and `|w| ≤ 13/14`, then `w` is `≥ 1/14` from every integer. -/
theorem far_from_int (w : ℝ) (h1 : (1 : ℝ) / 14 ≤ |w|) (h2 : |w| ≤ 13 / 14) :
    ∀ k : ℤ, (1 : ℝ) / 14 ≤ |w + k| := by
  intro k
  have hwle : w ≤ 13 / 14 := (abs_le.mp h2).2
  have hwge : -(13 / 14) ≤ w := (abs_le.mp h2).1
  rcases lt_trichotomy k 0 with hk | hk | hk
  · have hkR : (k : ℝ) ≤ -1 := by exact_mod_cast (by omega : k ≤ -1)
    rw [abs_of_nonpos (by linarith)]; linarith
  · subst hk; simpa using h1
  · have hkR : (1 : ℝ) ≤ (k : ℝ) := by exact_mod_cast (by omega : 1 ≤ k)
    rw [abs_of_pos (by linarith)]; linarith

/-- **THM-1010 (sharp descent recursion).**  A `μ`-lonely time `t₀` for the family
with `v⋆` removed lifts to a `1/14`-lonely time for the whole family, provided every
kept `v_i` obeys the sweep budget `v_i ≤ 14·v⋆·(μ − 1/14)`. -/
theorem descent_general
    (v : ι → ℤ) (vstar : ι) (t0 μ : ℝ)
    (hkept : ∀ i, i ≠ vstar → ∀ m : ℤ, μ ≤ |(v i : ℝ) * t0 - m|)
    (hpos : ∀ i, 0 < v i)
    (hbudget : ∀ i, i ≠ vstar → (v i : ℝ) ≤ 14 * (v vstar : ℝ) * (μ - 1 / 14)) :
    ∃ t : ℝ, Lonely 14 v t := by
  have hvspos : (0 : ℝ) < (v vstar : ℝ) := by exact_mod_cast hpos vstar
  have hvsne : (v vstar : ℝ) ≠ 0 := ne_of_gt hvspos
  set g : ℝ := (v vstar : ℝ) * t0 with hg
  set n0 : ℤ := round g with hn0
  set f : ℝ := g - n0 with hf
  have hf2 : |f| ≤ 1 / 2 := by rw [hf]; exact abs_sub_round g
  have hfle : f ≤ 1 / 2 := (abs_le.mp hf2).2
  have hfge : -(1 / 2) ≤ f := (abs_le.mp hf2).1
  have hgf : (v vstar : ℝ) * t0 = (n0 : ℝ) + f := by rw [← hg, hf]; ring
  -- shared small-kick bound: |v_i s| ≤ μ − 1/14 whenever |s| ≤ 1/(14 v⋆).
  have hsmall : ∀ (s : ℝ), |s| ≤ 1 / (14 * (v vstar : ℝ)) →
      ∀ i, i ≠ vstar → |(v i : ℝ) * s| ≤ μ - 1 / 14 := by
    intro s hs i hi
    have hvipos : (0 : ℝ) < (v i : ℝ) := by exact_mod_cast hpos i
    rw [abs_mul, abs_of_pos hvipos]
    calc (v i : ℝ) * |s| ≤ (v i : ℝ) * (1 / (14 * (v vstar : ℝ))) :=
            mul_le_mul_of_nonneg_left hs (le_of_lt hvipos)
      _ ≤ μ - 1 / 14 := by
            rw [mul_one_div, div_le_iff₀ (by positivity)]; nlinarith [hbudget i hi]
  -- shared tail: assemble Lonely from the vstar value and the kept-runner bound.
  have tail : ∀ (s w : ℝ),
      (v vstar : ℝ) * (t0 + s) = w + (n0 : ℝ) →
      (1 : ℝ) / 14 ≤ |w| → |w| ≤ 13 / 14 →
      (∀ i, i ≠ vstar → |(v i : ℝ) * s| ≤ μ - 1 / 14) →
      Lonely 14 v (t0 + s) := by
    intro s w hbig h1 h2 hsm i m
    by_cases hi : i = vstar
    · rw [hi]
      have hval : (v vstar : ℝ) * (t0 + s) - m = w + ((n0 - m : ℤ) : ℝ) := by
        rw [hbig]; push_cast; ring
      rw [hval]; exact far_from_int w h1 h2 (n0 - m)
    · have hkm := hkept i hi m
      have hsplit : (v i : ℝ) * (t0 + s) - m
          = ((v i : ℝ) * t0 - m) + (v i : ℝ) * s := by ring
      have hrt := abs_sub_abs_le_abs_sub ((v i : ℝ) * t0 - m) (-((v i : ℝ) * s))
      rw [abs_neg, sub_neg_eq_add, ← hsplit] at hrt
      linarith [hkm, hrt, hsm i hi]
  have hkickpos : (0 : ℝ) < 1 / (14 * (v vstar : ℝ)) := by positivity
  rcases le_or_gt 0 f with hf0 | hf0
  · refine ⟨t0 + 1 / (14 * (v vstar : ℝ)), ?_⟩
    have hbig : (v vstar : ℝ) * (t0 + 1 / (14 * (v vstar : ℝ))) = (f + 1 / 14) + (n0 : ℝ) := by
      have hmid : (v vstar : ℝ) * (1 / (14 * (v vstar : ℝ))) = 1 / 14 := by field_simp
      rw [mul_add, hmid, hgf]; ring
    have hsbound : |(1 : ℝ) / (14 * (v vstar : ℝ))| ≤ 1 / (14 * (v vstar : ℝ)) := by
      rw [abs_of_pos hkickpos]
    refine tail _ (f + 1 / 14) hbig ?_ ?_ (hsmall _ hsbound)
    · rw [abs_of_nonneg (by linarith [hf0])]; linarith [hf0]
    · rw [abs_of_nonneg (by linarith [hf0])]; linarith [hfle]
  · refine ⟨t0 + (-(1 / (14 * (v vstar : ℝ)))), ?_⟩
    have hbig : (v vstar : ℝ) * (t0 + (-(1 / (14 * (v vstar : ℝ))))) = (f - 1 / 14) + (n0 : ℝ) := by
      have hmid : (v vstar : ℝ) * (-(1 / (14 * (v vstar : ℝ)))) = -(1 / 14) := by field_simp
      rw [mul_add, hmid, hgf]; ring
    have hsbound : |(-(1 / (14 * (v vstar : ℝ))) : ℝ)| ≤ 1 / (14 * (v vstar : ℝ)) := by
      rw [abs_neg, abs_of_pos hkickpos]
    refine tail _ (f - 1 / 14) hbig ?_ ?_ (hsmall _ hsbound)
    · rw [abs_of_neg (by linarith [hf0])]; linarith [hf0]
    · rw [abs_of_neg (by linarith [hf0])]; linarith [hfge]

/-- **THM-1008 (descent floor, dominant regime).**  The `μ = 1/13` corollary: a
`1/13`-lonely time for `V∖{v⋆}` with `13·v_i ≤ v⋆` for every kept `v_i` (i.e.
`ρ = v_max/v_2nd ≥ 13`) yields a `1/14`-lonely time. -/
theorem descent_dominant
    (v : ι → ℤ) (vstar : ι) (t0 : ℝ)
    (hkept : ∀ i, i ≠ vstar → ∀ m : ℤ, (1 : ℝ) / 13 ≤ |(v i : ℝ) * t0 - m|)
    (hpos : ∀ i, 0 < v i)
    (hdom : ∀ i, i ≠ vstar → 13 * v i ≤ v vstar) :
    ∃ t : ℝ, Lonely 14 v t := by
  refine descent_general v vstar t0 (1 / 13) hkept hpos ?_
  intro i hi
  have hdomR : (13 : ℝ) * (v i : ℝ) ≤ (v vstar : ℝ) := by exact_mod_cast hdom i hi
  have hvipos : (0 : ℝ) < (v i : ℝ) := by exact_mod_cast hpos i
  nlinarith [hdomR, hvipos]

end LonelyRunner

#print axioms LonelyRunner.descent_general
#print axioms LonelyRunner.descent_dominant
