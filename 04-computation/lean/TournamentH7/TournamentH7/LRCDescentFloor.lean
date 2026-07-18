/-
  TournamentH7.LRCDescentFloor — THM-1008, the LRC(13)-descent floor, dominant regime
  (boxeph-2026-07-18-S83).

  Remove the largest speed `v⋆`; LRC(13) [SETTLED] gives a `1/13`-lonely time `t₀`
  for the other twelve.  In the DOMINANT regime `v⋆ ≥ 13·v_i` for every kept `v_i`
  (i.e. `ρ = v_max/v_2nd ≥ 13`), the minimal kick `s = ±1/(14 v⋆)` that lifts `v⋆`
  to the `1/14` band keeps every kept runner `≥ 1/13 − 1/182 = 1/14`, producing a
  `1/14`-lonely time for the whole family.

  Conditional on the LRC(13) witness (`hkept`), so kernel-pure and citation-free at
  the Lean level.  Companion to `sieve_frac` and `fill1_perturbation`.  Elementary:
  `round` + reverse triangle on ℝ/ℤ.
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

/-- **THM-1008 (descent floor, dominant regime).**  Given a `1/13`-lonely time `t₀`
for the family with `v⋆` removed, positive speeds, and `13·v_i ≤ v⋆` for every
kept `v_i`, the whole family has a `1/14`-lonely time. -/
theorem descent_dominant
    (v : ι → ℤ) (vstar : ι) (t0 : ℝ)
    (hkept : ∀ i, i ≠ vstar → ∀ m : ℤ, (1 : ℝ) / 13 ≤ |(v i : ℝ) * t0 - m|)
    (hpos : ∀ i, 0 < v i)
    (hdom : ∀ i, i ≠ vstar → 13 * v i ≤ v vstar) :
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
  -- shared small-kick bound: |v_i s| ≤ 1/182 whenever |s| ≤ 1/(14 v⋆).
  have hsmall : ∀ (s : ℝ), |s| ≤ 1 / (14 * (v vstar : ℝ)) →
      ∀ i, i ≠ vstar → |(v i : ℝ) * s| ≤ 1 / 182 := by
    intro s hs i hi
    have hvipos : (0 : ℝ) < (v i : ℝ) := by exact_mod_cast hpos i
    have hdomR : (13 : ℝ) * (v i : ℝ) ≤ (v vstar : ℝ) := by exact_mod_cast hdom i hi
    rw [abs_mul, abs_of_pos hvipos]
    calc (v i : ℝ) * |s| ≤ (v i : ℝ) * (1 / (14 * (v vstar : ℝ))) :=
            mul_le_mul_of_nonneg_left hs (le_of_lt hvipos)
      _ ≤ 1 / 182 := by
            rw [mul_one_div, div_le_iff₀ (by positivity)]; nlinarith [hdomR]
  -- shared tail: assemble Lonely from the vstar value and the kept-runner bound.
  have tail : ∀ (s w : ℝ),
      (v vstar : ℝ) * (t0 + s) = w + (n0 : ℝ) →
      (1 : ℝ) / 14 ≤ |w| → |w| ≤ 13 / 14 →
      (∀ i, i ≠ vstar → |(v i : ℝ) * s| ≤ 1 / 182) →
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
  -- kick size and its |s| bound
  have hkickpos : (0 : ℝ) < 1 / (14 * (v vstar : ℝ)) := by positivity
  rcases le_or_gt 0 f with hf0 | hf0
  · -- f ≥ 0: kick up.
    refine ⟨t0 + 1 / (14 * (v vstar : ℝ)), ?_⟩
    have hbig : (v vstar : ℝ) * (t0 + 1 / (14 * (v vstar : ℝ))) = (f + 1 / 14) + (n0 : ℝ) := by
      have hmid : (v vstar : ℝ) * (1 / (14 * (v vstar : ℝ))) = 1 / 14 := by field_simp
      rw [mul_add, hmid, hgf]; ring
    have hsbound : |(1 : ℝ) / (14 * (v vstar : ℝ))| ≤ 1 / (14 * (v vstar : ℝ)) := by
      rw [abs_of_pos hkickpos]
    refine tail _ (f + 1 / 14) hbig ?_ ?_ (hsmall _ hsbound)
    · rw [abs_of_nonneg (by linarith [hf0])]; linarith [hf0]
    · rw [abs_of_nonneg (by linarith [hf0])]; linarith [hfle]
  · -- f < 0: kick down.
    refine ⟨t0 + (-(1 / (14 * (v vstar : ℝ)))), ?_⟩
    have hbig : (v vstar : ℝ) * (t0 + (-(1 / (14 * (v vstar : ℝ))))) = (f - 1 / 14) + (n0 : ℝ) := by
      have hmid : (v vstar : ℝ) * (-(1 / (14 * (v vstar : ℝ)))) = -(1 / 14) := by field_simp
      rw [mul_add, hmid, hgf]; ring
    have hsbound : |(-(1 / (14 * (v vstar : ℝ))) : ℝ)| ≤ 1 / (14 * (v vstar : ℝ)) := by
      rw [abs_neg, abs_of_pos hkickpos]
    refine tail _ (f - 1 / 14) hbig ?_ ?_ (hsmall _ hsbound)
    · rw [abs_of_neg (by linarith [hf0])]; linarith [hf0]
    · rw [abs_of_neg (by linarith [hf0])]; linarith [hfge]

end LonelyRunner

#print axioms LonelyRunner.descent_dominant
