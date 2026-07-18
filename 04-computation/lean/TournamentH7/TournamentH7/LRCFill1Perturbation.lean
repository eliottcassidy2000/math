/-
  TournamentH7.LRCFill1Perturbation — THM-1003, the fill-1 perturbation base case
  (boxeph-2026-07-17-S82).

  The under-filled-circle clause of the resonance-fill crux, base case: if a
  single speed `v⋆` is divisible by `b` (`2 ≤ b < N`) and every other speed is
  NOT, then the minimal kick off the resonance center `a/b` that lifts `v⋆` to
  exactly `1/N`, namely `t = a/b + 1/(N·v⋆)`, is an `N`-lonely time — PROVIDED
  the body is dominated: `b·(v i + v⋆) ≤ N·v⋆` for every other speed `v i`
  (equivalently `b·B ≤ (N−b)·v⋆`, `B = body max`).

  Elementary: one reverse-triangle inequality on ℝ/ℤ. Companion to the sieve
  lemma (`LonelyRunner.sieve_frac`, the empty-circle f_b=0 witness).
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner

variable {ι : Type*}

/-- **THM-1003 (fill-1 perturbation, base case).**  One speed `v⋆` divisible by
`b`, all others not; the body dominated (`b·(v i + v⋆) ≤ N·v⋆`).  Then
`t = a/b + 1/(N·v⋆)` is `N`-lonely. -/
theorem fill1_perturbation
    (N b : ℕ) (a : ℤ) (v : ι → ℤ) (vstar : ι)
    (hN : 2 ≤ N) (hb2 : 2 ≤ b) (hbN : b < N)
    (hcop : IsCoprime (b : ℤ) a)
    (hstar_dvd : (b : ℤ) ∣ v vstar) (hstar_pos : 0 < v vstar)
    (hfill1 : ∀ i, i ≠ vstar → ¬ (b : ℤ) ∣ v i)
    (hbody_pos : ∀ i, i ≠ vstar → 0 < v i)
    (hdom : ∀ i, i ≠ vstar → (b : ℤ) * (v i + v vstar) ≤ (N : ℤ) * v vstar) :
    Lonely N v ((a : ℝ) / b + 1 / (N * (v vstar : ℝ))) := by
  have hN0 : (0 : ℝ) < N := by exact_mod_cast (by omega : 0 < N)
  have hb0 : (0 : ℝ) < b := by exact_mod_cast (by omega : 0 < b)
  have hvs0 : (0 : ℝ) < (v vstar : ℝ) := by exact_mod_cast hstar_pos
  have hNne : (N : ℝ) ≠ 0 := ne_of_gt hN0
  have hbne : (b : ℝ) ≠ 0 := ne_of_gt hb0
  have hvsne : (v vstar : ℝ) ≠ 0 := ne_of_gt hvs0
  have hinvN : (1 : ℝ) / N ≤ 1 / 2 :=
    one_div_le_one_div_of_le (by norm_num) (by exact_mod_cast hN)
  set t : ℝ := (a : ℝ) / b + 1 / (N * (v vstar : ℝ)) with ht
  intro i m
  by_cases hi : i = vstar
  · -- Stranded runner: v⋆·t = (integer) + 1/N, distance to ℤ is exactly 1/N.
    rw [hi]
    obtain ⟨w, hw⟩ := hstar_dvd
    have hwR : (v vstar : ℝ) = (b : ℝ) * (w : ℝ) := by rw [hw]; push_cast; ring
    have hwne : (w : ℝ) ≠ 0 := by
      have hbw : (b : ℝ) * (w : ℝ) ≠ 0 := by rw [← hwR]; exact hvsne
      exact right_ne_zero_of_mul hbw
    have hvt : (v vstar : ℝ) * t - m = ((w * a - m : ℤ) : ℝ) + 1 / N := by
      rw [ht, hwR]; push_cast; field_simp; ring
    rw [hvt]
    set k : ℤ := w * a - m with hk
    rcases le_or_gt 0 k with hk0 | hk0
    · have hkR : (0 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk0
      have h1N : (0 : ℝ) ≤ 1 / N := by positivity
      rw [abs_of_nonneg (by linarith)]
      linarith
    · have hkR : (k : ℝ) ≤ -1 := by exact_mod_cast (by omega : k ≤ -1)
      have hneg : (k : ℝ) + 1 / N < 0 := by linarith [hinvN]
      rw [abs_of_neg hneg]
      linarith [hinvN, hkR]
  · -- Body runner: reverse-triangle off the resonance center.
    have hcopvi : ¬ ((b : ℤ) ∣ v i * a) := fun h => hfill1 i hi (hcop.dvd_of_dvd_mul_right h)
    have hvi_pos : (0 : ℝ) < (v i : ℝ) := by exact_mod_cast hbody_pos i hi
    -- distance of v i · a/b to ℤ is ≥ 1/b at the integer m.
    have hdist : (1 : ℝ) / b ≤ |(v i : ℝ) * ((a : ℝ) / b) - m| := by
      have hne : (v i * a - m * b) ≠ 0 := by
        intro h; exact hcopvi ⟨m, by linarith [sub_eq_zero.mp h]⟩
      have h1 : (1 : ℝ) ≤ |((v i * a - m * b : ℤ) : ℝ)| := by
        rw [← Int.cast_abs]; exact_mod_cast Int.one_le_abs hne
      have key : (v i : ℝ) * ((a : ℝ) / b) - m = ((v i * a - m * b : ℤ) : ℝ) / b := by
        field_simp; push_cast; ring
      rw [key, abs_div, abs_of_pos hb0]
      gcongr
    -- the kick s' = v i /(N v⋆) ≤ 1/b − 1/N via the domination hypothesis.
    have hkick : (v i : ℝ) / (N * (v vstar : ℝ)) ≤ 1 / b - 1 / N := by
      have hdR : (b : ℝ) * ((v i : ℝ) + (v vstar : ℝ)) ≤ (N : ℝ) * (v vstar : ℝ) := by
        exact_mod_cast hdom i hi
      have e : 1 / (b : ℝ) - 1 / N - (v i : ℝ) / (N * (v vstar : ℝ))
          = ((N : ℝ) * v vstar - b * ((v i : ℝ) + v vstar)) / (b * (N * v vstar)) := by
        field_simp; ring
      have hden : (0 : ℝ) < (b : ℝ) * (N * v vstar) := by positivity
      have hnn : (0 : ℝ) ≤ 1 / (b : ℝ) - 1 / N - (v i : ℝ) / (N * (v vstar : ℝ)) := by
        rw [e]; exact div_nonneg (by linarith [hdR]) (le_of_lt hden)
      linarith [hnn]
    -- assemble via reverse triangle inequality.
    have hsplit : (v i : ℝ) * t - m
        = ((v i : ℝ) * ((a : ℝ) / b) - m) + (v i : ℝ) / (N * (v vstar : ℝ)) := by
      rw [ht]; ring
    have hs'pos : (0 : ℝ) ≤ (v i : ℝ) / (N * (v vstar : ℝ)) :=
      div_nonneg (le_of_lt hvi_pos) (by positivity)
    have hrt := abs_sub_abs_le_abs_sub ((v i : ℝ) * ((a : ℝ) / b) - m)
      (- ((v i : ℝ) / (N * (v vstar : ℝ))))
    rw [abs_neg, sub_neg_eq_add, ← hsplit, abs_of_nonneg hs'pos] at hrt
    linarith [hrt, hdist, hkick]

end LonelyRunner

#print axioms LonelyRunner.fill1_perturbation
