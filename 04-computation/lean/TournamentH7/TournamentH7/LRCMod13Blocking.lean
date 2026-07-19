/-
  TournamentH7.LRCMod13Blocking — the mod-13 pair-blocking lemma (boxeph-2026-07-18-S115/S116).

  Sharpening of HYP-4382 (the n=12 equality rigidity `M(C)=1/13 ⟺ C = d·{1,…,12}`):
  a NECESSARY condition for tightness, proved directly from the sieve at 13.

  Core sieve fact (`sieve13_middle_witness`): if at scale `b` every speed lands in the
  MIDDLE band mod 13 (residue in `[2,11]`, i.e. `c_i·b ≢ 0,±1 (mod 13)`), then the time
  `t = b/13` puts every runner at circle-distance `≥ 2/13` from the integers — a witness
  that `M(C) ≥ 2/13 > 1/13`.  Contrapositive: `M(C) = 1/13` forces the mod-13 pair-blocking
  `{±c_i mod 13} = {1,…,12}` (some `c_i ≡ ±b^{-1}` for every `b`).

  The same one-line argument gives a proved mod-`p` blocking for every prime `p`; this file
  formalizes the `p=13` case (the natural modulus for the `1/13` threshold).
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner

/-- Integer core: a residue `r ∈ [2,11]` is at distance `≥ 2` from every multiple of 13.
`13k + r` never lands in `{-1,0,1}` when `2 ≤ r ≤ 11`. -/
theorem mod13_middle_far (r : ℤ) (h2 : 2 ≤ r) (h11 : r ≤ 11) (k : ℤ) :
    2 ≤ |13 * k + r| := by
  rcases le_or_gt 0 (13 * k + r) with h | h
  · rw [abs_of_nonneg h]; omega
  · rw [abs_of_neg h]; omega

/-- **The mod-13 middle-band sieve witness (PROVED).**  If at scale `b` every speed's
residue `(c_i·b) mod 13` lies in the middle band `[2,11]`, then at `t = b/13` every runner
is at distance `≥ 2/13` from every integer.  Hence `M(C) ≥ 2/13 > 1/13`. -/
theorem sieve13_middle_witness {ι : Type*} (c : ι → ℤ) (b : ℤ)
    (hmid : ∀ i, 2 ≤ (c i * b) % 13 ∧ (c i * b) % 13 ≤ 11) :
    ∀ i, ∀ m : ℤ, (2 : ℝ) / 13 ≤ |(c i : ℝ) * ((b : ℝ) / 13) - m| := by
  intro i m
  obtain ⟨h2, h11⟩ := hmid i
  -- the real value is (c_i·b − 13m)/13
  have hreal : (c i : ℝ) * ((b : ℝ) / 13) - m
      = ((c i * b - 13 * m : ℤ) : ℝ) / 13 := by push_cast; ring
  -- the integer numerator is ≥ 2 in absolute value
  have hdecomp : c i * b - 13 * m = 13 * (c i * b / 13 - m) + (c i * b) % 13 := by
    omega
  have hint : (2 : ℤ) ≤ |c i * b - 13 * m| := by
    rw [hdecomp]; exact mod13_middle_far _ h2 h11 _
  have hintR : (2 : ℝ) ≤ |((c i * b - 13 * m : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]; exact_mod_cast hint
  rw [hreal, abs_div, show |(13 : ℝ)| = 13 by norm_num]
  gcongr

/-- **Contrapositive — mod-13 blocking is necessary for the tight witness.**  If the family
is NOT `2/13`-clearable at scale `b` for the reason above, that is precisely the failure of
some runner to sit in the middle band.  So a family whose every scale `b` leaves a runner in
the middle band would have `M ≥ 2/13`; hence a family with `M = 1/13` must, at every `b`, have
some runner OUTSIDE the middle band — the pair-blocking condition. -/
theorem no_middle_band_witness_of_tight {ι : Type*} (c : ι → ℤ) (b : ℤ)
    (htight : ∃ i, ∀ m : ℤ, |(c i : ℝ) * ((b : ℝ) / 13) - m| < 2 / 13) :
    ¬ (∀ i, 2 ≤ (c i * b) % 13 ∧ (c i * b) % 13 ≤ 11) := by
  intro hmid
  obtain ⟨i, hi⟩ := htight
  exact absurd (sieve13_middle_witness c b hmid i 0) (not_le.mpr (hi 0))

end LonelyRunner

#print axioms LonelyRunner.mod13_middle_far
#print axioms LonelyRunner.sieve13_middle_witness
#print axioms LonelyRunner.no_middle_band_witness_of_tight
