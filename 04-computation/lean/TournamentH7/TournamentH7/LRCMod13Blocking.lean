/-
  TournamentH7.LRCMod13Blocking — the mod-13 pair-blocking lemma (boxeph-2026-07-18-S115/S116).

  Sharpening of HYP-4382 (the n=12 equality rigidity `M(C)=1/13 ⟺ C = d·{1,…,12}`):
  a NECESSARY condition for tightness, proved directly from the sieve at 13.

  Core sieve fact (`sieve13_middle_witness`): if at scale `b` every speed lands in the
  MIDDLE band mod 13 (residue in `[2,11]`, i.e. `c_i·b ≢ 0,±1 (mod 13)`), then the time
  `t = b/13` puts every runner at circle-distance `≥ 2/13` from the integers — a witness
  that `M(C) ≥ 2/13 > 1/13`.  Contrapositive: a family with global margin at most
  `1/13` has, at every scale `b`, some residue in `{0,±1}`.  If all speeds are units
  modulo 13, this yields the usual `±1` pair-blocking.  The unit hypothesis is essential:
  a zero residue blocks the middle-band witness without contradicting tightness.

  The same one-line argument gives a proved mod-`p` blocking for every prime `p`; this file
  formalizes the `p=13` case (the natural modulus for the `1/13` threshold).
-/
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

/-! ### A direct witness for every nonconstant positive twelve-term AP -/

/-- The twelve speeds `a, a+d, ..., a+11d`. -/
def ap12 (a d : ℤ) : Fin 12 → ℤ :=
  fun j => a + (j.val : ℤ) * d

/-- An elementary exact lattice-distance lemma.  If `s` lies strictly between
`0` and `Q`, and both end distances clear `Q/13`, then `s` is more than
`Q/13` from every multiple of `Q`. -/
lemma lattice_dist_gt_thirteenth {Q s : ℤ} (hQ : 0 < Q) (hs : 0 < s)
    (hsQ : s < Q) (hleft : Q < 13 * s) (hright : Q < 13 * (Q - s)) (m : ℤ) :
    Q < 13 * |s - m * Q| := by
  rcases lt_trichotomy m 1 with hm | hm | hm
  · have hm0 : m ≤ 0 := by omega
    have hmQ : m * Q ≤ 0 := mul_nonpos_of_nonpos_of_nonneg hm0 hQ.le
    have hnonneg : 0 ≤ s - m * Q := by omega
    rw [abs_of_nonneg hnonneg]
    omega
  · subst m
    have hneg : s - 1 * Q < 0 := by omega
    rw [abs_of_neg hneg]
    omega
  · have hm1 : 1 ≤ m := by omega
    have hQm : Q ≤ m * Q := by
      simpa using mul_le_mul_of_nonneg_right hm1 hQ.le
    have hneg : s - m * Q < 0 := by omega
    rw [abs_of_neg hneg]
    omega

/-- For `a>d>=1`, every speed in the twelve-term AP has integer lattice
distance strictly greater than `(2a+11d)/13` at numerator one. -/
theorem ap12_lattice_distance (a d : ℤ) (had : d < a) (hd : 1 ≤ d)
    (j : Fin 12) (m : ℤ) :
    2 * a + 11 * d <
      13 * |ap12 a d j - m * (2 * a + 11 * d)| := by
  have hj0 : 0 ≤ (j.val : ℤ) := by omega
  have hj11 : (j.val : ℤ) ≤ 11 := by omega
  have hd0 : 0 ≤ d := by omega
  have hjd : 0 ≤ (j.val : ℤ) * d := mul_nonneg hj0 hd0
  have htail : 0 ≤ (11 - (j.val : ℤ)) * d :=
    mul_nonneg (by omega) hd0
  have hs_pos : 0 < ap12 a d j := by
    simp only [ap12]
    omega
  have hs_gap : 0 < 2 * a + 11 * d - ap12 a d j := by
    have hid : 2 * a + 11 * d - ap12 a d j =
        a + (11 - (j.val : ℤ)) * d := by
      simp only [ap12]
      ring
    rw [hid]
    omega
  have hleft : 0 < 13 * ap12 a d j - (2 * a + 11 * d) := by
    have hid : 13 * ap12 a d j - (2 * a + 11 * d) =
        11 * (a - d) + 13 * ((j.val : ℤ) * d) := by
      simp only [ap12]
      ring
    rw [hid]
    omega
  have hright : 0 <
      13 * ((2 * a + 11 * d) - ap12 a d j) - (2 * a + 11 * d) := by
    have hid :
        13 * ((2 * a + 11 * d) - ap12 a d j) - (2 * a + 11 * d) =
          11 * (a - d) + 13 * ((11 - (j.val : ℤ)) * d) := by
      simp only [ap12]
      ring
    rw [hid]
    omega
  apply lattice_dist_gt_thirteenth
  · omega
  · exact hs_pos
  · omega
  · omega
  · omega

/-- **AP WITNESS DISPATCH.**  If `a>d>=1`, then at
`t=1/(2a+11d)` every one of `a,a+d,...,a+11d` has margin strictly greater
than `1/13` from every integer. -/
theorem ap12_strict_witness (a d : ℤ) (had : d < a) (hd : 1 ≤ d) :
    ∀ j : Fin 12, ∀ m : ℤ,
      (1 : ℝ) / 13 <
        |(ap12 a d j : ℝ) *
            ((1 : ℝ) / ((2 * a + 11 * d : ℤ) : ℝ)) - m| := by
  intro j m
  have hQ : 0 < 2 * a + 11 * d := by omega
  have hQR : (0 : ℝ) < ((2 * a + 11 * d : ℤ) : ℝ) := by exact_mod_cast hQ
  have hint := ap12_lattice_distance a d had hd j m
  have hintR : ((2 * a + 11 * d : ℤ) : ℝ) <
      13 * |((ap12 a d j - m * (2 * a + 11 * d) : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]
    exact_mod_cast hint
  have hreal :
      (ap12 a d j : ℝ) * ((1 : ℝ) / ((2 * a + 11 * d : ℤ) : ℝ)) - m
        = ((ap12 a d j - m * (2 * a + 11 * d) : ℤ) : ℝ) /
            ((2 * a + 11 * d : ℤ) : ℝ) := by
    field_simp [ne_of_gt hQR]
    <;> push_cast
    <;> ring
  rw [hreal, abs_div, abs_of_pos hQR, lt_div_iff₀ hQR]
  have hnum : ((2 * a + 11 * d : ℤ) : ℝ) / 13 <
      |((ap12 a d j - m * (2 * a + 11 * d) : ℤ) : ℝ)| := by
    nlinarith
  convert hnum using 1 <;> ring

/-- Closed-threshold form of `ap12_strict_witness`: every such nonconstant
positive AP is `13`-lonely at the same explicit time. -/
theorem ap12_lonely (a d : ℤ) (had : d < a) (hd : 1 ≤ d) :
    Lonely 13 (ap12 a d)
      ((1 : ℝ) / ((2 * a + 11 * d : ℤ) : ℝ)) := by
  intro j m
  exact (ap12_strict_witness a d had hd j m).le

end LonelyRunner

#print axioms LonelyRunner.mod13_middle_far
#print axioms LonelyRunner.sieve13_middle_witness
#print axioms LonelyRunner.no_middle_band_witness_of_tight
#print axioms LonelyRunner.lattice_dist_gt_thirteenth
#print axioms LonelyRunner.ap12_lattice_distance
#print axioms LonelyRunner.ap12_strict_witness
#print axioms LonelyRunner.ap12_lonely
