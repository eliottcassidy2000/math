import Mathlib

/-!
# Static `(4,4,8,8)` scalar matching on the reduced eight-branch circle

The finite pattern census has two geometric shapes up to dihedral symmetry:
two opposite-parity antipodal q4 fibers, while the two q8 rows are either the
two step-one edges or the two step-three edges partitioning the complement.
Selecting from each q8 row the endpoint with the parity of either q4 fiber
gives the complementary code `-2*c₄+c₈a+c₈b = 0 (mod 8)`.

After signing the two primitive q8 numerators to their common residue type
`t ∈ {1,3} (mod 8)` and the chosen q4 numerator to `t (mod 4)`, three bad
row inequalities yield a scalar wall of radius `3/14`, sharper than `2/7`.

The natural vertices here are branch residue classes, not runners.  This
quotient preserves the exact parallel-partition predicate but destroys phase
chronology; in particular, the separate dynamic proof that at least one of
the two q4-derived frequencies is nonzero is intentionally not claimed here.
-/

namespace LonelyRunner
namespace LRCPairTowerGapTwoFrequency

open LonelyRunner
open scoped Classical

noncomputable section

/-- Three-row scalar frequency attached to one q4 row and the two q8 rows. -/
def fourEightEightPhaseFrequency (a₄ a₈a a₈b : ℤ) : ℤ :=
  (-2 * a₄ + a₈a + a₈b) / 8

/-- Inside one parity class, avoiding a fixed mod-four fiber means shifting
by exactly two modulo four. -/
theorem four_dvd_sub_sub_two_of_same_parity
    (fiber point : ℤ)
    (hparity : point ≡ fiber [ZMOD 2])
    (houtside : ¬ point ≡ fiber [ZMOD 4]) :
    (4 : ℤ) ∣ point - fiber - 2 := by
  change point % 2 = fiber % 2 at hparity
  change ¬ point % 4 = fiber % 4 at houtside
  have hpointParity4 : point % 4 % 2 = point % 2 :=
    Int.emod_emod_of_dvd point (by norm_num : (2 : ℤ) ∣ 4)
  have hfiberParity4 : fiber % 4 % 2 = fiber % 2 :=
    Int.emod_emod_of_dvd fiber (by norm_num : (2 : ℤ) ∣ 4)
  have hpoint0 : (0 : ℤ) ≤ point % 4 :=
    Int.emod_nonneg point (by norm_num)
  have hpoint4 : point % 4 < 4 :=
    Int.emod_lt_of_pos point (by norm_num)
  have hfiber0 : (0 : ℤ) ≤ fiber % 4 :=
    Int.emod_nonneg fiber (by norm_num)
  have hfiber4 : fiber % 4 < 4 :=
    Int.emod_lt_of_pos fiber (by norm_num)
  have hresidue : point % 4 - fiber % 4 = 2 ∨
      point % 4 - fiber % 4 = -2 := by
    omega
  have hpointDecomp := Int.emod_add_mul_ediv point 4
  have hfiberDecomp := Int.emod_add_mul_ediv fiber 4
  rcases hresidue with hresidue | hresidue
  · refine ⟨point / 4 - fiber / 4, ?_⟩
    omega
  · refine ⟨point / 4 - fiber / 4 - 1, ?_⟩
    omega

/-- The two residues of a fixed parity outside one mod-four fiber are
complementary modulo eight. -/
theorem eight_dvd_complementary_parity_sum
    (fiber first second : ℤ)
    (hfirstParity : first ≡ fiber [ZMOD 2])
    (hsecondParity : second ≡ fiber [ZMOD 2])
    (hfirstOutside : ¬ first ≡ fiber [ZMOD 4])
    (hsecondOutside : ¬ second ≡ fiber [ZMOD 4])
    (hdistinct : ¬ first ≡ second [ZMOD 8]) :
    (8 : ℤ) ∣ -2 * fiber + first + second := by
  obtain ⟨firstShift, hfirstShift⟩ :=
    four_dvd_sub_sub_two_of_same_parity fiber first
      hfirstParity hfirstOutside
  obtain ⟨secondShift, hsecondShift⟩ :=
    four_dvd_sub_sub_two_of_same_parity fiber second
      hsecondParity hsecondOutside
  have hshiftDistinct : ¬ (2 : ℤ) ∣ firstShift - secondShift := by
    intro heven
    obtain ⟨halfShift, hhalfShift⟩ := heven
    apply hdistinct
    rw [Int.modEq_iff_dvd]
    refine ⟨-halfShift, ?_⟩
    omega
  have hshiftOdd : (firstShift - secondShift) % 2 = 1 := by
    rcases Int.emod_two_eq_zero_or_one (firstShift - secondShift) with hzero | hone
    · exact False.elim (hshiftDistinct (Int.dvd_iff_emod_eq_zero.mpr hzero))
    · exact hone
  have hshiftDecomp := Int.emod_add_mul_ediv (firstShift - secondShift) 2
  have hsumEven : (2 : ℤ) ∣ 1 + firstShift + secondShift := by
    refine ⟨secondShift + (firstShift - secondShift) / 2 + 1, ?_⟩
    omega
  obtain ⟨halfSum, hhalfSum⟩ := hsumEven
  refine ⟨halfSum, ?_⟩
  omega

/-- **Gap-two scalar matching.**  If one normalized q4 numerator and the two
q8 numerators have a common signed odd residue type, and their selected row
offsets form the complementary-parity code, then their three bad inequalities
force the combined integer frequency into radius `3/14`. -/
theorem frequency_bad_of_fourEightEight_matching
    (a₄ a₈a a₈b residue c₄ c₈a c₈b : ℤ) (u : ℝ)
    (hres₄ : (4 : ℤ) ∣ a₄ - residue)
    (hres₈a : (8 : ℤ) ∣ a₈a - residue)
    (hres₈b : (8 : ℤ) ∣ a₈b - residue)
    (hcode : (8 : ℤ) ∣ -2 * c₄ + c₈a + c₈b)
    (hbad₄ : ∃ n : ℤ,
      |(a₄ : ℝ) * ((u + (c₄ : ℝ)) / 4) - n| < 1 / 14)
    (hbad₈a : ∃ n : ℤ,
      |(a₈a : ℝ) * ((u + (c₈a : ℝ)) / 8) - n| < 1 / 14)
    (hbad₈b : ∃ n : ℤ,
      |(a₈b : ℝ) * ((u + (c₈b : ℝ)) / 8) - n| < 1 / 14) :
    ∃ n : ℤ,
      |(fourEightEightPhaseFrequency a₄ a₈a a₈b : ℝ) * u - n| <
        3 / 14 := by
  obtain ⟨n₄, hn₄⟩ := hbad₄
  obtain ⟨n₈a, hn₈a⟩ := hbad₈a
  obtain ⟨n₈b, hn₈b⟩ := hbad₈b
  obtain ⟨k₄, hk₄⟩ := hres₄
  obtain ⟨k₈a, hk₈a⟩ := hres₈a
  obtain ⟨k₈b, hk₈b⟩ := hres₈b
  obtain ⟨kcode, hkcode⟩ := hcode
  have hfrequencyDvd :
      (8 : ℤ) ∣ -2 * a₄ + a₈a + a₈b := by
    refine ⟨-k₄ + k₈a + k₈b, ?_⟩
    omega
  have hoffsetDvd :
      (8 : ℤ) ∣ -2 * a₄ * c₄ + a₈a * c₈a + a₈b * c₈b := by
    refine ⟨residue * kcode - k₄ * c₄ +
      k₈a * c₈a + k₈b * c₈b, ?_⟩
    have ha₄ : a₄ = 4 * k₄ + residue := by omega
    have ha₈a : a₈a = 8 * k₈a + residue := by omega
    have ha₈b : a₈b = 8 * k₈b + residue := by omega
    have hcode' : -2 * c₄ + c₈a + c₈b = 8 * kcode := hkcode
    rw [ha₄, ha₈a, ha₈b]
    calc
      -2 * (4 * k₄ + residue) * c₄ +
          (8 * k₈a + residue) * c₈a +
          (8 * k₈b + residue) * c₈b =
        8 * (-k₄ * c₄ + k₈a * c₈a + k₈b * c₈b) +
          residue * (-2 * c₄ + c₈a + c₈b) := by ring
      _ = 8 * (residue * kcode - k₄ * c₄ +
          k₈a * c₈a + k₈b * c₈b) := by rw [hcode']; ring
  let frequency : ℤ := (-2 * a₄ + a₈a + a₈b) / 8
  let offset : ℤ :=
    (-2 * a₄ * c₄ + a₈a * c₈a + a₈b * c₈b) / 8
  have hfrequency : frequency * 8 = -2 * a₄ + a₈a + a₈b :=
    Int.ediv_mul_cancel hfrequencyDvd
  have hoffset : offset * 8 =
      -2 * a₄ * c₄ + a₈a * c₈a + a₈b * c₈b :=
    Int.ediv_mul_cancel hoffsetDvd
  refine ⟨n₈a + n₈b - n₄ - offset, ?_⟩
  have heq :
      (frequency : ℝ) * u -
          ((n₈a + n₈b - n₄ - offset : ℤ) : ℝ) =
        -((a₄ : ℝ) * ((u + (c₄ : ℝ)) / 4) - n₄) +
          ((a₈a : ℝ) * ((u + (c₈a : ℝ)) / 8) - n₈a) +
          ((a₈b : ℝ) * ((u + (c₈b : ℝ)) / 8) - n₈b) := by
    push_cast
    have hfrequencyReal : (frequency : ℝ) * 8 =
        -2 * (a₄ : ℝ) + a₈a + a₈b := by exact_mod_cast hfrequency
    have hoffsetReal : (offset : ℝ) * 8 =
        -2 * (a₄ : ℝ) * c₄ + a₈a * c₈a + a₈b * c₈b := by
      exact_mod_cast hoffset
    have hfrequencyDiv : (frequency : ℝ) =
        (-2 * (a₄ : ℝ) + a₈a + a₈b) / 8 := by linarith
    have hoffsetDiv : (offset : ℝ) =
        (-2 * (a₄ : ℝ) * c₄ + a₈a * c₈a + a₈b * c₈b) / 8 := by
      linarith
    rw [hfrequencyDiv, hoffsetDiv]
    ring
  change |(frequency : ℝ) * u - _| < _
  rw [heq]
  calc
    |-((a₄ : ℝ) * ((u + (c₄ : ℝ)) / 4) - n₄) +
        ((a₈a : ℝ) * ((u + (c₈a : ℝ)) / 8) - n₈a) +
        ((a₈b : ℝ) * ((u + (c₈b : ℝ)) / 8) - n₈b)| ≤
      |(a₄ : ℝ) * ((u + (c₄ : ℝ)) / 4) - n₄| +
        |(a₈a : ℝ) * ((u + (c₈a : ℝ)) / 8) - n₈a| +
        |(a₈b : ℝ) * ((u + (c₈b : ℝ)) / 8) - n₈b| := by
          nlinarith [
            abs_add_le
              (-((a₄ : ℝ) * ((u + (c₄ : ℝ)) / 4) - n₄))
              ((a₈a : ℝ) * ((u + (c₈a : ℝ)) / 8) - n₈a),
            abs_add_le
              (-((a₄ : ℝ) * ((u + (c₄ : ℝ)) / 4) - n₄) +
                ((a₈a : ℝ) * ((u + (c₈a : ℝ)) / 8) - n₈a))
              ((a₈b : ℝ) * ((u + (c₈b : ℝ)) / 8) - n₈b),
            abs_neg ((a₄ : ℝ) * ((u + (c₄ : ℝ)) / 4) - n₄)]
    _ < 3 / 14 := by linarith

#print axioms eight_dvd_complementary_parity_sum
#print axioms four_dvd_sub_sub_two_of_same_parity
#print axioms frequency_bad_of_fourEightEight_matching

end
end LRCPairTowerGapTwoFrequency
end LonelyRunner
