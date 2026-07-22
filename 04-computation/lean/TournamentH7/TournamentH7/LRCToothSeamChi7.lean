/-
  TournamentH7.LRCToothSeamChi7 -- arithmetic kernel for THM-1156.

  For the open radius-1/14 teeth

      ((14m-1)/(14a), (14m+1)/(14a)),
      ((14n-1)/(14b), (14n+1)/(14b)),

  the oriented seam from the right end of the first tooth to the left end of
  the second has numerator

      N = 14(an-bm) - (a+b).

  This file proves the endpoint identity, the exact Bezout criterion for an
  abutment, the gcd quantum on the positive (uncovered-gap) side, the
  primitive mod-14 complement law, and the open-endpoint third-support
  implication.  The negative sign of N is deliberately not identified with
  an intersection length: nested teeth have a containment plateau, so such a
  geometric use needs an additional facing-endpoint/non-containment premise.

  Kernel-pure: no `sorry`, no `admit`, and no `native_decide`.
-/
import Mathlib

namespace LonelyRunner
namespace ToothSeamChi7

/-- Left endpoint of the `m`-th open radius-`1/14` tooth of speed `a`. -/
def toothLeft (a m : ℤ) : ℚ := (14 * m - 1) / (14 * a)

/-- Right endpoint of the `m`-th open radius-`1/14` tooth of speed `a`. -/
def toothRight (a m : ℤ) : ℚ := (14 * m + 1) / (14 * a)

/-- Cleared numerator of the oriented seam `toothLeft b n - toothRight a m`. -/
def seamNumerator (a b m n : ℤ) : ℤ := 14 * (a * n - b * m) - (a + b)

/-- The exact cleared endpoint identity. -/
theorem seam_endpoint_identity (a b m n : ℤ) (ha : a ≠ 0) (hb : b ≠ 0) :
    toothLeft b n - toothRight a m =
      (seamNumerator a b m n : ℚ) / (14 * a * b) := by
  simp only [toothLeft, toothRight, seamNumerator]
  push_cast
  field_simp
  ring

/-- Exact two-tooth abutment is a single linear Diophantine condition. -/
theorem exists_exact_abutment_iff (a b : ℤ) :
    (∃ m n : ℤ, seamNumerator a b m n = 0) ↔
      (14 * (Int.gcd a b : ℤ)) ∣ a + b := by
  constructor
  · rintro ⟨m, n, hzero⟩
    have heq : 14 * (a * n - b * m) = a + b := sub_eq_zero.mp hzero
    have hga : (Int.gcd a b : ℤ) ∣ a := Int.gcd_dvd_left a b
    have hgb : (Int.gcd a b : ℤ) ∣ b := Int.gcd_dvd_right a b
    have hlin : (Int.gcd a b : ℤ) ∣ a * n - b * m :=
      dvd_sub (dvd_mul_of_dvd_left hga n) (dvd_mul_of_dvd_left hgb m)
    obtain ⟨c, hc⟩ := hlin
    refine ⟨c, ?_⟩
    rw [← heq, hc]
    ring
  · rintro ⟨c, hc⟩
    let A : ℤ := Int.gcdA a b
    let B : ℤ := Int.gcdB a b
    have hbez : (Int.gcd a b : ℤ) = a * A + b * B := by
      exact Int.gcd_eq_gcd_ab a b
    refine ⟨-B * c, A * c, ?_⟩
    simp only [seamNumerator]
    rw [hc, hbez]
    ring

/-- Every seam numerator is congruent to `-(a+b)` modulo `14*gcd(a,b)`. -/
theorem seam_congruence (a b m n : ℤ) :
    (14 * (Int.gcd a b : ℤ)) ∣ seamNumerator a b m n + (a + b) := by
  have hga : (Int.gcd a b : ℤ) ∣ a := Int.gcd_dvd_left a b
  have hgb : (Int.gcd a b : ℤ) ∣ b := Int.gcd_dvd_right a b
  have hlin : (Int.gcd a b : ℤ) ∣ a * n - b * m :=
    dvd_sub (dvd_mul_of_dvd_left hga n) (dvd_mul_of_dvd_left hgb m)
  obtain ⟨c, hc⟩ := hlin
  refine ⟨c, ?_⟩
  simp only [seamNumerator]
  rw [hc]
  ring

/-- Every multiple of the gcd occurs as a determinant `a*n-b*m`.  This is
the attainability half of the directed seam quantum. -/
theorem exists_determinant_eq_gcd_mul (a b z : ℤ) :
    ∃ m n : ℤ, a * n - b * m = (Int.gcd a b : ℤ) * z := by
  let A : ℤ := Int.gcdA a b
  let B : ℤ := Int.gcdB a b
  have hbez : (Int.gcd a b : ℤ) = a * A + b * B := Int.gcd_eq_gcd_ab a b
  refine ⟨-B * z, A * z, ?_⟩
  rw [hbez]
  ring

/-- Division-free floor form of the one-sided strict-penetration quantum.
If `q=floor((s-1)/M)` is supplied through its two defining inequalities and
`r=s-Mq`, then every strict lattice penetration `Mz<s` has depth at least
`r`. -/
theorem strict_penetration_floor_quantum (M s q z r : ℤ) (hM : 0 < M)
    (hqlo : M * q ≤ s - 1) (hqup : s - 1 < M * (q + 1))
    (hr : r = s - M * q) (hz : M * z < s) :
    r ≤ s - M * z := by
  have hzint : M * z ≤ s - 1 := by omega
  have hzlt : M * z < M * (q + 1) := hzint.trans_lt hqup
  have hzq : z ≤ q := by
    have : z < q + 1 := (Int.mul_lt_mul_left hM).mp hzlt
    omega
  have hmul : M * z ≤ M * q := mul_le_mul_of_nonneg_left hzq hM.le
  rw [hr]
  linarith

/-- The floor quantum is attained by an actual pair of tooth indices. -/
theorem exists_seam_at_floor_quantum (a b q r : ℤ)
    (hr : r = a + b - 14 * (Int.gcd a b : ℤ) * q) :
    ∃ m n : ℤ, seamNumerator a b m n = -r := by
  obtain ⟨m, n, hdet⟩ := exists_determinant_eq_gcd_mul a b q
  refine ⟨m, n, ?_⟩
  simp only [seamNumerator]
  rw [hdet, hr]
  ring

/-- For an abutment-compatible pair, all oriented seam numerators lie on the
`14*gcd` lattice. -/
theorem seam_quantized_of_compatible (a b m n : ℤ)
    (hcompat : (14 * (Int.gcd a b : ℤ)) ∣ a + b) :
    (14 * (Int.gcd a b : ℤ)) ∣ seamNumerator a b m n := by
  obtain ⟨c, hc⟩ := seam_congruence a b m n
  obtain ⟨d, hd⟩ := hcompat
  refine ⟨c - d, ?_⟩
  linear_combination hc - hd

/-- One-sided strict seam quantum.  For positive speeds, a positive seam in
an abutment-compatible projective pair has cleared numerator at least
`14*gcd(a,b)`.  Dividing by `14ab` gives the geometric lower bound
`gcd(a,b)/(ab)`. -/
theorem positive_seam_quantum (a b m n : ℤ) (ha : 0 < a)
    (hcompat : (14 * (Int.gcd a b : ℤ)) ∣ a + b)
    (hpos : 0 < seamNumerator a b m n) :
    14 * (Int.gcd a b : ℤ) ≤ seamNumerator a b m n := by
  have hgpos : 0 < (Int.gcd a b : ℤ) := by
    exact_mod_cast Int.gcd_pos_of_ne_zero_left b ha.ne'
  obtain ⟨c, hc⟩ := seam_quantized_of_compatible a b m n hcompat
  have hcpos : 0 < c := by
    rw [hc] at hpos
    nlinarith
  rw [hc]
  nlinarith

/-- The rational form of the positive seam quantum. -/
theorem positive_endpoint_gap_quantum (a b m n : ℤ) (ha : 0 < a) (hb : 0 < b)
    (hcompat : (14 * (Int.gcd a b : ℤ)) ∣ a + b)
    (hgap : toothRight a m < toothLeft b n) :
    ((Int.gcd a b : ℤ) : ℚ) / (a * b) ≤
      toothLeft b n - toothRight a m := by
  have ha0 : a ≠ 0 := ne_of_gt ha
  have hb0 : b ≠ 0 := ne_of_gt hb
  have habpos : (0 : ℚ) < 14 * a * b := by positivity
  have hnumpos : 0 < seamNumerator a b m n := by
    have hdiff : 0 < toothLeft b n - toothRight a m := sub_pos.mpr hgap
    rw [seam_endpoint_identity a b m n ha0 hb0] at hdiff
    have hnumQ : (0 : ℚ) < (seamNumerator a b m n : ℚ) := by
      rcases div_pos_iff.mp hdiff with hsame | hsame
      · exact hsame.1
      · exact False.elim ((not_lt_of_ge habpos.le) hsame.2)
    exact_mod_cast hnumQ
  have hquant := positive_seam_quantum a b m n ha hcompat hnumpos
  rw [seam_endpoint_identity a b m n ha0 hb0]
  calc
    ((Int.gcd a b : ℤ) : ℚ) / (a * b) =
        (14 * ((Int.gcd a b : ℤ) : ℚ)) / (14 * a * b) := by
          field_simp
    _ ≤ (seamNumerator a b m n : ℚ) / (14 * a * b) := by
      apply (div_le_div_iff_of_pos_right habpos).2
      exact_mod_cast hquant

/-- In primitive coordinates, exact compatibility is complementarity mod 14. -/
theorem primitive_mod14_complement {A B : ℕ} (hcompat : 14 ∣ A + B) :
    (A : ZMod 14) + B = 0 := by
  have hzero : ((A + B : ℕ) : ZMod 14) = 0 := by
    rw [CharP.cast_eq_zero_iff (ZMod 14) 14]
    exact hcompat
  simpa using hzero

/-- Coprime complementary residues modulo seven are both nonzero. -/
theorem primitive_mod7_nonzero {A B : ℕ} (hcop : Nat.Coprime A B)
    (hseven : 7 ∣ A + B) :
    (A : ZMod 7) ≠ 0 ∧ (B : ZMod 7) ≠ 0 := by
  have hnotA : ¬ 7 ∣ A := by
    intro hA
    have hB : 7 ∣ B := by omega
    exact (Nat.not_coprime_of_dvd_of_dvd (by norm_num : 1 < 7) hA hB) hcop
  have hnotB : ¬ 7 ∣ B := by
    intro hB
    have hA : 7 ∣ A := by omega
    exact (Nat.not_coprime_of_dvd_of_dvd (by norm_num : 1 < 7) hA hB) hcop
  constructor
  · rw [Ne, CharP.cast_eq_zero_iff (ZMod 7) 7]
    exact hnotA
  · rw [Ne, CharP.cast_eq_zero_iff (ZMod 7) 7]
    exact hnotB

local instance : Fact (Nat.Prime 7) := ⟨by norm_num⟩

/-- The canonical quadratic sign modulo seven. -/
def chi7 (x : ℤ) : ℤ := legendreSym 7 x

/-- Negation flips the nonzero quadratic sign modulo seven. -/
theorem chi7_neg (x : ℤ) : chi7 (-x) = -chi7 x := by
  rw [chi7, chi7, legendreSym.at_neg (p := 7) (by norm_num)]
  rw [ZMod.χ₄_nat_three_mod_four (by norm_num : 7 % 4 = 3)]
  ring

/-- The primitive complementary pair crosses the canonical `chi7`
bipartition. -/
theorem primitive_chi7_bipartition {A B : ℕ} (hcop : Nat.Coprime A B)
    (hseven : 7 ∣ A + B) :
    chi7 (B : ℤ) = -chi7 (A : ℤ) := by
  have _hnonzero := primitive_mod7_nonzero hcop hseven
  obtain ⟨k, hk⟩ := hseven
  have hB : (B : ℤ) = -(A : ℤ) + 7 * (k : ℤ) := by
    have hkZ : (A : ℤ) + (B : ℤ) = 7 * (k : ℤ) := by
      exact_mod_cast hk
    omega
  have hmod : (B : ℤ) % 7 = (-(A : ℤ)) % 7 := by
    rw [hB]
    simp [Int.add_emod, Int.mul_emod]
  calc
    chi7 (B : ℤ) = chi7 (-(A : ℤ)) := by
      rw [chi7, chi7]
      calc
        legendreSym 7 (B : ℤ) = legendreSym 7 ((B : ℤ) % 7) := legendreSym.mod 7 _
        _ = legendreSym 7 ((-(A : ℤ)) % 7) := by rw [hmod]
        _ = legendreSym 7 (-(A : ℤ)) := (legendreSym.mod 7 _).symm
    _ = -chi7 (A : ℤ) := chi7_neg (A : ℤ)

/-- Membership in one selected open tooth. -/
def InOpenTooth (a m : ℤ) (t : ℚ) : Prop := toothLeft a m < t ∧ t < toothRight a m

/-- An exact abutment point belongs to neither of the two selected open teeth. -/
theorem exact_abutment_omits_seam (a b m n : ℤ)
    (habut : toothRight a m = toothLeft b n) :
    ¬ InOpenTooth a m (toothRight a m) ∧
      ¬ InOpenTooth b n (toothRight a m) := by
  constructor
  · intro h
    exact (lt_irrefl _ h.2)
  · intro h
    rw [habut] at h
    exact (lt_irrefl _ h.1)

/-- Qualitative open-cover consequence: if a family covers an exact seam
point, its owner cannot be either of the two abutting teeth.  This is the
formal third-support obligation; the extra support may be a tooth of a third
comb in the LRC application. -/
theorem exact_abutment_needs_third_support {I : Type*} (a b m n : ℤ)
    (owner : I → ℤ × ℤ) (i j : I)
    (hi : owner i = (a, m)) (hj : owner j = (b, n))
    (habut : toothRight a m = toothLeft b n)
    (hcover : ∃ k : I, InOpenTooth (owner k).1 (owner k).2 (toothRight a m)) :
    ∃ k : I, k ≠ i ∧ k ≠ j ∧
      InOpenTooth (owner k).1 (owner k).2 (toothRight a m) := by
  obtain ⟨k, hk⟩ := hcover
  refine ⟨k, ?_, ?_, hk⟩
  · intro hki
    subst k
    rw [hi] at hk
    exact (exact_abutment_omits_seam a b m n habut).1 hk
  · intro hkj
    subst k
    rw [hj] at hk
    exact (exact_abutment_omits_seam a b m n habut).2 hk

end ToothSeamChi7
end LonelyRunner

#print axioms LonelyRunner.ToothSeamChi7.seam_endpoint_identity
#print axioms LonelyRunner.ToothSeamChi7.exists_exact_abutment_iff
#print axioms LonelyRunner.ToothSeamChi7.strict_penetration_floor_quantum
#print axioms LonelyRunner.ToothSeamChi7.exists_seam_at_floor_quantum
#print axioms LonelyRunner.ToothSeamChi7.positive_endpoint_gap_quantum
#print axioms LonelyRunner.ToothSeamChi7.primitive_chi7_bipartition
#print axioms LonelyRunner.ToothSeamChi7.exact_abutment_needs_third_support
