import TournamentH7.LRCPairTowerGapTwo
import TournamentH7.LRCPairTowerGapTwoFrequency

/-!
# Dynamic producer for the `(4,4,8,8)` parallel-class wall

An exact four-row failure partition forces each saturated q8 row to meet both
parity sheets.  Choosing the endpoint in the parity sheet of either q4 row
produces the complementary mod-eight offset code and the three normalized bad
witnesses required by `frequency_bad_of_fourEightEight_matching`.

The two q4-derived scalar frequencies are then considered jointly.  Their q4
numerators cannot coincide, so both frequencies cannot vanish.  Either one
crosses the sharp LRC(9) escape threshold, or every surviving frequency is a
concrete small signed relation.

Tournament-analysis audit: vertices are branch residues on the parallel-class
circle.  The binary relation is same-parity q8 compatibility, and the two q4
choices form the tie path.  Quotienting to runners destroys both the branch
offset code and the phasewise bad witnesses.  The challenged assumption is
that an arbitrary primitive signing already has a common q8 residue type.  The
parallel rectangle below proves instead that the two types agree only after a
possible sign flip; no common signing is assumed from denominator data alone.
-/

namespace LonelyRunner
namespace LRCPairTowerGapTwoProducer

open LonelyRunner LRCPairTowerGapTwo LRCPairTowerGapTwoFrequency
open scoped Classical

noncomputable section

def normalizedBadBranches (numerator denominator g : ℤ) (u : ℝ) : Finset ℤ :=
  (Finset.Ico (0 : ℤ) g).filter fun branch =>
    ∃ integer : ℤ,
      |(numerator : ℝ) *
        ((u + (branch : ℝ)) / (denominator : ℝ)) - integer| < 1 / 14

structure SignedReducedNormalForm
    (δ g denominator numerator : ℤ) where
  divisor : ℤ
  sign : ℤ
  divisor_ne : divisor ≠ 0
  delta_eq : δ = divisor * (sign * numerator)
  modulus_eq : g = divisor * denominator
  sign_eq : sign = 1 ∨ sign = -1
  coprime : IsCoprime numerator denominator

theorem detunedBadBranches_eq_normalizedBadBranches
    (δ g divisor denominator sign numerator : ℤ) (u : ℝ)
    (hdivisor : divisor ≠ 0) (hdenominator : denominator ≠ 0)
    (hδ : δ = divisor * (sign * numerator))
    (hg : g = divisor * denominator)
    (hsign : sign = 1 ∨ sign = -1) :
    detunedBadBranches δ g u =
      normalizedBadBranches numerator denominator g u := by
  ext branch
  simp only [detunedBadBranches, normalizedBadBranches, Finset.mem_filter]
  constructor
  · rintro ⟨hbranch, integer, hinteger⟩
    refine ⟨hbranch, ?_⟩
    rcases hsign with rfl | rfl
    · refine ⟨integer, ?_⟩
      have heq :
          (δ : ℝ) * ((u + (branch : ℝ)) / (g : ℝ)) - integer =
            (numerator : ℝ) *
              ((u + (branch : ℝ)) / (denominator : ℝ)) - integer := by
        rw [hδ, hg]
        push_cast
        field_simp [show (divisor : ℝ) ≠ 0 by exact_mod_cast hdivisor,
          show (denominator : ℝ) ≠ 0 by exact_mod_cast hdenominator]
      rwa [heq] at hinteger
    · refine ⟨-integer, ?_⟩
      have heq :
          (δ : ℝ) * ((u + (branch : ℝ)) / (g : ℝ)) - integer =
            -((numerator : ℝ) *
              ((u + (branch : ℝ)) / (denominator : ℝ)) - (-integer : ℤ)) := by
        rw [hδ, hg]
        push_cast
        field_simp [show (divisor : ℝ) ≠ 0 by exact_mod_cast hdivisor,
          show (denominator : ℝ) ≠ 0 by exact_mod_cast hdenominator]
        ring
      rw [heq, abs_neg] at hinteger
      exact hinteger
  · rintro ⟨hbranch, integer, hinteger⟩
    refine ⟨hbranch, ?_⟩
    rcases hsign with rfl | rfl
    · refine ⟨integer, ?_⟩
      have heq :
          (δ : ℝ) * ((u + (branch : ℝ)) / (g : ℝ)) - integer =
            (numerator : ℝ) *
              ((u + (branch : ℝ)) / (denominator : ℝ)) - integer := by
        rw [hδ, hg]
        push_cast
        field_simp [show (divisor : ℝ) ≠ 0 by exact_mod_cast hdivisor,
          show (denominator : ℝ) ≠ 0 by exact_mod_cast hdenominator]
      rwa [heq]
    · refine ⟨-integer, ?_⟩
      have heq :
          (δ : ℝ) * ((u + (branch : ℝ)) / (g : ℝ)) - (-integer : ℤ) =
            -((numerator : ℝ) *
              ((u + (branch : ℝ)) / (denominator : ℝ)) - integer) := by
        rw [hδ, hg]
        push_cast
        field_simp [show (divisor : ℝ) ≠ 0 by exact_mod_cast hdivisor,
          show (denominator : ℝ) ≠ 0 by exact_mod_cast hdenominator]
        ring
      rw [heq, abs_neg]
      exact hinteger

theorem normalizedBadBranches_mem_of_modEq
    (numerator denominator g : ℤ) (u : ℝ)
    (hdenominator : denominator ≠ 0)
    {branch branch' : ℤ}
    (hbranch : branch ∈ normalizedBadBranches numerator denominator g u)
    (hbranch'Ico : branch' ∈ Finset.Ico (0 : ℤ) g)
    (hmod : branch ≡ branch' [ZMOD denominator]) :
    branch' ∈ normalizedBadBranches numerator denominator g u := by
  obtain ⟨integer, hinteger⟩ := (Finset.mem_filter.mp hbranch).2
  obtain ⟨shift, hshift⟩ := hmod.dvd
  rw [normalizedBadBranches, Finset.mem_filter]
  refine ⟨hbranch'Ico, integer + numerator * shift, ?_⟩
  have hbranchEq : branch' = branch + denominator * shift := by omega
  have hbranchEqReal : (branch' : ℝ) =
      (branch : ℝ) + (denominator : ℝ) * (shift : ℝ) := by
    exact_mod_cast hbranchEq
  have heq :
      (numerator : ℝ) *
          ((u + (branch' : ℝ)) / (denominator : ℝ)) -
          ((integer + numerator * shift : ℤ) : ℝ) =
        (numerator : ℝ) *
          ((u + (branch : ℝ)) / (denominator : ℝ)) - integer := by
    rw [hbranchEqReal]
    push_cast
    field_simp [show (denominator : ℝ) ≠ 0 by exact_mod_cast hdenominator]
    ring
  rwa [heq]

theorem normalized_bad_witness_of_detuned_mem
    (δ g denominator numerator : ℤ) (u : ℝ)
    (hnorm : SignedReducedNormalForm δ g denominator numerator)
    (hdenominator : denominator ≠ 0) {branch : ℤ}
    (hbranch : branch ∈ detunedBadBranches δ g u) :
    ∃ integer : ℤ,
      |(numerator : ℝ) *
          ((u + (branch : ℝ)) / (denominator : ℝ)) - integer| < 1 / 14 := by
  have hrowEq := detunedBadBranches_eq_normalizedBadBranches
    δ g hnorm.divisor denominator hnorm.sign numerator u
      hnorm.divisor_ne hdenominator hnorm.delta_eq hnorm.modulus_eq
      hnorm.sign_eq
  have hnormalized : branch ∈
      normalizedBadBranches numerator denominator g u := by
    rw [← hrowEq]
    exact hbranch
  exact (Finset.mem_filter.mp hnormalized).2

theorem normalizedQEight_pair_modEq_of_same_parity
    (numerator g : ℤ) (u : ℝ)
    (hcoprime : IsCoprime numerator 8)
    {first second : ℤ}
    (hfirst : first ∈ normalizedBadBranches numerator 8 g u)
    (hsecond : second ∈ normalizedBadBranches numerator 8 g u)
    (hparity : first ≡ second [ZMOD 2]) :
    first ≡ second [ZMOD 8] := by
  obtain ⟨integerFirst, hnearFirst⟩ := (Finset.mem_filter.mp hfirst).2
  obtain ⟨integerSecond, hnearSecond⟩ := (Finset.mem_filter.mp hsecond).2
  let defect : ℤ :=
    numerator * (first - second) - 8 * (integerFirst - integerSecond)
  have hdefectReal :
      (defect : ℝ) / 8 =
        ((numerator : ℝ) * ((u + (first : ℝ)) / 8) - integerFirst) -
        ((numerator : ℝ) * ((u + (second : ℝ)) / 8) - integerSecond) := by
    dsimp [defect]
    push_cast
    ring
  have hdefectAbs : |(defect : ℝ)| < 2 := by
    have hquotient : |(defect : ℝ) / 8| < 1 / 7 := by
      rw [hdefectReal]
      calc
        |_ - _| ≤
            |(numerator : ℝ) * ((u + (first : ℝ)) / 8) - integerFirst| +
            |(numerator : ℝ) * ((u + (second : ℝ)) / 8) - integerSecond| :=
          abs_sub _ _
        _ < 1 / 14 + 1 / 14 := add_lt_add hnearFirst hnearSecond
        _ = 1 / 7 := by norm_num
    rw [abs_div, abs_of_pos (by norm_num : (0 : ℝ) < 8)] at hquotient
    nlinarith [hquotient]
  have hdefectInt : |defect| < 2 := by exact_mod_cast hdefectAbs
  have hdiffEven : (2 : ℤ) ∣ first - second := by
    have := hparity.dvd
    simpa only [neg_sub] using dvd_neg.mpr this
  have hdefectEven : (2 : ℤ) ∣ defect := by
    apply dvd_sub
    · exact dvd_mul_of_dvd_right hdiffEven numerator
    · exact dvd_mul_of_dvd_left (by norm_num : (2 : ℤ) ∣ 8) _
  obtain ⟨halfDefect, hhalfDefect⟩ := hdefectEven
  have hdefectBounds := abs_lt.mp hdefectInt
  have hdefectZero : defect = 0 := by omega
  have h8dvdProduct : (8 : ℤ) ∣ numerator * (first - second) := by
    refine ⟨integerFirst - integerSecond, ?_⟩
    dsimp [defect] at hdefectZero
    omega
  have h8dvdDiff : (8 : ℤ) ∣ first - second :=
    hcoprime.symm.dvd_of_dvd_mul_left h8dvdProduct
  rw [Int.modEq_iff_dvd]
  simpa only [neg_sub] using dvd_neg.mpr h8dvdDiff

theorem Ico_zero_filter_modEq_card_of_dvd_local
    (g modulus representative : ℤ) (hg : 0 ≤ g)
    (hmodulus : 0 < modulus) (hmodulusDvd : modulus ∣ g) :
    ({x ∈ Finset.Ico (0 : ℤ) g |
      x ≡ representative [ZMOD modulus]}).card = (g / modulus).toNat := by
  rw [← Int.ofNat_inj]
  rw [Int.Ico_filter_modEq_card (r := modulus) 0 g hmodulus representative]
  obtain ⟨factor, rfl⟩ := hmodulusDvd
  have hfactor : 0 ≤ factor := by nlinarith
  push_cast
  have hrational :
      (((modulus : ℚ) * (factor : ℚ) - (representative : ℚ)) /
          (modulus : ℚ)) =
        (-(representative : ℚ) / (modulus : ℚ)) + (factor : ℚ) := by
    have hmodulusQ : (modulus : ℚ) ≠ 0 := by exact_mod_cast hmodulus.ne'
    field_simp
    ring
  rw [hrational, Int.ceil_add_intCast]
  have hzero :
      (((0 : ℚ) - (representative : ℚ)) / (modulus : ℚ)) =
        (-(representative : ℚ) / (modulus : ℚ)) := by ring
  rw [hzero]
  have hdiv : modulus * factor / modulus = factor :=
    Int.mul_ediv_cancel_left factor hmodulus.ne'
  rw [hdiv]
  omega

theorem modEq_two_of_not_modEq_local {first second parity : ℤ}
    (hfirst : ¬ first ≡ parity [ZMOD 2])
    (hsecond : ¬ second ≡ parity [ZMOD 2]) :
    first ≡ second [ZMOD 2] := by
  change first % 2 = second % 2
  change ¬ first % 2 = parity % 2 at hfirst
  change ¬ second % 2 = parity % 2 at hsecond
  rcases Int.emod_two_eq_zero_or_one first with hfirst0 | hfirst1 <;>
    rcases Int.emod_two_eq_zero_or_one second with hsecond0 | hsecond1 <;>
      rcases Int.emod_two_eq_zero_or_one parity with hparity0 | hparity1 <;>
        omega

/-- A saturated q8 row meets either parity class. -/
theorem exists_qEight_badBranch_modEq_two
    (δ g divisor sign numerator : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq8 : g / (Int.gcd δ g : ℤ) = 8)
    (hdivisor : divisor ≠ 0)
    (hδ : δ = divisor * (sign * numerator))
    (hgFactor : g = divisor * 8)
    (hsign : sign = 1 ∨ sign = -1)
    (hcoprime : IsCoprime numerator 8)
    (hcard : 4 * (detunedBadBranches δ g u).card = g.toNat)
    (parity : ℤ) :
    ∃ branch ∈ detunedBadBranches δ g u,
      branch ≡ parity [ZMOD 2] := by
  have hg0 : (0 : ℤ) < g := by omega
  have hg1 : (1 : ℤ) ≤ g := by omega
  have hgNat : 0 < g.toNat := by omega
  have hrow : (detunedBadBranches δ g u).Nonempty := by
    rw [← Finset.card_pos]
    omega
  obtain ⟨base, hbase⟩ := hrow
  have hrowEq := detunedBadBranches_eq_normalizedBadBranches
    δ g divisor 8 sign numerator u hdivisor (by norm_num) hδ hgFactor hsign
  by_contra hnone
  push Not at hnone
  have hsameParity : ∀ branch ∈ detunedBadBranches δ g u,
      branch ≡ base [ZMOD 2] := by
    intro branch hbranch
    exact modEq_two_of_not_modEq_local
      (hnone branch hbranch) (hnone base hbase)
  have hsubset : detunedBadBranches δ g u ⊆
      {branch ∈ Finset.Ico (0 : ℤ) g | branch ≡ base [ZMOD 8]} := by
    intro branch hbranch
    rw [Finset.mem_filter]
    refine ⟨detunedBadBranches_subset_Ico δ g u hbranch,
      ?_⟩
    have hbranchNorm : branch ∈ normalizedBadBranches numerator 8 g u := by
      rw [← hrowEq]
      exact hbranch
    have hbaseNorm : base ∈ normalizedBadBranches numerator 8 g u := by
      rw [← hrowEq]
      exact hbase
    exact normalizedQEight_pair_modEq_of_same_parity
      numerator g u hcoprime hbranchNorm hbaseNorm
        (hsameParity branch hbranch)
  have h8dvd : (8 : ℤ) ∣ g := by
    have hfactor := Int.mul_ediv_cancel' (Int.gcd_dvd_right δ g)
    rw [hq8] at hfactor
    exact ⟨(Int.gcd δ g : ℤ), by simpa [mul_comm] using hfactor.symm⟩
  have hfilterCard := Ico_zero_filter_modEq_card_of_dvd_local
    g 8 base (by omega) (by norm_num) h8dvd
  have hcardLe := Finset.card_le_card hsubset
  rw [hfilterCard] at hcardLe
  have hfactorZ := Int.mul_ediv_cancel' (Int.gcd_dvd_right δ g)
  rw [hq8] at hfactorZ
  have hdiv : g / 8 = (Int.gcd δ g : ℤ) := by
    calc
      g / 8 = ((Int.gcd δ g : ℤ) * 8) / 8 :=
        congrArg (fun value : ℤ => value / 8) hfactorZ.symm
      _ = (Int.gcd δ g : ℤ) := by omega
  rw [hdiv] at hcardLe
  have hdpos : 0 < Int.gcd δ g := by
    rw [Int.gcd_pos_iff]
    right
    omega
  have hfactorNat : 8 * Int.gcd δ g = g.toNat := by omega
  omega

/-- Exact residue matching extracted from one failing q4488 partition, for
both choices of q4 row. -/
theorem two_matching_codes_of_four_four_eight_eight_failure
    (δ₄a δ₄b δ₈a δ₈b g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (hq8a : g / (Int.gcd δ₈a g : ℤ) = 8)
    (hq8b : g / (Int.gcd δ₈b g : ℤ) = 8)
    (a₄a a₄b a₈a a₈b : ℤ)
    (hnorm4a : SignedReducedNormalForm δ₄a g 4 a₄a)
    (hnorm4b : SignedReducedNormalForm δ₄b g 4 a₄b)
    (hnorm8a : SignedReducedNormalForm δ₈a g 8 a₈a)
    (hnorm8b : SignedReducedNormalForm δ₈b g 8 a₈b)
    (hfailure : ¬ HasFourDetunedGoodBranch δ₄a δ₄b δ₈a δ₈b g u) :
    (∃ c₄a c₈a c₈b : ℤ,
      c₄a ∈ detunedBadBranches δ₄a g u ∧
      c₈a ∈ detunedBadBranches δ₈a g u ∧
      c₈b ∈ detunedBadBranches δ₈b g u ∧
      c₈a ≡ c₄a [ZMOD 2] ∧
      c₈b ≡ c₄a [ZMOD 2] ∧
      (8 : ℤ) ∣ -2 * c₄a + c₈a + c₈b) ∧
    (∃ c₄b c₈a c₈b : ℤ,
      c₄b ∈ detunedBadBranches δ₄b g u ∧
      c₈a ∈ detunedBadBranches δ₈a g u ∧
      c₈b ∈ detunedBadBranches δ₈b g u ∧
      c₈a ≡ c₄b [ZMOD 2] ∧
      c₈b ≡ c₄b [ZMOD 2] ∧
      (8 : ℤ) ∣ -2 * c₄b + c₈a + c₈b) := by
  obtain ⟨hpartition, hcard4a, hcard4b, hcard8a, hcard8b⟩ :=
    four_four_eight_eight_failure_parallel_partition
      δ₄a δ₄b δ₈a δ₈b g u hg hq4a hq4b hq8a hq8b hfailure
  have hgNat : 0 < g.toNat := by omega
  have hrow4a : (detunedBadBranches δ₄a g u).Nonempty := by
    rw [← Finset.card_pos]
    omega
  have hrow4b : (detunedBadBranches δ₄b g u).Nonempty := by
    rw [← Finset.card_pos]
    omega
  obtain ⟨c₄a, hc₄a⟩ := hrow4a
  obtain ⟨c₄b, hc₄b⟩ := hrow4b
  have hrowEq4a := detunedBadBranches_eq_normalizedBadBranches
    δ₄a g hnorm4a.divisor 4 hnorm4a.sign a₄a u
      hnorm4a.divisor_ne (by norm_num) hnorm4a.delta_eq
      hnorm4a.modulus_eq hnorm4a.sign_eq
  have hrowEq4b := detunedBadBranches_eq_normalizedBadBranches
    δ₄b g hnorm4b.divisor 4 hnorm4b.sign a₄b u
      hnorm4b.divisor_ne (by norm_num) hnorm4b.delta_eq
      hnorm4b.modulus_eq hnorm4b.sign_eq
  have hrowEq8a := detunedBadBranches_eq_normalizedBadBranches
    δ₈a g hnorm8a.divisor 8 hnorm8a.sign a₈a u
      hnorm8a.divisor_ne (by norm_num) hnorm8a.delta_eq
      hnorm8a.modulus_eq hnorm8a.sign_eq
  have hrowEq8b := detunedBadBranches_eq_normalizedBadBranches
    δ₈b g hnorm8b.divisor 8 hnorm8b.sign a₈b u
      hnorm8b.divisor_ne (by norm_num) hnorm8b.delta_eq
      hnorm8b.modulus_eq hnorm8b.sign_eq
  have buildMatching
      (δ₄ a₄ : ℤ) (hrowEq4 : detunedBadBranches δ₄ g u =
        normalizedBadBranches a₄ 4 g u)
      (c₄ : ℤ) (hc₄ : c₄ ∈ detunedBadBranches δ₄ g u)
      (hdisjoint4a : Disjoint (detunedBadBranches δ₄ g u)
        (detunedBadBranches δ₈a g u))
      (hdisjoint4b : Disjoint (detunedBadBranches δ₄ g u)
        (detunedBadBranches δ₈b g u)) :
      ∃ c₈a c₈b : ℤ,
        c₈a ∈ detunedBadBranches δ₈a g u ∧
        c₈b ∈ detunedBadBranches δ₈b g u ∧
        c₈a ≡ c₄ [ZMOD 2] ∧
        c₈b ≡ c₄ [ZMOD 2] ∧
        (8 : ℤ) ∣ -2 * c₄ + c₈a + c₈b := by
    obtain ⟨c₈a, hc₈a, hpar8a⟩ :=
      exists_qEight_badBranch_modEq_two
        δ₈a g hnorm8a.divisor hnorm8a.sign a₈a u hg hq8a
          hnorm8a.divisor_ne hnorm8a.delta_eq hnorm8a.modulus_eq
          hnorm8a.sign_eq hnorm8a.coprime hcard8a c₄
    obtain ⟨c₈b, hc₈b, hpar8b⟩ :=
      exists_qEight_badBranch_modEq_two
        δ₈b g hnorm8b.divisor hnorm8b.sign a₈b u hg hq8b
          hnorm8b.divisor_ne hnorm8b.delta_eq hnorm8b.modulus_eq
          hnorm8b.sign_eq hnorm8b.coprime hcard8b c₄
    have hout8a : ¬ c₈a ≡ c₄ [ZMOD 4] := by
      intro hmod
      have hc₈a4 : c₈a ∈ detunedBadBranches δ₄ g u := by
        rw [hrowEq4]
        apply normalizedBadBranches_mem_of_modEq
          a₄ 4 g u (by norm_num)
        · rw [← hrowEq4]
          exact hc₄
        · exact detunedBadBranches_subset_Ico δ₈a g u hc₈a
        · exact hmod.symm
      exact (Finset.disjoint_left.mp hdisjoint4a) hc₈a4 hc₈a
    have hout8b : ¬ c₈b ≡ c₄ [ZMOD 4] := by
      intro hmod
      have hc₈b4 : c₈b ∈ detunedBadBranches δ₄ g u := by
        rw [hrowEq4]
        apply normalizedBadBranches_mem_of_modEq
          a₄ 4 g u (by norm_num)
        · rw [← hrowEq4]
          exact hc₄
        · exact detunedBadBranches_subset_Ico δ₈b g u hc₈b
        · exact hmod.symm
      exact (Finset.disjoint_left.mp hdisjoint4b) hc₈b4 hc₈b
    have hdistinct : ¬ c₈a ≡ c₈b [ZMOD 8] := by
      intro hmod
      have hc₈b8a : c₈b ∈ detunedBadBranches δ₈a g u := by
        rw [hrowEq8a]
        apply normalizedBadBranches_mem_of_modEq
          a₈a 8 g u (by norm_num)
        · rw [← hrowEq8a]
          exact hc₈a
        · exact detunedBadBranches_subset_Ico δ₈b g u hc₈b
        · exact hmod
      exact (Finset.disjoint_left.mp hpartition.disjoint₃₄)
        hc₈b8a hc₈b
    exact ⟨c₈a, c₈b, hc₈a, hc₈b, hpar8a, hpar8b,
      eight_dvd_complementary_parity_sum c₄ c₈a c₈b
        hpar8a hpar8b hout8a hout8b hdistinct⟩
  obtain ⟨c₈aa, c₈ba, hc₈aa, hc₈ba, hpar8aa, hpar8ba, hcodea⟩ :=
    buildMatching δ₄a a₄a hrowEq4a c₄a hc₄a
      hpartition.disjoint₁₃ hpartition.disjoint₁₄
  obtain ⟨c₈ab, c₈bb, hc₈ab, hc₈bb, hpar8ab, hpar8bb, hcodeb⟩ :=
    buildMatching δ₄b a₄b hrowEq4b c₄b hc₄b
      hpartition.disjoint₂₃ hpartition.disjoint₂₄
  exact ⟨⟨c₄a, c₈aa, c₈ba, hc₄a, hc₈aa, hc₈ba,
      hpar8aa, hpar8ba, hcodea⟩,
    ⟨c₄b, c₈ab, c₈bb, hc₄b, hc₈ab, hc₈bb,
      hpar8ab, hpar8bb, hcodeb⟩⟩

/-- The exact q4488 failure partition supplies the complete normalized
matching-wall input for either q4 row. -/
theorem two_matching_walls_of_four_four_eight_eight_failure
    (δ₄a δ₄b δ₈a δ₈b g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (hq8a : g / (Int.gcd δ₈a g : ℤ) = 8)
    (hq8b : g / (Int.gcd δ₈b g : ℤ) = 8)
    (a₄a a₄b a₈a a₈b : ℤ)
    (hnorm4a : SignedReducedNormalForm δ₄a g 4 a₄a)
    (hnorm4b : SignedReducedNormalForm δ₄b g 4 a₄b)
    (hnorm8a : SignedReducedNormalForm δ₈a g 8 a₈a)
    (hnorm8b : SignedReducedNormalForm δ₈b g 8 a₈b)
    (hfailure : ¬ HasFourDetunedGoodBranch δ₄a δ₄b δ₈a δ₈b g u) :
    (∃ c₄a c₈a c₈b : ℤ,
      (8 : ℤ) ∣ -2 * c₄a + c₈a + c₈b ∧
      (∃ n : ℤ,
        |(a₄a : ℝ) * ((u + (c₄a : ℝ)) / 4) - n| < 1 / 14) ∧
      (∃ n : ℤ,
        |(a₈a : ℝ) * ((u + (c₈a : ℝ)) / 8) - n| < 1 / 14) ∧
      (∃ n : ℤ,
        |(a₈b : ℝ) * ((u + (c₈b : ℝ)) / 8) - n| < 1 / 14)) ∧
    (∃ c₄b c₈a c₈b : ℤ,
      (8 : ℤ) ∣ -2 * c₄b + c₈a + c₈b ∧
      (∃ n : ℤ,
        |(a₄b : ℝ) * ((u + (c₄b : ℝ)) / 4) - n| < 1 / 14) ∧
      (∃ n : ℤ,
        |(a₈a : ℝ) * ((u + (c₈a : ℝ)) / 8) - n| < 1 / 14) ∧
      (∃ n : ℤ,
        |(a₈b : ℝ) * ((u + (c₈b : ℝ)) / 8) - n| < 1 / 14)) := by
  obtain ⟨⟨c₄a, c₈aa, c₈ba, hc₄a, hc₈aa, hc₈ba, -, -, hcodea⟩,
      ⟨c₄b, c₈ab, c₈bb, hc₄b, hc₈ab, hc₈bb, -, -, hcodeb⟩⟩ :=
    two_matching_codes_of_four_four_eight_eight_failure
      δ₄a δ₄b δ₈a δ₈b g u hg hq4a hq4b hq8a hq8b
        a₄a a₄b a₈a a₈b hnorm4a hnorm4b hnorm8a hnorm8b hfailure
  refine ⟨⟨c₄a, c₈aa, c₈ba, hcodea, ?_, ?_, ?_⟩,
    ⟨c₄b, c₈ab, c₈bb, hcodeb, ?_, ?_, ?_⟩⟩
  · exact normalized_bad_witness_of_detuned_mem
      δ₄a g 4 a₄a u hnorm4a (by norm_num) hc₄a
  · exact normalized_bad_witness_of_detuned_mem
      δ₈a g 8 a₈a u hnorm8a (by norm_num) hc₈aa
  · exact normalized_bad_witness_of_detuned_mem
      δ₈b g 8 a₈b u hnorm8b (by norm_num) hc₈ba
  · exact normalized_bad_witness_of_detuned_mem
      δ₄b g 4 a₄b u hnorm4b (by norm_num) hc₄b
  · exact normalized_bad_witness_of_detuned_mem
      δ₈a g 8 a₈a u hnorm8a (by norm_num) hc₈ab
  · exact normalized_bad_witness_of_detuned_mem
      δ₈b g 8 a₈b u hnorm8b (by norm_num) hc₈bb

/-- Disjoint saturated q4 rows cannot have the same signed primitive
numerator: equality would make their normalized bad rows identical. -/
theorem normalized_qFour_numerators_ne_of_failure
    (δ₄a δ₄b δ₈a δ₈b g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (hq8a : g / (Int.gcd δ₈a g : ℤ) = 8)
    (hq8b : g / (Int.gcd δ₈b g : ℤ) = 8)
    (a₄a a₄b : ℤ)
    (hnorm4a : SignedReducedNormalForm δ₄a g 4 a₄a)
    (hnorm4b : SignedReducedNormalForm δ₄b g 4 a₄b)
    (hfailure : ¬ HasFourDetunedGoodBranch δ₄a δ₄b δ₈a δ₈b g u) :
    a₄a ≠ a₄b := by
  intro hequal
  obtain ⟨hpartition, hcard4a, -, -, -⟩ :=
    four_four_eight_eight_failure_parallel_partition
      δ₄a δ₄b δ₈a δ₈b g u hg hq4a hq4b hq8a hq8b hfailure
  have hrow4a : (detunedBadBranches δ₄a g u).Nonempty := by
    rw [← Finset.card_pos]
    have hgNat : 0 < g.toNat := by omega
    omega
  obtain ⟨branch, hbranch4a⟩ := hrow4a
  have hrowEq4a := detunedBadBranches_eq_normalizedBadBranches
    δ₄a g hnorm4a.divisor 4 hnorm4a.sign a₄a u
      hnorm4a.divisor_ne (by norm_num) hnorm4a.delta_eq
      hnorm4a.modulus_eq hnorm4a.sign_eq
  have hrowEq4b := detunedBadBranches_eq_normalizedBadBranches
    δ₄b g hnorm4b.divisor 4 hnorm4b.sign a₄b u
      hnorm4b.divisor_ne (by norm_num) hnorm4b.delta_eq
      hnorm4b.modulus_eq hnorm4b.sign_eq
  have hrowsEqual : detunedBadBranches δ₄a g u =
      detunedBadBranches δ₄b g u := by
    rw [hrowEq4a, hrowEq4b, hequal]
  have hbranch4b : branch ∈ detunedBadBranches δ₄b g u := by
    rw [← hrowsEqual]
    exact hbranch4a
  exact (Finset.disjoint_left.mp hpartition.disjoint₁₂)
    hbranch4a hbranch4b

/-- Common signed residue data makes at least one of the two q4-derived
q4488 frequencies nonzero when the q4 numerators differ. -/
theorem one_fourEightEightPhaseFrequency_ne_zero
    (a₄a a₄b a₈a a₈b residue : ℤ)
    (hres₄a : (4 : ℤ) ∣ a₄a - residue)
    (hres₄b : (4 : ℤ) ∣ a₄b - residue)
    (hres₈a : (8 : ℤ) ∣ a₈a - residue)
    (hres₈b : (8 : ℤ) ∣ a₈b - residue)
    (hdistinct : a₄a ≠ a₄b) :
    fourEightEightPhaseFrequency a₄a a₈a a₈b ≠ 0 ∨
      fourEightEightPhaseFrequency a₄b a₈a a₈b ≠ 0 := by
  have hdivA : (8 : ℤ) ∣ -2 * a₄a + a₈a + a₈b := by
    obtain ⟨k₄a, hk₄a⟩ := hres₄a
    obtain ⟨k₈a, hk₈a⟩ := hres₈a
    obtain ⟨k₈b, hk₈b⟩ := hres₈b
    refine ⟨-k₄a + k₈a + k₈b, ?_⟩
    omega
  have hdivB : (8 : ℤ) ∣ -2 * a₄b + a₈a + a₈b := by
    obtain ⟨k₄b, hk₄b⟩ := hres₄b
    obtain ⟨k₈a, hk₈a⟩ := hres₈a
    obtain ⟨k₈b, hk₈b⟩ := hres₈b
    refine ⟨-k₄b + k₈a + k₈b, ?_⟩
    omega
  have hmulA := Int.ediv_mul_cancel hdivA
  have hmulB := Int.ediv_mul_cancel hdivB
  change fourEightEightPhaseFrequency a₄a a₈a a₈b * 8 =
    -2 * a₄a + a₈a + a₈b at hmulA
  change fourEightEightPhaseFrequency a₄b a₈a a₈b * 8 =
    -2 * a₄b + a₈a + a₈b at hmulB
  by_contra hzero
  rw [not_or] at hzero
  push Not at hzero
  simp only [hzero.1, hzero.2, zero_mul] at hmulA hmulB
  exact hdistinct (by omega)

/-- Joint large/small split for the two q4-derived q4488 frequencies.  If
neither frequency crosses the sharp LRC(9) escape threshold, each is either
an exact signed three-term relation or a nonzero relation below threshold,
and the two exact relations cannot occur simultaneously. -/
theorem two_fourEightEight_frequencies_large_or_relation_small
    (a₄a a₄b a₈a a₈b residue : ℤ) (B : ℝ)
    (hB0 : 0 < B)
    (hres₄a : (4 : ℤ) ∣ a₄a - residue)
    (hres₄b : (4 : ℤ) ∣ a₄b - residue)
    (hres₈a : (8 : ℤ) ∣ a₈a - residue)
    (hres₈b : (8 : ℤ) ∣ a₈b - residue)
    (hdistinct : a₄a ≠ a₄b) :
    (fourEightEightPhaseFrequency a₄a a₈a a₈b ≠ 0 ∧
      15 * B ≤ 2 *
        |(fourEightEightPhaseFrequency a₄a a₈a a₈b : ℝ)|) ∨
    (fourEightEightPhaseFrequency a₄b a₈a a₈b ≠ 0 ∧
      15 * B ≤ 2 *
        |(fourEightEightPhaseFrequency a₄b a₈a a₈b : ℝ)|) ∨
    (((-2 * a₄a + a₈a + a₈b = 0) ∨
        (fourEightEightPhaseFrequency a₄a a₈a a₈b ≠ 0 ∧
          2 * |(fourEightEightPhaseFrequency a₄a a₈a a₈b : ℝ)| <
            15 * B)) ∧
      ((-2 * a₄b + a₈a + a₈b = 0) ∨
        (fourEightEightPhaseFrequency a₄b a₈a a₈b ≠ 0 ∧
          2 * |(fourEightEightPhaseFrequency a₄b a₈a a₈b : ℝ)| <
            15 * B)) ∧
      (fourEightEightPhaseFrequency a₄a a₈a a₈b ≠ 0 ∨
        fourEightEightPhaseFrequency a₄b a₈a a₈b ≠ 0)) := by
  have hnonzero := one_fourEightEightPhaseFrequency_ne_zero
    a₄a a₄b a₈a a₈b residue hres₄a hres₄b hres₈a hres₈b hdistinct
  have hdivA : (8 : ℤ) ∣ -2 * a₄a + a₈a + a₈b := by
    obtain ⟨k₄a, hk₄a⟩ := hres₄a
    obtain ⟨k₈a, hk₈a⟩ := hres₈a
    obtain ⟨k₈b, hk₈b⟩ := hres₈b
    refine ⟨-k₄a + k₈a + k₈b, by omega⟩
  have hdivB : (8 : ℤ) ∣ -2 * a₄b + a₈a + a₈b := by
    obtain ⟨k₄b, hk₄b⟩ := hres₄b
    obtain ⟨k₈a, hk₈a⟩ := hres₈a
    obtain ⟨k₈b, hk₈b⟩ := hres₈b
    refine ⟨-k₄b + k₈a + k₈b, by omega⟩
  have hmulA := Int.ediv_mul_cancel hdivA
  have hmulB := Int.ediv_mul_cancel hdivB
  change fourEightEightPhaseFrequency a₄a a₈a a₈b * 8 =
    -2 * a₄a + a₈a + a₈b at hmulA
  change fourEightEightPhaseFrequency a₄b a₈a a₈b * 8 =
    -2 * a₄b + a₈a + a₈b at hmulB
  by_cases hlargeA : 15 * B ≤ 2 *
      |(fourEightEightPhaseFrequency a₄a a₈a a₈b : ℝ)|
  · have hneA : fourEightEightPhaseFrequency a₄a a₈a a₈b ≠ 0 := by
      intro hzeroA
      rw [hzeroA] at hlargeA
      norm_num at hlargeA
      linarith
    exact Or.inl ⟨hneA, hlargeA⟩
  · have hsmallA : 2 *
        |(fourEightEightPhaseFrequency a₄a a₈a a₈b : ℝ)| < 15 * B :=
      lt_of_not_ge hlargeA
    by_cases hlargeB : 15 * B ≤ 2 *
        |(fourEightEightPhaseFrequency a₄b a₈a a₈b : ℝ)|
    · have hneB : fourEightEightPhaseFrequency a₄b a₈a a₈b ≠ 0 := by
        intro hzeroB
        rw [hzeroB] at hlargeB
        norm_num at hlargeB
        linarith
      exact Or.inr (Or.inl ⟨hneB, hlargeB⟩)
    · have hsmallB : 2 *
          |(fourEightEightPhaseFrequency a₄b a₈a a₈b : ℝ)| < 15 * B :=
        lt_of_not_ge hlargeB
      refine Or.inr (Or.inr ⟨?_, ?_, hnonzero⟩)
      · by_cases hzeroA : fourEightEightPhaseFrequency a₄a a₈a a₈b = 0
        · exact Or.inl (by omega)
        · exact Or.inr ⟨hzeroA, hsmallA⟩
      · by_cases hzeroB : fourEightEightPhaseFrequency a₄b a₈a a₈b = 0
        · exact Or.inl (by omega)
        · exact Or.inr ⟨hzeroB, hsmallB⟩

theorem normalizedQEight_pair_unit_relation_of_opposite_parity
    (numerator g : ℤ) (u : ℝ)
    (hcoprime : IsCoprime numerator 8)
    {first second : ℤ}
    (hfirst : first ∈ normalizedBadBranches numerator 8 g u)
    (hsecond : second ∈ normalizedBadBranches numerator 8 g u)
    (hparity : ¬ first ≡ second [ZMOD 2]) :
    (8 : ℤ) ∣ numerator * (first - second) - 1 ∨
      (8 : ℤ) ∣ numerator * (first - second) + 1 := by
  obtain ⟨integerFirst, hnearFirst⟩ := (Finset.mem_filter.mp hfirst).2
  obtain ⟨integerSecond, hnearSecond⟩ := (Finset.mem_filter.mp hsecond).2
  let defect : ℤ :=
    numerator * (first - second) - 8 * (integerFirst - integerSecond)
  have hdefectReal :
      (defect : ℝ) / 8 =
        ((numerator : ℝ) * ((u + (first : ℝ)) / 8) - integerFirst) -
        ((numerator : ℝ) * ((u + (second : ℝ)) / 8) - integerSecond) := by
    dsimp [defect]
    push_cast
    ring
  have hdefectAbs : |(defect : ℝ)| < 2 := by
    have hquotient : |(defect : ℝ) / 8| < 1 / 7 := by
      rw [hdefectReal]
      calc
        |_ - _| ≤
            |(numerator : ℝ) * ((u + (first : ℝ)) / 8) - integerFirst| +
            |(numerator : ℝ) * ((u + (second : ℝ)) / 8) - integerSecond| :=
          abs_sub _ _
        _ < 1 / 14 + 1 / 14 := add_lt_add hnearFirst hnearSecond
        _ = 1 / 7 := by norm_num
    rw [abs_div, abs_of_pos (by norm_num : (0 : ℝ) < 8)] at hquotient
    nlinarith [hquotient]
  have hdefectInt : |defect| < 2 := by exact_mod_cast hdefectAbs
  have hnumeratorOdd : Odd numerator := by
    apply Int.isCoprime_two_right.mp
    have hcoprime' : IsCoprime numerator ((4 : ℤ) * 2) := by
      norm_num
      exact hcoprime
    exact hcoprime'.of_mul_right_right
  have hdiffOdd : Odd (first - second) := by
    apply Int.not_even_iff_odd.mp
    intro hdiffEven
    apply hparity
    rw [Int.modEq_iff_dvd]
    have hdiffDvd : (2 : ℤ) ∣ first - second :=
      even_iff_two_dvd.mp hdiffEven
    simpa only [neg_sub] using dvd_neg.mpr hdiffDvd
  have hdefectOdd : Odd defect := by
    apply Odd.sub_even (hnumeratorOdd.mul hdiffOdd)
    exact ⟨4 * (integerFirst - integerSecond), by ring⟩
  obtain ⟨halfDefect, hhalfDefect⟩ := hdefectOdd
  have hdefectBounds := abs_lt.mp hdefectInt
  have hdefectCases : defect = 1 ∨ defect = -1 := by omega
  rcases hdefectCases with hdefectOne | hdefectNegOne
  · exact Or.inl ⟨integerFirst - integerSecond, by
      dsimp [defect] at hdefectOne
      omega⟩
  · exact Or.inr ⟨integerFirst - integerSecond, by
      dsimp [defect] at hdefectNegOne
      omega⟩

theorem qEight_numerators_modEq_or_neg_of_parallel_rectangle
    (numeratorA numeratorB g : ℤ) (u : ℝ)
    (hcoprimeA : IsCoprime numeratorA 8)
    (hcoprimeB : IsCoprime numeratorB 8)
    {firstA secondA firstB secondB : ℤ}
    (hfirstA : firstA ∈ normalizedBadBranches numeratorA 8 g u)
    (hsecondA : secondA ∈ normalizedBadBranches numeratorA 8 g u)
    (hfirstB : firstB ∈ normalizedBadBranches numeratorB 8 g u)
    (hsecondB : secondB ∈ normalizedBadBranches numeratorB 8 g u)
    (hopposite : ¬ firstA ≡ secondA [ZMOD 2])
    (hparallel : firstA - secondA ≡ firstB - secondB [ZMOD 8]) :
    numeratorA ≡ numeratorB [ZMOD 8] ∨
      numeratorA ≡ -numeratorB [ZMOD 8] := by
  have hoppositeB : ¬ firstB ≡ secondB [ZMOD 2] := by
    intro hsameB
    have hdiffBEven : (2 : ℤ) ∣ firstB - secondB := by
      have := hsameB.dvd
      simpa only [neg_sub] using dvd_neg.mpr this
    have hdiffAEven : (2 : ℤ) ∣ firstA - secondA := by
      have hparallelTwo : firstA - secondA ≡ firstB - secondB [ZMOD 2] :=
        hparallel.of_dvd (by norm_num)
      have hparallelDvd := hparallelTwo.dvd
      have := dvd_sub hdiffBEven hparallelDvd
      have heq :
          firstB - secondB - ((firstB - secondB) - (firstA - secondA)) =
            firstA - secondA := by ring
      rw [heq] at this
      exact this
    apply hopposite
    rw [Int.modEq_iff_dvd]
    simpa only [neg_sub] using dvd_neg.mpr hdiffAEven
  have hunitA := normalizedQEight_pair_unit_relation_of_opposite_parity
    numeratorA g u hcoprimeA hfirstA hsecondA hopposite
  have hunitB0 := normalizedQEight_pair_unit_relation_of_opposite_parity
    numeratorB g u hcoprimeB hfirstB hsecondB hoppositeB
  have hunitB :
      (8 : ℤ) ∣ numeratorB * (firstA - secondA) - 1 ∨
      (8 : ℤ) ∣ numeratorB * (firstA - secondA) + 1 := by
    rcases hunitB0 with hplus | hminus
    · left
      rw [Int.modEq_iff_dvd] at hparallel
      obtain ⟨parallelShift, hparallelShift⟩ := hparallel
      obtain ⟨unitShift, hunitShift⟩ := hplus
      refine ⟨unitShift - numeratorB * parallelShift, ?_⟩
      calc
        numeratorB * (firstA - secondA) - 1 =
            (numeratorB * (firstB - secondB) - 1) -
              numeratorB * ((firstB - secondB) - (firstA - secondA)) := by ring
        _ = 8 * unitShift - numeratorB * (8 * parallelShift) := by
          rw [hunitShift, hparallelShift]
        _ = 8 * (unitShift - numeratorB * parallelShift) := by ring
    · right
      rw [Int.modEq_iff_dvd] at hparallel
      obtain ⟨parallelShift, hparallelShift⟩ := hparallel
      obtain ⟨unitShift, hunitShift⟩ := hminus
      refine ⟨unitShift - numeratorB * parallelShift, ?_⟩
      calc
        numeratorB * (firstA - secondA) + 1 =
            (numeratorB * (firstB - secondB) + 1) -
              numeratorB * ((firstB - secondB) - (firstA - secondA)) := by ring
        _ = 8 * unitShift - numeratorB * (8 * parallelShift) := by
          rw [hunitShift, hparallelShift]
        _ = 8 * (unitShift - numeratorB * parallelShift) := by ring
  have hcoprimeDifference : IsCoprime (firstA - secondA) 8 := by
    rcases hunitA with hplus | hminus
    · obtain ⟨unitShift, hunitShift⟩ := hplus
      refine ⟨numeratorA, -unitShift, ?_⟩
      calc
        numeratorA * (firstA - secondA) + -unitShift * 8 =
            (numeratorA * (firstA - secondA) - 1) -
              8 * unitShift + 1 := by ring
        _ = 1 := by rw [hunitShift]; ring
    · obtain ⟨unitShift, hunitShift⟩ := hminus
      refine ⟨-numeratorA, unitShift, ?_⟩
      calc
        -numeratorA * (firstA - secondA) + unitShift * 8 =
            -(numeratorA * (firstA - secondA) + 1) +
              8 * unitShift + 1 := by ring
        _ = 1 := by rw [hunitShift]; ring
  rcases hunitA with hplusA | hminusA <;>
    rcases hunitB with hplusB | hminusB
  · left
    rw [Int.modEq_iff_dvd]
    apply hcoprimeDifference.symm.dvd_of_dvd_mul_right
    obtain ⟨shiftA, hshiftA⟩ := hplusA
    obtain ⟨shiftB, hshiftB⟩ := hplusB
    refine ⟨shiftB - shiftA, ?_⟩
    calc
      (numeratorB - numeratorA) * (firstA - secondA) =
          (numeratorB * (firstA - secondA) - 1) -
            (numeratorA * (firstA - secondA) - 1) := by ring
      _ = 8 * shiftB - 8 * shiftA := by rw [hshiftA, hshiftB]
      _ = 8 * (shiftB - shiftA) := by ring
  · right
    rw [Int.modEq_iff_dvd]
    apply hcoprimeDifference.symm.dvd_of_dvd_mul_right
    obtain ⟨shiftA, hshiftA⟩ := hplusA
    obtain ⟨shiftB, hshiftB⟩ := hminusB
    refine ⟨-(shiftA + shiftB), ?_⟩
    calc
      (-numeratorB - numeratorA) * (firstA - secondA) =
          -((numeratorA * (firstA - secondA) - 1) +
            (numeratorB * (firstA - secondA) + 1)) := by ring
      _ = -(8 * shiftA + 8 * shiftB) := by rw [hshiftA, hshiftB]
      _ = 8 * (-(shiftA + shiftB)) := by ring
  · right
    rw [Int.modEq_iff_dvd]
    apply hcoprimeDifference.symm.dvd_of_dvd_mul_right
    obtain ⟨shiftA, hshiftA⟩ := hminusA
    obtain ⟨shiftB, hshiftB⟩ := hplusB
    refine ⟨-(shiftA + shiftB), ?_⟩
    calc
      (-numeratorB - numeratorA) * (firstA - secondA) =
          -((numeratorA * (firstA - secondA) + 1) +
            (numeratorB * (firstA - secondA) - 1)) := by ring
      _ = -(8 * shiftA + 8 * shiftB) := by rw [hshiftA, hshiftB]
      _ = 8 * (-(shiftA + shiftB)) := by ring
  · left
    rw [Int.modEq_iff_dvd]
    apply hcoprimeDifference.symm.dvd_of_dvd_mul_right
    obtain ⟨shiftA, hshiftA⟩ := hminusA
    obtain ⟨shiftB, hshiftB⟩ := hminusB
    refine ⟨shiftB - shiftA, ?_⟩
    calc
      (numeratorB - numeratorA) * (firstA - secondA) =
          (numeratorB * (firstA - secondA) + 1) -
            (numeratorA * (firstA - secondA) + 1) := by ring
      _ = 8 * shiftB - 8 * shiftA := by rw [hshiftA, hshiftB]
      _ = 8 * (shiftB - shiftA) := by ring

theorem modEq_four_of_same_parity_outside
    (fiber first second : ℤ)
    (hfirstParity : first ≡ fiber [ZMOD 2])
    (hsecondParity : second ≡ fiber [ZMOD 2])
    (hfirstOutside : ¬ first ≡ fiber [ZMOD 4])
    (hsecondOutside : ¬ second ≡ fiber [ZMOD 4]) :
    first ≡ second [ZMOD 4] := by
  obtain ⟨firstShift, hfirstShift⟩ :=
    four_dvd_sub_sub_two_of_same_parity
      fiber first hfirstParity hfirstOutside
  obtain ⟨secondShift, hsecondShift⟩ :=
    four_dvd_sub_sub_two_of_same_parity
      fiber second hsecondParity hsecondOutside
  rw [Int.modEq_iff_dvd]
  refine ⟨secondShift - firstShift, by omega⟩

theorem modEq_eight_add_four_of_modEq_four_not_modEq_eight
    (first second : ℤ)
    (hfour : first ≡ second [ZMOD 4])
    (height : ¬ first ≡ second [ZMOD 8]) :
    first ≡ second + 4 [ZMOD 8] := by
  obtain ⟨shift, hshift⟩ := hfour.dvd
  have hshiftNotEven : ¬ Even shift := by
    intro hshiftEven
    obtain ⟨halfShift, hhalfShift⟩ := hshiftEven
    apply height
    rw [Int.modEq_iff_dvd]
    refine ⟨halfShift, by omega⟩
  obtain ⟨halfShift, hhalfShift⟩ :=
    Int.not_even_iff_odd.mp hshiftNotEven
  rw [Int.modEq_iff_dvd]
  refine ⟨halfShift + 1, by omega⟩

/-- The parallel-class rectangle produced by an exact q4488 partition forces
the two primitive q8 numerators into one signed residue type modulo eight. -/
theorem qEight_numerators_modEq_or_neg_of_four_four_eight_eight_failure
    (δ₄a δ₄b δ₈a δ₈b g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (hq8a : g / (Int.gcd δ₈a g : ℤ) = 8)
    (hq8b : g / (Int.gcd δ₈b g : ℤ) = 8)
    (a₄a a₄b a₈a a₈b : ℤ)
    (hnorm4a : SignedReducedNormalForm δ₄a g 4 a₄a)
    (hnorm4b : SignedReducedNormalForm δ₄b g 4 a₄b)
    (hnorm8a : SignedReducedNormalForm δ₈a g 8 a₈a)
    (hnorm8b : SignedReducedNormalForm δ₈b g 8 a₈b)
    (hfailure : ¬ HasFourDetunedGoodBranch δ₄a δ₄b δ₈a δ₈b g u) :
    a₈a ≡ a₈b [ZMOD 8] ∨ a₈a ≡ -a₈b [ZMOD 8] := by
  obtain ⟨hpartition, -, -, -, -⟩ :=
    four_four_eight_eight_failure_parallel_partition
      δ₄a δ₄b δ₈a δ₈b g u hg hq4a hq4b hq8a hq8b hfailure
  obtain ⟨
      ⟨c₄a, c₈aa, c₈ba, hc₄a, hc₈aa, hc₈ba,
        hpar8aa, hpar8ba, -⟩,
      ⟨c₄b, c₈ab, c₈bb, hc₄b, hc₈ab, hc₈bb,
        hpar8ab, hpar8bb, -⟩⟩ :=
    two_matching_codes_of_four_four_eight_eight_failure
      δ₄a δ₄b δ₈a δ₈b g u hg hq4a hq4b hq8a hq8b
        a₄a a₄b a₈a a₈b hnorm4a hnorm4b hnorm8a hnorm8b hfailure
  have hrowEq4a := detunedBadBranches_eq_normalizedBadBranches
    δ₄a g hnorm4a.divisor 4 hnorm4a.sign a₄a u
      hnorm4a.divisor_ne (by norm_num) hnorm4a.delta_eq
      hnorm4a.modulus_eq hnorm4a.sign_eq
  have hrowEq4b := detunedBadBranches_eq_normalizedBadBranches
    δ₄b g hnorm4b.divisor 4 hnorm4b.sign a₄b u
      hnorm4b.divisor_ne (by norm_num) hnorm4b.delta_eq
      hnorm4b.modulus_eq hnorm4b.sign_eq
  have hrowEq8a := detunedBadBranches_eq_normalizedBadBranches
    δ₈a g hnorm8a.divisor 8 hnorm8a.sign a₈a u
      hnorm8a.divisor_ne (by norm_num) hnorm8a.delta_eq
      hnorm8a.modulus_eq hnorm8a.sign_eq
  have hrowEq8b := detunedBadBranches_eq_normalizedBadBranches
    δ₈b g hnorm8b.divisor 8 hnorm8b.sign a₈b u
      hnorm8b.divisor_ne (by norm_num) hnorm8b.delta_eq
      hnorm8b.modulus_eq hnorm8b.sign_eq
  have hc4Outside : ¬ c₄b ≡ c₄a [ZMOD 4] := by
    intro hmod
    have hc₄b4a : c₄b ∈ detunedBadBranches δ₄a g u := by
      rw [hrowEq4a]
      apply normalizedBadBranches_mem_of_modEq a₄a 4 g u (by norm_num)
      · rw [← hrowEq4a]
        exact hc₄a
      · exact detunedBadBranches_subset_Ico δ₄b g u hc₄b
      · exact hmod.symm
    exact (Finset.disjoint_left.mp hpartition.disjoint₁₂) hc₄b4a hc₄b
  have hc₈aaOutside4a : ¬ c₈aa ≡ c₄a [ZMOD 4] := by
    intro hmod
    have hc₈aa4a : c₈aa ∈ detunedBadBranches δ₄a g u := by
      rw [hrowEq4a]
      apply normalizedBadBranches_mem_of_modEq a₄a 4 g u (by norm_num)
      · rw [← hrowEq4a]
        exact hc₄a
      · exact detunedBadBranches_subset_Ico δ₈a g u hc₈aa
      · exact hmod.symm
    exact (Finset.disjoint_left.mp hpartition.disjoint₁₃) hc₈aa4a hc₈aa
  have hc4OppositeParity : ¬ c₄a ≡ c₄b [ZMOD 2] := by
    intro hsame
    have hc₈aaMod4c₄b : c₈aa ≡ c₄b [ZMOD 4] :=
      modEq_four_of_same_parity_outside c₄a c₈aa c₄b
        hpar8aa hsame.symm hc₈aaOutside4a hc4Outside
    have hc₈aa4b : c₈aa ∈ detunedBadBranches δ₄b g u := by
      rw [hrowEq4b]
      apply normalizedBadBranches_mem_of_modEq a₄b 4 g u (by norm_num)
      · rw [← hrowEq4b]
        exact hc₄b
      · exact detunedBadBranches_subset_Ico δ₈a g u hc₈aa
      · exact hc₈aaMod4c₄b.symm
    exact (Finset.disjoint_left.mp hpartition.disjoint₂₃) hc₈aa4b hc₈aa
  have outside8aa : ¬ c₈aa ≡ c₄a [ZMOD 4] := hc₈aaOutside4a
  have outside8ba : ¬ c₈ba ≡ c₄a [ZMOD 4] := by
    intro hmod
    have hc₈ba4a : c₈ba ∈ detunedBadBranches δ₄a g u := by
      rw [hrowEq4a]
      apply normalizedBadBranches_mem_of_modEq a₄a 4 g u (by norm_num)
      · rw [← hrowEq4a]
        exact hc₄a
      · exact detunedBadBranches_subset_Ico δ₈b g u hc₈ba
      · exact hmod.symm
    exact (Finset.disjoint_left.mp hpartition.disjoint₁₄) hc₈ba4a hc₈ba
  have outside8ab : ¬ c₈ab ≡ c₄b [ZMOD 4] := by
    intro hmod
    have hc₈ab4b : c₈ab ∈ detunedBadBranches δ₄b g u := by
      rw [hrowEq4b]
      apply normalizedBadBranches_mem_of_modEq a₄b 4 g u (by norm_num)
      · rw [← hrowEq4b]
        exact hc₄b
      · exact detunedBadBranches_subset_Ico δ₈a g u hc₈ab
      · exact hmod.symm
    exact (Finset.disjoint_left.mp hpartition.disjoint₂₃) hc₈ab4b hc₈ab
  have outside8bb : ¬ c₈bb ≡ c₄b [ZMOD 4] := by
    intro hmod
    have hc₈bb4b : c₈bb ∈ detunedBadBranches δ₄b g u := by
      rw [hrowEq4b]
      apply normalizedBadBranches_mem_of_modEq a₄b 4 g u (by norm_num)
      · rw [← hrowEq4b]
        exact hc₄b
      · exact detunedBadBranches_subset_Ico δ₈b g u hc₈bb
      · exact hmod.symm
    exact (Finset.disjoint_left.mp hpartition.disjoint₂₄) hc₈bb4b hc₈bb
  have hfourA : c₈aa ≡ c₈ba [ZMOD 4] :=
    modEq_four_of_same_parity_outside c₄a c₈aa c₈ba
      hpar8aa hpar8ba outside8aa outside8ba
  have hfourB : c₈ab ≡ c₈bb [ZMOD 4] :=
    modEq_four_of_same_parity_outside c₄b c₈ab c₈bb
      hpar8ab hpar8bb outside8ab outside8bb
  have heightA : ¬ c₈aa ≡ c₈ba [ZMOD 8] := by
    intro hmod
    have hc₈ba8a : c₈ba ∈ detunedBadBranches δ₈a g u := by
      rw [hrowEq8a]
      apply normalizedBadBranches_mem_of_modEq a₈a 8 g u (by norm_num)
      · rw [← hrowEq8a]
        exact hc₈aa
      · exact detunedBadBranches_subset_Ico δ₈b g u hc₈ba
      · exact hmod
    exact (Finset.disjoint_left.mp hpartition.disjoint₃₄) hc₈ba8a hc₈ba
  have heightB : ¬ c₈ab ≡ c₈bb [ZMOD 8] := by
    intro hmod
    have hc₈bb8a : c₈bb ∈ detunedBadBranches δ₈a g u := by
      rw [hrowEq8a]
      apply normalizedBadBranches_mem_of_modEq a₈a 8 g u (by norm_num)
      · rw [← hrowEq8a]
        exact hc₈ab
      · exact detunedBadBranches_subset_Ico δ₈b g u hc₈bb
      · exact hmod
    exact (Finset.disjoint_left.mp hpartition.disjoint₃₄) hc₈bb8a hc₈bb
  have hcrossA : c₈aa ≡ c₈ba + 4 [ZMOD 8] :=
    modEq_eight_add_four_of_modEq_four_not_modEq_eight
      c₈aa c₈ba hfourA heightA
  have hcrossB : c₈ab ≡ c₈bb + 4 [ZMOD 8] :=
    modEq_eight_add_four_of_modEq_four_not_modEq_eight
      c₈ab c₈bb hfourB heightB
  have hparallel : c₈aa - c₈ab ≡ c₈ba - c₈bb [ZMOD 8] := by
    have hsub := hcrossA.sub hcrossB
    convert hsub using 1
    all_goals ring
  have hopposite : ¬ c₈aa ≡ c₈ab [ZMOD 2] := by
    intro hsame
    exact hc4OppositeParity (hpar8aa.symm.trans (hsame.trans hpar8ab))
  have hc₈aaNorm : c₈aa ∈ normalizedBadBranches a₈a 8 g u := by
    rw [← hrowEq8a]
    exact hc₈aa
  have hc₈abNorm : c₈ab ∈ normalizedBadBranches a₈a 8 g u := by
    rw [← hrowEq8a]
    exact hc₈ab
  have hc₈baNorm : c₈ba ∈ normalizedBadBranches a₈b 8 g u := by
    rw [← hrowEq8b]
    exact hc₈ba
  have hc₈bbNorm : c₈bb ∈ normalizedBadBranches a₈b 8 g u := by
    rw [← hrowEq8b]
    exact hc₈bb
  exact qEight_numerators_modEq_or_neg_of_parallel_rectangle
    a₈a a₈b g u hnorm8a.coprime hnorm8b.coprime
      hc₈aaNorm hc₈abNorm hc₈baNorm hc₈bbNorm hopposite hparallel

def SignedReducedNormalForm.negNumerator
    {δ g denominator numerator : ℤ}
    (hnorm : SignedReducedNormalForm δ g denominator numerator) :
    SignedReducedNormalForm δ g denominator (-numerator) := by
  refine {
    divisor := hnorm.divisor
    sign := -hnorm.sign
    divisor_ne := hnorm.divisor_ne
    delta_eq := ?_
    modulus_eq := hnorm.modulus_eq
    sign_eq := ?_
    coprime := hnorm.coprime.neg_left }
  · calc
      δ = hnorm.divisor * (hnorm.sign * numerator) := hnorm.delta_eq
      _ = hnorm.divisor * (-hnorm.sign * -numerator) := by ring
  · rcases hnorm.sign_eq with hsign | hsign
    · right
      simp [hsign]
    · left
      simp [hsign]

theorem modEq_four_or_neg_modEq_four_of_coprime
    (numerator residue : ℤ)
    (hnumerator : IsCoprime numerator 4)
    (hresidue : IsCoprime residue 8) :
    numerator ≡ residue [ZMOD 4] ∨
      -numerator ≡ residue [ZMOD 4] := by
  have hnumeratorOdd : Odd numerator := by
    apply Int.isCoprime_two_right.mp
    have hcoprime' : IsCoprime numerator ((2 : ℤ) * 2) := by
      norm_num
      exact hnumerator
    exact hcoprime'.of_mul_right_right
  have hresidueOdd : Odd residue := by
    apply Int.isCoprime_two_right.mp
    have hcoprime' : IsCoprime residue ((4 : ℤ) * 2) := by
      norm_num
      exact hresidue
    exact hcoprime'.of_mul_right_right
  have hnumeratorParity : numerator % 4 % 2 = 1 := by
    rw [Int.emod_emod_of_dvd numerator (by norm_num : (2 : ℤ) ∣ 4)]
    exact Int.odd_iff.mp hnumeratorOdd
  have hresidueParity : residue % 4 % 2 = 1 := by
    rw [Int.emod_emod_of_dvd residue (by norm_num : (2 : ℤ) ∣ 4)]
    exact Int.odd_iff.mp hresidueOdd
  have hnumerator0 : (0 : ℤ) ≤ numerator % 4 :=
    Int.emod_nonneg numerator (by norm_num)
  have hnumerator4 : numerator % 4 < 4 :=
    Int.emod_lt_of_pos numerator (by norm_num)
  have hresidue0 : (0 : ℤ) ≤ residue % 4 :=
    Int.emod_nonneg residue (by norm_num)
  have hresidue4 : residue % 4 < 4 :=
    Int.emod_lt_of_pos residue (by norm_num)
  have hnumeratorCases : numerator % 4 = 1 ∨ numerator % 4 = 3 := by
    have hmod2nonneg : (0 : ℤ) ≤ numerator % 4 % 2 :=
      Int.emod_nonneg _ (by norm_num)
    have hmod2lt : numerator % 4 % 2 < 2 :=
      Int.emod_lt_of_pos _ (by norm_num)
    omega
  have hresidueCases : residue % 4 = 1 ∨ residue % 4 = 3 := by
    have hmod2nonneg : (0 : ℤ) ≤ residue % 4 % 2 :=
      Int.emod_nonneg _ (by norm_num)
    have hmod2lt : residue % 4 % 2 < 2 :=
      Int.emod_lt_of_pos _ (by norm_num)
    omega
  have hnumeratorDecomp := Int.emod_add_mul_ediv numerator 4
  have hresidueDecomp := Int.emod_add_mul_ediv residue 4
  rcases hnumeratorCases with hnumeratorOne | hnumeratorThree <;>
    rcases hresidueCases with hresidueOne | hresidueThree
  · left
    change numerator % 4 = residue % 4
    omega
  · right
    rw [Int.modEq_iff_dvd]
    refine ⟨residue / 4 + numerator / 4 + 1, by omega⟩
  · right
    rw [Int.modEq_iff_dvd]
    refine ⟨residue / 4 + numerator / 4 + 1, by omega⟩
  · left
    change numerator % 4 = residue % 4
    omega

/-- Any primitive normal forms can be signed, after one failing q4488 phase,
so that all four numerators share the residue type needed by the scalar wall.
The two signed q4 numerators remain distinct. -/
theorem exists_common_signed_normal_forms_of_failure
    (δ₄a δ₄b δ₈a δ₈b g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (hq8a : g / (Int.gcd δ₈a g : ℤ) = 8)
    (hq8b : g / (Int.gcd δ₈b g : ℤ) = 8)
    (a₄a a₄b a₈a a₈b : ℤ)
    (hnorm4a : SignedReducedNormalForm δ₄a g 4 a₄a)
    (hnorm4b : SignedReducedNormalForm δ₄b g 4 a₄b)
    (hnorm8a : SignedReducedNormalForm δ₈a g 8 a₈a)
    (hnorm8b : SignedReducedNormalForm δ₈b g 8 a₈b)
    (hfailure : ¬ HasFourDetunedGoodBranch δ₄a δ₄b δ₈a δ₈b g u) :
    ∃ signed₄a signed₄b signed₈a signed₈b residue : ℤ,
      ∃ _hnorm4a' : SignedReducedNormalForm δ₄a g 4 signed₄a,
      ∃ _hnorm4b' : SignedReducedNormalForm δ₄b g 4 signed₄b,
      ∃ _hnorm8a' : SignedReducedNormalForm δ₈a g 8 signed₈a,
      ∃ _hnorm8b' : SignedReducedNormalForm δ₈b g 8 signed₈b,
      (4 : ℤ) ∣ signed₄a - residue ∧
      (4 : ℤ) ∣ signed₄b - residue ∧
      (8 : ℤ) ∣ signed₈a - residue ∧
      (8 : ℤ) ∣ signed₈b - residue ∧
      signed₄a ≠ signed₄b := by
  have hcommon :=
    qEight_numerators_modEq_or_neg_of_four_four_eight_eight_failure
      δ₄a δ₄b δ₈a δ₈b g u hg hq4a hq4b hq8a hq8b
        a₄a a₄b a₈a a₈b hnorm4a hnorm4b hnorm8a hnorm8b hfailure
  obtain ⟨signed₈b, hnorm8b', hres8b⟩ :
      ∃ signed₈b : ℤ,
        ∃ hnorm8b' : SignedReducedNormalForm δ₈b g 8 signed₈b,
        (8 : ℤ) ∣ signed₈b - a₈a := by
    rcases hcommon with hsame | hopposite
    · exact ⟨a₈b, hnorm8b, hsame.dvd⟩
    · exact ⟨-a₈b, hnorm8b.negNumerator, hopposite.dvd⟩
  have alignFour
      (δ numerator : ℤ)
      (hnorm : SignedReducedNormalForm δ g 4 numerator) :
      ∃ signed : ℤ,
        ∃ hnorm' : SignedReducedNormalForm δ g 4 signed,
        (4 : ℤ) ∣ signed - a₈a := by
    rcases modEq_four_or_neg_modEq_four_of_coprime
        numerator a₈a hnorm.coprime hnorm8a.coprime with
      hsame | hopposite
    · refine ⟨numerator, hnorm, ?_⟩
      simpa only [neg_sub] using dvd_neg.mpr hsame.dvd
    · refine ⟨-numerator, hnorm.negNumerator, ?_⟩
      simpa only [neg_sub] using dvd_neg.mpr hopposite.dvd
  obtain ⟨signed₄a, hnorm4a', hres4a⟩ := alignFour δ₄a a₄a hnorm4a
  obtain ⟨signed₄b, hnorm4b', hres4b⟩ := alignFour δ₄b a₄b hnorm4b
  have hdistinct : signed₄a ≠ signed₄b :=
    normalized_qFour_numerators_ne_of_failure
      δ₄a δ₄b δ₈a δ₈b g u hg hq4a hq4b hq8a hq8b
        signed₄a signed₄b hnorm4a' hnorm4b' hfailure
  exact ⟨signed₄a, signed₄b, a₈a, signed₈b, a₈a,
    hnorm4a', hnorm4b', hnorm8a, hnorm8b',
    hres4a, hres4b, by simp, hres8b, hdistinct⟩

def canonicalReducedNumerator (δ g : ℤ) : ℤ :=
  δ / (Int.gcd δ g : ℤ)

def canonicalReducedNormalForm
    (δ g denominator : ℤ) (hg : 0 < g)
    (hdenominator : g / (Int.gcd δ g : ℤ) = denominator) :
    SignedReducedNormalForm δ g denominator
      (canonicalReducedNumerator δ g) := by
  let divisor : ℤ := (Int.gcd δ g : ℤ)
  have hdivisorPos : 0 < divisor := by
    dsimp [divisor]
    have hne : Int.gcd δ g ≠ 0 := by
      intro hzero
      rw [Int.gcd_eq_zero_iff] at hzero
      omega
    exact_mod_cast Nat.pos_of_ne_zero hne
  have hdivisorNe : divisor ≠ 0 := ne_of_gt hdivisorPos
  have hdivisorDelta : divisor ∣ δ := by
    dsimp [divisor]
    exact Int.gcd_dvd_left δ g
  have hdivisorG : divisor ∣ g := by
    dsimp [divisor]
    exact Int.gcd_dvd_right δ g
  have hdeltaFactor :
      δ = divisor * canonicalReducedNumerator δ g := by
    dsimp [canonicalReducedNumerator, divisor]
    exact (Int.mul_ediv_cancel' (Int.gcd_dvd_left δ g)).symm
  have hgFactor : g = divisor * denominator := by
    calc
      g = divisor * (g / divisor) := (Int.mul_ediv_cancel' hdivisorG).symm
      _ = divisor * denominator := by
        dsimp [divisor]
        rw [hdenominator]
  refine {
    divisor := divisor
    sign := 1
    divisor_ne := hdivisorNe
    delta_eq := ?_
    modulus_eq := hgFactor
    sign_eq := Or.inl rfl
    coprime := ?_ }
  · simpa using hdeltaFactor
  · let bezoutA : ℤ := Int.gcdA δ g
    let bezoutB : ℤ := Int.gcdB δ g
    have hbezout : divisor = δ * bezoutA + g * bezoutB := by
      dsimp [divisor, bezoutA, bezoutB]
      exact Int.gcd_eq_gcd_ab δ g
    refine ⟨bezoutA, bezoutB, ?_⟩
    apply mul_left_cancel₀ hdivisorNe
    calc
      divisor *
          (bezoutA * canonicalReducedNumerator δ g +
            bezoutB * denominator) =
          (divisor * canonicalReducedNumerator δ g) * bezoutA +
            (divisor * denominator) * bezoutB := by ring
      _ = δ * bezoutA + g * bezoutB := by
        rw [← hdeltaFactor, ← hgFactor]
      _ = divisor := hbezout.symm
      _ = divisor * 1 := by ring

/-- Denominator data alone supplies the phase-independent primitive forms;
one failing phase then signs all four into the common q4488 residue type. -/
theorem exists_common_signed_normal_forms_of_failure_from_denominators
    (δ₄a δ₄b δ₈a δ₈b g : ℤ) (u : ℝ) (hg : 2 ≤ g)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (hq8a : g / (Int.gcd δ₈a g : ℤ) = 8)
    (hq8b : g / (Int.gcd δ₈b g : ℤ) = 8)
    (hfailure : ¬ HasFourDetunedGoodBranch δ₄a δ₄b δ₈a δ₈b g u) :
    ∃ signed₄a signed₄b signed₈a signed₈b residue : ℤ,
      ∃ _hnorm4a : SignedReducedNormalForm δ₄a g 4 signed₄a,
      ∃ _hnorm4b : SignedReducedNormalForm δ₄b g 4 signed₄b,
      ∃ _hnorm8a : SignedReducedNormalForm δ₈a g 8 signed₈a,
      ∃ _hnorm8b : SignedReducedNormalForm δ₈b g 8 signed₈b,
      (4 : ℤ) ∣ signed₄a - residue ∧
      (4 : ℤ) ∣ signed₄b - residue ∧
      (8 : ℤ) ∣ signed₈a - residue ∧
      (8 : ℤ) ∣ signed₈b - residue ∧
      signed₄a ≠ signed₄b := by
  exact exists_common_signed_normal_forms_of_failure
    δ₄a δ₄b δ₈a δ₈b g u hg hq4a hq4b hq8a hq8b
      (canonicalReducedNumerator δ₄a g)
      (canonicalReducedNumerator δ₄b g)
      (canonicalReducedNumerator δ₈a g)
      (canonicalReducedNumerator δ₈b g)
      (canonicalReducedNormalForm δ₄a g 4 (by omega) hq4a)
      (canonicalReducedNormalForm δ₄b g 4 (by omega) hq4b)
      (canonicalReducedNormalForm δ₈a g 8 (by omega) hq8a)
      (canonicalReducedNormalForm δ₈b g 8 (by omega) hq8b)
      hfailure

def FourEightEightMatchingAt
    (a₄ a₈a a₈b : ℤ) (u : ℝ) : Prop :=
  ∃ c₄ c₈a c₈b : ℤ,
    (8 : ℤ) ∣ -2 * c₄ + c₈a + c₈b ∧
    (∃ n : ℤ,
      |(a₄ : ℝ) * ((u + (c₄ : ℝ)) / 4) - n| < 1 / 14) ∧
    (∃ n : ℤ,
      |(a₈a : ℝ) * ((u + (c₈a : ℝ)) / 8) - n| < 1 / 14) ∧
    (∃ n : ℤ,
      |(a₈b : ℝ) * ((u + (c₈b : ℝ)) / 8) - n| < 1 / 14)

def TwoFourEightEightFrequencyEscapeOrRelation
    (a₄a a₄b a₈a a₈b : ℤ) (B : ℝ) : Prop :=
  (fourEightEightPhaseFrequency a₄a a₈a a₈b ≠ 0 ∧
    15 * B ≤ 2 *
      |(fourEightEightPhaseFrequency a₄a a₈a a₈b : ℝ)|) ∨
  (fourEightEightPhaseFrequency a₄b a₈a a₈b ≠ 0 ∧
    15 * B ≤ 2 *
      |(fourEightEightPhaseFrequency a₄b a₈a a₈b : ℝ)|) ∨
  (((-2 * a₄a + a₈a + a₈b = 0) ∨
      (fourEightEightPhaseFrequency a₄a a₈a a₈b ≠ 0 ∧
        2 * |(fourEightEightPhaseFrequency a₄a a₈a a₈b : ℝ)| <
          15 * B)) ∧
    ((-2 * a₄b + a₈a + a₈b = 0) ∨
      (fourEightEightPhaseFrequency a₄b a₈a a₈b ≠ 0 ∧
        2 * |(fourEightEightPhaseFrequency a₄b a₈a a₈b : ℝ)| <
          15 * B)) ∧
    (fourEightEightPhaseFrequency a₄a a₈a a₈b ≠ 0 ∨
      fourEightEightPhaseFrequency a₄b a₈a a₈b ≠ 0))

/-- **Unconditional q4488 producer.**  A single failing phase fixes common
signed primitive numerators using only the four denominator equations.  Those
same numerators produce both matching walls at every failing phase.  Jointly,
one q4-derived frequency crosses the sharp escape threshold, or the exact
residual consists of two zero/small signed relations with at least one
nonzero. -/
theorem exists_unconditional_q4488_failure_producer
    (δ₄a δ₄b δ₈a δ₈b g : ℤ) (anchor : ℝ) (hg : 2 ≤ g)
    (hq4a : g / (Int.gcd δ₄a g : ℤ) = 4)
    (hq4b : g / (Int.gcd δ₄b g : ℤ) = 4)
    (hq8a : g / (Int.gcd δ₈a g : ℤ) = 8)
    (hq8b : g / (Int.gcd δ₈b g : ℤ) = 8)
    (B : ℝ) (hB0 : 0 < B)
    (hfailure : ¬
      HasFourDetunedGoodBranch δ₄a δ₄b δ₈a δ₈b g anchor) :
    ∃ a₄a a₄b a₈a a₈b residue : ℤ,
      ∃ _hnorm4a : SignedReducedNormalForm δ₄a g 4 a₄a,
      ∃ _hnorm4b : SignedReducedNormalForm δ₄b g 4 a₄b,
      ∃ _hnorm8a : SignedReducedNormalForm δ₈a g 8 a₈a,
      ∃ _hnorm8b : SignedReducedNormalForm δ₈b g 8 a₈b,
      (4 : ℤ) ∣ a₄a - residue ∧
      (4 : ℤ) ∣ a₄b - residue ∧
      (8 : ℤ) ∣ a₈a - residue ∧
      (8 : ℤ) ∣ a₈b - residue ∧
      a₄a ≠ a₄b ∧
      (∀ u : ℝ,
        ¬ HasFourDetunedGoodBranch δ₄a δ₄b δ₈a δ₈b g u →
        FourEightEightMatchingAt a₄a a₈a a₈b u ∧
          FourEightEightMatchingAt a₄b a₈a a₈b u) ∧
      TwoFourEightEightFrequencyEscapeOrRelation a₄a a₄b a₈a a₈b B := by
  obtain ⟨a₄a, a₄b, a₈a, a₈b, residue,
      hnorm4a, hnorm4b, hnorm8a, hnorm8b,
      hres4a, hres4b, hres8a, hres8b, hdistinct⟩ :=
    exists_common_signed_normal_forms_of_failure_from_denominators
      δ₄a δ₄b δ₈a δ₈b g anchor hg hq4a hq4b hq8a hq8b hfailure
  refine ⟨a₄a, a₄b, a₈a, a₈b, residue,
    hnorm4a, hnorm4b, hnorm8a, hnorm8b,
    hres4a, hres4b, hres8a, hres8b, hdistinct, ?_, ?_⟩
  · intro u hfail
    simpa [FourEightEightMatchingAt] using
      (two_matching_walls_of_four_four_eight_eight_failure
        δ₄a δ₄b δ₈a δ₈b g u hg hq4a hq4b hq8a hq8b
          a₄a a₄b a₈a a₈b hnorm4a hnorm4b hnorm8a hnorm8b hfail)
  · simpa [TwoFourEightEightFrequencyEscapeOrRelation] using
      (two_fourEightEight_frequencies_large_or_relation_small
        a₄a a₄b a₈a a₈b residue B hB0
          hres4a hres4b hres8a hres8b hdistinct)

#print axioms exists_qEight_badBranch_modEq_two
#print axioms two_matching_codes_of_four_four_eight_eight_failure
#print axioms two_matching_walls_of_four_four_eight_eight_failure
#print axioms normalized_qFour_numerators_ne_of_failure
#print axioms two_fourEightEight_frequencies_large_or_relation_small
#print axioms normalizedQEight_pair_unit_relation_of_opposite_parity
#print axioms qEight_numerators_modEq_or_neg_of_parallel_rectangle
#print axioms qEight_numerators_modEq_or_neg_of_four_four_eight_eight_failure
#print axioms exists_common_signed_normal_forms_of_failure
#print axioms canonicalReducedNormalForm
#print axioms exists_common_signed_normal_forms_of_failure_from_denominators
#print axioms exists_unconditional_q4488_failure_producer

end
end LRCPairTowerGapTwoProducer
end LonelyRunner
