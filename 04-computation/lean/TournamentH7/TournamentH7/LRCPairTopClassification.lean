import TournamentH7.LRCPairCovarianceKernel
import TournamentH7.LRCWeightedRatioLayer

/-!
Exact classification of the two highest negative Bernoulli pair superlevels.
After primitive gcd reduction, the strict ratio-11 superlevel has colors
`12` and `13`, while the strict ratio-12 superlevel has only color `13`.
Each fixed multiplicative color is a path forest on distinct magnitudes, giving
the concrete THM-954 count caps `24` and `12`.
-/

namespace LonelyRunner
namespace LRCPairTopClassification

open LRCPairCovarianceKernel
open LRCWeightedRatioLayer

def primitiveKernel (first second : ℕ) : ℚ :=
  (bernoulliResidue14 ((first : ℤ) - second) -
      bernoulliResidue14 ((first : ℤ) + second)) /
    ((first : ℚ) * second)

def primitiveDeficit (first second : ℕ) : ℚ :=
  max 0 (-primitiveKernel first second)

theorem pairDeficit_eq_primitiveDeficit (first second : ℤ) :
    pairDeficit first second =
      primitiveDeficit (reducedFirst first second) (reducedSecond first second) := by
  rfl

/-- Equal nonzero magnitudes have zero negative pair weight.  This does not
permit quotienting multiplicities in the path-count theorem: two repeated
magnitude classes can still create many cross edges of one ratio color. -/
theorem pairDeficit_eq_zero_of_natAbs_eq
    (first second : ℤ) (hfirst : first ≠ 0)
    (hequal : first.natAbs = second.natAbs) :
    pairDeficit first second = 0 := by
  have habsPos : 0 < first.natAbs := Int.natAbs_pos.mpr hfirst
  unfold pairDeficit pairKernel reducedFirst reducedSecond
  rw [← hequal]
  simp only [Nat.gcd_self, Nat.div_self habsPos]
  norm_num [bernoulliResidue14]

theorem bernoulli_difference_abs_le (first second : ℤ) :
    |bernoulliResidue14 (first - second) -
      bernoulliResidue14 (first + second)| ≤ (12 / 49 : ℚ) := by
  let lower := (first - second) % 14
  let upper := (first + second) % 14
  have hlower0 : 0 ≤ lower := Int.emod_nonneg _ (by norm_num)
  have hlower14 : lower < 14 := Int.emod_lt_of_pos _ (by norm_num)
  have hupper0 : 0 ≤ upper := Int.emod_nonneg _ (by norm_num)
  have hupper14 : upper < 14 := Int.emod_lt_of_pos _ (by norm_num)
  have hlowerDecomp :
      14 * ((first - second) / 14) + lower = first - second := by
    simpa [lower] using Int.mul_ediv_add_emod (first - second) 14
  have hupperDecomp :
      14 * ((first + second) / 14) + upper = first + second := by
    simpa [upper] using Int.mul_ediv_add_emod (first + second) 14
  unfold bernoulliResidue14
  change abs (
    ((lower : ℚ) ^ 2 / 196 - (lower : ℚ) / 14 + 1 / 6) -
    ((upper : ℚ) ^ 2 / 196 - (upper : ℚ) / 14 + 1 / 6)) ≤ 12 / 49
  interval_cases lower <;> interval_cases upper <;> norm_num <;> omega

theorem primitiveDeficit_le_envelope (first second : ℕ)
    (hfirst : 0 < first) (hsecond : 0 < second) :
    primitiveDeficit first second ≤
      (12 / 49 : ℚ) / ((first : ℚ) * second) := by
  have hdenominator : (0 : ℚ) < (first : ℚ) * second := by positivity
  have hnumerator := bernoulli_difference_abs_le (first : ℤ) (second : ℤ)
  unfold primitiveDeficit primitiveKernel
  rw [max_le_iff]
  constructor
  · positivity
  · rw [show -((bernoulliResidue14 ((first : ℤ) - second) -
          bernoulliResidue14 ((first : ℤ) + second)) /
          ((first : ℚ) * second)) =
        (-(bernoulliResidue14 ((first : ℤ) - second) -
          bernoulliResidue14 ((first : ℤ) + second))) /
          ((first : ℚ) * second) by ring]
    apply (div_le_div_iff_of_pos_right hdenominator).2
    calc
      -(bernoulliResidue14 ((first : ℤ) - second) -
          bernoulliResidue14 ((first : ℤ) + second))
          ≤ |bernoulliResidue14 ((first : ℤ) - second) -
              bernoulliResidue14 ((first : ℤ) + second)| := neg_le_abs _
      _ ≤ 12 / 49 := hnumerator

theorem product_lt_thirty_three_of_ratio11_lt
    (first second : ℕ) (hfirst : 0 < first) (hsecond : 0 < second)
    (habove : ratio11Level < primitiveDeficit first second) :
    first * second < 33 := by
  have hdenominator : (0 : ℚ) < (first : ℚ) * second := by positivity
  have hbound := habove.trans_le
    (primitiveDeficit_le_envelope first second hfirst hsecond)
  have hscaled := (lt_div_iff₀ hdenominator).mp hbound
  have hcast : ((first * second : ℕ) : ℚ) < 33 := by
    push_cast
    norm_num [ratio11Level] at hscaled ⊢
    nlinarith
  exact_mod_cast hcast

set_option maxHeartbeats 2000000 in
theorem primitive_top_colors
    (first second : ℕ) (hfirst : 0 < first) (hsecond : 0 < second)
    (hcoprime : Nat.Coprime first second)
    (habove : ratio11Level < primitiveDeficit first second) :
    (first = 1 ∧ second = 12) ∨ (first = 12 ∧ second = 1) ∨
      (first = 1 ∧ second = 13) ∨ (first = 13 ∧ second = 1) := by
  have hproduct := product_lt_thirty_three_of_ratio11_lt first second
    hfirst hsecond habove
  have hfirstLt : first < 33 := by
    calc
      first = first * 1 := by omega
      _ ≤ first * second := Nat.mul_le_mul_left first hsecond
      _ < 33 := hproduct
  have hsecondLt : second < 33 := by
    calc
      second = 1 * second := by omega
      _ ≤ first * second := Nat.mul_le_mul_right second hfirst
      _ < 33 := hproduct
  interval_cases first <;> interval_cases second <;>
    norm_num [primitiveDeficit, primitiveKernel, bernoulliResidue14,
      ratio11Level] at *

theorem primitive_one_top_color
    (first second : ℕ) (hfirst : 0 < first) (hsecond : 0 < second)
    (hcoprime : Nat.Coprime first second)
    (habove : ratio12Level < primitiveDeficit first second) :
    (first = 1 ∧ second = 13) ∨ (first = 13 ∧ second = 1) := by
  have hcut : ratio11Level < ratio12Level := by
    norm_num [ratio11Level, ratio12Level]
  rcases primitive_top_colors first second hfirst hsecond hcoprime
      (hcut.trans habove) with h12 | h12 | h13 | h13
  · rcases h12 with ⟨rfl, rfl⟩
    norm_num [primitiveDeficit, primitiveKernel, bernoulliResidue14,
      ratio12Level] at habove
  · rcases h12 with ⟨rfl, rfl⟩
    norm_num [primitiveDeficit, primitiveKernel, bernoulliResidue14,
      ratio12Level] at habove
  · exact Or.inl h13
  · exact Or.inr h13

theorem reducedFirst_pos (first second : ℤ) (hfirst : first ≠ 0) :
    0 < reducedFirst first second := by
  have habs : 0 < first.natAbs := Int.natAbs_pos.mpr hfirst
  have hgcdPos : 0 < Nat.gcd first.natAbs second.natAbs :=
    Nat.gcd_pos_of_pos_left second.natAbs habs
  exact Nat.div_pos (Nat.le_of_dvd habs (Nat.gcd_dvd_left _ _)) hgcdPos

theorem reducedSecond_pos (first second : ℤ) (hsecond : second ≠ 0) :
    0 < reducedSecond first second := by
  have habs : 0 < second.natAbs := Int.natAbs_pos.mpr hsecond
  have hgcdPos : 0 < Nat.gcd first.natAbs second.natAbs :=
    Nat.gcd_pos_of_pos_right first.natAbs habs
  exact Nat.div_pos (Nat.le_of_dvd habs (Nat.gcd_dvd_right _ _)) hgcdPos

theorem reduced_coprime (first second : ℤ) (hfirst : first ≠ 0) :
    Nat.Coprime (reducedFirst first second) (reducedSecond first second) := by
  unfold reducedFirst reducedSecond
  apply Nat.coprime_div_gcd_div_gcd
  exact Nat.gcd_pos_of_pos_left second.natAbs (Int.natAbs_pos.mpr hfirst)

theorem pairDeficit_top_colors
    (first second : ℤ) (hfirst : first ≠ 0) (hsecond : second ≠ 0)
    (habove : ratio11Level < pairDeficit first second) :
    second.natAbs = 12 * first.natAbs ∨
      first.natAbs = 12 * second.natAbs ∨
      second.natAbs = 13 * first.natAbs ∨
      first.natAbs = 13 * second.natAbs := by
  let common := Nat.gcd first.natAbs second.natAbs
  have hfirstFactor : reducedFirst first second * common = first.natAbs := by
    exact Nat.div_mul_cancel (Nat.gcd_dvd_left first.natAbs second.natAbs)
  have hsecondFactor : reducedSecond first second * common = second.natAbs := by
    exact Nat.div_mul_cancel (Nat.gcd_dvd_right first.natAbs second.natAbs)
  rw [pairDeficit_eq_primitiveDeficit] at habove
  rcases primitive_top_colors (reducedFirst first second)
      (reducedSecond first second) (reducedFirst_pos first second hfirst)
      (reducedSecond_pos first second hsecond)
      (reduced_coprime first second hfirst) habove with
    h12 | h12 | h13 | h13
  · left
    rcases h12 with ⟨hfirstReduced, hsecondReduced⟩
    simp only [hfirstReduced, hsecondReduced] at hfirstFactor hsecondFactor
    omega
  · right; left
    rcases h12 with ⟨hfirstReduced, hsecondReduced⟩
    simp only [hfirstReduced, hsecondReduced] at hfirstFactor hsecondFactor
    omega
  · right; right; left
    rcases h13 with ⟨hfirstReduced, hsecondReduced⟩
    simp only [hfirstReduced, hsecondReduced] at hfirstFactor hsecondFactor
    omega
  · right; right; right
    rcases h13 with ⟨hfirstReduced, hsecondReduced⟩
    simp only [hfirstReduced, hsecondReduced] at hfirstFactor hsecondFactor
    omega

theorem pairDeficit_one_top_color
    (first second : ℤ) (hfirst : first ≠ 0) (hsecond : second ≠ 0)
    (habove : ratio12Level < pairDeficit first second) :
    second.natAbs = 13 * first.natAbs ∨
      first.natAbs = 13 * second.natAbs := by
  let common := Nat.gcd first.natAbs second.natAbs
  have hfirstFactor : reducedFirst first second * common = first.natAbs := by
    exact Nat.div_mul_cancel (Nat.gcd_dvd_left first.natAbs second.natAbs)
  have hsecondFactor : reducedSecond first second * common = second.natAbs := by
    exact Nat.div_mul_cancel (Nat.gcd_dvd_right first.natAbs second.natAbs)
  rw [pairDeficit_eq_primitiveDeficit] at habove
  rcases primitive_one_top_color (reducedFirst first second)
      (reducedSecond first second) (reducedFirst_pos first second hfirst)
      (reducedSecond_pos first second hsecond)
      (reduced_coprime first second hfirst) habove with h13 | h13
  · left
    rcases h13 with ⟨hfirstReduced, hsecondReduced⟩
    simp only [hfirstReduced, hsecondReduced] at hfirstFactor hsecondFactor
    omega
  · right
    rcases h13 with ⟨hfirstReduced, hsecondReduced⟩
    simp only [hfirstReduced, hsecondReduced] at hfirstFactor hsecondFactor
    omega

def ratioPairs (magnitude : Fin 13 → ℤ) (ratio : ℤ) :
    Finset (Fin 13 × Fin 13) :=
  Finset.univ.filter fun pair =>
    magnitude pair.2 = ratio * magnitude pair.1

theorem ratioPairs_card_le_twelve
    (magnitude : Fin 13 → ℤ) (ratio : ℤ)
    (hpositive : ∀ index, 0 < magnitude index)
    (hratio : 1 < ratio)
    (hinjective : Function.Injective magnitude) :
    (ratioPairs magnitude ratio).card ≤ 12 := by
  obtain ⟨minimum, -, hminimum⟩ :=
    Finset.exists_min_image (Finset.univ : Finset (Fin 13)) magnitude
      ⟨0, Finset.mem_univ 0⟩
  have htarget : ∀ pair ∈ ratioPairs magnitude ratio,
      pair.2 ∈ Finset.univ.erase minimum := by
    intro pair hpair
    rw [ratioPairs, Finset.mem_filter] at hpair
    rw [Finset.mem_erase]
    refine ⟨?_, Finset.mem_univ _⟩
    intro heq
    have hle := hminimum pair.1 (Finset.mem_univ pair.1)
    have hstrict : magnitude pair.1 < ratio * magnitude pair.1 := by
      have hproduct : 0 < (ratio - 1) * magnitude pair.1 :=
        mul_pos (sub_pos.mpr hratio) (hpositive pair.1)
      nlinarith
    rw [heq] at hpair
    linarith
  have hinjectiveOn : ∀ first ∈ ratioPairs magnitude ratio,
      ∀ second ∈ ratioPairs magnitude ratio,
      first.2 = second.2 → first = second := by
    intro first hfirstPair second hsecondPair htargetEq
    rw [ratioPairs, Finset.mem_filter] at hfirstPair hsecondPair
    have hsourcesScaled : ratio * magnitude first.1 =
        ratio * magnitude second.1 := by
      rw [← hfirstPair.2, ← hsecondPair.2, htargetEq]
    have hratioNe : ratio ≠ 0 := by omega
    have hsources : magnitude first.1 = magnitude second.1 :=
      mul_left_cancel₀ hratioNe hsourcesScaled
    exact Prod.ext (hinjective hsources) htargetEq
  calc
    (ratioPairs magnitude ratio).card
        ≤ (Finset.univ.erase minimum).card :=
      Finset.card_le_card_of_injOn (fun pair => pair.2) htarget
        hinjectiveOn
    _ = 12 := by
      rw [Finset.card_erase_of_mem (Finset.mem_univ minimum)]
      simp

def ratioSupports (magnitude : Fin 13 → ℤ) (ratio : ℤ) :
    Finset PairSupport :=
  (ratioPairs magnitude ratio).image fun pair => {pair.1, pair.2}

theorem ratioSupports_card_le_twelve
    (magnitude : Fin 13 → ℤ) (ratio : ℤ)
    (hpositive : ∀ index, 0 < magnitude index)
    (hratio : 1 < ratio)
    (hinjective : Function.Injective magnitude) :
    (ratioSupports magnitude ratio).card ≤ 12 :=
  Finset.card_image_le.trans
    (ratioPairs_card_le_twelve magnitude ratio hpositive hratio hinjective)

theorem pair_support_mem_ratioSupports_of_natAbs_relation
    (v : Fin 13 → ℤ) (first second : Fin 13) (ratio : ℕ)
    (hrelation : (v second).natAbs = ratio * (v first).natAbs) :
    {first, second} ∈ ratioSupports (fun index => |v index|) ratio := by
  rw [ratioSupports, Finset.mem_image]
  refine ⟨(first, second), ?_, rfl⟩
  rw [ratioPairs, Finset.mem_filter]
  refine ⟨Finset.mem_univ _, ?_⟩
  have hrelationInt : ((v second).natAbs : ℤ) =
      (ratio : ℤ) * (v first).natAbs := by
    exact_mod_cast hrelation
  simpa only [Int.natCast_natAbs] using hrelationInt

theorem pair_support_mem_ratioSupports_of_reverse_natAbs_relation
    (v : Fin 13 → ℤ) (first second : Fin 13) (ratio : ℕ)
    (hrelation : (v first).natAbs = ratio * (v second).natAbs) :
    {first, second} ∈ ratioSupports (fun index => |v index|) ratio := by
  simpa [Finset.pair_comm] using
    pair_support_mem_ratioSupports_of_natAbs_relation v second first ratio hrelation

theorem ratio12_superlevel_subset_ratio13
    (v : Fin 13 → ℤ) (weight : PairSupport → ℚ)
    (hv : ∀ index, v index ≠ 0)
    (hweight : ∀ first second, first ≠ second →
      weight {first, second} = pairDeficit (v first) (v second)) :
    ∀ support ∈ pairSupports,
      ratio12Level < weight support →
        support ∈ ratioSupports (fun index => |v index|) 13 := by
  intro support hsupport habove
  rw [pairSupports, Finset.mem_powersetCard] at hsupport
  obtain ⟨_, hcard⟩ := hsupport
  obtain ⟨first, second, hne, rfl⟩ := Finset.card_eq_two.mp hcard
  rw [hweight first second hne] at habove
  rcases pairDeficit_one_top_color (v first) (v second) (hv first)
      (hv second) habove with hforward | hreverse
  · exact pair_support_mem_ratioSupports_of_natAbs_relation
      v first second 13 hforward
  · exact pair_support_mem_ratioSupports_of_reverse_natAbs_relation
      v first second 13 hreverse

theorem ratio11_superlevel_subset_ratio12_union_ratio13
    (v : Fin 13 → ℤ) (weight : PairSupport → ℚ)
    (hv : ∀ index, v index ≠ 0)
    (hweight : ∀ first second, first ≠ second →
      weight {first, second} = pairDeficit (v first) (v second)) :
    ∀ support ∈ pairSupports,
      ratio11Level < weight support →
        support ∈ ratioSupports (fun index => |v index|) 12 ∪
          ratioSupports (fun index => |v index|) 13 := by
  intro support hsupport habove
  rw [pairSupports, Finset.mem_powersetCard] at hsupport
  obtain ⟨_, hcard⟩ := hsupport
  obtain ⟨first, second, hne, rfl⟩ := Finset.card_eq_two.mp hcard
  rw [hweight first second hne] at habove
  rcases pairDeficit_top_colors (v first) (v second) (hv first)
      (hv second) habove with h12 | h12 | h13 | h13
  · rw [Finset.mem_union]
    exact Or.inl (pair_support_mem_ratioSupports_of_natAbs_relation
      v first second 12 h12)
  · rw [Finset.mem_union]
    exact Or.inl (pair_support_mem_ratioSupports_of_reverse_natAbs_relation
      v first second 12 h12)
  · rw [Finset.mem_union]
    exact Or.inr (pair_support_mem_ratioSupports_of_natAbs_relation
      v first second 13 h13)
  · rw [Finset.mem_union]
    exact Or.inr (pair_support_mem_ratioSupports_of_reverse_natAbs_relation
      v first second 13 h13)

/-- The two concrete top-layer count fields required by THM-954. -/
theorem pairDeficit_top_path_caps
    (v : Fin 13 → ℤ) (weight : PairSupport → ℚ)
    (hv : ∀ index, v index ≠ 0)
    (hdistinct : ∀ first second, first ≠ second →
      |v first| ≠ |v second|)
    (hweight : ∀ first second, first ≠ second →
      weight {first, second} = pairDeficit (v first) (v second)) :
    countAbove pairSupports weight ratio11Level ≤ 24 ∧
      countAbove pairSupports weight ratio12Level ≤ 12 := by
  let magnitude : Fin 13 → ℤ := fun index => |v index|
  have hpositive : ∀ index, 0 < magnitude index := fun index =>
    abs_pos.mpr (hv index)
  have hinjective : Function.Injective magnitude := by
    intro first second hequal
    by_contra hne
    exact hdistinct first second hne hequal
  have hone := ratio12_superlevel_subset_ratio13 v weight hv hweight
  have htwo := ratio11_superlevel_subset_ratio12_union_ratio13 v weight hv hweight
  constructor
  · calc
      countAbove pairSupports weight ratio11Level
          ≤ (ratioSupports magnitude 12 ∪ ratioSupports magnitude 13).card := by
        apply Finset.card_le_card
        intro support hsupport
        rw [Finset.mem_filter] at hsupport
        exact htwo support hsupport.1 hsupport.2
      _ ≤ (ratioSupports magnitude 12).card +
            (ratioSupports magnitude 13).card := Finset.card_union_le _ _
      _ ≤ 12 + 12 := Nat.add_le_add
        (ratioSupports_card_le_twelve magnitude 12 hpositive (by norm_num)
          hinjective)
        (ratioSupports_card_le_twelve magnitude 13 hpositive (by norm_num)
          hinjective)
      _ = 24 := by norm_num
  · apply (Finset.card_le_card ?_).trans
      (ratioSupports_card_le_twelve magnitude 13 hpositive (by norm_num)
        hinjective)
    intro support hsupport
    rw [Finset.mem_filter] at hsupport
    exact hone support hsupport.1 hsupport.2
#print axioms bernoulli_difference_abs_le
#print axioms pairDeficit_eq_zero_of_natAbs_eq
#print axioms primitiveDeficit_le_envelope
#print axioms product_lt_thirty_three_of_ratio11_lt
#print axioms primitive_top_colors
#print axioms primitive_one_top_color
#print axioms pairDeficit_top_colors
#print axioms pairDeficit_one_top_color
#print axioms ratioPairs_card_le_twelve
#print axioms ratio12_superlevel_subset_ratio13
#print axioms ratio11_superlevel_subset_ratio12_union_ratio13
#print axioms pairDeficit_top_path_caps

end LRCPairTopClassification
end LonelyRunner
