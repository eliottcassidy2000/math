import TournamentH7.LRCPairRatioQuotient
import TournamentH7.LRCWeightedRatioLayer

/-!
# Finite primitive-ratio covers for the THM-954 thresholds

The Bernoulli residue numerator is uniformly bounded by `12/49`.  After
dividing a nonzero speed pair by its gcd, every primitive pair above a positive
threshold therefore lies in a finite product box.  This file packages that
argument uniformly for all seven THM-954 thresholds and supplies a decidable
finite quotient graph for replaying clique certificates.
-/

namespace LonelyRunner
namespace LRCPairRatioFiniteCover

open Finset SimpleGraph
open LRCPairCovarianceKernel LRCPairRatioQuotient LRCWeightedRatioLayer

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
    exact (neg_le_abs _).trans hnumerator

def primitivePairCandidates (threshold : ℚ) (cap : ℕ) :
    Finset (ℕ × ℕ) :=
  (Finset.range (cap + 1)).biUnion fun first =>
    ((Finset.range (cap / first + 1)).filter fun second =>
      0 < first ∧ 0 < second ∧ Nat.Coprime first second ∧
        threshold < primitiveDeficit first second).image fun second =>
          (first, second)

def finiteAllowedRatios (threshold : ℚ) (cap : ℕ) : Finset ℚ :=
  (primitivePairCandidates threshold cap).image fun pair =>
    (pair.1 : ℚ) / pair.2

theorem product_lt_cap_succ
    (threshold : ℚ) (cap first second : ℕ)
    (hthreshold : 0 < threshold)
    (hcap : (12 / 49 : ℚ) ≤
      threshold * (((cap + 1 : ℕ) : ℚ)))
    (hfirst : 0 < first) (hsecond : 0 < second)
    (habove : threshold < primitiveDeficit first second) :
    first * second < cap + 1 := by
  have hdenominator : (0 : ℚ) < (first : ℚ) * second := by positivity
  have hbound := habove.trans_le
    (primitiveDeficit_le_envelope first second hfirst hsecond)
  have hscaled := (lt_div_iff₀ hdenominator).mp hbound
  have hmul : threshold * ((first : ℚ) * second) <
      threshold * ((cap + 1 : ℕ) : ℚ) := hscaled.trans_le hcap
  have hdenominatorLt : (first : ℚ) * second < ((cap + 1 : ℕ) : ℚ) :=
    (mul_lt_mul_iff_right₀ hthreshold).mp hmul
  have hcast : ((first * second : ℕ) : ℚ) < ((cap + 1 : ℕ) : ℚ) := by
    simpa only [Nat.cast_mul] using hdenominatorLt
  exact_mod_cast hcast

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

theorem speed_ratio_eq_reduced_ratio
    (first second : ℤ) (hfirst : first ≠ 0) :
    ((first.natAbs : ℚ) / second.natAbs) =
      (reducedFirst first second : ℚ) / reducedSecond first second := by
  let common := Nat.gcd first.natAbs second.natAbs
  have hcommon : 0 < common :=
    Nat.gcd_pos_of_pos_left second.natAbs (Int.natAbs_pos.mpr hfirst)
  have hfirstFactor : reducedFirst first second * common = first.natAbs :=
    Nat.div_mul_cancel (Nat.gcd_dvd_left first.natAbs second.natAbs)
  have hsecondFactor : reducedSecond first second * common = second.natAbs :=
    Nat.div_mul_cancel (Nat.gcd_dvd_right first.natAbs second.natAbs)
  rw [← hfirstFactor, ← hsecondFactor]
  push_cast
  exact mul_div_mul_right _ _ (by exact_mod_cast hcommon.ne')

/-- One exact envelope proof supplies a finite allowed-ratio cover at any
positive threshold once `cap + 1` clears `12 / (49 * threshold)`. -/
theorem finiteAllowedRatios_cover
    (threshold : ℚ) (cap : ℕ)
    (hthreshold : 0 < threshold)
    (hcap : (12 / 49 : ℚ) ≤
      threshold * (((cap + 1 : ℕ) : ℚ))) :
    ∀ ratio, ratioAllowed threshold ratio →
      ratio ∈ finiteAllowedRatios threshold cap := by
  intro ratio hallowed
  obtain ⟨first, second, hfirst, hsecond, rfl, habove⟩ := hallowed
  let reducedFirst' := reducedFirst first second
  let reducedSecond' := reducedSecond first second
  have hfirstPos : 0 < reducedFirst' := reducedFirst_pos first second hfirst
  have hsecondPos : 0 < reducedSecond' := reducedSecond_pos first second hsecond
  have hcoprime : Nat.Coprime reducedFirst' reducedSecond' :=
    reduced_coprime first second hfirst
  rw [pairDeficit_eq_primitiveDeficit] at habove
  have hproduct := product_lt_cap_succ threshold cap reducedFirst' reducedSecond'
    hthreshold hcap hfirstPos hsecondPos habove
  have hfirstLt : reducedFirst' < cap + 1 := by
    calc
      reducedFirst' = reducedFirst' * 1 := by omega
      _ ≤ reducedFirst' * reducedSecond' :=
        Nat.mul_le_mul_left reducedFirst' hsecondPos
      _ < cap + 1 := hproduct
  have hsecondLeDiv : reducedSecond' ≤ cap / reducedFirst' := by
    apply (Nat.le_div_iff_mul_le hfirstPos).2
    have hproductLe : reducedFirst' * reducedSecond' ≤ cap := by omega
    simpa [Nat.mul_comm] using hproductLe
  have hsecondLtHyperbola : reducedSecond' < cap / reducedFirst' + 1 := by
    omega
  rw [speed_ratio_eq_reduced_ratio first second hfirst]
  rw [finiteAllowedRatios, Finset.mem_image]
  refine ⟨(reducedFirst', reducedSecond'), ?_, rfl⟩
  rw [primitivePairCandidates, Finset.mem_biUnion]
  refine ⟨reducedFirst', Finset.mem_range.mpr hfirstLt, ?_⟩
  rw [Finset.mem_image]
  refine ⟨reducedSecond', ?_, rfl⟩
  rw [Finset.mem_filter]
  exact ⟨Finset.mem_range.mpr hsecondLtHyperbola,
    hfirstPos, hsecondPos, hcoprime, habove⟩

def tau3FiniteRatios : Finset ℚ := finiteAllowedRatios tau3 53
def tau4FiniteRatios : Finset ℚ := finiteAllowedRatios tau4 129
def tau5FiniteRatios : Finset ℚ := finiteAllowedRatios tau5 431
def tau6FiniteRatios : Finset ℚ := finiteAllowedRatios tau6 767
def tau7FiniteRatios : Finset ℚ := finiteAllowedRatios tau7 1843
def tau8FiniteRatios : Finset ℚ := finiteAllowedRatios tau8 2447
def tau9FiniteRatios : Finset ℚ := finiteAllowedRatios tau9 3455

theorem tier_envelope_caps :
    (12 / 49 : ℚ) ≤ tau3 * (53 + 1) ∧
    (12 / 49 : ℚ) ≤ tau4 * (129 + 1) ∧
    (12 / 49 : ℚ) ≤ tau5 * (431 + 1) ∧
    (12 / 49 : ℚ) ≤ tau6 * (767 + 1) ∧
    (12 / 49 : ℚ) ≤ tau7 * (1843 + 1) ∧
    (12 / 49 : ℚ) ≤ tau8 * (2447 + 1) ∧
    (12 / 49 : ℚ) ≤ tau9 * (3455 + 1) := by
  norm_num [tau3, tau4, tau5, tau6, tau7, tau8, tau9]

set_option maxRecDepth 100000 in
theorem tau3FiniteRatios_cover :
    ∀ ratio, ratioAllowed tau3 ratio → ratio ∈ tau3FiniteRatios :=
  finiteAllowedRatios_cover tau3 53
    (by norm_num [tau3]) (by norm_num [tau3])

set_option maxRecDepth 100000 in
theorem tau4FiniteRatios_cover :
    ∀ ratio, ratioAllowed tau4 ratio → ratio ∈ tau4FiniteRatios :=
  finiteAllowedRatios_cover tau4 129
    (by norm_num [tau4]) (by norm_num [tau4])

set_option maxRecDepth 100000 in
theorem tau5FiniteRatios_cover :
    ∀ ratio, ratioAllowed tau5 ratio → ratio ∈ tau5FiniteRatios :=
  finiteAllowedRatios_cover tau5 431
    (by norm_num [tau5]) (by norm_num [tau5])

set_option maxRecDepth 100000 in
theorem tau6FiniteRatios_cover :
    ∀ ratio, ratioAllowed tau6 ratio → ratio ∈ tau6FiniteRatios :=
  finiteAllowedRatios_cover tau6 767
    (by norm_num [tau6]) (by norm_num [tau6])

set_option maxRecDepth 100000 in
theorem tau7FiniteRatios_cover :
    ∀ ratio, ratioAllowed tau7 ratio → ratio ∈ tau7FiniteRatios :=
  finiteAllowedRatios_cover tau7 1843
    (by norm_num [tau7]) (by norm_num [tau7])

set_option maxRecDepth 100000 in
theorem tau8FiniteRatios_cover :
    ∀ ratio, ratioAllowed tau8 ratio → ratio ∈ tau8FiniteRatios :=
  finiteAllowedRatios_cover tau8 2447
    (by norm_num [tau8]) (by norm_num [tau8])

set_option maxRecDepth 100000 in
theorem tau9FiniteRatios_cover :
    ∀ ratio, ratioAllowed tau9 ratio → ratio ∈ tau9FiniteRatios :=
  finiteAllowedRatios_cover tau9 3455
    (by norm_num [tau9]) (by norm_num [tau9])

/-- Decidable quotient supergraph induced by a finite ratio cover. -/
def finiteCoverGraph (vertices : Finset ℚ) : SimpleGraph ℚ where
  Adj first second := first ≠ second ∧
    (first / second ∈ vertices ∨ second / first ∈ vertices)
  symm first second hadj := by
    rcases hadj with ⟨hne, hquotient | hquotient⟩
    · exact ⟨hne.symm, Or.inr hquotient⟩
    · exact ⟨hne.symm, Or.inl hquotient⟩
  loopless := ⟨fun _ hadj => hadj.1 rfl⟩

instance finiteCoverGraph_decidableAdj (vertices : Finset ℚ) :
    DecidableRel (finiteCoverGraph vertices).Adj := by
  intro first second
  change Decidable (first ≠ second ∧
    (first / second ∈ vertices ∨ second / first ∈ vertices))
  infer_instance

theorem allowedRatioGraph_le_finiteCoverGraph
    (threshold : ℚ) (vertices : Finset ℚ)
    (hcover : ∀ ratio, ratioAllowed threshold ratio → ratio ∈ vertices) :
    allowedRatioGraph threshold ≤ finiteCoverGraph vertices := by
  intro first second hadj
  rcases hadj.2.2.2 with hquotient | hquotient
  · exact ⟨hadj.1, Or.inl (hcover _ hquotient)⟩
  · exact ⟨hadj.1, Or.inr (hcover _ hquotient)⟩

/-- Soundness of an exhaustive powerset replay on a finite vertex list. -/
theorem cliqueFreeOn_of_powerset_check
    (graph : SimpleGraph ℚ) (vertices : Finset ℚ) (cliqueSize : ℕ)
    (hcheck : ∀ clique ∈ vertices.powersetCard cliqueSize,
      ¬ graph.IsClique clique) :
    graph.CliqueFreeOn (vertices : Set ℚ) cliqueSize := by
  intro clique hsubset hclique
  apply hcheck clique
  rw [Finset.mem_powersetCard]
  exact ⟨Finset.coe_subset.mp hsubset, hclique.card_eq⟩
  exact hclique.isClique

theorem finiteCoverGraph_cliqueFreeOn_two_of_quotient_free
    (vertices : Finset ℚ)
    (hfree : (vertices : Set ℚ).Pairwise fun first second =>
      first / second ∉ vertices ∧ second / first ∉ vertices) :
    (finiteCoverGraph vertices).CliqueFreeOn (vertices : Set ℚ) 2 := by
  intro clique hsubset hclique
  obtain ⟨first, second, hne, hpair⟩ :=
    Finset.card_eq_two.mp hclique.card_eq
  subst clique
  have hfirst : first ∈ vertices := hsubset (by simp)
  have hsecond : second ∈ vertices := hsubset (by simp)
  have hadj : (finiteCoverGraph vertices).Adj first second :=
    hclique.isClique (by simp) (by simp) hne
  rcases hadj.2 with hquotient | hquotient
  · exact (hfree hfirst hsecond hne).1 hquotient
  · exact (hfree hfirst hsecond hne).2 hquotient

/-- The fourteen exact tau3 ratios, retained here as a small graph replay
certificate independently of the heavier primitive-pair enumeration. -/
def tau3CertificateRatios : Finset ℚ :=
  {1 / 10, 1 / 11, 1 / 12, 1 / 13, 2 / 11, 3 / 10, 3 / 11,
    10 / 1, 10 / 3, 11 / 1, 11 / 2, 11 / 3, 12 / 1, 13 / 1}

set_option maxHeartbeats 3000000 in
theorem tau3CertificateRatios_quotient_free :
    (tau3CertificateRatios : Set ℚ).Pairwise fun first second =>
      first / second ∉ tau3CertificateRatios ∧
        second / first ∉ tau3CertificateRatios := by
  norm_num [tau3CertificateRatios, Set.Pairwise]

theorem tau3FiniteCoverGraph_cliqueFreeOn_two :
    (finiteCoverGraph tau3CertificateRatios).CliqueFreeOn
      (tau3CertificateRatios : Set ℚ) 2 :=
  finiteCoverGraph_cliqueFreeOn_two_of_quotient_free
    tau3CertificateRatios tau3CertificateRatios_quotient_free

#print axioms finiteAllowedRatios_cover
#print axioms tier_envelope_caps
#print axioms cliqueFreeOn_of_powerset_check
#print axioms tau3FiniteCoverGraph_cliqueFreeOn_two

end LRCPairRatioFiniteCover
end LonelyRunner
