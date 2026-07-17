import TournamentH7.LRCPairRatioQuotient
import TournamentH7.LRCPairTopClassification

/-!
# The first direct THM-954 ratio certificate

At threshold `tau3 = 2/441`, the Bernoulli envelope reduces the primitive
ratio search to `a*b < 54`.  Exact kernel enumeration leaves fourteen oriented
ratios.  No quotient of two distinct listed ratios is listed, so the anchored
ratio graph is `K₂`-free and the runner threshold graph is `K₃`-free.
-/

namespace LonelyRunner
namespace LRCPairRatioTau3

open Finset SimpleGraph
open LRCB5ContinuumFloor LRCWeightedRatioLayer
open LRCPairCovarianceKernel LRCPairTopClassification
open LRCPairRatioQuotient

def tau3Vertices : Finset ℚ :=
  {1 / 10, 1 / 11, 1 / 12, 1 / 13, 2 / 11, 3 / 10, 3 / 11,
    10 / 1, 10 / 3, 11 / 1, 11 / 2, 11 / 3, 12 / 1, 13 / 1}

theorem product_lt_fifty_four_of_tau3_lt
    (first second : ℕ) (hfirst : 0 < first) (hsecond : 0 < second)
    (habove : tau3 < primitiveDeficit first second) :
    first * second < 54 := by
  have hdenominator : (0 : ℚ) < (first : ℚ) * second := by positivity
  have hbound := habove.trans_le
    (primitiveDeficit_le_envelope first second hfirst hsecond)
  have hscaled := (lt_div_iff₀ hdenominator).mp hbound
  have hcast : ((first * second : ℕ) : ℚ) < 54 := by
    push_cast
    norm_num [tau3] at hscaled ⊢
    nlinarith
  exact_mod_cast hcast

set_option maxHeartbeats 8000000 in
theorem primitive_tau3_mem
    (first second : ℕ) (hfirst : 0 < first) (hsecond : 0 < second)
    (hcoprime : Nat.Coprime first second)
    (habove : tau3 < primitiveDeficit first second) :
    (first : ℚ) / second ∈ tau3Vertices := by
  have hproduct := product_lt_fifty_four_of_tau3_lt first second
    hfirst hsecond habove
  have hfirstLt : first < 54 := by
    calc
      first = first * 1 := by omega
      _ ≤ first * second := Nat.mul_le_mul_left first hsecond
      _ < 54 := hproduct
  have hsecondLt : second < 54 := by
    calc
      second = 1 * second := by omega
      _ ≤ first * second := Nat.mul_le_mul_right second hfirst
      _ < 54 := hproduct
  interval_cases first <;> interval_cases second <;>
    norm_num [primitiveDeficit, primitiveKernel, bernoulliResidue14,
      tau3Vertices, tau3] at *

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

theorem tau3_ratio_cover
    (first second : ℤ) (hfirst : first ≠ 0) (hsecond : second ≠ 0)
    (habove : tau3 < pairDeficit first second) :
    ((first.natAbs : ℚ) / second.natAbs) ∈ tau3Vertices := by
  rw [pairDeficit_eq_primitiveDeficit] at habove
  rw [speed_ratio_eq_reduced_ratio first second hfirst]
  exact primitive_tau3_mem _ _ (reducedFirst_pos first second hfirst)
    (reducedSecond_pos first second hsecond)
    (reduced_coprime first second hfirst) habove

theorem ratioAllowed_tau3_mem (ratio : ℚ)
    (hallowed : ratioAllowed tau3 ratio) :
    ratio ∈ tau3Vertices := by
  obtain ⟨first, second, hfirst, hsecond, rfl, habove⟩ := hallowed
  exact tau3_ratio_cover first second hfirst hsecond habove

set_option maxHeartbeats 3000000 in
theorem tau3Vertices_quotient_free :
    (tau3Vertices : Set ℚ).Pairwise fun first second =>
      first / second ∉ tau3Vertices ∧ second / first ∉ tau3Vertices := by
  norm_num [tau3Vertices, Set.Pairwise]

theorem allowedRatioGraph_tau3_cliqueFreeOn_two :
    (allowedRatioGraph tau3).CliqueFreeOn (tau3Vertices : Set ℚ) 2 := by
  intro clique hsubset hclique
  obtain ⟨first, second, hne, hpair⟩ :=
    Finset.card_eq_two.mp hclique.card_eq
  subst clique
  have hfirst : first ∈ tau3Vertices := hsubset (by simp)
  have hsecond : second ∈ tau3Vertices := hsubset (by simp)
  have hadj : (allowedRatioGraph tau3).Adj first second :=
    hclique.isClique (by simp) (by simp) hne
  rcases hadj.2.2.2 with hquotient | hquotient
  · exact (tau3Vertices_quotient_free hfirst hsecond hne).1
      (ratioAllowed_tau3_mem _ hquotient)
  · exact (tau3Vertices_quotient_free hfirst hsecond hne).2
      (ratioAllowed_tau3_mem _ hquotient)

theorem allowedRatioGraph_tau3_cliqueFree_two :
    (allowedRatioGraph tau3).CliqueFree 2 :=
  allowedRatioGraph_cliqueFree_of_finite_cover tau3 tau3Vertices 2
    (by omega) ratioAllowed_tau3_mem
    allowedRatioGraph_tau3_cliqueFreeOn_two

/-- The actual `free_tau3` producer required by the THM-954 continuum
certificate. -/
theorem thresholdGraph_tau3_cliqueFree_three
    (v : Fin 13 → ℤ)
    (weight : LRCB5ContinuumFloor.PairSupport → ℚ)
    (hv : ∀ index, v index ≠ 0)
    (hdistinct : ∀ first second, first ≠ second →
      |v first| ≠ |v second|)
    (hweight : ∀ first second, first ≠ second →
      weight {first, second} = pairDeficit (v first) (v second)) :
    (thresholdGraph weight tau3).CliqueFree 3 := by
  simpa using thresholdGraph_cliqueFree_of_allowedRatioGraph
    v weight tau3 2 (by omega) hv hdistinct hweight
      allowedRatioGraph_tau3_cliqueFree_two

#print axioms primitive_tau3_mem
#print axioms tau3_ratio_cover
#print axioms tau3Vertices_quotient_free
#print axioms allowedRatioGraph_tau3_cliqueFree_two
#print axioms thresholdGraph_tau3_cliqueFree_three

end LRCPairRatioTau3
end LonelyRunner
