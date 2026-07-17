import TournamentH7.LRCB5PairGridBridge
import TournamentH7.LRCWeightedRatioLayer
import Mathlib.Combinatorics.SimpleGraph.Extremal.Turan

/-!
Formal adapter from termwise pair-covariance minorants and strict threshold
clique certificates to the conservative continuum pair floor of THM-954.
Mathlib's Turán theorem converts the seven clique-free obligations into the
exact edge-count caps consumed by the weighted ratio-layer arithmetic.
-/

namespace LonelyRunner
namespace LRCB5ContinuumFloor

open Finset
open SimpleGraph
open LRCB5CleanModulus
open LRCB5NormalizedBridge
open LRCB5PairGridBridge
open LRCPairRatioLayerArithmetic
open LRCWeightedRatioLayer

noncomputable section

abbrev PairSupport := Finset (Fin 13)

def pairSupports : Finset PairSupport :=
  (Finset.univ : Finset (Fin 13)).powersetCard 2

/-- The strict threshold graph on the thirteen runner indices. -/
def thresholdGraph (weight : PairSupport → ℚ) (threshold : ℚ) :
    SimpleGraph (Fin 13) where
  Adj first second := first ≠ second ∧ threshold < weight {first, second}
  symm first second hadj :=
    ⟨hadj.1.symm, by simpa [Finset.pair_comm] using hadj.2⟩
  loopless := ⟨fun _vertex hadj => hadj.1 rfl⟩

instance thresholdGraph_decidableAdj
    (weight : PairSupport → ℚ) (threshold : ℚ) :
    DecidableRel (thresholdGraph weight threshold).Adj := by
  intro first second
  change Decidable (first ≠ second ∧ threshold < weight {first, second})
  infer_instance

/-- Every filtered two-support is represented by an edge of the corresponding
strict threshold graph.  Equality of the image is enough for cardinal bounds;
injectivity of `Sym2.toFinset` is deliberately unnecessary here. -/
theorem thresholdGraph_edge_image
    (weight : PairSupport → ℚ) (threshold : ℚ) :
    (thresholdGraph weight threshold).edgeFinset.image Sym2.toFinset =
      pairSupports.filter (fun support => threshold < weight support) := by
  classical
  ext support
  constructor
  · intro hsupport
    rw [Finset.mem_image] at hsupport
    obtain ⟨edge, hedge, rfl⟩ := hsupport
    rw [Finset.mem_filter]
    constructor
    · rw [pairSupports, Finset.mem_powersetCard]
      exact ⟨by simp, (thresholdGraph weight threshold).card_toFinset_mem_edgeFinset
        ⟨edge, hedge⟩⟩
    · induction edge using Sym2.inductionOn with
      | _ first second =>
          rw [SimpleGraph.mem_edgeFinset, SimpleGraph.mem_edgeSet] at hedge
          change first ≠ second ∧ threshold < weight {first, second} at hedge
          simpa [Sym2.toFinset_mk_eq] using hedge.2
  · intro hsupport
    rw [Finset.mem_filter, pairSupports, Finset.mem_powersetCard] at hsupport
    obtain ⟨_, hcard⟩ := hsupport.1
    obtain ⟨first, second, hne, rfl⟩ := Finset.card_eq_two.mp hcard
    rw [Finset.mem_image]
    refine ⟨s(first, second), ?_, Sym2.toFinset_mk_eq⟩
    rw [SimpleGraph.mem_edgeFinset, SimpleGraph.mem_edgeSet]
    change first ≠ second ∧ threshold < weight {first, second}
    exact ⟨hne, hsupport.2⟩

theorem countAbove_le_thresholdGraph_edges
    (weight : PairSupport → ℚ) (threshold : ℚ) :
    countAbove pairSupports weight threshold ≤
      (thresholdGraph weight threshold).edgeFinset.card := by
  classical
  unfold countAbove
  rw [← thresholdGraph_edge_image]
  exact Finset.card_image_le

/-- Minimal interface between a concrete pair-covariance calculation and the
abstract weighted ratio-layer arithmetic.  A producer may use the exact
Bernoulli formula, but the consumer only needs a rational termwise minorant
and the nine threshold-cardinality bounds. -/
structure ContinuumTierCertificate (v : Fin 13 → ℤ) where
  weight : PairSupport → ℚ
  correlation_lower : ∀ support ∈ pairSupports,
    -(weight support : ℝ) ≤
      pairContinuumCorrelation v support - 1 / 49
  weight_le_tau2 : ∀ support ∈ pairSupports, weight support ≤ tau2
  count_tau9 : countAbove pairSupports weight tau9 ≤ 73
  count_tau8 : countAbove pairSupports weight tau8 ≤ 72
  count_tau7 : countAbove pairSupports weight tau7 ≤ 70
  count_tau6 : countAbove pairSupports weight tau6 ≤ 67
  count_tau5 : countAbove pairSupports weight tau5 ≤ 63
  count_tau4 : countAbove pairSupports weight tau4 ≤ 56
  count_tau3 : countAbove pairSupports weight tau3 ≤ 42
  count_ratio11 : countAbove pairSupports weight ratio11Level ≤ 24
  count_ratio12 : countAbove pairSupports weight ratio12Level ≤ 12

/-- Producer-facing version: the seven finite ratio DAGs need only establish
clique-freeness of the strict threshold graphs.  Mathlib's exact Turán theorem
then supplies all seven numerical edge caps. -/
structure ContinuumCliqueCertificate (v : Fin 13 → ℤ) where
  weight : PairSupport → ℚ
  correlation_lower : ∀ support ∈ pairSupports,
    -(weight support : ℝ) ≤
      pairContinuumCorrelation v support - 1 / 49
  weight_le_tau2 : ∀ support ∈ pairSupports, weight support ≤ tau2
  free_tau9 : (thresholdGraph weight tau9).CliqueFree 9
  free_tau8 : (thresholdGraph weight tau8).CliqueFree 8
  free_tau7 : (thresholdGraph weight tau7).CliqueFree 7
  free_tau6 : (thresholdGraph weight tau6).CliqueFree 6
  free_tau5 : (thresholdGraph weight tau5).CliqueFree 5
  free_tau4 : (thresholdGraph weight tau4).CliqueFree 4
  free_tau3 : (thresholdGraph weight tau3).CliqueFree 3
  count_ratio11 : countAbove pairSupports weight ratio11Level ≤ 24
  count_ratio12 : countAbove pairSupports weight ratio12Level ≤ 12

def ContinuumCliqueCertificate.toTierCertificate
    {v : Fin 13 → ℤ} (certificate : ContinuumCliqueCertificate v) :
    ContinuumTierCertificate v where
  weight := certificate.weight
  correlation_lower := certificate.correlation_lower
  weight_le_tau2 := certificate.weight_le_tau2
  count_tau9 := (countAbove_le_thresholdGraph_edges certificate.weight tau9).trans
    (by simpa using (certificate.free_tau9.card_edgeFinset_le (r := 8)))
  count_tau8 := (countAbove_le_thresholdGraph_edges certificate.weight tau8).trans
    (by simpa using (certificate.free_tau8.card_edgeFinset_le (r := 7)))
  count_tau7 := (countAbove_le_thresholdGraph_edges certificate.weight tau7).trans
    (by simpa using (certificate.free_tau7.card_edgeFinset_le (r := 6)))
  count_tau6 := (countAbove_le_thresholdGraph_edges certificate.weight tau6).trans
    (by simpa using (certificate.free_tau6.card_edgeFinset_le (r := 5)))
  count_tau5 := (countAbove_le_thresholdGraph_edges certificate.weight tau5).trans
    (by simpa using (certificate.free_tau5.card_edgeFinset_le (r := 4)))
  count_tau4 := (countAbove_le_thresholdGraph_edges certificate.weight tau4).trans
    (by simpa using (certificate.free_tau4.card_edgeFinset_le (r := 3)))
  count_tau3 := (countAbove_le_thresholdGraph_edges certificate.weight tau3).trans
    (by simpa using (certificate.free_tau3.card_edgeFinset_le (r := 2)))
  count_ratio11 := certificate.count_ratio11
  count_ratio12 := certificate.count_ratio12

theorem pairSupports_card : pairSupports.card = 78 := by
  exact pair_support_card

/-- The exact formal join requested by `LRCB5PairGridBridge`: any rational
pair-weight certificate satisfying the existing layer premises supplies its
continuum floor. -/
theorem continuumMass2_ge_neg_path_bound
    (v : Fin 13 → ℤ) (certificate : ContinuumTierCertificate v) :
    -(negativePairTierBoundPathOnly : ℝ) ≤ continuumMass2 v := by
  have hweightQ :
      ∑ support ∈ pairSupports, certificate.weight support ≤
        negativePairTierBoundPathOnly :=
    weighted_pair_sum_le_path_bound pairSupports certificate.weight
      certificate.weight_le_tau2 (by rw [pairSupports_card])
      certificate.count_tau9 certificate.count_tau8 certificate.count_tau7
      certificate.count_tau6 certificate.count_tau5 certificate.count_tau4
      certificate.count_tau3 certificate.count_ratio11
      certificate.count_ratio12
  have hweightR :
      ∑ support ∈ pairSupports, (certificate.weight support : ℝ) ≤
        (negativePairTierBoundPathOnly : ℝ) := by
    exact_mod_cast hweightQ
  calc
    -(negativePairTierBoundPathOnly : ℝ) ≤
        -(∑ support ∈ pairSupports, (certificate.weight support : ℝ)) := by
      linarith
    _ = ∑ support ∈ pairSupports, -(certificate.weight support : ℝ) := by
      simp
    _ ≤ ∑ support ∈ pairSupports,
        (pairContinuumCorrelation v support - 1 / 49) := by
      exact Finset.sum_le_sum certificate.correlation_lower
    _ = continuumMass2 v := by
      rfl

theorem continuumMass2_ge_neg_path_bound_of_clique_certificate
    (v : Fin 13 → ℤ) (certificate : ContinuumCliqueCertificate v) :
    -(negativePairTierBoundPathOnly : ℝ) ≤ continuumMass2 v :=
  continuumMass2_ge_neg_path_bound v certificate.toTierCertificate

/-- Consumer-ready form for the strict clean-modulus pair socket. -/
theorem normalizedMass2_clean_534_gt_neg_target_of_tier_certificate
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (certificate : ContinuumTierCertificate v)
    (hraw : RawPairGridDiscrepancyAt v (cleanModulus v 534)) :
    -(13 / 30 : ℚ) < normalizedMass2 v (cleanModulus v 534) := by
  exact normalizedMass2_clean_534_gt_neg_target_of_continuum_floor
    v hv (continuumMass2_ge_neg_path_bound v certificate) hraw

theorem normalizedMass2_clean_534_gt_neg_target_of_clique_certificate
    (v : Fin 13 → ℤ) (hv : ∀ i, v i ≠ 0)
    (certificate : ContinuumCliqueCertificate v)
    (hraw : RawPairGridDiscrepancyAt v (cleanModulus v 534)) :
    -(13 / 30 : ℚ) < normalizedMass2 v (cleanModulus v 534) := by
  exact normalizedMass2_clean_534_gt_neg_target_of_continuum_floor
    v hv (continuumMass2_ge_neg_path_bound_of_clique_certificate v certificate) hraw

#print axioms continuumMass2_ge_neg_path_bound
#print axioms continuumMass2_ge_neg_path_bound_of_clique_certificate
#print axioms normalizedMass2_clean_534_gt_neg_target_of_tier_certificate
#print axioms normalizedMass2_clean_534_gt_neg_target_of_clique_certificate

end
end LRCB5ContinuumFloor
end LonelyRunner

