import TournamentH7.LRCAnchoredCliqueTransfer
import TournamentH7.LRCB5ContinuumFloor
import TournamentH7.LRCPairCovarianceKernel

/-!
Concrete anchored quotient for the seven THM-954 pair-ratio layers.  An
anchor in a runner clique sends every other runner to its absolute-speed ratio
with the anchor.  Distinct absolute speeds make this map injective, and exact
field cancellation sends runner adjacency to quotient adjacency.

This isolates the remaining finite task: for each strict layer threshold,
certify clique-freeness of `allowedRatioGraph` from the primitive-ratio DAG.
-/

namespace LonelyRunner
namespace LRCPairRatioQuotient

open Finset SimpleGraph
open LRCAnchoredCliqueTransfer
open LRCB5ContinuumFloor
open LRCPairCovarianceKernel

/-- A rational ratio is allowed at a strict threshold when some nonzero speed
pair realizing it has pair deficit strictly above the threshold. -/
def ratioAllowed (threshold ratio : ℚ) : Prop :=
  ∃ first second : ℤ,
    first ≠ 0 ∧ second ≠ 0 ∧
    ratio = ((first.natAbs : ℚ) / (second.natAbs : ℚ)) ∧
    threshold < pairDeficit first second

/-- The infinite presentation of the finite strict primitive-ratio graph.
Disallowed rationals are isolated. -/
def allowedRatioGraph (threshold : ℚ) : SimpleGraph ℚ where
  Adj first second :=
    first ≠ second ∧ ratioAllowed threshold first ∧
      ratioAllowed threshold second ∧
      (ratioAllowed threshold (first / second) ∨
        ratioAllowed threshold (second / first))
  symm first second hadj := by
    rcases hadj with ⟨hne, hfirst, hsecond, hquotient | hquotient⟩
    · exact ⟨hne.symm, hsecond, hfirst, Or.inr hquotient⟩
    · exact ⟨hne.symm, hsecond, hfirst, Or.inl hquotient⟩
  loopless := ⟨fun _ hadj => hadj.1 rfl⟩

def anchoredRatio (v : Fin 13 → ℤ) (anchor index : Fin 13) : ℚ :=
  (v index).natAbs / (v anchor).natAbs

theorem anchoredRatio_allowed
    (v : Fin 13 → ℤ) (weight : LRCB5ContinuumFloor.PairSupport → ℚ)
    (threshold : ℚ)
    (hv : ∀ index, v index ≠ 0)
    (hweight : ∀ first second, first ≠ second →
      weight {first, second} = pairDeficit (v first) (v second))
    {clique : Finset (Fin 13)} {anchor index : Fin 13}
    (hclique : (thresholdGraph weight threshold).IsClique clique)
    (hanchor : anchor ∈ clique) (hindex : index ∈ clique.erase anchor) :
    ratioAllowed threshold (anchoredRatio v anchor index) := by
  have hindexClique : index ∈ clique := (Finset.mem_erase.mp hindex).2
  have hne : index ≠ anchor := (Finset.mem_erase.mp hindex).1
  have hadj := hclique hindexClique hanchor hne
  change index ≠ anchor ∧ threshold < weight {index, anchor} at hadj
  refine ⟨v index, v anchor, hv index, hv anchor, rfl, ?_⟩
  simpa [hweight index anchor hne] using hadj.2

theorem anchoredRatio_injectiveOn
    (v : Fin 13 → ℤ) (anchor : Fin 13)
    (hv : ∀ index, v index ≠ 0)
    (hdistinct : ∀ first second, first ≠ second →
      |v first| ≠ |v second|)
    (indices : Finset (Fin 13)) :
    Set.InjOn (anchoredRatio v anchor) indices := by
  intro first _ second _ hequal
  have hdenominator : ((v anchor).natAbs : ℚ) ≠ 0 := by
    exact_mod_cast (Int.natAbs_ne_zero.mpr (hv anchor))
  have habsCast : ((v first).natAbs : ℚ) = (v second).natAbs :=
    (div_left_inj' hdenominator).mp hequal
  have habsNat : (v first).natAbs = (v second).natAbs := by
    exact_mod_cast habsCast
  by_contra hne
  apply hdistinct first second hne
  have habsInt : ((v first).natAbs : ℤ) = (v second).natAbs := by
    exact_mod_cast habsNat
  simpa only [Int.natCast_natAbs] using habsInt

theorem anchoredRatio_quotient_allowed
    (v : Fin 13 → ℤ) (weight : LRCB5ContinuumFloor.PairSupport → ℚ)
    (threshold : ℚ)
    (hv : ∀ index, v index ≠ 0)
    (hweight : ∀ first second, first ≠ second →
      weight {first, second} = pairDeficit (v first) (v second))
    {clique : Finset (Fin 13)} {anchor first second : Fin 13}
    (hclique : (thresholdGraph weight threshold).IsClique clique)
    (hfirst : first ∈ clique.erase anchor)
    (hsecond : second ∈ clique.erase anchor) (hne : first ≠ second) :
    ratioAllowed threshold
      (anchoredRatio v anchor first / anchoredRatio v anchor second) := by
  have hfirstClique : first ∈ clique := (Finset.mem_erase.mp hfirst).2
  have hsecondClique : second ∈ clique := (Finset.mem_erase.mp hsecond).2
  have hadj := hclique hfirstClique hsecondClique hne
  change first ≠ second ∧ threshold < weight {first, second} at hadj
  refine ⟨v first, v second, hv first, hv second, ?_, ?_⟩
  · unfold anchoredRatio
    exact div_div_div_cancel_right₀
      (by exact_mod_cast (Int.natAbs_ne_zero.mpr (hv anchor))) _ _
  · simpa [hweight first second hne] using hadj.2

/-- Parameterized producer bridge for all seven strict thresholds.  A finite
certificate only has to prove clique-freeness of `allowedRatioGraph`; anchored
arithmetic then produces the runner threshold certificate. -/
theorem thresholdGraph_cliqueFree_of_allowedRatioGraph
    (v : Fin 13 → ℤ)
    (weight : LRCB5ContinuumFloor.PairSupport → ℚ)
    (threshold : ℚ) (cliqueSize : ℕ) (hcliqueSize : 0 < cliqueSize)
    (hv : ∀ index, v index ≠ 0)
    (hdistinct : ∀ first second, first ≠ second →
      |v first| ≠ |v second|)
    (hweight : ∀ first second, first ≠ second →
      weight {first, second} = pairDeficit (v first) (v second))
    (hratioFree : (allowedRatioGraph threshold).CliqueFree cliqueSize) :
    (thresholdGraph weight threshold).CliqueFree (cliqueSize + 1) := by
  apply cliqueFree_succ_of_anchored_quotient
    (thresholdGraph weight threshold) (allowedRatioGraph threshold)
    cliqueSize hcliqueSize hratioFree
  intro clique anchor hclique hanchor
  refine ⟨anchoredRatio v anchor,
    anchoredRatio_injectiveOn v anchor hv hdistinct (clique.erase anchor), ?_⟩
  intro first hfirst second hsecond hne
  refine ⟨?_,
    anchoredRatio_allowed v weight threshold hv hweight hclique.isClique
      hanchor hfirst,
    anchoredRatio_allowed v weight threshold hv hweight hclique.isClique
      hanchor hsecond,
    Or.inl (anchoredRatio_quotient_allowed v weight threshold hv hweight
      hclique.isClique hfirst hsecond hne)⟩
  intro hequal
  exact hne (anchoredRatio_injectiveOn v anchor hv hdistinct
    (clique.erase anchor) hfirst hsecond hequal)

/-- A finite list covering every allowed rational carries every nontrivial
clique of the infinite presentation. -/
theorem allowedRatioGraph_cliqueFree_of_finite_cover
    (threshold : ℚ) (vertices : Finset ℚ) (cliqueSize : ℕ)
    (hcliqueSize : 2 ≤ cliqueSize)
    (hcover : ∀ ratio, ratioAllowed threshold ratio → ratio ∈ vertices)
    (hfinite : (allowedRatioGraph threshold).CliqueFreeOn
      (vertices : Set ℚ) cliqueSize) :
    (allowedRatioGraph threshold).CliqueFree cliqueSize := by
  intro clique hclique
  apply hfinite ?_ hclique
  intro ratio hratio
  have hratioFinset : ratio ∈ clique := hratio
  have heraseCard : (clique.erase ratio).card = cliqueSize - 1 := by
    rw [Finset.card_erase_of_mem hratioFinset, hclique.card_eq]
  have herasePositive : 0 < (clique.erase ratio).card := by
    rw [heraseCard]
    omega
  obtain ⟨other, hother⟩ := Finset.card_pos.mp herasePositive
  have hotherClique : other ∈ clique := (Finset.mem_erase.mp hother).2
  have hne : ratio ≠ other := (Finset.mem_erase.mp hother).1.symm
  have hadj := hclique.isClique hratioFinset hotherClique hne
  exact hcover ratio hadj.2.1

#print axioms anchoredRatio_allowed
#print axioms anchoredRatio_injectiveOn
#print axioms anchoredRatio_quotient_allowed
#print axioms thresholdGraph_cliqueFree_of_allowedRatioGraph
#print axioms allowedRatioGraph_cliqueFree_of_finite_cover

end LRCPairRatioQuotient
end LonelyRunner
