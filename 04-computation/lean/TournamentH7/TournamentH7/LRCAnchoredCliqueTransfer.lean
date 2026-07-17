import Mathlib.Combinatorics.SimpleGraph.Clique

/-!
Generic anchored-clique transfer used by the THM-954 ratio graphs.  Erasing an
anchor from a runner clique and injecting the remaining vertices into a ratio
clique lifts a `k`-clique exclusion on the ratio graph to a `k+1` exclusion on
the runner graph.
-/

namespace LonelyRunner
namespace LRCAnchoredCliqueTransfer

open Finset SimpleGraph

/-- If every clique can be anchored and its remaining vertices inject into a
clique of a quotient graph, one quotient clique exclusion lifts by one. -/
theorem cliqueFree_succ_of_anchored_quotient
    {Runner Ratio : Type*} [DecidableEq Runner] [DecidableEq Ratio]
    (runnerGraph : SimpleGraph Runner) (ratioGraph : SimpleGraph Ratio)
    (cliqueSize : ℕ) (hpositive : 0 < cliqueSize)
    (hratioFree : ratioGraph.CliqueFree cliqueSize)
    (hquotient : ∀ (clique : Finset Runner) (anchor : Runner),
      runnerGraph.IsNClique (cliqueSize + 1) clique → anchor ∈ clique →
      ∃ quotient : Runner → Ratio,
        Set.InjOn quotient (clique.erase anchor) ∧
        ∀ first ∈ clique.erase anchor, ∀ second ∈ clique.erase anchor,
          first ≠ second → ratioGraph.Adj (quotient first) (quotient second)) :
    runnerGraph.CliqueFree (cliqueSize + 1) := by
  intro clique hclique
  have hcardPositive : 0 < clique.card := by
    rw [hclique.card_eq]
    omega
  obtain ⟨anchor, hanchor⟩ := Finset.card_pos.mp hcardPositive
  obtain ⟨quotient, hinjective, hadjacent⟩ :=
    hquotient clique anchor hclique hanchor
  let quotientClique : Finset Ratio :=
    (clique.erase anchor).image quotient
  apply hratioFree quotientClique
  constructor
  · intro first hfirst second hsecond hne
    have hfirstFinset : first ∈ quotientClique := hfirst
    have hsecondFinset : second ∈ quotientClique := hsecond
    simp only [quotientClique, Finset.mem_image] at hfirstFinset hsecondFinset
    obtain ⟨firstSource, hfirstSource, rfl⟩ := hfirstFinset
    obtain ⟨secondSource, hsecondSource, rfl⟩ := hsecondFinset
    apply hadjacent firstSource hfirstSource secondSource hsecondSource
    intro hequal
    subst secondSource
    exact hne rfl
  · simp only [quotientClique]
    rw [Finset.card_image_iff.mpr hinjective,
      Finset.card_erase_of_mem hanchor, hclique.card_eq]
    omega

#print axioms cliqueFree_succ_of_anchored_quotient

end LRCAnchoredCliqueTransfer
end LonelyRunner
