import TournamentH7.LRCPairRatioTau5Replay5

/-!
# The tau5 local-neighborhood replay

Seventeen kernel shards check that every reduced cross-product edge maps to the
explicit sparse table.  The table has independent common neighborhoods, hence
no four-clique.  This file lifts that finite result to the `tau5` allowed-ratio
graph and supplies the `free_tau5` endpoint for THM-954.
-/

namespace LonelyRunner
namespace LRCPairRatioTau5

open Finset SimpleGraph
open LRCB5ContinuumFloor LRCPairCovarianceKernel LRCWeightedRatioLayer
open LRCPairRatioQuotient LRCPairRatioFiniteCover

set_option maxRecDepth 100000

theorem edgeTableCheckBlock_true (block : Fin 17) :
    edgeTableCheckBlock block = true := by
  fin_cases block
  · exact edgeTableCheckBlock_zero_true
  · exact edgeTableCheckBlock_1_true
  · exact edgeTableCheckBlock_2_true
  · exact edgeTableCheckBlock_3_true
  · exact edgeTableCheckBlock_4_true
  · exact edgeTableCheckBlock_5_true
  · exact edgeTableCheckBlock_6_true
  · exact edgeTableCheckBlock_7_true
  · exact edgeTableCheckBlock_8_true
  · exact edgeTableCheckBlock_9_true
  · exact edgeTableCheckBlock_10_true
  · exact edgeTableCheckBlock_11_true
  · exact edgeTableCheckBlock_12_true
  · exact edgeTableCheckBlock_13_true
  · exact edgeTableCheckBlock_14_true
  · exact edgeTableCheckBlock_15_true
  · exact edgeTableCheckBlock_16_true

theorem localTableCheckBlock_true (block : Fin 17) :
    localTableCheckBlock block = true := by
  fin_cases block
  · exact localTableCheckBlock_zero_true
  · exact localTableCheckBlock_1_true
  · exact localTableCheckBlock_2_true
  · exact localTableCheckBlock_3_true
  · exact localTableCheckBlock_4_true
  · exact localTableCheckBlock_5_true
  · exact localTableCheckBlock_6_true
  · exact localTableCheckBlock_7_true
  · exact localTableCheckBlock_8_true
  · exact localTableCheckBlock_9_true
  · exact localTableCheckBlock_10_true
  · exact localTableCheckBlock_11_true
  · exact localTableCheckBlock_12_true
  · exact localTableCheckBlock_13_true
  · exact localTableCheckBlock_14_true
  · exact localTableCheckBlock_15_true
  · exact localTableCheckBlock_16_true

theorem mem_indexBlock (index : Fin 272) :
    index ∈ indexBlock (indexBlockIndex index) := by
  decide +kernel +revert

theorem indexGraph_le_tableGraph : indexGraph ≤ tableGraph := by
  intro first second hadj
  have hfirst := (List.all_eq_true.mp
      (edgeTableCheckBlock_true (indexBlockIndex first)))
    first (mem_indexBlock first)
  have hsecond := (List.all_eq_true.mp hfirst)
    second (List.mem_finRange second)
  rw [if_pos hadj] at hsecond
  exact of_decide_eq_true hsecond

/-- Every table-graph edge has an independent common neighborhood. -/
theorem tableGraph_locally_triangleFree :
    ∀ first second, tableGraph.Adj first second →
      ∀ third, tableGraph.Adj first third →
        tableGraph.Adj second third →
          ∀ fourth, tableGraph.Adj first fourth →
            tableGraph.Adj second fourth →
              ¬tableGraph.Adj third fourth := by
  intro first second hfirstSecond third hfirstThird hsecondThird
    fourth hfirstFourth hsecondFourth hthirdFourth
  have hfirst := (List.all_eq_true.mp
      (localTableCheckBlock_true (indexBlockIndex first)))
    first (mem_indexBlock first)
  have hsecond := (List.all_eq_true.mp hfirst)
    second (List.mem_finRange second)
  rw [if_pos hfirstSecond] at hsecond
  have hthird := (List.all_eq_true.mp hsecond)
    third (List.mem_finRange third)
  rw [if_pos hfirstThird, if_pos hsecondThird] at hthird
  have hfourth := (List.all_eq_true.mp hthird)
    fourth (List.mem_finRange fourth)
  rw [if_pos hfirstFourth, if_pos hsecondFourth] at hfourth
  exact (of_decide_eq_true hfourth) hthirdFourth

theorem indexGraph_no_four
    (first second third fourth : Fin 272)
    (hfirstSecond : indexGraph.Adj first second)
    (hfirstThird : indexGraph.Adj first third)
    (hfirstFourth : indexGraph.Adj first fourth)
    (hsecondThird : indexGraph.Adj second third)
    (hsecondFourth : indexGraph.Adj second fourth)
    (hthirdFourth : indexGraph.Adj third fourth) : False := by
  exact (tableGraph_locally_triangleFree first second
    (indexGraph_le_tableGraph hfirstSecond) third
    (indexGraph_le_tableGraph hfirstThird)
    (indexGraph_le_tableGraph hsecondThird) fourth
    (indexGraph_le_tableGraph hfirstFourth)
    (indexGraph_le_tableGraph hsecondFourth))
    (indexGraph_le_tableGraph hthirdFourth)

abbrev graph : SimpleGraph ℚ := finiteCoverGraph tau5Vertices

theorem index_edge_maps
    (first second : Fin 272)
    (hadj : graph.Adj (ratio first) (ratio second)) :
    indexGraph.Adj first second := by
  have hne : first ≠ second := by
    intro hequal
    subst second
    exact hadj.1 rfl
  refine ⟨hne, ?_⟩
  rcases hadj.2 with hquotient | hquotient
  · apply Or.inl
    apply reducedCross_above_of_tau5_mem
    rw [tau5FiniteRatios_eq_tau5Vertices]
    exact hquotient
  · apply Or.inr
    apply reducedCross_above_of_tau5_mem
    rw [tau5FiniteRatios_eq_tau5Vertices]
    exact hquotient

theorem finiteCoverGraph_tau5_cliqueFreeOn_four :
    graph.CliqueFreeOn (tau5Vertices : Set ℚ) 4 := by
  intro clique hsubset hclique
  obtain ⟨first, second, third, fourth, hfirstSecond, hfirstThird,
      hfirstFourth, hsecondThird, hsecondFourth, hthirdFourth,
      hcliqueEq⟩ := Finset.card_eq_four.mp hclique.card_eq
  subst clique
  have hfirst : first ∈ tau5Vertices := hsubset (by simp)
  have hsecond : second ∈ tau5Vertices := hsubset (by simp)
  have hthird : third ∈ tau5Vertices := hsubset (by simp)
  have hfourth : fourth ∈ tau5Vertices := hsubset (by simp)
  rw [tau5Vertices, Finset.mem_image] at hfirst hsecond hthird hfourth
  obtain ⟨firstIndex, _, rfl⟩ := hfirst
  obtain ⟨secondIndex, _, rfl⟩ := hsecond
  obtain ⟨thirdIndex, _, rfl⟩ := hthird
  obtain ⟨fourthIndex, _, rfl⟩ := hfourth
  apply indexGraph_no_four firstIndex secondIndex thirdIndex fourthIndex
  · apply index_edge_maps
    exact hclique.isClique (by simp) (by simp) hfirstSecond
  · apply index_edge_maps
    exact hclique.isClique (by simp) (by simp) hfirstThird
  · apply index_edge_maps
    exact hclique.isClique (by simp) (by simp) hfirstFourth
  · apply index_edge_maps
    exact hclique.isClique (by simp) (by simp) hsecondThird
  · apply index_edge_maps
    exact hclique.isClique (by simp) (by simp) hsecondFourth
  · apply index_edge_maps
    exact hclique.isClique (by simp) (by simp) hthirdFourth

theorem ratioAllowed_tau5_mem (candidate : ℚ)
    (hallowed : ratioAllowed tau5 candidate) :
    candidate ∈ tau5Vertices := by
  rw [← tau5FiniteRatios_eq_tau5Vertices]
  exact tau5FiniteRatios_cover candidate hallowed

theorem allowedRatioGraph_tau5_le_graph :
    allowedRatioGraph tau5 ≤ graph :=
  allowedRatioGraph_le_finiteCoverGraph tau5 tau5Vertices
    ratioAllowed_tau5_mem

theorem allowedRatioGraph_tau5_cliqueFreeOn_four :
    (allowedRatioGraph tau5).CliqueFreeOn (tau5Vertices : Set ℚ) 4 :=
  SimpleGraph.CliqueFreeOn.anti
    (allowedRatioGraph tau5) graph
    allowedRatioGraph_tau5_le_graph
    finiteCoverGraph_tau5_cliqueFreeOn_four

theorem allowedRatioGraph_tau5_cliqueFree_four :
    (allowedRatioGraph tau5).CliqueFree 4 :=
  allowedRatioGraph_cliqueFree_of_finite_cover tau5 tau5Vertices 4
    (by omega) ratioAllowed_tau5_mem
    allowedRatioGraph_tau5_cliqueFreeOn_four

/-- The `free_tau5` producer required by the THM-954 continuum certificate. -/
theorem thresholdGraph_tau5_cliqueFree_five
    (v : Fin 13 → ℤ)
    (weight : LRCB5ContinuumFloor.PairSupport → ℚ)
    (hv : ∀ index, v index ≠ 0)
    (hdistinct : ∀ first second, first ≠ second →
      |v first| ≠ |v second|)
    (hweight : ∀ first second, first ≠ second →
      weight {first, second} = pairDeficit (v first) (v second)) :
    (thresholdGraph weight tau5).CliqueFree 5 := by
  simpa using thresholdGraph_cliqueFree_of_allowedRatioGraph
    v weight tau5 4 (by omega) hv hdistinct hweight
      allowedRatioGraph_tau5_cliqueFree_four

#print axioms primitiveDeficit_above_of_mem_finiteAllowedRatios
#print axioms ratio_div_eq_reducedCross
#print axioms reducedCross_above_of_tau5_mem
#print axioms edgeTableCheckBlock_zero_true
#print axioms edgeTableCheckBlock_true
#print axioms indexGraph_le_tableGraph
#print axioms tableGraph_locally_triangleFree
#print axioms quotient_graph_fingerprint
#print axioms tau5FiniteRatios_eq_tau5Vertices
#print axioms finiteCoverGraph_tau5_cliqueFreeOn_four
#print axioms allowedRatioGraph_tau5_cliqueFree_four
#print axioms thresholdGraph_tau5_cliqueFree_five

end LRCPairRatioTau5
end LonelyRunner
