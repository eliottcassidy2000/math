import TournamentH7.LRCPairRatioFiniteCover

/-!
# The tau4 parallel-class circle certificate

The exact primitive-ratio cover at `tau4` has sixty vertices.  Its quotient
graph maps to the Wagner graph: eight parallel classes arranged on a circle,
with allowed class displacements `±1` and the antipodal displacement `4`.
Since the Wagner graph is triangle-free, the ratio graph is `K₃`-free and the
runner threshold graph is `K₄`-free.

The finite replay uses natural cross-products rather than rational reduction.
For memory stability its Boolean table is split into a `6 × 6` grid of ten by
ten index blocks.  Every cell is checked by kernel reduction; this file uses no
`native_decide` or native SAT oracle.

The quotient is not `K₂,₂`-free: the sixty-vertex graph contains eighteen
four-cycles.  Even the parity-cross core (class masses `18 × 20`) has forty-four
edges and four four-cycles.  Thus a raw Zarankiewicz cap would destroy faithful
information.  The correct certificate is Wagner `K₃`-freeness: its `C₈`
bipartite spine has eight edges (below `z(4,4;2,2) = 9`), while four antipodal
same-parity matching edges supply the necessary non-bipartite correction.
-/

namespace LonelyRunner
namespace LRCPairRatioTau4

open Finset SimpleGraph
open LRCB5ContinuumFloor LRCPairCovarianceKernel LRCWeightedRatioLayer
open LRCPairRatioQuotient LRCPairRatioFiniteCover

def numerator : Fin 60 → ℕ :=
  ![1, 1, 1, 1, 1, 1, 1, 1, 2, 1,
    2, 1, 1, 1, 1, 3, 3, 1, 2, 2,
    2, 3, 3, 3, 4, 3, 4, 5, 5, 5,
    8, 9, 11, 9, 8, 11, 10, 11, 13, 9,
    11, 13, 8, 25, 26, 9, 10, 11, 12, 25,
    13, 39, 24, 25, 26, 27, 39, 40, 41, 55]

def denominator : Fin 60 → ℕ :=
  ![55, 41, 40, 39, 27, 26, 25, 24, 39, 13,
    25, 12, 11, 10, 9, 26, 25, 8, 13, 11,
    9, 13, 11, 10, 11, 8, 9, 11, 9, 8,
    5, 5, 5, 4, 3, 4, 3, 3, 3, 2,
    2, 2, 1, 3, 3, 1, 1, 1, 1, 2,
    1, 2, 1, 1, 1, 1, 1, 1, 1, 1]

def ratio (index : Fin 60) : ℚ :=
  (numerator index : ℚ) / denominator index

def tau4Vertices : Finset ℚ := Finset.univ.image ratio

/-- Parallel-class labels around the eight-class circle. -/
def color : Fin 60 → Fin 8 :=
  ![0, 0, 5, 1, 2, 0, 1, 1, 4, 0,
    0, 6, 0, 7, 5, 0, 0, 0, 0, 7,
    5, 5, 0, 1, 0, 2, 3, 1, 0, 4,
    0, 0, 0, 7, 6, 0, 2, 0, 1, 1,
    6, 0, 0, 0, 4, 1, 3, 0, 2, 0,
    0, 0, 2, 0, 0, 3, 0, 1, 0, 7]

def isIsolated : Fin 60 → Bool :=
  ![false, true, false, false, false, true, false, false, false, true,
    true, false, true, false, false, false, true, true, true, false,
    false, false, true, false, true, false, false, false, true, false,
    false, true, false, false, false, true, false, true, false, false,
    false, true, true, true, false, false, false, true, false, true,
    true, false, false, false, true, false, false, false, true, false]

def isolatedIndices : Finset (Fin 60) :=
  Finset.univ.filter fun index => isIsolated index

def activeFiber (fiber : Fin 8) : Finset (Fin 60) :=
  Finset.univ.filter fun index => !isIsolated index && color index == fiber

def wagnerAdj (first second : Fin 8) : Bool :=
  ((first.val + 1) % 8 == second.val) ||
  ((second.val + 1) % 8 == first.val) ||
  ((first.val + 4) % 8 == second.val)

def wagnerGraph : SimpleGraph (Fin 8) where
  Adj first second := wagnerAdj first second = true
  symm := by
    intro first second hadj
    fin_cases first <;> fin_cases second <;> simp_all [wagnerAdj]
  loopless := ⟨by
    intro vertex hadj
    fin_cases vertex <;> simp_all [wagnerAdj]⟩

instance : DecidableRel wagnerGraph.Adj := by
  intro first second
  change Decidable (wagnerAdj first second = true)
  infer_instance

abbrev graph : SimpleGraph ℚ := finiteCoverGraph tau4Vertices

theorem numerator_pos (index : Fin 60) : 0 < numerator index := by
  fin_cases index <;> decide

theorem denominator_pos (index : Fin 60) : 0 < denominator index := by
  fin_cases index <;> decide

theorem quotient_eq_ratio_iff_cross (first second quotient : Fin 60) :
    ratio first / ratio second = ratio quotient ↔
      numerator first * denominator second * denominator quotient =
        denominator first * numerator second * numerator quotient := by
  have hnfirst : (numerator first : ℚ) ≠ 0 := by
    exact_mod_cast (numerator_pos first).ne'
  have hnsecond : (numerator second : ℚ) ≠ 0 := by
    exact_mod_cast (numerator_pos second).ne'
  have hdfirst : (denominator first : ℚ) ≠ 0 := by
    exact_mod_cast (denominator_pos first).ne'
  have hdsecond : (denominator second : ℚ) ≠ 0 := by
    exact_mod_cast (denominator_pos second).ne'
  have hdquotient : (denominator quotient : ℚ) ≠ 0 := by
    exact_mod_cast (denominator_pos quotient).ne'
  unfold ratio
  constructor <;> intro h
  · field_simp [hnfirst, hnsecond, hdfirst, hdsecond, hdquotient] at h
    exact_mod_cast h
  · field_simp [hnfirst, hnsecond, hdfirst, hdsecond, hdquotient]
    exact_mod_cast h

def crossEntry (first second quotient : Fin 60) : Bool :=
  if numerator first * denominator second * denominator quotient =
      denominator first * numerator second * numerator quotient then
    wagnerAdj (color first) (color second)
  else true

def indexBlock : Fin 6 → List (Fin 60) :=
  ![[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
    [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
    [30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
    [40, 41, 42, 43, 44, 45, 46, 47, 48, 49],
    [50, 51, 52, 53, 54, 55, 56, 57, 58, 59]]

def indexBlockIndex (index : Fin 60) : Fin 6 :=
  ⟨index.val / 10, by omega⟩

def crossCheckBlock (firstBlock secondBlock : Fin 6) : Bool :=
  (indexBlock firstBlock).all fun first =>
    (indexBlock secondBlock).all fun second =>
      (List.finRange 60).all fun quotient =>
        crossEntry first second quotient

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_zero_zero_true :
    crossCheckBlock 0 0 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_zero_one_true :
    crossCheckBlock 0 1 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_zero_two_true :
    crossCheckBlock 0 2 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_zero_three_true :
    crossCheckBlock 0 3 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_zero_four_true :
    crossCheckBlock 0 4 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_zero_five_true :
    crossCheckBlock 0 5 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_one_zero_true :
    crossCheckBlock 1 0 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_one_one_true :
    crossCheckBlock 1 1 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_one_two_true :
    crossCheckBlock 1 2 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_one_three_true :
    crossCheckBlock 1 3 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_one_four_true :
    crossCheckBlock 1 4 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_one_five_true :
    crossCheckBlock 1 5 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_two_zero_true :
    crossCheckBlock 2 0 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_two_one_true :
    crossCheckBlock 2 1 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_two_two_true :
    crossCheckBlock 2 2 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_two_three_true :
    crossCheckBlock 2 3 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_two_four_true :
    crossCheckBlock 2 4 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_two_five_true :
    crossCheckBlock 2 5 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_three_zero_true :
    crossCheckBlock 3 0 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_three_one_true :
    crossCheckBlock 3 1 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_three_two_true :
    crossCheckBlock 3 2 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_three_three_true :
    crossCheckBlock 3 3 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_three_four_true :
    crossCheckBlock 3 4 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_three_five_true :
    crossCheckBlock 3 5 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_four_zero_true :
    crossCheckBlock 4 0 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_four_one_true :
    crossCheckBlock 4 1 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_four_two_true :
    crossCheckBlock 4 2 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_four_three_true :
    crossCheckBlock 4 3 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_four_four_true :
    crossCheckBlock 4 4 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_four_five_true :
    crossCheckBlock 4 5 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_five_zero_true :
    crossCheckBlock 5 0 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_five_one_true :
    crossCheckBlock 5 1 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_five_two_true :
    crossCheckBlock 5 2 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_five_three_true :
    crossCheckBlock 5 3 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_five_four_true :
    crossCheckBlock 5 4 = true := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 10000000 in
theorem crossCheckBlock_five_five_true :
    crossCheckBlock 5 5 = true := by
  decide +kernel

theorem crossCheckBlock_true (firstBlock secondBlock : Fin 6) :
    crossCheckBlock firstBlock secondBlock = true := by
  fin_cases firstBlock <;> fin_cases secondBlock
  · exact crossCheckBlock_zero_zero_true
  · exact crossCheckBlock_zero_one_true
  · exact crossCheckBlock_zero_two_true
  · exact crossCheckBlock_zero_three_true
  · exact crossCheckBlock_zero_four_true
  · exact crossCheckBlock_zero_five_true
  · exact crossCheckBlock_one_zero_true
  · exact crossCheckBlock_one_one_true
  · exact crossCheckBlock_one_two_true
  · exact crossCheckBlock_one_three_true
  · exact crossCheckBlock_one_four_true
  · exact crossCheckBlock_one_five_true
  · exact crossCheckBlock_two_zero_true
  · exact crossCheckBlock_two_one_true
  · exact crossCheckBlock_two_two_true
  · exact crossCheckBlock_two_three_true
  · exact crossCheckBlock_two_four_true
  · exact crossCheckBlock_two_five_true
  · exact crossCheckBlock_three_zero_true
  · exact crossCheckBlock_three_one_true
  · exact crossCheckBlock_three_two_true
  · exact crossCheckBlock_three_three_true
  · exact crossCheckBlock_three_four_true
  · exact crossCheckBlock_three_five_true
  · exact crossCheckBlock_four_zero_true
  · exact crossCheckBlock_four_one_true
  · exact crossCheckBlock_four_two_true
  · exact crossCheckBlock_four_three_true
  · exact crossCheckBlock_four_four_true
  · exact crossCheckBlock_four_five_true
  · exact crossCheckBlock_five_zero_true
  · exact crossCheckBlock_five_one_true
  · exact crossCheckBlock_five_two_true
  · exact crossCheckBlock_five_three_true
  · exact crossCheckBlock_five_four_true
  · exact crossCheckBlock_five_five_true

theorem mem_indexBlock (index : Fin 60) :
    index ∈ indexBlock (indexBlockIndex index) := by
  fin_cases index <;> decide

theorem cross_edge_maps
    (first second quotient : Fin 60)
    (hcross : numerator first * denominator second * denominator quotient =
      denominator first * numerator second * numerator quotient) :
    wagnerGraph.Adj (color first) (color second) := by
  have hfirst := (List.all_eq_true.mp
      (crossCheckBlock_true (indexBlockIndex first) (indexBlockIndex second)))
    first (mem_indexBlock first)
  have hsecond := (List.all_eq_true.mp hfirst)
    second (mem_indexBlock second)
  have hquotient := (List.all_eq_true.mp hsecond)
    quotient (List.mem_finRange quotient)
  simpa [crossEntry, hcross, wagnerGraph] using hquotient

theorem index_edge_maps
    (first second : Fin 60)
    (hadj : graph.Adj (ratio first) (ratio second)) :
    wagnerGraph.Adj (color first) (color second) := by
  rcases hadj.2 with hquotient | hquotient
  · rw [tau4Vertices, Finset.mem_image] at hquotient
    obtain ⟨quotient, _, hequal⟩ := hquotient
    apply cross_edge_maps first second quotient
    exact (quotient_eq_ratio_iff_cross first second quotient).mp hequal.symm
  · rw [tau4Vertices, Finset.mem_image] at hquotient
    obtain ⟨quotient, _, hequal⟩ := hquotient
    apply wagnerGraph.symm
    apply cross_edge_maps second first quotient
    exact (quotient_eq_ratio_iff_cross second first quotient).mp hequal.symm

set_option maxRecDepth 100000 in
theorem wagner_no_triangle :
    ∀ first second third : Fin 8,
      wagnerGraph.Adj first second →
      wagnerGraph.Adj first third →
      wagnerGraph.Adj second third → False := by
  decide

theorem finiteCoverGraph_tau4_cliqueFreeOn_three :
    graph.CliqueFreeOn (tau4Vertices : Set ℚ) 3 := by
  intro clique hsubset hclique
  obtain ⟨first, second, third, hfirstSecond, hfirstThird,
      hsecondThird, hcliqueEq⟩ := Finset.card_eq_three.mp hclique.card_eq
  subst clique
  have hfirst : first ∈ tau4Vertices := hsubset (by simp)
  have hsecond : second ∈ tau4Vertices := hsubset (by simp)
  have hthird : third ∈ tau4Vertices := hsubset (by simp)
  rw [tau4Vertices, Finset.mem_image] at hfirst hsecond hthird
  obtain ⟨firstIndex, _, rfl⟩ := hfirst
  obtain ⟨secondIndex, _, rfl⟩ := hsecond
  obtain ⟨thirdIndex, _, rfl⟩ := hthird
  apply wagner_no_triangle (color firstIndex) (color secondIndex)
    (color thirdIndex)
  · apply index_edge_maps
    exact hclique.isClique (by simp) (by simp) hfirstSecond
  · apply index_edge_maps
    exact hclique.isClique (by simp) (by simp) hfirstThird
  · apply index_edge_maps
    exact hclique.isClique (by simp) (by simp) hsecondThird

set_option maxRecDepth 100000 in
theorem parallel_class_fingerprint :
    isolatedIndices.card = 22 ∧
      (activeFiber 0).card = 7 ∧
      (activeFiber 1).card = 9 ∧
      (activeFiber 2).card = 5 ∧
      (activeFiber 3).card = 3 ∧
      (activeFiber 4).card = 3 ∧
      (activeFiber 5).card = 4 ∧
      (activeFiber 6).card = 3 ∧
      (activeFiber 7).card = 4 := by
  decide

set_option maxRecDepth 100000 in
set_option maxHeartbeats 20000000 in
theorem tau4FiniteRatios_eq_tau4Vertices :
    tau4FiniteRatios = tau4Vertices := by
  decide +kernel

set_option maxRecDepth 100000 in
theorem ratioAllowed_tau4_mem (candidate : ℚ)
    (hallowed : ratioAllowed tau4 candidate) :
    candidate ∈ tau4Vertices := by
  rw [← tau4FiniteRatios_eq_tau4Vertices]
  exact tau4FiniteRatios_cover candidate hallowed

theorem allowedRatioGraph_tau4_le_graph :
    allowedRatioGraph tau4 ≤ graph :=
  allowedRatioGraph_le_finiteCoverGraph tau4 tau4Vertices
    ratioAllowed_tau4_mem

theorem allowedRatioGraph_tau4_cliqueFreeOn_three :
    (allowedRatioGraph tau4).CliqueFreeOn (tau4Vertices : Set ℚ) 3 :=
  SimpleGraph.CliqueFreeOn.anti
    (allowedRatioGraph tau4) graph
    allowedRatioGraph_tau4_le_graph
    finiteCoverGraph_tau4_cliqueFreeOn_three

theorem allowedRatioGraph_tau4_cliqueFree_three :
    (allowedRatioGraph tau4).CliqueFree 3 :=
  allowedRatioGraph_cliqueFree_of_finite_cover tau4 tau4Vertices 3
    (by omega) ratioAllowed_tau4_mem
    allowedRatioGraph_tau4_cliqueFreeOn_three

/-- The `free_tau4` producer required by the THM-954 continuum certificate. -/
theorem thresholdGraph_tau4_cliqueFree_four
    (v : Fin 13 → ℤ)
    (weight : LRCB5ContinuumFloor.PairSupport → ℚ)
    (hv : ∀ index, v index ≠ 0)
    (hdistinct : ∀ first second, first ≠ second →
      |v first| ≠ |v second|)
    (hweight : ∀ first second, first ≠ second →
      weight {first, second} = pairDeficit (v first) (v second)) :
    (thresholdGraph weight tau4).CliqueFree 4 := by
  simpa using thresholdGraph_cliqueFree_of_allowedRatioGraph
    v weight tau4 3 (by omega) hv hdistinct hweight
      allowedRatioGraph_tau4_cliqueFree_three

#print axioms crossCheckBlock_true
#print axioms cross_edge_maps
#print axioms finiteCoverGraph_tau4_cliqueFreeOn_three
#print axioms parallel_class_fingerprint
#print axioms tau4FiniteRatios_eq_tau4Vertices
#print axioms ratioAllowed_tau4_mem
#print axioms allowedRatioGraph_tau4_cliqueFree_three
#print axioms thresholdGraph_tau4_cliqueFree_four

end LRCPairRatioTau4
end LonelyRunner
