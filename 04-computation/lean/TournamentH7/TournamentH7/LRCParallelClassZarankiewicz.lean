/-
  Exact three-row Zarankiewicz ledger for the parallel-class circle. It shows
  that the selected-witness residuals attain z(3,g;2,1)=g, so static incidence
  has no slack; cyclic phase or wall-event information is indispensable.
  No \`sorry\`; no \`native_decide\`.
-/

import TournamentH7.LRCDetunedOverlap

namespace LonelyRunner
namespace LRC14Grand

open scoped Classical

/-- Three-row incidence ledger: excess edge count above the occupied right
vertices is paid by pair codegrees. -/
theorem three_row_edges_le_union_add_pair_codegrees
    {α : Type*} [DecidableEq α] (first second third : Finset α) :
    first.card + second.card + third.card ≤
      (first ∪ second ∪ third).card +
        (first ∩ second).card + (first ∩ third).card +
        (second ∩ third).card := by
  have hfirstSecond := Finset.card_union_add_card_inter first second
  have hstep := Finset.card_union_add_card_inter (first ∪ second) third
  have hinterEq :
      (first ∪ second) ∩ third =
        (first ∩ third) ∪ (second ∩ third) := by
    ext element
    simp only [Finset.mem_inter, Finset.mem_union]
    tauto
  have hinterBound :
      ((first ∪ second) ∩ third).card ≤
        (first ∩ third).card + (second ∩ third).card := by
    rw [hinterEq]
    exact Finset.card_union_le _ _
  omega


/-- Exact three-row inclusion-exclusion, retaining the triple codegree. -/
theorem three_row_edges_add_triple_eq_union_add_pair_codegrees
    {α : Type*} [DecidableEq α] (first second third : Finset α) :
    first.card + second.card + third.card +
        (first ∩ second ∩ third).card =
      (first ∪ second ∪ third).card +
        (first ∩ second).card + (first ∩ third).card +
        (second ∩ third).card := by
  have hfirstSecond := Finset.card_union_add_card_inter first second
  have hstep := Finset.card_union_add_card_inter (first ∪ second) third
  have hinterEq :
      (first ∪ second) ∩ third =
        (first ∩ third) ∪ (second ∩ third) := by
    ext element
    simp only [Finset.mem_inter, Finset.mem_union]
    tauto
  have hpairStep := Finset.card_union_add_card_inter
    (first ∩ third) (second ∩ third)
  have htripleEq :
      (first ∩ third) ∩ (second ∩ third) = first ∩ second ∩ third := by
    ext element
    simp only [Finset.mem_inter]
    tauto
  rw [hinterEq] at hstep
  rw [htripleEq] at hpairStep
  omega

/-- Pair codegrees themselves form a three-row family whose only repeated
vertices are the triple incidences. -/
theorem pair_codegrees_le_union_add_two_triples
    {α : Type*} [DecidableEq α] (first second third : Finset α) :
    (first ∩ second).card + (first ∩ third).card +
        (second ∩ third).card ≤
      (first ∪ second ∪ third).card +
        2 * (first ∩ second ∩ third).card := by
  let pair12 := first ∩ second
  let pair13 := first ∩ third
  let pair23 := second ∩ third
  let triple := first ∩ second ∩ third
  have h1213 : pair12 ∩ pair13 = triple := by
    ext element
    simp [pair12, pair13, triple, and_left_comm]
  have hunion23 : (pair12 ∪ pair13) ∩ pair23 = triple := by
    ext element
    simp [pair12, pair13, pair23, triple]
    tauto
  have hfirst := Finset.card_union_add_card_inter pair12 pair13
  have hsecond := Finset.card_union_add_card_inter (pair12 ∪ pair13) pair23
  rw [h1213] at hfirst
  rw [hunion23] at hsecond
  have hsubset : pair12 ∪ pair13 ∪ pair23 ⊆ first ∪ second ∪ third := by
    intro element helement
    simp only [pair12, pair13, pair23, Finset.mem_union, Finset.mem_inter] at helement ⊢
    tauto
  have hcard := Finset.card_le_card hsubset
  dsimp [pair12, pair13, pair23, triple] at hfirst hsecond hcard ⊢
  omega

/-- The middle Zarankiewicz envelope: twice the edge count is at most three
times the carrier plus the total pair-codegree budget. -/
theorem two_mul_three_row_edges_le_three_mul_carrier_add_pair_codegrees
    {α : Type*} [DecidableEq α]
    (branches first second third : Finset α)
    (hfirst : first ⊆ branches) (hsecond : second ⊆ branches)
    (hthird : third ⊆ branches) :
    2 * (first.card + second.card + third.card) ≤
      3 * branches.card + (first ∩ second).card +
        (first ∩ third).card + (second ∩ third).card := by
  have hunionSubset : first ∪ second ∪ third ⊆ branches := by
    intro element helement
    simp only [Finset.mem_union] at helement
    rcases helement with (helement | helement) | helement
    · exact hfirst helement
    · exact hsecond helement
    · exact hthird helement
  have hunionCard := Finset.card_le_card hunionSubset
  have hexact :=
    three_row_edges_add_triple_eq_union_add_pair_codegrees first second third
  have hpairs := pair_codegrees_le_union_add_two_triples first second third
  omega

/-- Full three-left-vertex upper envelope for `K_{2,t}`-free incidence. -/
theorem three_row_zarankiewicz_full_upper
    {α : Type*} [DecidableEq α]
    (branches first second third : Finset α)
    (hfirst : first ⊆ branches) (hsecond : second ⊆ branches)
    (hthird : third ⊆ branches) (r : ℕ)
    (h12 : (first ∩ second).card ≤ r)
    (h13 : (first ∩ third).card ≤ r)
    (h23 : (second ∩ third).card ≤ r) :
    first.card + second.card + third.card ≤ 3 * branches.card ∧
    2 * (first.card + second.card + third.card) ≤
      3 * (branches.card + r) ∧
    first.card + second.card + third.card ≤ branches.card + 3 * r := by
  have hcardFirst := Finset.card_le_card hfirst
  have hcardSecond := Finset.card_le_card hsecond
  have hcardThird := Finset.card_le_card hthird
  have hmiddle :=
    two_mul_three_row_edges_le_three_mul_carrier_add_pair_codegrees
      branches first second third hfirst hsecond hthird
  have hunionSubset : first ∪ second ∪ third ⊆ branches := by
    intro element helement
    simp only [Finset.mem_union] at helement
    rcases helement with (helement | helement) | helement
    · exact hfirst helement
    · exact hsecond helement
    · exact hthird helement
  have hunionCard := Finset.card_le_card hunionSubset
  have hedge := three_row_edges_le_union_add_pair_codegrees first second third
  omega

/-- Zarankiewicz upper bound for three left vertices when every left pair has
codegree at most `r`. -/
theorem three_row_zarankiewicz_bound
    {α : Type*} [DecidableEq α]
    (branches first second third : Finset α)
    (hfirst : first ⊆ branches) (hsecond : second ⊆ branches)
    (hthird : third ⊆ branches) (r : ℕ)
    (h12 : (first ∩ second).card ≤ r)
    (h13 : (first ∩ third).card ≤ r)
    (h23 : (second ∩ third).card ≤ r) :
    first.card + second.card + third.card ≤ branches.card + 3 * r := by
  have hunionSubset : first ∪ second ∪ third ⊆ branches := by
    intro element helement
    simp only [Finset.mem_union] at helement
    rcases helement with (helement | helement) | helement
    · exact hfirst helement
    · exact hsecond helement
    · exact hthird helement
  have hunionCard := Finset.card_le_card hunionSubset
  have hedge := three_row_edges_le_union_add_pair_codegrees first second third
  omega

/-- At `t=1`, the relevant Zarankiewicz value is exactly the carrier size:
three pairwise-disjoint rows have at most one incidence per right vertex. -/
theorem three_row_zarankiewicz_t_one
    {α : Type*} [DecidableEq α]
    (branches first second third : Finset α)
    (hfirst : first ⊆ branches) (hsecond : second ⊆ branches)
    (hthird : third ⊆ branches)
    (h12 : Disjoint first second) (h13 : Disjoint first third)
    (h23 : Disjoint second third) :
    first.card + second.card + third.card ≤ branches.card := by
  apply three_row_zarankiewicz_bound branches first second third
      hfirst hsecond hthird 0
  · rw [Finset.disjoint_iff_inter_eq_empty.mp h12]
    simp
  · rw [Finset.disjoint_iff_inter_eq_empty.mp h13]
    simp
  · rw [Finset.disjoint_iff_inter_eq_empty.mp h23]
    simp

/-- Equality in the `t=1` bound is the exact parallel-class partition: every
right vertex has degree one. -/
theorem three_row_zarankiewicz_t_one_equality_iff_cover
    {α : Type*} [DecidableEq α]
    (branches first second third : Finset α)
    (hfirst : first ⊆ branches) (hsecond : second ⊆ branches)
    (hthird : third ⊆ branches)
    (h12 : Disjoint first second) (h13 : Disjoint first third)
    (h23 : Disjoint second third) :
    first.card + second.card + third.card = branches.card ↔
      first ∪ second ∪ third = branches := by
  have hpair : (first ∪ second).card = first.card + second.card := by
    exact Finset.card_union_of_disjoint h12
  have hUnionThird : Disjoint (first ∪ second) third := by
    rw [Finset.disjoint_left]
    intro element helement hthirdElement
    simp only [Finset.mem_union] at helement
    rcases helement with hfirstElement | hsecondElement
    · exact (Finset.disjoint_left.mp h13) hfirstElement hthirdElement
    · exact (Finset.disjoint_left.mp h23) hsecondElement hthirdElement
  have htotal :
      (first ∪ second ∪ third).card =
        first.card + second.card + third.card := by
    rw [Finset.card_union_of_disjoint hUnionThird, hpair]
  have hunionSubset : first ∪ second ∪ third ⊆ branches := by
    intro element helement
    simp only [Finset.mem_union] at helement
    rcases helement with (helement | helement) | helement
    · exact hfirst helement
    · exact hsecond helement
    · exact hthird helement
  constructor
  · intro heq
    exact Finset.eq_of_subset_of_card_le hunionSubset (by rw [htotal, heq])
  · intro heq
    rw [← htotal, heq]


/-! ## Axiom audit -/

#print axioms three_row_zarankiewicz_bound
#print axioms three_row_edges_add_triple_eq_union_add_pair_codegrees
#print axioms pair_codegrees_le_union_add_two_triples
#print axioms two_mul_three_row_edges_le_three_mul_carrier_add_pair_codegrees
#print axioms three_row_zarankiewicz_full_upper
#print axioms three_row_zarankiewicz_t_one
#print axioms three_row_zarankiewicz_t_one_equality_iff_cover

end LRC14Grand
end LonelyRunner
