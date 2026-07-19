import Mathlib

/-!
# The common-period beat-hole arithmetic core (THM-1216)

For a reduced period `Q`, the strict radius-`1/14` danger mask has

`A(Q) = 2 * ceil(Q/14) - 1`

residues.  Five masks are relevant to a fast-fast difference beat: one for
the defining pair and four for the complementary fast speeds.  Every mask
contains residue zero.  Their union therefore has cardinality at most
`5*A(Q)-4`, and

`B(Q) = 5*A(Q)-3`

consecutive distinct residues force a point outside the union.  The key
arithmetic fact is `B(Q) <= Q` for every `Q >= 2`.

This file kernel-checks that arithmetic and the five-finset common-point
ledger.  The geometric bridge from a slow-gap interval to a consecutive
integer block, and the gcd reduction of each speed to a unit mask modulo
`Q`, remain explicit provider hypotheses in `commonPeriod_bridge_impossible`.
-/

namespace LRC14.CommonPeriodBeatHole

open scoped Classical

/-- Natural-number ceiling of `Q/14`. -/
def ceilDiv14 (Q : ℕ) : ℕ :=
  (Q + 13) / 14

/-- Number of strict dangerous residues in one reduced period. -/
def windowCount (Q : ℕ) : ℕ :=
  2 * ceilDiv14 Q - 1

/-- Uniform five-mask escape threshold. -/
def blockThreshold (Q : ℕ) : ℕ :=
  5 * windowCount Q - 3

/-- Class-sensitive threshold when exactly `r` distinct masks occur. -/
def classThreshold (Q r : ℕ) : ℕ :=
  2 + r * (windowCount Q - 1)

theorem one_le_ceilDiv14 (Q : ℕ) (hQ : 1 ≤ Q) :
    1 ≤ ceilDiv14 Q := by
  unfold ceilDiv14
  omega

theorem one_le_windowCount (Q : ℕ) (hQ : 1 ≤ Q) :
    1 ≤ windowCount Q := by
  unfold windowCount
  have := one_le_ceilDiv14 Q hQ
  omega

/-- The class-sensitive formula specializes to the uniform five-mask one. -/
theorem classThreshold_five (Q : ℕ) (hQ : 1 ≤ Q) :
    classThreshold Q 5 = blockThreshold Q := by
  unfold classThreshold blockThreshold
  have := one_le_windowCount Q hQ
  omega

/-- The forced common-zero saving makes the five-mask threshold fit inside
one reduced period for every nontrivial period. -/
theorem blockThreshold_le_period (Q : ℕ) (hQ : 2 ≤ Q) :
    blockThreshold Q ≤ Q := by
  unfold blockThreshold windowCount ceilDiv14
  omega

/-- At reduced period fourteen the exact threshold is two, not six. -/
theorem blockThreshold_fourteen : blockThreshold 14 = 2 := by
  decide

/-- Every nontrivial reduced period at most fourteen has singleton danger
masks, so the class-sensitive threshold is exactly two regardless of how
many mask labels are presented. -/
theorem classThreshold_eq_two_of_small_period
    (Q r : ℕ) (hQ : 1 ≤ Q) (hsmall : Q ≤ 14) :
    classThreshold Q r = 2 := by
  unfold classThreshold windowCount ceilDiv14
  omega

/-- One common point saves four units from the naive sum of five row sizes. -/
theorem five_union_card_add_four_le_sum_cards
    {α : Type*} [DecidableEq α]
    (x : α) (s0 s1 s2 s3 s4 : Finset α)
    (hx0 : x ∈ s0) (hx1 : x ∈ s1) (hx2 : x ∈ s2)
    (hx3 : x ∈ s3) (hx4 : x ∈ s4) :
    (s0 ∪ s1 ∪ s2 ∪ s3 ∪ s4).card + 4 ≤
      s0.card + s1.card + s2.card + s3.card + s4.card := by
  have h01 := Finset.card_union_add_card_inter s0 s1
  have h012 := Finset.card_union_add_card_inter (s0 ∪ s1) s2
  have h0123 := Finset.card_union_add_card_inter (s0 ∪ s1 ∪ s2) s3
  have h01234 :=
    Finset.card_union_add_card_inter (s0 ∪ s1 ∪ s2 ∪ s3) s4
  have hi01 : 1 ≤ (s0 ∩ s1).card := by
    apply Finset.card_pos.mpr
    exact ⟨x, by simp [hx0, hx1]⟩
  have hi012 : 1 ≤ ((s0 ∪ s1) ∩ s2).card := by
    apply Finset.card_pos.mpr
    exact ⟨x, by simp [hx0, hx2]⟩
  have hi0123 : 1 ≤ ((s0 ∪ s1 ∪ s2) ∩ s3).card := by
    apply Finset.card_pos.mpr
    exact ⟨x, by simp [hx0, hx3]⟩
  have hi01234 : 1 ≤ ((s0 ∪ s1 ∪ s2 ∪ s3) ∩ s4).card := by
    apply Finset.card_pos.mpr
    exact ⟨x, by simp [hx0, hx4]⟩
  omega

/-- Five common-zero masks of individual size at most `A` have union size
at most `5*A-4`, written subtraction-free. -/
theorem five_union_card_add_four_le_five_mul
    {α : Type*} [DecidableEq α]
    (x : α) (s0 s1 s2 s3 s4 : Finset α) (A : ℕ)
    (hx0 : x ∈ s0) (hx1 : x ∈ s1) (hx2 : x ∈ s2)
    (hx3 : x ∈ s3) (hx4 : x ∈ s4)
    (h0 : s0.card ≤ A) (h1 : s1.card ≤ A) (h2 : s2.card ≤ A)
    (h3 : s3.card ≤ A) (h4 : s4.card ≤ A) :
    (s0 ∪ s1 ∪ s2 ∪ s3 ∪ s4).card + 4 ≤ 5 * A := by
  have hsaving := five_union_card_add_four_le_sum_cards
    x s0 s1 s2 s3 s4 hx0 hx1 hx2 hx3 hx4
  omega

/-- Class-sensitive arithmetic bridge.  Once a provider has compressed the
five labelled masks to `r` distinct classes and proved the common-point union
bound `1 + r*(A-1)`, `B_r=2+r*(A-1)` block residues cannot all be covered.
This theorem deliberately separates that mask-nerve provider from the finite
cardinality consumer. -/
theorem classThreshold_bridge_impossible
    {α : Type*} [DecidableEq α]
    (Q r : ℕ) (block maskUnion : Finset α)
    (hblock : classThreshold Q r ≤ block.card)
    (hunion : maskUnion.card ≤ 1 + r * (windowCount Q - 1))
    (hcover : block ⊆ maskUnion) : False := by
  have hcard := Finset.card_le_card hcover
  unfold classThreshold at hblock
  omega

/-- Typed finite bridge: a block of at least `B(Q)` distinct residues cannot
be contained in five common-zero masks, each of size at most `A(Q)`. -/
theorem commonPeriod_bridge_impossible
    {α : Type*} [DecidableEq α]
    (Q : ℕ) (hQ : 2 ≤ Q) (x : α)
    (block s0 s1 s2 s3 s4 : Finset α)
    (hx0 : x ∈ s0) (hx1 : x ∈ s1) (hx2 : x ∈ s2)
    (hx3 : x ∈ s3) (hx4 : x ∈ s4)
    (hblock : blockThreshold Q ≤ block.card)
    (h0 : s0.card ≤ windowCount Q) (h1 : s1.card ≤ windowCount Q)
    (h2 : s2.card ≤ windowCount Q) (h3 : s3.card ≤ windowCount Q)
    (h4 : s4.card ≤ windowCount Q)
    (hcover : block ⊆ s0 ∪ s1 ∪ s2 ∪ s3 ∪ s4) : False := by
  have hunion := five_union_card_add_four_le_five_mul
    x s0 s1 s2 s3 s4 (windowCount Q)
    hx0 hx1 hx2 hx3 hx4 h0 h1 h2 h3 h4
  have hcard := Finset.card_le_card hcover
  have hA := one_le_windowCount Q (by omega)
  unfold blockThreshold at hblock
  omega

#print axioms blockThreshold_le_period
#print axioms blockThreshold_fourteen
#print axioms classThreshold_eq_two_of_small_period
#print axioms five_union_card_add_four_le_sum_cards
#print axioms five_union_card_add_four_le_five_mul
#print axioms classThreshold_bridge_impossible
#print axioms commonPeriod_bridge_impossible

end LRC14.CommonPeriodBeatHole
