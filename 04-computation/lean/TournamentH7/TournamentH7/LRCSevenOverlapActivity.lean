/-
  TournamentH7.LRCSevenOverlapActivity

  The first activity-weighted bridge from THM-947's colored overlap graph to
  THM-950's multiplier census.  For a fixed ordered triple, bad multipliers
  split exactly into rank-one (zero base color) and genuinely colored events.
  Every colored event above the dense pair pays at least three units of spoke
  mass; a nonzero lower triangle pays at least five.  Two unit spokes route an
  event back to the rank-one side.

  Tournament-analysis audit: the proof-relevant vertices here are multiplier
  events, not runners.  Each vertex carries a colored runner graph, with
  `overlapDet` as the pair observable and zero/nonzero as the switch.  The
  eventwise partition and summed color mass are preserved.  An uncolored
  runner tournament destroys both coefficient magnitude and activity
  multiplicity, so it cannot feed the census.  The challenged assumption is
  that one static runner graph should organize all phases; the faithful object
  is instead a fibered movie of colored graphs over the multiplier set.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCSevenOverlapDenseCore

namespace LonelyRunner
namespace LRC14Grand

open LonelyRunner
open LRC14Concrete
open scoped Classical

/-- Multipliers at which one fixed ordered runner triple is simultaneously
bad. -/
def badTripleMultipliers (v : Fin 13 → ℤ) (q : ℕ)
    (σ : Equiv.Perm (Fin 13)) (left right top : Fin 13) : Finset ℕ :=
  (Finset.Ioo 0 q).filter fun p =>
    ¬ inBand v q p (σ left) ∧
    ¬ inBand v q p (σ right) ∧
    ¬ inBand v q p (σ top)

/-- Rank-one events: the two lower witnesses are exactly proportional to the
two lower speeds. -/
def alignedBadTripleMultipliers (v : Fin 13 → ℤ) (q : ℕ)
    (σ : Equiv.Perm (Fin 13)) (left right top : Fin 13) : Finset ℕ :=
  (badTripleMultipliers v q σ left right top).filter fun p =>
    overlapDet v q p (σ left) (σ right) = 0

/-- Genuinely colored events: the lower pair carries a nonzero determinant. -/
def coloredBadTripleMultipliers (v : Fin 13 → ℤ) (q : ℕ)
    (σ : Equiv.Perm (Fin 13)) (left right top : Fin 13) : Finset ℕ :=
  (badTripleMultipliers v q σ left right top).filter fun p =>
    overlapDet v q p (σ left) (σ right) ≠ 0

/-- The rank-one/colored dichotomy is an exact partition of the fixed-triple
activity. -/
theorem badTripleMultipliers_eq_aligned_union_colored
    (v : Fin 13 → ℤ) (q : ℕ) (σ : Equiv.Perm (Fin 13))
    (left right top : Fin 13) :
    badTripleMultipliers v q σ left right top =
      alignedBadTripleMultipliers v q σ left right top ∪
        coloredBadTripleMultipliers v q σ left right top := by
  ext p
  by_cases hcolor : overlapDet v q p (σ left) (σ right) = 0 <;>
    simp [alignedBadTripleMultipliers, coloredBadTripleMultipliers, hcolor]

theorem alignedBadTripleMultipliers_disjoint_colored
    (v : Fin 13 → ℤ) (q : ℕ) (σ : Equiv.Perm (Fin 13))
    (left right top : Fin 13) :
    Disjoint (alignedBadTripleMultipliers v q σ left right top)
      (coloredBadTripleMultipliers v q σ left right top) := by
  rw [Finset.disjoint_left]
  intro p haligned hcolored
  have hzero := (Finset.mem_filter.mp haligned).2
  have hne := (Finset.mem_filter.mp hcolored).2
  exact hne hzero

/-- Exact activity ledger for a fixed triple. -/
theorem badTripleMultipliers_card_eq_aligned_add_colored
    (v : Fin 13 → ℤ) (q : ℕ) (σ : Equiv.Perm (Fin 13))
    (left right top : Fin 13) :
    (badTripleMultipliers v q σ left right top).card =
      (alignedBadTripleMultipliers v q σ left right top).card +
        (coloredBadTripleMultipliers v q σ left right top).card := by
  rw [badTripleMultipliers_eq_aligned_union_colored,
    Finset.card_union_of_disjoint
      (alignedBadTripleMultipliers_disjoint_colored
        v q σ left right top)]

/-- **Activity-weighted three-unit charge.**  Summing the pointwise Plücker
constraint over all genuinely colored bad events retains phase multiplicity:
each event contributes at least three units on the two high spokes. -/
theorem three_mul_coloredBadTriple_card_le_spokeMass
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q : ℕ) (σ : Equiv.Perm (Fin 13))
    (hmono : Monotone (fun index => |v (σ index)|))
    (dense : Fin 12)
    (hladder : ∀ r : Fin 12, dense < r →
      3 * |v (σ r.castSucc)| ≤ |v (σ r.succ)|)
    (left right top : Fin 13)
    (hleft : (left : ℕ) < (top : ℕ))
    (hright : (right : ℕ) < (top : ℕ))
    (htopIndex : (dense : ℕ) + 2 ≤ (top : ℕ)) :
    (3 : ℤ) * (coloredBadTripleMultipliers
      v q σ left right top).card ≤
      ∑ p ∈ coloredBadTripleMultipliers v q σ left right top,
        (|overlapDet v q p (σ right) (σ top)| +
          |overlapDet v q p (σ left) (σ top)|) := by
  calc
    (3 : ℤ) * (coloredBadTripleMultipliers
        v q σ left right top).card =
        ∑ _p ∈ coloredBadTripleMultipliers v q σ left right top,
          (3 : ℤ) := by simp [mul_comm]
    _ ≤ _ := by
      apply Finset.sum_le_sum
      intro p hp
      have hbase : overlapDet v q p (σ left) (σ right) ≠ 0 :=
        (Finset.mem_filter.mp hp).2
      exact overlapDet_incident_mass_ge_three_of_high_base_ne_zero
        v hv q p σ hmono dense hladder left right top
          hleft hright htopIndex hbase

/-- Four-runner activity whose three lower edges all carry nonzero colors. -/
def coloredBadTriangleMultipliers (v : Fin 13 → ℤ) (q : ℕ)
    (σ : Equiv.Perm (Fin 13)) (a b c top : Fin 13) : Finset ℕ :=
  (Finset.Ioo 0 q).filter fun p =>
    ¬ inBand v q p (σ a) ∧
    ¬ inBand v q p (σ b) ∧
    ¬ inBand v q p (σ c) ∧
    ¬ inBand v q p (σ top) ∧
    overlapDet v q p (σ a) (σ b) ≠ 0 ∧
    overlapDet v q p (σ a) (σ c) ≠ 0 ∧
    overlapDet v q p (σ b) (σ c) ≠ 0

/-- **Activity-weighted five-unit charge.**  Every active nonzero lower
triangle contributes at least five determinant units on its three spokes. -/
theorem five_mul_coloredBadTriangle_card_le_spokeMass
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q : ℕ) (σ : Equiv.Perm (Fin 13))
    (hmono : Monotone (fun index => |v (σ index)|))
    (dense : Fin 12)
    (hladder : ∀ r : Fin 12, dense < r →
      3 * |v (σ r.castSucc)| ≤ |v (σ r.succ)|)
    (a b c top : Fin 13)
    (ha : (a : ℕ) < (top : ℕ))
    (hb : (b : ℕ) < (top : ℕ))
    (hc : (c : ℕ) < (top : ℕ))
    (htopIndex : (dense : ℕ) + 2 ≤ (top : ℕ)) :
    (5 : ℤ) * (coloredBadTriangleMultipliers
      v q σ a b c top).card ≤
      ∑ p ∈ coloredBadTriangleMultipliers v q σ a b c top,
        (|overlapDet v q p (σ a) (σ top)| +
          |overlapDet v q p (σ b) (σ top)| +
          |overlapDet v q p (σ c) (σ top)|) := by
  calc
    (5 : ℤ) * (coloredBadTriangleMultipliers
        v q σ a b c top).card =
        ∑ _p ∈ coloredBadTriangleMultipliers v q σ a b c top,
          (5 : ℤ) := by simp [mul_comm]
    _ ≤ _ := by
      apply Finset.sum_le_sum
      intro p hp
      obtain ⟨_haBad, _hbBad, _hcBad, _htopBad, hab, hac, hbc⟩ :=
        (Finset.mem_filter.mp hp).2
      exact overlapDet_three_spoke_mass_ge_five_of_high_nonzero_triangle
        v hv q p σ hmono dense hladder a b c top ha hb hc htopIndex
          hab hac hbc

/-- Bad events whose two high spokes both have unit color. -/
def unitSpokeBadTripleMultipliers (v : Fin 13 → ℤ) (q : ℕ)
    (σ : Equiv.Perm (Fin 13)) (left right top : Fin 13) : Finset ℕ :=
  (badTripleMultipliers v q σ left right top).filter fun p =>
    |overlapDet v q p (σ left) (σ top)| ≤ 1 ∧
    |overlapDet v q p (σ right) (σ top)| ≤ 1

/-- Unit-spoke activity is forced into the rank-one side of the exact event
partition. -/
theorem unitSpokeBadTripleMultipliers_subset_aligned
    (v : Fin 13 → ℤ) (hv : ∀ index, v index ≠ 0)
    (q : ℕ) (σ : Equiv.Perm (Fin 13))
    (hmono : Monotone (fun index => |v (σ index)|))
    (dense : Fin 12)
    (hladder : ∀ r : Fin 12, dense < r →
      3 * |v (σ r.castSucc)| ≤ |v (σ r.succ)|)
    (left right top : Fin 13)
    (hleft : (left : ℕ) < (top : ℕ))
    (hright : (right : ℕ) < (top : ℕ))
    (htopIndex : (dense : ℕ) + 2 ≤ (top : ℕ)) :
  unitSpokeBadTripleMultipliers v q σ left right top ⊆
      alignedBadTripleMultipliers v q σ left right top := by
  intro p hp
  unfold unitSpokeBadTripleMultipliers at hp
  unfold alignedBadTripleMultipliers
  rw [Finset.mem_filter] at hp ⊢
  refine ⟨hp.1, ?_⟩
  exact overlapDet_base_eq_zero_of_two_incident_unit
    v hv q p σ hmono dense hladder left right top hleft hright htopIndex
      hp.2.1 hp.2.2

/-! ## Axiom audit -/

#print axioms badTripleMultipliers_eq_aligned_union_colored
#print axioms badTripleMultipliers_card_eq_aligned_add_colored
#print axioms three_mul_coloredBadTriple_card_le_spokeMass
#print axioms five_mul_coloredBadTriangle_card_le_spokeMass
#print axioms unitSpokeBadTripleMultipliers_subset_aligned

end LRC14Grand
end LonelyRunner
