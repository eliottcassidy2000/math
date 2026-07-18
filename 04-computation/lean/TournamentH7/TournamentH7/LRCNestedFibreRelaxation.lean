/-
  TournamentH7.LRCNestedFibreRelaxation -- the generic finite relaxation
  behind THM-994's scale-twenty-seven, THM-1072's scale-twenty-eight,
  THM-1090's complementary scale-thirty, and THM-1096's scale-thirty-two
  fibre obstructions.

  Fix an anchor union `Q`.  Every nonanchor mask contributes only its part
  outside `Q`; overlaps among those deviations can only reduce the literal
  union.  Consequently

    |Q ∪ ⋃ i, M i| ≤ |Q| + ∑ i, |M i \ Q|.

  Replacing each chosen deviation by its maximum over an independent finite
  choice space is therefore a sound upper relaxation.  No periodicity premise
  is needed for this finite-set inequality; in THM-994 periodicity is what
  makes the anchor bank small enough to enumerate.  The concrete finite
  banks remain external certificates; this file kernel-checks the generic
  relaxation they consume.

  Kernel-pure: no `sorry`, no `native_decide`, and no new axioms.
-/
import Mathlib

namespace LonelyRunner
namespace LRC14
namespace NestedFibreRelaxation

open Finset

variable {provider choice point : Type*}
variable [DecidableEq point]

/-- The union supplied by a finite set of providers at fixed choices. -/
def providerUnion (providers : Finset provider)
    (mask : provider → choice → Finset point) (chosen : provider → choice) :
    Finset point :=
  providers.biUnion fun i => mask i (chosen i)

/-- The masks of the designated anchor providers. -/
def anchorUnion (anchors : Finset provider)
    (mask : provider → choice → Finset point) (chosen : provider → choice) :
    Finset point :=
  providerUnion anchors mask chosen

/-- Adding masks to an anchor union costs at most the sum of their deviations
outside the anchor.  This form is independent of how `Q` was constructed. -/
theorem card_union_biUnion_le_anchor_add_deviations
    (Q : Finset point) (providers : Finset provider) (mask : provider → Finset point) :
    (Q ∪ providers.biUnion mask).card ≤
      Q.card + ∑ i ∈ providers, (mask i \ Q).card := by
  have hunion : Q ∪ providers.biUnion mask =
      Q ∪ providers.biUnion (fun i => mask i \ Q) := by
    ext x
    by_cases hx : x ∈ Q <;> simp [hx]
  rw [hunion]
  calc
    (Q ∪ providers.biUnion fun i => mask i \ Q).card
        ≤ Q.card + (providers.biUnion fun i => mask i \ Q).card :=
          Finset.card_union_le _ _
    _ ≤ Q.card + ∑ i ∈ providers, (mask i \ Q).card :=
      Nat.add_le_add_left Finset.card_biUnion_le Q.card

variable [DecidableEq provider]

/-- Exact decomposition of a total provider union into its anchor union and
the remaining providers. -/
theorem providerUnion_eq_anchor_union_remainder
    (providers anchors : Finset provider)
    (mask : provider → choice → Finset point) (chosen : provider → choice)
    (hanchors : anchors ⊆ providers) :
    providerUnion providers mask chosen =
      anchorUnion anchors mask chosen ∪
        providerUnion (providers \ anchors) mask chosen := by
  ext x
  constructor
  · intro hx
    obtain ⟨i, hi, hxi⟩ := Finset.mem_biUnion.mp hx
    by_cases hia : i ∈ anchors
    · exact Finset.mem_union_left _
        (Finset.mem_biUnion.mpr ⟨i, hia, hxi⟩)
    · exact Finset.mem_union_right _
        (Finset.mem_biUnion.mpr ⟨i, Finset.mem_sdiff.mpr ⟨hi, hia⟩, hxi⟩)
  · intro hx
    rcases Finset.mem_union.mp hx with hx | hx
    · obtain ⟨i, hi, hxi⟩ := Finset.mem_biUnion.mp hx
      exact Finset.mem_biUnion.mpr ⟨i, hanchors hi, hxi⟩
    · obtain ⟨i, hi, hxi⟩ := Finset.mem_biUnion.mp hx
      exact Finset.mem_biUnion.mpr ⟨i, (Finset.mem_sdiff.mp hi).1, hxi⟩

/-- **Nested-fibre anchor relaxation, literal form.**  The total literal
union is bounded by the anchor union plus the sum of nonanchor deviations. -/
theorem card_providerUnion_le_anchor_add_deviations
    (providers anchors : Finset provider)
    (mask : provider → choice → Finset point) (chosen : provider → choice)
    (hanchors : anchors ⊆ providers) :
    (providerUnion providers mask chosen).card ≤
      (anchorUnion anchors mask chosen).card +
        ∑ i ∈ providers \ anchors,
          (mask i (chosen i) \ anchorUnion anchors mask chosen).card := by
  rw [providerUnion_eq_anchor_union_remainder providers anchors mask chosen hanchors]
  exact card_union_biUnion_le_anchor_add_deviations
    (anchorUnion anchors mask chosen) (providers \ anchors)
    (fun i => mask i (chosen i))

/-- Any pointwise bounds on nonanchor deviations can replace the literal
deviation terms in the relaxation. -/
theorem card_providerUnion_le_anchor_add_pointwiseBounds
    (providers anchors : Finset provider)
    (mask : provider → choice → Finset point) (chosen : provider → choice)
    (bound : provider → ℕ) (hanchors : anchors ⊆ providers)
    (hbound : ∀ i ∈ providers \ anchors,
      (mask i (chosen i) \ anchorUnion anchors mask chosen).card ≤ bound i) :
    (providerUnion providers mask chosen).card ≤
      (anchorUnion anchors mask chosen).card +
        ∑ i ∈ providers \ anchors, bound i := by
  refine (card_providerUnion_le_anchor_add_deviations
    providers anchors mask chosen hanchors).trans ?_
  exact Nat.add_le_add_left (Finset.sum_le_sum hbound)
    (anchorUnion anchors mask chosen).card

/-- The independent-choice deviation maximum for one provider.  `Finset.sup`
returns zero on an empty choice space, so consumers should supply membership
of every actually chosen nonanchor option. -/
def deviationMax (Q : Finset point) (options : provider → Finset choice)
    (mask : provider → choice → Finset point) (i : provider) : ℕ :=
  (options i).sup fun e => (mask i e \ Q).card

/-- **The THM-994/1072/1090/1096 pointwise-maximum relaxation.**  Each nonanchor provider may
independently take its best option.  This forgets compatibility and overlap,
so it is a sound upper bound on every literal simultaneous choice. -/
theorem card_providerUnion_le_anchor_add_pointwiseMax
    (providers anchors : Finset provider)
    (options : provider → Finset choice)
    (mask : provider → choice → Finset point) (chosen : provider → choice)
    (hanchors : anchors ⊆ providers)
    (hchosen : ∀ i ∈ providers \ anchors, chosen i ∈ options i) :
    (providerUnion providers mask chosen).card ≤
      (anchorUnion anchors mask chosen).card +
        ∑ i ∈ providers \ anchors,
          deviationMax (anchorUnion anchors mask chosen) options mask i := by
  apply card_providerUnion_le_anchor_add_pointwiseBounds
    providers anchors mask chosen
      (deviationMax (anchorUnion anchors mask chosen) options mask) hanchors
  intro i hi
  exact Finset.le_sup (f := fun e =>
    (mask i e \ anchorUnion anchors mask chosen).card) (hchosen i hi)

/-! ## Axiom audit -/

#print axioms card_union_biUnion_le_anchor_add_deviations
#print axioms providerUnion_eq_anchor_union_remainder
#print axioms card_providerUnion_le_anchor_add_deviations
#print axioms card_providerUnion_le_anchor_add_pointwiseBounds
#print axioms card_providerUnion_le_anchor_add_pointwiseMax

end NestedFibreRelaxation
end LRC14
end LonelyRunner
