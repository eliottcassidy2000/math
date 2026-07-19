/-
  THM-1258: kernel-pure consumers for the complementary Z/4 and Z/9
  scale-thirty-six fibre certificate.

  The external exact bank proves that every scalar-surviving row has an owner
  at which at least one independent-choice relaxation is strictly below 36.
  This file checks both logical steps used by that certificate:

    * an anchor/deviation ceiling below 36 prevents full 36-sheet cover;
    * if literal cover implies liveness in both gauges at every owner, one
      certified failed gauge at one owner prevents all-owner cover.

  The 206,725,596 scalar contexts remain in the companion Python certificate.
  No global n=12 sporadic-branch conclusion is asserted here.
-/
import TournamentH7.LRCNestedFibreRelaxation

namespace LonelyRunner
namespace LRC14
namespace ScaleThirtySixComplementaryFibre

open Finset
open NestedFibreRelaxation

variable {provider choice : Type*}
variable [DecidableEq provider]

/-- A concrete independent-choice anchor relaxation below 36 rules out the
literal full sheet set.  This applies unchanged to the Z/4 and Z/9 anchors. -/
theorem strict_relaxation_blocks_full_cover
    (providers anchors : Finset provider)
    (options : provider → Finset choice)
    (mask : provider → choice → Finset (Fin 36))
    (chosen : provider → choice)
    (hanchors : anchors ⊆ providers)
    (hchosen : ∀ i ∈ providers \ anchors, chosen i ∈ options i)
    (hceiling :
      (anchorUnion anchors mask chosen).card +
          ∑ i ∈ providers \ anchors,
            deviationMax (anchorUnion anchors mask chosen) options mask i < 36) :
    providerUnion providers mask chosen ≠ univ := by
  have hcard : (providerUnion providers mask chosen).card < 36 :=
    (card_providerUnion_le_anchor_add_pointwiseMax
      providers anchors options mask chosen hanchors hchosen).trans_lt hceiling
  intro hfull
  rw [hfull] at hcard
  norm_num at hcard

variable {row owner : Type*}

/-- The exact paired-flag terminality statement.  A literal cover must be
live in every sound upper-bound gauge.  The finite certificate supplies one
owner where at least one of the two gauges is below 36, so all-owner cover is
impossible. -/
theorem complementary_failed_owner_blocks_all_owner_cover
    (covers : row → owner → Prop)
    (z4Bound z9Bound : row → owner → ℕ)
    (hcover4 : ∀ r o, covers r o → 36 ≤ z4Bound r o)
    (hcover9 : ∀ r o, covers r o → 36 ≤ z9Bound r o)
    (hcertificate : ∀ r, ∃ o,
      z4Bound r o < 36 ∨ z9Bound r o < 36) :
    ∀ r, ¬ ∀ o, covers r o := by
  intro r hall
  obtain ⟨o, h4 | h9⟩ := hcertificate r
  · exact (Nat.not_le_of_lt h4) (hcover4 r o (hall o))
  · exact (Nat.not_le_of_lt h9) (hcover9 r o (hall o))

#print axioms strict_relaxation_blocks_full_cover
#print axioms complementary_failed_owner_blocks_all_owner_cover

end ScaleThirtySixComplementaryFibre
end LRC14
end LonelyRunner
