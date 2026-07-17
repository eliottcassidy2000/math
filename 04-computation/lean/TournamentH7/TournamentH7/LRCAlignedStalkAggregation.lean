/-
  TournamentH7.LRCAlignedStalkAggregation

  The incidence-transposition bridge between the rooted six-face carrier and
  the fixed-stalk arithmetic of `LRCAlignedStalkGluing`.  At a fixed root,
  summing complete zero-color events first over multipliers or first over
  lower six-faces gives the same activity.  The fixed-stalk resonance bound
  can therefore be summed while retaining the true multiplier/face
  multiplicity.  This is an exact accounting identity; it does not itself
  bound the resulting global mass or eliminate reuse.

  The proof-relevant carrier is the bipartite incidence graph

      multiplier event  --  fixed rooted six-face.

  A static runner tournament forgets both sides of this incidence and cannot
  express the Fubini identity.  The zero-color switch is retained on each
  incidence; a lexicographic order on the fixed faces is a tie Hamiltonian
  path.

  No `sorry`; no `native_decide`.
-/

import TournamentH7.LRCSevenStalkReuseBudget
import TournamentH7.LRCAlignedStalkGluing

namespace LonelyRunner
namespace LRCAlignedStalkAggregation

open Finset
open LRC14Concrete
open scoped BigOperators Classical

/-- The fixed lower-six-face carrier attached to one root. -/
def rootedSixFaces (root : Fin 13) : Finset (Finset (Fin 13)) :=
  (Finset.univ.erase root).powersetCard 6

/-- All rooted seven-stalk events attached to a fixed root. -/
def rootedFaceActivity
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13) : ℕ :=
  ∑ p ∈ Finset.Ioo 0 q,
    ((rootedSixFaces root).filter fun face =>
      ∀ i ∈ insert root face, ¬ inBand v q p i).card

/-- Rooted events on which at least one root spoke has nonzero color. -/
def coloredRootSpokeFaceActivity
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13) : ℕ :=
  ∑ p ∈ Finset.Ioo 0 q,
    ((rootedSixFaces root).filter fun face =>
      (∀ i ∈ insert root face, ¬ inBand v q p i) ∧
      ∃ i ∈ face, overlapDet v q p i root ≠ 0).card

/-- Root-spoke mass restricted to the active rooted-face carrier. -/
def activeRootedFaceSpokeMass
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13) : ℕ :=
  ∑ p ∈ Finset.Ioo 0 q,
    ∑ face ∈ (rootedSixFaces root).filter (fun face =>
      ∀ i ∈ insert root face, ¬ inBand v q p i),
      ∑ i ∈ face,
        LRCSevenStalkReuseBudget.rootSpokeMass v q p root i

/-- Complete zero-color rooted-stalk activity, counted on the multiplier-first
side of the incidence graph. -/
def zeroColorRootedFaceActivity
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13) : ℕ :=
  ∑ p ∈ Finset.Ioo 0 q,
    ((rootedSixFaces root).filter fun face =>
      (∀ i ∈ insert root face, ¬ inBand v q p i) ∧
      (∀ i ∈ insert root face, ∀ j ∈ insert root face,
        overlapDet v q p i j = 0)).card

/-- **Exact incidence transposition.**  Multiplier-first complete zero-color
activity is exactly the sum of the fixed-stalk multiplier fibers. -/
theorem zeroColorRootedFaceActivity_eq_sum_stalkFibers
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13) :
    zeroColorRootedFaceActivity v q root =
      ∑ face ∈ rootedSixFaces root,
        (zeroColorStalkMultipliers v q (insert root face)).card := by
  unfold zeroColorRootedFaceActivity
  calc
    (∑ p ∈ Finset.Ioo 0 q,
        ((rootedSixFaces root).filter fun face =>
          (∀ i ∈ insert root face, ¬ inBand v q p i) ∧
          (∀ i ∈ insert root face, ∀ j ∈ insert root face,
            overlapDet v q p i j = 0)).card) =
        ∑ p ∈ Finset.Ioo 0 q,
          ∑ face ∈ rootedSixFaces root,
            if ((∀ i ∈ insert root face, ¬ inBand v q p i) ∧
              (∀ i ∈ insert root face, ∀ j ∈ insert root face,
                overlapDet v q p i j = 0)) then (1 : ℕ) else 0 := by
      apply Finset.sum_congr rfl
      intro p _hp
      rw [Finset.sum_boole]
      norm_cast
    _ = ∑ face ∈ rootedSixFaces root,
          ∑ p ∈ Finset.Ioo 0 q,
            if ((∀ i ∈ insert root face, ¬ inBand v q p i) ∧
              (∀ i ∈ insert root face, ∀ j ∈ insert root face,
                overlapDet v q p i j = 0)) then (1 : ℕ) else 0 := by
      rw [Finset.sum_comm]
    _ = ∑ face ∈ rootedSixFaces root,
          (zeroColorStalkMultipliers v q (insert root face)).card := by
      apply Finset.sum_congr rfl
      intro face _hface
      rw [Finset.sum_boole]
      rfl

/-- The same rooted activity expressed by the cheaper anchor-star certificate:
only the six colors incident to the root are retained. -/
def anchoredZeroColorRootedFaceActivity
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13) : ℕ :=
  ∑ p ∈ Finset.Ioo 0 q,
    ((rootedSixFaces root).filter fun face =>
      (∀ i ∈ insert root face, ¬ inBand v q p i) ∧
      (∀ i ∈ insert root face, overlapDet v q p i root = 0)).card

/-- Exact incidence transposition for the anchor-star carrier. -/
theorem anchoredZeroColorRootedFaceActivity_eq_sum_stalkFibers
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13) :
    anchoredZeroColorRootedFaceActivity v q root =
      ∑ face ∈ rootedSixFaces root,
        (anchoredZeroColorStalkMultipliers
          v q (insert root face) root).card := by
  unfold anchoredZeroColorRootedFaceActivity
  calc
    (∑ p ∈ Finset.Ioo 0 q,
        ((rootedSixFaces root).filter fun face =>
          (∀ i ∈ insert root face, ¬ inBand v q p i) ∧
          (∀ i ∈ insert root face,
            overlapDet v q p i root = 0)).card) =
        ∑ p ∈ Finset.Ioo 0 q,
          ∑ face ∈ rootedSixFaces root,
            if ((∀ i ∈ insert root face, ¬ inBand v q p i) ∧
              (∀ i ∈ insert root face,
                overlapDet v q p i root = 0)) then (1 : ℕ) else 0 := by
      apply Finset.sum_congr rfl
      intro p _hp
      rw [Finset.sum_boole]
      norm_cast
    _ = ∑ face ∈ rootedSixFaces root,
          ∑ p ∈ Finset.Ioo 0 q,
            if ((∀ i ∈ insert root face, ¬ inBand v q p i) ∧
              (∀ i ∈ insert root face,
                overlapDet v q p i root = 0)) then (1 : ℕ) else 0 := by
      rw [Finset.sum_comm]
    _ = ∑ face ∈ rootedSixFaces root,
          (anchoredZeroColorStalkMultipliers
            v q (insert root face) root).card := by
      apply Finset.sum_congr rfl
      intro face _hface
      rw [Finset.sum_boole]
      rfl

/-- Exact event partition on the rooted-face carrier: every active face has
either an all-zero root star or a nonzero root spoke, and never both. -/
theorem rootedFaceActivity_eq_anchored_add_colored
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13) :
    rootedFaceActivity v q root =
      anchoredZeroColorRootedFaceActivity v q root +
        coloredRootSpokeFaceActivity v q root := by
  unfold rootedFaceActivity anchoredZeroColorRootedFaceActivity
    coloredRootSpokeFaceActivity
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro p _hp
  let active := (rootedSixFaces root).filter fun face =>
    ∀ i ∈ insert root face, ¬ inBand v q p i
  let aligned := (rootedSixFaces root).filter fun face =>
    (∀ i ∈ insert root face, ¬ inBand v q p i) ∧
    (∀ i ∈ insert root face, overlapDet v q p i root = 0)
  let colored := (rootedSixFaces root).filter fun face =>
    (∀ i ∈ insert root face, ¬ inBand v q p i) ∧
    ∃ i ∈ face, overlapDet v q p i root ≠ 0
  have hunion : active = aligned ∪ colored := by
    ext face
    constructor
    · intro hactive
      have hface := Finset.mem_filter.mp hactive
      by_cases hzero : ∀ i ∈ insert root face,
          overlapDet v q p i root = 0
      · exact Finset.mem_union.mpr <| Or.inl <|
          Finset.mem_filter.mpr ⟨hface.1, hface.2, hzero⟩
      · apply Finset.mem_union.mpr
        apply Or.inr
        refine Finset.mem_filter.mpr ⟨hface.1, hface.2, ?_⟩
        by_contra hnone
        push Not at hnone
        apply hzero
        intro i hi
        rcases Finset.mem_insert.mp hi with rfl | hiface
        · simp [overlapDet]
        · exact hnone i hiface
    · intro hunionMem
      rcases Finset.mem_union.mp hunionMem with haligned | hcolored
      · exact Finset.mem_filter.mpr
          ⟨(Finset.mem_filter.mp haligned).1,
            (Finset.mem_filter.mp haligned).2.1⟩
      · exact Finset.mem_filter.mpr
          ⟨(Finset.mem_filter.mp hcolored).1,
            (Finset.mem_filter.mp hcolored).2.1⟩
  have hdisjoint : Disjoint aligned colored := by
    rw [Finset.disjoint_left]
    intro face haligned hcolored
    obtain ⟨_hface, _hbad, hzero⟩ := Finset.mem_filter.mp haligned
    obtain ⟨_hface', _hbad', i, hi, hne⟩ :=
      Finset.mem_filter.mp hcolored
    exact hne (hzero i (Finset.mem_insert_of_mem hi))
  change active.card = aligned.card + colored.card
  rw [hunion, Finset.card_union_of_disjoint hdisjoint]

/-- Every event in the colored complement contributes at least one unit of
root-spoke mass.  This universal one-unit statement is deliberately separate
from the stronger three-unit ladder charge. -/
theorem coloredRootSpokeFaceActivity_le_activeSpokeMass
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13) :
    coloredRootSpokeFaceActivity v q root ≤
      activeRootedFaceSpokeMass v q root := by
  unfold coloredRootSpokeFaceActivity activeRootedFaceSpokeMass
  apply Finset.sum_le_sum
  intro p _hp
  let active := (rootedSixFaces root).filter fun face =>
    ∀ i ∈ insert root face, ¬ inBand v q p i
  let colored := (rootedSixFaces root).filter fun face =>
    (∀ i ∈ insert root face, ¬ inBand v q p i) ∧
    ∃ i ∈ face, overlapDet v q p i root ≠ 0
  change colored.card ≤
    ∑ face ∈ active,
      ∑ i ∈ face,
        LRCSevenStalkReuseBudget.rootSpokeMass v q p root i
  calc
    colored.card = ∑ _face ∈ colored, 1 := by simp
    _ ≤ ∑ face ∈ colored,
        ∑ i ∈ face,
          LRCSevenStalkReuseBudget.rootSpokeMass v q p root i := by
      apply Finset.sum_le_sum
      intro face hface
      obtain ⟨_faceMem, _hbad, i, hi, hne⟩ := Finset.mem_filter.mp hface
      have hone : 1 ≤
          LRCSevenStalkReuseBudget.rootSpokeMass v q p root i := by
        unfold LRCSevenStalkReuseBudget.rootSpokeMass
        exact Int.natAbs_pos.mpr hne
      exact hone.trans
        (Finset.single_le_sum (fun j _ => Nat.zero_le _) hi)
    _ ≤ ∑ face ∈ active,
        ∑ i ∈ face,
          LRCSevenStalkReuseBudget.rootSpokeMass v q p root i := by
      apply Finset.sum_le_sum_of_subset_of_nonneg
      · intro face hface
        obtain ⟨hfaceCarrier, hbad, _hcolor⟩ := Finset.mem_filter.mp hface
        exact Finset.mem_filter.mpr ⟨hfaceCarrier, hbad⟩
      · intro face _hactive _hnotColored
        exact Nat.zero_le _

/-- When the root is bad, the active fixed faces are exactly the six-faces of
the bad-lower-vertex set used by the reuse ledger. -/
theorem activeRootedFaces_eq_badLower
    (v : Fin 13 → ℤ) (q p : ℕ) (root : Fin 13)
    (hbadroot : ¬ inBand v q p root) :
    (rootedSixFaces root).filter (fun face =>
      ∀ i ∈ insert root face, ¬ inBand v q p i) =
      (LRCSevenStalkReuseBudget.badLowerVertices
        v q p root).powersetCard 6 := by
  ext face
  constructor
  · intro hface
    obtain ⟨hcarrier, hallBad⟩ := Finset.mem_filter.mp hface
    have hcarrierData := Finset.mem_powersetCard.mp hcarrier
    apply Finset.mem_powersetCard.mpr
    refine ⟨?_, hcarrierData.2⟩
    intro i hi
    unfold LRCSevenStalkReuseBudget.badLowerVertices
    exact Finset.mem_filter.mpr
      ⟨hcarrierData.1 hi, hallBad i (Finset.mem_insert_of_mem hi)⟩
  · intro hface
    have hbadData := Finset.mem_powersetCard.mp hface
    have hcarrier : face ∈ rootedSixFaces root := by
      unfold rootedSixFaces
      apply Finset.mem_powersetCard.mpr
      refine ⟨?_, hbadData.2⟩
      intro i hi
      exact (Finset.mem_filter.mp (hbadData.1 hi)).1
    apply Finset.mem_filter.mpr
    refine ⟨hcarrier, ?_⟩
    intro i hi
    rcases Finset.mem_insert.mp hi with rfl | hiface
    · exact hbadroot
    · exact (Finset.mem_filter.mp (hbadData.1 hiface)).2

/-- If the root is in band, there is no active rooted face attached to it. -/
theorem activeRootedFaces_eq_empty
    (v : Fin 13 → ℤ) (q p : ℕ) (root : Fin 13)
    (hgoodroot : inBand v q p root) :
    (rootedSixFaces root).filter (fun face =>
      ∀ i ∈ insert root face, ¬ inBand v q p i) = ∅ := by
  apply Finset.eq_empty_iff_forall_notMem.mpr
  intro face hface
  have hallBad := (Finset.mem_filter.mp hface).2
  exact (hallBad root (Finset.mem_insert_self root face)) hgoodroot

/-- **Exact reuse composition.**  The active rooted-face spoke mass is the
sum of the imported exact face-reuse factors, with a zero contribution when
the selected root is not bad. -/
theorem activeRootedFaceSpokeMass_eq_reuseTransport
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13) :
    activeRootedFaceSpokeMass v q root =
      ∑ p ∈ Finset.Ioo 0 q,
        if ¬ inBand v q p root then
          ((LRCSevenStalkReuseBudget.badLowerVertices
              v q p root).card - 1).choose 5 *
            ∑ i ∈ LRCSevenStalkReuseBudget.badLowerVertices
                v q p root,
              LRCSevenStalkReuseBudget.rootSpokeMass v q p root i
        else 0 := by
  unfold activeRootedFaceSpokeMass
  apply Finset.sum_congr rfl
  intro p _hp
  by_cases hbadroot : ¬ inBand v q p root
  · rw [if_pos hbadroot, activeRootedFaces_eq_badLower
      v q p root hbadroot]
    exact LRCSevenStalkReuseBudget.badLower_sixFace_spokeMass_eq
      v q p root
  · have hgoodroot : inBand v q p root := not_not.mp hbadroot
    rw [if_neg hbadroot, activeRootedFaces_eq_empty
      v q p root hgoodroot]
    simp

/-- A nonzero root speed makes the gcd of every inserted rooted face strictly
positive.  Thus positivity is not an extra aggregation hypothesis. -/
theorem insertedRootFace_gcd_pos
    (v : Fin 13 → ℤ) (root : Fin 13) (face : Finset (Fin 13))
    (hvroot : v root ≠ 0) :
    0 < (insert root face).gcd v := by
  have hne : (insert root face).gcd v ≠ 0 :=
    Finset.gcd_ne_zero_iff.mpr ⟨root, Finset.mem_insert_self root face, hvroot⟩
  have hnonneg : 0 ≤ (insert root face).gcd v :=
    Int.nonneg_of_normalize_eq_self Finset.normalize_gcd
  omega

/-- On an inserted rooted face the anchor-star certificate is exact, not a
relaxation: the root lies in the stalk and its nonzero speed completes every
missing zero-color edge. -/
theorem anchoredZeroColorRootedFaceActivity_eq_zeroColor
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13) (hvroot : v root ≠ 0) :
    anchoredZeroColorRootedFaceActivity v q root =
      zeroColorRootedFaceActivity v q root := by
  rw [anchoredZeroColorRootedFaceActivity_eq_sum_stalkFibers,
    zeroColorRootedFaceActivity_eq_sum_stalkFibers]
  apply Finset.sum_congr rfl
  intro face _hface
  congr 1
  apply Finset.Subset.antisymm
  · exact anchoredZeroColorStalkMultipliers_subset
      v q (insert root face) root hvroot
  · intro p hp
    unfold zeroColorStalkMultipliers at hp
    unfold anchoredZeroColorStalkMultipliers
    obtain ⟨hpWindow, hbad, hzero⟩ := Finset.mem_filter.mp hp
    exact Finset.mem_filter.mpr
      ⟨hpWindow, hbad, fun i hi =>
        hzero i hi root (Finset.mem_insert_self root face)⟩

/-- **Summed aligned-stalk budget.**  Once the top-window inequality is known
for every rooted face, all complete zero-color multiplier activity is bounded
by the explicit sum of the corresponding gcd-resonance budgets.  Overlapping
stalks are counted with their true multiplicity on both sides. -/
theorem zeroColorRootedFaceActivity_le_gcdBudget
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13)
    (hq : 0 < q) (hvroot : 0 < v root)
    (hwindow : ∀ face ∈ rootedSixFaces root,
      (insert root face).gcd v * (q : ℤ) ≤ 14 * v root) :
    zeroColorRootedFaceActivity v q root ≤
      ∑ face ∈ rootedSixFaces root,
        (Nat.gcd ((insert root face).gcd v).natAbs q - 1) := by
  rw [zeroColorRootedFaceActivity_eq_sum_stalkFibers]
  apply Finset.sum_le_sum
  intro face hface
  exact zeroColorStalkMultipliers_card_le
    v q (insert root face) hq
      (insertedRootFace_gcd_pos v root face (ne_of_gt hvroot))
      root (Finset.mem_insert_self root face) hvroot (hwindow face hface)

/-- **Anchor-star summed aligned budget.**  This is the form that composes
directly with the root-spoke reuse ledger: vanishing root-spoke mass is paid by
the gcd budget, while its complement has a nonzero root-spoke charge. -/
theorem anchoredZeroColorRootedFaceActivity_le_gcdBudget
    (v : Fin 13 → ℤ) (q : ℕ) (root : Fin 13)
    (hq : 0 < q) (hvroot : 0 < v root)
    (hwindow : ∀ face ∈ rootedSixFaces root,
      (insert root face).gcd v * (q : ℤ) ≤ 14 * v root) :
    anchoredZeroColorRootedFaceActivity v q root ≤
      ∑ face ∈ rootedSixFaces root,
        (Nat.gcd ((insert root face).gcd v).natAbs q - 1) := by
  rw [anchoredZeroColorRootedFaceActivity_eq_sum_stalkFibers]
  apply Finset.sum_le_sum
  intro face hface
  exact anchoredZeroColorStalkMultipliers_card_le
    v q (insert root face) hq
      (insertedRootFace_gcd_pos v root face (ne_of_gt hvroot))
      root root (ne_of_gt hvroot) (Finset.mem_insert_self root face)
      hvroot (hwindow face hface)

/-! ## Axiom audit -/

#print axioms zeroColorRootedFaceActivity_eq_sum_stalkFibers
#print axioms anchoredZeroColorRootedFaceActivity_eq_sum_stalkFibers
#print axioms insertedRootFace_gcd_pos
#print axioms anchoredZeroColorRootedFaceActivity_eq_zeroColor
#print axioms zeroColorRootedFaceActivity_le_gcdBudget
#print axioms anchoredZeroColorRootedFaceActivity_le_gcdBudget

end LRCAlignedStalkAggregation
end LonelyRunner
