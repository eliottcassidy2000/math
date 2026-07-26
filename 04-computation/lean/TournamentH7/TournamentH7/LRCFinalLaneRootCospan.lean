/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: codex (LRC multi-agent project, 2026-07-26)
-/
import TournamentH7.LRCCommensuration

/-!
# Finite thirteen-root and Haar kernels for THM-2426

This module formalizes three dependency-light pieces of the audited paper proof
`THM-2426-compositional-thirteen-root-final-septimal-lane-exclusion`:

* multiplication by a unit and translation give a permutation of the thirteen
  roots;
* deleting one root from a thirteen-root fibre and imposing a seven-root phase
  condition leaves at least six roots, and the bound is sharp;
* the exact rational mass arithmetic used after root disintegration, together
  with Haar invariance under every nonzero integer dilation of the circle.

The module deliberately does **not** identify these finite sets with the
paper's danger combs, prove that endpoint-excluded fibres have the advertised
root patterns, perform conditional Haar disintegration, or consume the
THM-2391 scalar-cover hypotheses.  It is therefore a formalized kernel of
THM-2426, not a Lean proof of the final-lane exclusion or of LRC(14).
-/

open MeasureTheory MeasureTheory.Measure Set
open scoped ENNReal

namespace LonelyRunner
namespace LRCFinalLaneRootCospan

/-! ## Affine permutations of the thirteen roots -/

/-- Multiplication by a unit followed by translation is an equivalence of
the thirteen residue roots. -/
def affineRootEquiv (u : (ZMod 13)ˣ) (a : ZMod 13) : ZMod 13 ≃ ZMod 13 where
  toFun x := (u : ZMod 13) * x + a
  invFun y := (↑(u⁻¹) : ZMod 13) * (y - a)
  left_inv x := by
    simp
  right_inv y := by
    simp

/-- The affine image of a finite root set. -/
def affineRootSet (u : (ZMod 13)ˣ) (a : ZMod 13)
    (roots : Finset (ZMod 13)) : Finset (ZMod 13) :=
  roots.map (affineRootEquiv u a).toEmbedding

@[simp]
theorem card_affineRootSet (u : (ZMod 13)ˣ) (a : ZMod 13)
    (roots : Finset (ZMod 13)) :
    (affineRootSet u a roots).card = roots.card := by
  simp [affineRootSet]

/-- A fixed six-residue interval model.  Every generic six-root phase band is
an affine image of this model. -/
def innerSix : Finset (ZMod 13) := {0, 1, 2, 3, 4, 5}

/-- The complementary seven-residue phase model. -/
def outerSeven : Finset (ZMod 13) := Finset.univ \ innerSix

theorem innerSix_card : innerSix.card = 6 := by
  decide

theorem outerSeven_card : outerSeven.card = 7 := by
  decide

/-- Any unit-affine copy of the phase-good block has seven roots. -/
def phaseGoodRoots (h : (ZMod 13)ˣ) (a : ZMod 13) : Finset (ZMod 13) :=
  affineRootSet h a outerSeven

@[simp]
theorem phaseGoodRoots_card (h : (ZMod 13)ˣ) (a : ZMod 13) :
    (phaseGoodRoots h a).card = 7 := by
  simp [phaseGoodRoots, outerSeven_card]

/-- A unit-affine one-root phase section, used by the `M=1` branch. -/
def onePhaseRoot (h : (ZMod 13)ˣ) (a : ZMod 13) : Finset (ZMod 13) :=
  affineRootSet h a {0}

@[simp]
theorem onePhaseRoot_card (h : (ZMod 13)ˣ) (a : ZMod 13) :
    (onePhaseRoot h a).card = 1 := by
  simp [onePhaseRoot]

/-- The transition-root model after the unique narrow-danger root has been
deleted from the full thirteen-root fibre. -/
def transitionRoots (u : (ZMod 13)ˣ) (a : ZMod 13) : Finset (ZMod 13) :=
  Finset.univ.erase (affineRootEquiv u a 0)

@[simp]
theorem transitionRoots_card (u : (ZMod 13)ˣ) (a : ZMod 13) :
    (transitionRoots u a).card = 12 := by
  simp [transitionRoots]

/-! ## The finite cospan count -/

/-- Inclusion-exclusion on a thirteen-element universe: a twelve-root set and
a seven-root set have at least six common roots. -/
theorem card_inter_ge_six_of_card_twelve_seven
    {α : Type*} [Fintype α] [DecidableEq α]
    (hcard : Fintype.card α = 13) (transition phaseGood : Finset α)
    (htransition : transition.card = 12) (hphase : phaseGood.card = 7) :
    6 ≤ (transition ∩ phaseGood).card := by
  have hunion : (transition ∪ phaseGood).card ≤ 13 := by
    calc
      (transition ∪ phaseGood).card ≤ (Finset.univ : Finset α).card :=
        Finset.card_le_card (Finset.subset_univ _)
      _ = Fintype.card α := Finset.card_univ
      _ = 13 := hcard
  have hinc := Finset.card_union_add_card_inter transition phaseGood
  omega

/-- The precise finite root-cospan kernel of the `M=2` count in THM-2426:
the unique narrow-danger deletion and any seven-root phase-good block leave
at least six useful roots. -/
theorem affine_thirteen_root_cospan_floor
    (u h : (ZMod 13)ˣ) (narrowShift phaseShift : ZMod 13) :
    6 ≤
      ((transitionRoots u narrowShift) ∩
        phaseGoodRoots h phaseShift).card := by
  apply card_inter_ge_six_of_card_twelve_seven
      (α := ZMod 13) (by norm_num)
  · exact transitionRoots_card u narrowShift
  · exact phaseGoodRoots_card h phaseShift

/-- Six is the best possible conclusion from the two marginal root counts
alone.  The witness deletes root `6` from the fixed outer-seven block. -/
theorem thirteen_root_cospan_floor_sharp :
    ∃ transition phaseGood : Finset (ZMod 13),
      transition.card = 12 ∧
      phaseGood.card = 7 ∧
      (transition ∩ phaseGood).card = 6 := by
  refine ⟨Finset.univ.erase 6, outerSeven, ?_, outerSeven_card, ?_⟩
  · decide
  · decide

/-! ## Exact mass arithmetic -/

/-- The different-septimal-depth safe mass followed by the six-of-thirteen
root floor. -/
theorem m2_section_mass_rat :
    (6 : ℚ) / 13 * (6 / 49) = 36 / 637 := by
  norm_num

/-- The final one-of-seven phase landing in the `M=2` branch. -/
theorem m2_landed_mass_rat :
    ((6 : ℚ) / 13 * (6 / 49)) / 7 = 36 / 4459 := by
  norm_num

/-- The exact positive remainder in the `M=1` branch after the sharp
narrow-comb half-overlap loss. -/
theorem m1_residual_mass_rat :
    (6 : ℚ) / 637 - 1 / 182 = 5 / 1274 := by
  norm_num

theorem m1_residual_mass_rat_pos :
    0 < (6 : ℚ) / 637 - 1 / 182 := by
  norm_num

/-- Real-valued form suited to a future `Measure.toReal` consumer. -/
theorem m1_residual_mass_real :
    (6 : ℝ) / 637 - 1 / 182 = 5 / 1274 ∧
      0 < (6 : ℝ) / 637 - 1 / 182 := by
  constructor <;> norm_num

/-! ## Haar pullback under integer dilation -/

/-- Every nonzero integer dilation of the unit circle preserves normalized
Haar measure.  This is the formal pullback used when the paper scales out a
common integer factor. -/
theorem volume_preimage_integerDilation {k : ℤ} (hk : k ≠ 0)
    (s : Set UnitAddCircle) (hs : NullMeasurableSet s volume) :
    volume ((fun x : UnitAddCircle => k • x) ⁻¹' s) = volume s :=
  (measurePreserving_zsmul volume hk).measure_preimage hs

/-- The concrete seven-dilation pullback used in the `M=2` branch. -/
theorem volume_preimage_sevenDilation
    (s : Set UnitAddCircle) (hs : NullMeasurableSet s volume) :
    volume ((fun x : UnitAddCircle => (7 : ℤ) • x) ⁻¹' s) = volume s := by
  exact volume_preimage_integerDilation (by norm_num) s hs

#print axioms affine_thirteen_root_cospan_floor
#print axioms thirteen_root_cospan_floor_sharp
#print axioms m2_landed_mass_rat
#print axioms m1_residual_mass_real
#print axioms volume_preimage_sevenDilation

end LRCFinalLaneRootCospan
end LonelyRunner
