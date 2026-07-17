/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: codex-c9-lean (LRC multi-agent project, 2026-07-17)
-/
import Mathlib

/-!
# The scale-nine terminal owner-nerve obstruction

This file formalizes only the final combinatorial implication in THM-969,
equations (17)--(20):

* in the all-order-nine bank, the pair-intersection nerve of the six owner
  obligations is the three-edge matching `3K₂` (cycle-antipodal pairs);
* in the mixed bank, the nerve induced by the four order-nine owners is the
  two-edge matching `2K₂` (the two sign pairs).

In either case one nonedge is a disjoint pair of obligations, so the total
owner intersection is empty.  We package this implication for arbitrary
obligation predicates and also verify the two finite graph tables with
ordinary kernel `decide`.

The graph tables are faithful encodings of the *reported pair nerves*, not of
the literal CRT mask bank.  In particular, this file does **not** formalize or
recheck THM-969's preceding reduction of 482,294,736 raw contexts to 76
owner-local contexts, its mask cardinalities, or its pair-intersection counts.
There is no `sorry`, no `native_decide`, and no added axiom.
-/

namespace LonelyRunner
namespace ScaleNineOwnerNerve

/-- A pair-intersection nerve says exactly which two distinct obligations
have a common witness.  Loops are deliberately ignored. -/
def ExactPairNerve {I α : Type*}
    (obligation : I → α → Prop) (edge : I → I → Bool) : Prop :=
  ∀ i j, i ≠ j →
    ((∃ x, obligation i x ∧ obligation j x) ↔ edge i j = true)

/-- The exact pair nerve induced on a selected subfamily of obligations. -/
def ExactInducedPairNerve {I J α : Type*}
    (obligation : I → α → Prop) (select : J → I)
    (edge : J → J → Bool) : Prop :=
  ExactPairNerve (fun j x ↦ obligation (select j) x) edge

/-! ## The all-order-nine `3K₂` nerve -/

/-- The six vertices are cyclic owner positions. -/
abbrev AllD9Owner := Fin 6

/-- Equation (18): two distinct owner vertices meet exactly when they are
antipodal on the signed six-cycle. -/
def threeKTwoEdge (i j : AllD9Owner) : Bool :=
  ((i.val + 3) % 6 == j.val)

/-- The 36-entry adjacency table is symmetric. -/
theorem threeKTwoEdge_symmetric :
    ∀ i j, threeKTwoEdge i j = threeKTwoEdge j i := by
  decide

/-- The table has no loops. -/
theorem threeKTwoEdge_irreflexive :
    ∀ i, threeKTwoEdge i i = false := by
  decide

/-- Every vertex has exactly one neighbor, as required for `3K₂`. -/
theorem threeKTwoEdge_degree_one :
    ∀ i,
      (Finset.univ.filter (fun j : AllD9Owner ↦
        threeKTwoEdge i j = true)).card = 1 := by
  decide

/-- There are exactly three unordered edges in the all-order-nine table. -/
theorem threeKTwoEdge_unordered_count :
    (Finset.univ.filter (fun p : AllD9Owner × AllD9Owner ↦
      p.1.val < p.2.val ∧ threeKTwoEdge p.1 p.2 = true)).card = 3 := by
  decide

/-- The matching has no triangle. -/
theorem threeKTwoEdge_triangle_free :
    ∀ i j k : AllD9Owner,
      i ≠ j → i ≠ k → j ≠ k →
      ¬ (threeKTwoEdge i j = true ∧
         threeKTwoEdge j k = true ∧
         threeKTwoEdge k i = true) := by
  decide

/-- An exact `3K₂` pair nerve cannot have a witness in three distinct owner
obligations. -/
theorem threeKTwo_no_threefold_intersection { α : Type* }
    (obligation : AllD9Owner → α → Prop)
    (hnerve : ExactPairNerve obligation threeKTwoEdge)
    (i j k : AllD9Owner) (hij : i ≠ j) (hik : i ≠ k) (hjk : j ≠ k) :
    ¬ ∃ x, obligation i x ∧ obligation j x ∧ obligation k x := by
  rintro ⟨x, hi, hj, hk⟩
  have heij : threeKTwoEdge i j = true :=
    (hnerve i j hij).mp ⟨x, hi, hj⟩
  have hejk : threeKTwoEdge j k = true :=
    (hnerve j k hjk).mp ⟨x, hj, hk⟩
  have heki : threeKTwoEdge k i = true :=
    (hnerve k i (Ne.symm hik)).mp ⟨x, hk, hi⟩
  exact threeKTwoEdge_triangle_free i j k hij hik hjk
    ⟨heij, hejk, heki⟩

/-- Therefore the intersection of all six all-order-nine owner obligations is
empty. -/
theorem threeKTwo_no_total_owner_intersection { α : Type* }
    (obligation : AllD9Owner → α → Prop)
    (hnerve : ExactPairNerve obligation threeKTwoEdge) :
    ¬ ∃ x, ∀ i, obligation i x := by
  rintro ⟨x, hx⟩
  exact threeKTwo_no_threefold_intersection obligation hnerve
    0 1 2 (by decide) (by decide) (by decide)
    ⟨x, hx 0, hx 1, hx 2⟩

/-! ## The mixed-bank `2K₂` induced sub-nerve -/

/-- The four order-nine owners are indexed in the order
`(2a, -2a, 3a, -3a)`. -/
abbrev MixedD9Owner := Fin 4

/-- Equation (20): the two sign-pairs meet internally, and every cross-pair
is disjoint. -/
def twoKTwoEdge (i j : MixedD9Owner) : Bool :=
  (i.val / 2 == j.val / 2) && (i.val != j.val)

/-- The 16-entry adjacency table is symmetric. -/
theorem twoKTwoEdge_symmetric :
    ∀ i j, twoKTwoEdge i j = twoKTwoEdge j i := by
  decide

/-- The table has no loops. -/
theorem twoKTwoEdge_irreflexive :
    ∀ i, twoKTwoEdge i i = false := by
  decide

/-- Every vertex has exactly one neighbor, as required for `2K₂`. -/
theorem twoKTwoEdge_degree_one :
    ∀ i,
      (Finset.univ.filter (fun j : MixedD9Owner ↦
        twoKTwoEdge i j = true)).card = 1 := by
  decide

/-- There are exactly two unordered edges in the mixed order-nine table. -/
theorem twoKTwoEdge_unordered_count :
    (Finset.univ.filter (fun p : MixedD9Owner × MixedD9Owner ↦
      p.1.val < p.2.val ∧ twoKTwoEdge p.1 p.2 = true)).card = 2 := by
  decide

/-- The four-element table is exactly the two sign-pairs
`{0,1}` and `{2,3}` (with both orientations). -/
theorem twoKTwoEdge_iff_sign_pair :
    ∀ i j : MixedD9Owner,
      twoKTwoEdge i j = true ↔
        (i = 0 ∧ j = 1) ∨ (i = 1 ∧ j = 0) ∨
        (i = 2 ∧ j = 3) ∨ (i = 3 ∧ j = 2) := by
  decide

/-- The cross-pair represented by vertices zero and two is a nonedge. -/
theorem twoKTwoEdge_zero_two : twoKTwoEdge 0 2 = false := by
  decide

/-- An exact `2K₂` pair nerve on four obligations has empty fourfold
intersection. -/
theorem twoKTwo_no_fourfold_intersection { α : Type* }
    (obligation : MixedD9Owner → α → Prop)
    (hnerve : ExactPairNerve obligation twoKTwoEdge) :
    ¬ ∃ x, ∀ i, obligation i x := by
  rintro ⟨x, hx⟩
  have hedge : twoKTwoEdge 0 2 = true :=
    (hnerve 0 2 (by decide)).mp ⟨x, hx 0, hx 2⟩
  simp [twoKTwoEdge_zero_two] at hedge

/-- If four selected obligations among a larger owner family have the mixed
`2K₂` induced nerve, then the total intersection of the larger family is
empty.  Injectivity records that the selection really consists of four
distinct owners, although exactness of the reported nerve already supplies
the contradiction used in the proof. -/
theorem twoKTwo_subnerve_no_total_owner_intersection
    {I α : Type*}
    (obligation : I → α → Prop) (select : MixedD9Owner → I)
    (_hselect : Function.Injective select)
    (hnerve : ExactInducedPairNerve obligation select twoKTwoEdge) :
    ¬ ∃ x, ∀ i, obligation i x := by
  rintro ⟨x, hx⟩
  exact twoKTwo_no_fourfold_intersection
    (fun j y ↦ obligation (select j) y) hnerve
    ⟨x, fun j ↦ hx (select j)⟩

#print axioms threeKTwoEdge_unordered_count
#print axioms threeKTwo_no_threefold_intersection
#print axioms threeKTwo_no_total_owner_intersection
#print axioms twoKTwoEdge_unordered_count
#print axioms twoKTwo_no_fourfold_intersection
#print axioms twoKTwo_subnerve_no_total_owner_intersection

end ScaleNineOwnerNerve
end LonelyRunner
