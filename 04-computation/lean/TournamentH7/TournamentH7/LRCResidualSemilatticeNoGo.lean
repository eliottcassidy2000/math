import Mathlib.Algebra.Group.Basic
import Mathlib.Data.Finset.Basic

/-!
# The residual semilattice cannot retain an ordered central clutch

Literal deletion by danger labels depends only on the finite set of labels
already used.  Its universal carrier is therefore `Finset Label` with union:
commutative, associative, and idempotent.

This module records three losses that are easy to obscure in a geometric
presentation:

* swapping two deletion labels changes no residual state;
* repeating a deletion label changes no residual state, even after any later
  unmarked continuation;
* a union-multiplicative invariant with values in a group is necessarily
  trivial.

Consequently an invariant which distinguishes two ordered quaternionic lifts
cannot descend through the unmarked union/residual carrier.  No assertion
about the existence of a marked physical lift is made here.
-/

namespace LonelyRunner
namespace ResidualSemilattice

variable {Label State GroupValue LiftValue : Type*}

/-- An unmarked residual statistic cannot see the order of two inserted
labels because both words have the same finite support. -/
theorem insert_order_blind [DecidableEq Label]
    (residual : Finset Label → State) (base : Finset Label)
    (x y : Label) :
    residual (insert y (insert x base)) =
      residual (insert x (insert y base)) := by
  rw [Finset.insert_comm]

/-- Repeating a deletion label is invisible to every statistic which factors
through finite support. -/
theorem insert_idempotent [DecidableEq Label]
    (residual : Finset Label → State) (base : Finset Label)
    (x : Label) :
    residual (insert x (insert x base)) = residual (insert x base) := by
  simp

/-- Appending any further unmarked set of labels cannot recover the order
which finite support has already forgotten. -/
theorem continuation_order_blind [DecidableEq Label]
    (residual : Finset Label → State) (base continuation : Finset Label)
    (x y : Label) :
    residual (continuation ∪ insert y (insert x base)) =
      residual (continuation ∪ insert x (insert y base)) := by
  rw [Finset.insert_comm]

/-- No function on the unmarked residual carrier can assign distinct values
to the two orders of the same pair of labels. -/
theorem no_order_sensitive_lift [DecidableEq Label]
    (lift : Finset Label → LiftValue) (base : Finset Label)
    (x y : Label) (xy yx : LiftValue) (hne : xy ≠ yx) :
    ¬ (lift (insert y (insert x base)) = xy ∧
        lift (insert x (insert y base)) = yx) := by
  intro h
  apply hne
  calc
    xy = lift (insert y (insert x base)) := h.1.symm
    _ = lift (insert x (insert y base)) := by rw [Finset.insert_comm]
    _ = yx := h.2

/-- Every union-multiplicative group-valued invariant of finite support is
trivial.  The load-bearing point is idempotence:
`value(S) = value(S ∪ S) = value(S)^2`. -/
theorem union_multiplicative_group_invariant_trivial
    [DecidableEq Label] [Group GroupValue]
    (value : Finset Label → GroupValue)
    (map_union :
      ∀ left right, value (left ∪ right) = value left * value right) :
    ∀ support, value support = 1 := by
  intro support
  have hidempotent :
      value support = value support * value support := by
    simpa using map_union support support
  have hcancel : (1 : GroupValue) = value support := by
    calc
      1 = (value support)⁻¹ * value support := by simp
      _ = (value support)⁻¹ * (value support * value support) := by
        rw [← hidempotent]
      _ = value support := by simp
  exact hcancel.symm

end ResidualSemilattice
end LonelyRunner

#print axioms LonelyRunner.ResidualSemilattice.insert_order_blind
#print axioms LonelyRunner.ResidualSemilattice.insert_idempotent
#print axioms LonelyRunner.ResidualSemilattice.continuation_order_blind
#print axioms LonelyRunner.ResidualSemilattice.no_order_sensitive_lift
#print axioms LonelyRunner.ResidualSemilattice.union_multiplicative_group_invariant_trivial
