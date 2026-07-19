/-
  LRCA12Chipwalk.lean -- arithmetic kernel for THM-1143.

  The geometric theorem identifies a shallow full-residue packet with grouped
  edge slides on thirteen sheets.  This module formalizes the reusable finite
  state layer: an edge slide is an A_12 root, roots have zero total mass,
  simultaneous groups preserve mass and are order-independent, and every
  prefix starting from the eleven-chip star stays on the affine mass-eleven
  hyperplane.  The analytic floor/danger-edge equivalence and the finite
  h <= 2 certificate remain external to this module.
-/
import Mathlib

open scoped BigOperators

namespace LRC14.A12Chipwalk

abbrev Sheet := Fin 13
abbrev State := Sheet → ℤ

/-- Unit mass at one sheet. -/
def delta (a : Sheet) : State := fun j => if j = a then 1 else 0

/-- Total chip mass. -/
def total (e : State) : ℤ := ∑ j, e j

/-- The root that transports one chip from `source` to `target`. -/
def root (source target : Sheet) : State := fun j => delta target j - delta source j

/-- Incidence vector of a two-endpoint sheet edge. -/
def edgeIncidence (left right : Sheet) : State := fun j => delta left j + delta right j

/-- Apply one root transport. -/
def applyRoot (e : State) (source target : Sheet) : State :=
  fun j => e j + root source target j

/-- Apply one simultaneous wall group.  Prefix safety is tested only after
    this whole fold, never at an artificial ordering inside the tie. -/
def applyGroup (e : State) (moves : List (Sheet × Sheet)) : State :=
  moves.foldl (fun state move => applyRoot state move.1 move.2) e

/-- Apply a chronological list of simultaneous wall groups. -/
def runFrom (e : State) (groups : List (List (Sheet × Sheet))) : State :=
  groups.foldl applyGroup e

/-- The eleven-chip star at the initial sheet. -/
def initial (rootSheet : Sheet) : State :=
  fun j => if j = rootSheet then 11 else 0

/-- The exact grouped-prefix ballot predicate. -/
def PrefixSafe (rootSheet : Sheet) (groups : List (List (Sheet × Sheet))) : Prop :=
  ∀ n, n ≤ groups.length →
    ∀ j, 0 ≤ runFrom (initial rootSheet) (groups.take n) j

theorem total_delta (a : Sheet) : total (delta a) = 1 := by
  simp [total, delta]

theorem total_root (source target : Sheet) : total (root source target) = 0 := by
  change (∑ j, (delta target j - delta source j)) = 0
  rw [Finset.sum_sub_distrib]
  change total (delta target) - total (delta source) = 0
  rw [total_delta, total_delta]
  norm_num

/-- Sliding `{source,pivot}` to `{pivot,target}` cancels the persistent
    endpoint and leaves exactly the root `delta_target-delta_source`. -/
theorem edge_slide_eq_root (source pivot target : Sheet) :
    (fun j => edgeIncidence pivot target j - edgeIncidence source pivot j) =
      root source target := by
  funext j
  simp [edgeIncidence, root]
  ring

theorem total_applyRoot (e : State) (source target : Sheet) :
    total (applyRoot e source target) = total e := by
  change (∑ j, (e j + root source target j)) = ∑ j, e j
  rw [Finset.sum_add_distrib]
  change total e + total (root source target) = total e
  rw [total_root, add_zero]

/-- The result of a tied pair does not depend on its artificial order. -/
theorem applyRoot_comm (e : State) (a b c d : Sheet) :
    applyRoot (applyRoot e a b) c d = applyRoot (applyRoot e c d) a b := by
  funext j
  simp [applyRoot]
  ring

theorem total_applyGroup (e : State) (moves : List (Sheet × Sheet)) :
    total (applyGroup e moves) = total e := by
  induction moves generalizing e with
  | nil => rfl
  | cons move moves ih =>
      change total (applyGroup (applyRoot e move.1 move.2) moves) = total e
      calc
        _ = total (applyRoot e move.1 move.2) := ih _
        _ = total e := total_applyRoot e move.1 move.2

theorem total_runFrom (e : State) (groups : List (List (Sheet × Sheet))) :
    total (runFrom e groups) = total e := by
  induction groups generalizing e with
  | nil => rfl
  | cons group groups ih =>
      simpa [runFrom, total_applyGroup] using ih (applyGroup e group)

theorem total_initial (rootSheet : Sheet) : total (initial rootSheet) = 11 := by
  simp [total, initial]

/-- Every grouped prefix, safe or torn, lies on the affine mass-eleven
    hyperplane.  Prefix safety adds coordinatewise nonnegativity. -/
theorem prefix_total_eleven (rootSheet : Sheet)
    (groups : List (List (Sheet × Sheet))) (n : ℕ) :
    total (runFrom (initial rootSheet) (groups.take n)) = 11 := by
  rw [total_runFrom, total_initial]

theorem prefixSafe_initial (rootSheet : Sheet)
    (groups : List (List (Sheet × Sheet)))
    (h : PrefixSafe rootSheet groups) : ∀ j, 0 ≤ initial rootSheet j := by
  simpa [PrefixSafe, runFrom] using h 0 (Nat.zero_le _)

#print axioms edge_slide_eq_root
#print axioms total_applyGroup
#print axioms prefix_total_eleven

end LRC14.A12Chipwalk
