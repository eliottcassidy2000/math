/-
  TournamentH7.RootSigns — type-A root-sign atoms

  This module records a small project-axiom-free bridge from tournament arcs
  to the positive roots of type A.  An arc i→j is represented by the signed root
  e_i - e_j in the integer root lattice `Fin n → ℤ`.

  The main tiny theorem is `root_cycle_sum`: every directed triangle carries
  the minimal type-A root relation

      (e_i - e_j) + (e_j - e_k) + (e_k - e_i) = 0.

  The slightly more general theorem is `walkRootSum_append_single`: the root
  sum along any vertex walk telescopes to the root from its first vertex to
  its last vertex.  Closed walks therefore carry zero total root.

  This is the formal version of the representation-lens slogan:
  "a 3-cycle is the first root relation."
-/

import Mathlib.Algebra.Group.Pi.Basic
import Mathlib.Data.Fin.Basic
import Mathlib.Data.Int.Basic

namespace TypeA

/-- The integer root lattice for type `A_{n-1}`, represented in the ambient
    coordinate lattice. -/
abbrev RootSpace (n : ℕ) := Fin n → ℤ

/-- The signed type-A root `e_i - e_j`.  The definition works uniformly even
    when `i = j`, in which case the root is zero. -/
def root {n : ℕ} (i j : Fin n) : RootSpace n :=
  fun k => (if k = i then (1 : ℤ) else 0) - (if k = j then (1 : ℤ) else 0)

/-- The degenerate root `e_i - e_i` is zero. -/
@[simp] theorem root_self {n : ℕ} (i : Fin n) : root i i = 0 := by
  funext k
  simp [root]

/-- Reversing an oriented root negates it. -/
theorem root_swap {n : ℕ} (i j : Fin n) : root j i = -root i j := by
  funext k
  simp [root]
  omega

/-- A two-edge backtrack has zero total root. -/
@[simp] theorem root_backtrack_sum {n : ℕ} (i j : Fin n) :
    root i j + root j i = 0 := by
  funext k
  simp [root]
  omega

/-! ### Telescoping calculus -/

/-- Consecutive type-A roots telescope. -/
@[simp] theorem root_add_root {n : ℕ} (i j k : Fin n) :
    root i j + root j k = root i k := by
  funext x
  simp [root]
  omega

/-- The only zero signed root is a degenerate root. -/
@[simp] theorem root_eq_zero_iff {n : ℕ} (i j : Fin n) :
    root i j = 0 ↔ i = j := by
  constructor
  · intro h
    by_contra hne
    have hcoord := congr_fun h i
    simp [root, hne] at hcoord
  · intro h
    subst j
    simp

/-- The root contributions around any directed triangle telescope to zero.

    This is the type-A root-lattice form of the directed 3-cycle relation:
    `i → j → k → i` carries `(e_i-e_j)+(e_j-e_k)+(e_k-e_i)=0`. -/
@[simp] theorem root_cycle_sum {n : ℕ} (i j k : Fin n) :
    root i j + root j k + root k i = 0 := by
  rw [root_add_root]
  exact root_backtrack_sum i k

/-- Root sum over consecutive edges in a finite vertex walk. -/
def walkRootSum {n : ℕ} : List (Fin n) -> RootSpace n
  | [] => 0
  | [_] => 0
  | a :: b :: rest => root a b + walkRootSum (b :: rest)

@[simp] theorem walkRootSum_nil {n : ℕ} :
    walkRootSum ([] : List (Fin n)) = 0 := rfl

@[simp] theorem walkRootSum_single {n : ℕ} (a : Fin n) :
    walkRootSum [a] = 0 := rfl

@[simp] theorem walkRootSum_cons_cons {n : ℕ} (a b : Fin n) (rest : List (Fin n)) :
    walkRootSum (a :: b :: rest) = root a b + walkRootSum (b :: rest) := rfl

/-- The root sum of any finite walk telescopes to the root from first vertex
    to last vertex. -/
@[simp] theorem walkRootSum_append_single {n : ℕ} (a b : Fin n) :
    ∀ middle : List (Fin n), walkRootSum (a :: middle ++ [b]) = root a b := by
  intro middle
  induction middle generalizing a with
  | nil =>
      funext x
      simp [walkRootSum]
  | cons c rest ih =>
      change root a c + walkRootSum (c :: rest ++ [b]) = root a b
      rw [ih c]
      exact root_add_root a c b

/-- Closed walks have zero total type-A root. -/
@[simp] theorem walkRootSum_closed {n : ℕ} (a : Fin n) (middle : List (Fin n)) :
    walkRootSum (a :: middle ++ [a]) = 0 := by
  simpa using walkRootSum_append_single a a middle

end TypeA
