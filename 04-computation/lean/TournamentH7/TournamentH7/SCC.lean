/-
  TournamentH7.SCC — Strongly connected components and Hamilton paths

  Defines:
   • `Tournament.Reaches T u v`     : there is a directed walk u → v in T
                                       (reflexive-transitive closure of beats).
   • `Tournament.IsSCC T S`         : `S ⊆ Fin n` is a strongly connected
                                       component of T.
   • `Tournament.IsHamiltonianPath` : an ordering of all n vertices with all
                                       consecutive arcs present in T.
   • `Tournament.H`                 : the Hamiltonian path count H(T)
                                       = #{σ : Perm (Fin n) | …}.
-/

import TournamentH7.Basic
import Mathlib.Data.Finset.Card
import Mathlib.Data.Fintype.Pi
import Mathlib.Data.Fintype.Perm

namespace Tournament

variable {n : ℕ}

/-- Reachability: `Reaches T u v` says there's a directed walk u → v in T. -/
inductive Reaches (T : Tournament n) : Fin n → Fin n → Prop
  | refl  (v : Fin n) : Reaches T v v
  | step  {u v w : Fin n} : T.arc u v = true → Reaches T v w → Reaches T u w

/-- A subset `S ⊆ Fin n` is a strongly-connected component of T if every
    pair of vertices in S is mutually reachable, and S is maximal with this
    property. -/
def IsSCC (T : Tournament n) (S : Finset (Fin n)) : Prop :=
  S.Nonempty ∧
  (∀ u ∈ S, ∀ v ∈ S, Reaches T u v) ∧
  (∀ T' : Finset (Fin n), S ⊆ T' →
    (∀ u ∈ T', ∀ v ∈ T', Reaches T u v) → T' = S)

/-- A Hamiltonian path in T: a bijection `σ : Fin n ≃ Fin n` such that
    every consecutive pair `(σ i, σ (i+1))` is an arc of T. -/
def IsHamiltonianPath (T : Tournament n) (σ : Equiv.Perm (Fin n)) : Prop :=
  ∀ (i : Fin n) (h : i.val + 1 < n),
    T.arc (σ i) (σ ⟨i.val + 1, h⟩) = true

/-- The set of Hamiltonian paths in T as a `Finset`. -/
def hamiltonianPaths (T : Tournament n) [DecidableEq (Fin n)] :
    Finset (Equiv.Perm (Fin n)) :=
  Finset.univ.filter (fun σ => ∀ (i : Fin n) (h : i.val + 1 < n),
    T.arc (σ i) (σ ⟨i.val + 1, h⟩) = true)

/-- `H(T)` = Hamiltonian path count. Defined classically (noncomputable) via
    `Fintype.card` of the satisfying subtype, because membership in
    `IsHamiltonianPath` is propositional. -/
noncomputable def H (T : Tournament n) : ℕ := by
  classical
  exact (Finset.univ.filter
    (fun σ : Equiv.Perm (Fin n) =>
      ∀ (i : Fin n) (h : i.val + 1 < n),
        T.arc (σ i) (σ ⟨i.val + 1, h⟩) = true)).card

end Tournament
