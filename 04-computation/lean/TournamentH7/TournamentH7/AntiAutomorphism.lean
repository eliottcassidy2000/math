/-
  TournamentH7.AntiAutomorphism — Anti-automorphisms and Hamiltonian
                                   path counts

  ─── What this module provides ─────────────────────────────────────────
  An *anti-automorphism* of a tournament T is a permutation
  `φ : Fin n ≃ Fin n` such that
        T.arc u v = true   ⟺   T.arc (φ v) (φ u) = true.
  Equivalently, φ relabels T into its op (reverse) tournament.

  ─── The abstract anti-palindrome theorem ──────────────────────────────
  *If T admits an anti-automorphism φ, then the number of Hamiltonian
  paths of T starting at v equals the number ending at φ v.*

  This is the core abstract content behind THM-316 (project-novel,
  opus-2026-05-22-S3): the all-0 staircase tournament `T_k` on 2k
  vertices admits the anti-automorphism v ↦ 2k − 1 − v, hence
        epStart(v, T_k)  =  epEnd(2k − 1 − v, T_k).

  ─── Lean status ───────────────────────────────────────────────────────
  We define `IsAntiAutomorphism`, `epStart`, `epEnd`, and state the
  abstract anti-palindrome theorem as an axiom (its proof is a clean
  bijection on the satisfying `Equiv.Perm` finset, deferred to a future
  session).
-/

import TournamentH7.SCC
import TournamentH7.GridReflection
import Mathlib.Data.Fintype.Card
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

open Finset

namespace Tournament

variable {n : ℕ}

/-! ### Anti-automorphisms -/

/-- An *anti-automorphism* of a tournament T is a permutation that
    relabels T into its op tournament. -/
def IsAntiAutomorphism (T : Tournament n) (φ : Equiv.Perm (Fin n)) : Prop :=
  ∀ u v : Fin n, T.arc u v = true ↔ T.arc (φ v) (φ u) = true

/-! ### Endpoint counts -/

/-- `epStart T v` = number of Hamiltonian paths of T starting at `v`.
    A Hamiltonian path on `n ≥ 1` vertices is encoded as a permutation
    `σ : Fin n ≃ Fin n` (the listing of vertices in order). -/
noncomputable def epStart (T : Tournament n) (hn : 0 < n) (v : Fin n) : ℕ := by
  classical
  exact (Finset.univ.filter
    (fun σ : Equiv.Perm (Fin n) =>
      IsHamiltonianPath T σ ∧
      σ ⟨0, hn⟩ = v)).card

/-- `epEnd T v` = number of Hamiltonian paths of T ending at `v`. -/
noncomputable def epEnd (T : Tournament n) (hn : 0 < n) (v : Fin n) : ℕ := by
  classical
  exact (Finset.univ.filter
    (fun σ : Equiv.Perm (Fin n) =>
      IsHamiltonianPath T σ ∧
      σ ⟨n - 1, by omega⟩ = v)).card

/-! ### The abstract anti-palindrome theorem -/

/-- **Axiom (abstract anti-palindrome).**

    If T admits an anti-automorphism φ, then the number of Hamiltonian
    paths starting at v equals the number ending at φ v.

    Proof sketch (informal): the map
        σ : HP(T) starting at v   ↦   φ ∘ σ ∘ reversal : HP(T) ending at φ v
    where reversal(i) = n−1−i.  Using the anti-automorphism property,
    if σ encodes the path v_0 → v_1 → … → v_{n−1}, then
    `φ(v_{n−1}) → φ(v_{n−2}) → … → φ(v_0)` is a HP of T ending at φ(v_0)
    = φ(v). The bijection is an involution by `φ ∘ φ ≈ id` and reversal
    being self-inverse. -/
axiom abstract_anti_palindrome
    (T : Tournament n) (hn : 0 < n) (φ : Equiv.Perm (Fin n))
    (hφ : IsAntiAutomorphism T φ) (v : Fin n) :
    epStart T hn v = epEnd T hn (φ v)

/-! ### Sum form -/

/-- Total: summing `epStart` over all starting vertices recovers H(T). -/
axiom epStart_sum_eq_H (T : Tournament n) (hn : 0 < n) :
    ∑ v : Fin n, epStart T hn v = H T

axiom epEnd_sum_eq_H (T : Tournament n) (hn : 0 < n) :
    ∑ v : Fin n, epEnd T hn v = H T

end Tournament
