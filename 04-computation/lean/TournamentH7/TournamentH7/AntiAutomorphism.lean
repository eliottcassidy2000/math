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
  We define `IsAntiAutomorphism`, `epStart`, `epEnd`, and prove the
  abstract anti-palindrome theorem by the explicit bijection
  `σ ↦ φ * σ * vertexReversal n`.
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

private lemma vertexReversal_succ_val (i : Fin n) (h : i.val + 1 < n) :
    ((vertexReversal n ⟨i.val + 1, h⟩).val + 1) =
      (vertexReversal n i).val := by
  simp [vertexReversal_apply]
  omega

private lemma vertexReversal_succ_fin (i : Fin n) (h : i.val + 1 < n) :
    (⟨(vertexReversal n ⟨i.val + 1, h⟩).val + 1,
        by
          rw [vertexReversal_succ_val i h]
          exact (vertexReversal n i).is_lt⟩ : Fin n) =
      vertexReversal n i := by
  apply Fin.ext
  exact vertexReversal_succ_val i h

private lemma vertexReversal_zero (hn : 0 < n) :
    vertexReversal n ⟨0, hn⟩ = ⟨n - 1, by omega⟩ := by
  apply Fin.ext
  simp [vertexReversal_apply]

private lemma vertexReversal_last (hn : 0 < n) :
    vertexReversal n ⟨n - 1, by omega⟩ = ⟨0, hn⟩ := by
  apply Fin.ext
  simp [vertexReversal_apply]

/-- **Theorem (abstract anti-palindrome).**

    If T admits an anti-automorphism φ, then the number of Hamiltonian
    paths starting at v equals the number ending at φ v.

    Proof: the map
        σ : HP(T) starting at v   ↦   φ ∘ σ ∘ reversal : HP(T) ending at φ v
    where reversal(i) = n−1−i.  Using the anti-automorphism property,
    if σ encodes the path v_0 → v_1 → … → v_{n−1}, then
    `φ(v_{n−1}) → φ(v_{n−2}) → … → φ(v_0)` is a HP of T ending at φ(v_0)
    = φ(v). The inverse map is `τ ↦ φ.symm * τ * reversal`. -/
theorem abstract_anti_palindrome
    (T : Tournament n) (hn : 0 < n) (φ : Equiv.Perm (Fin n))
    (hφ : IsAntiAutomorphism T φ) (v : Fin n) :
    epStart T hn v = epEnd T hn (φ v) := by
  classical
  unfold epStart epEnd
  apply Finset.card_bij
    (fun σ (_ : σ ∈ _) => φ * σ * vertexReversal n)
  · intro σ hσ
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hσ ⊢
    rcases hσ with ⟨hpath, hstart⟩
    constructor
    · intro i hi
      simp only [Equiv.Perm.mul_apply]
      have hprev := hpath (vertexReversal n ⟨i.val + 1, hi⟩)
        (by
          rw [vertexReversal_succ_val i hi]
          exact (vertexReversal n i).is_lt)
      rw [vertexReversal_succ_fin i hi] at hprev
      exact (hφ (σ (vertexReversal n ⟨i.val + 1, hi⟩))
        (σ (vertexReversal n i))).mp hprev
    · simp only [Equiv.Perm.mul_apply]
      rw [vertexReversal_last hn, hstart]
  · intro σ hσ τ hτ h_eq
    ext x
    have hx := congrArg
      (fun ρ : Equiv.Perm (Fin n) => φ.symm (ρ (vertexReversal n x))) h_eq
    exact congrArg Fin.val (by
      simpa [Equiv.Perm.mul_apply, vertexReversal_apply_twice] using hx)
  · intro τ hτ
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hτ
    rcases hτ with ⟨hpath, hend⟩
    refine ⟨φ.symm * τ * vertexReversal n, ?_, ?_⟩
    · simp only [Finset.mem_filter, Finset.mem_univ, true_and]
      constructor
      · intro i hi
        simp only [Equiv.Perm.mul_apply]
        have hprev := hpath (vertexReversal n ⟨i.val + 1, hi⟩)
          (by
            rw [vertexReversal_succ_val i hi]
            exact (vertexReversal n i).is_lt)
        rw [vertexReversal_succ_fin i hi] at hprev
        have hprev' :
            T.arc (φ (φ.symm (τ (vertexReversal n ⟨i.val + 1, hi⟩))))
              (φ (φ.symm (τ (vertexReversal n i)))) = true := by
          simpa using hprev
        exact (hφ (φ.symm (τ (vertexReversal n i)))
          (φ.symm (τ (vertexReversal n ⟨i.val + 1, hi⟩)))).mpr hprev'
      · simp only [Equiv.Perm.mul_apply]
        rw [vertexReversal_zero hn, hend]
        simp
    · ext x
      simp [Equiv.Perm.mul_apply, vertexReversal_apply_twice]

/-! ### Sum form -/

/-- Total: summing `epStart` over all starting vertices recovers H(T). -/
theorem epStart_sum_eq_H (T : Tournament n) (hn : 0 < n) :
    ∑ v : Fin n, epStart T hn v = H T := by
  classical
  unfold epStart H IsHamiltonianPath
  symm
  rw [Finset.card_eq_sum_card_fiberwise
    (s := Finset.univ.filter
      (fun σ : Equiv.Perm (Fin n) =>
        ∀ (i : Fin n) (h : i.val + 1 < n),
          T.arc (σ i) (σ ⟨i.val + 1, h⟩) = true))
    (t := (Finset.univ : Finset (Fin n)))
    (f := fun σ : Equiv.Perm (Fin n) => σ ⟨0, hn⟩)]
  · apply Finset.sum_congr rfl
    intro v _hv
    congr 1
    ext σ
    simp [and_comm]
  · intro σ hσ
    simp

theorem epEnd_sum_eq_H (T : Tournament n) (hn : 0 < n) :
    ∑ v : Fin n, epEnd T hn v = H T := by
  classical
  unfold epEnd H IsHamiltonianPath
  symm
  rw [Finset.card_eq_sum_card_fiberwise
    (s := Finset.univ.filter
      (fun σ : Equiv.Perm (Fin n) =>
        ∀ (i : Fin n) (h : i.val + 1 < n),
          T.arc (σ i) (σ ⟨i.val + 1, h⟩) = true))
    (t := (Finset.univ : Finset (Fin n)))
    (f := fun σ : Equiv.Perm (Fin n) => σ ⟨n - 1, by omega⟩)]
  · apply Finset.sum_congr rfl
    intro v _hv
    congr 1
    ext σ
    simp [and_comm]
  · intro σ hσ
    simp

end Tournament
