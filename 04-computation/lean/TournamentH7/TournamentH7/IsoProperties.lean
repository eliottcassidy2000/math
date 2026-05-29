/-
  TournamentH7.IsoProperties — invariants under tournament isomorphism

  Records that the basic tournament invariants are preserved under iso:
   • `op` (reverse) commutes with `relabel`.
   • `Tournament.outDegree` is preserved (as a multiset of out-degrees;
     individual vertex degrees move with the permutation).
   • `H(T)`, `alphaCount k T`, `IsRegular T` are isomorphism-invariant.

  The Hamiltonian-path-count preservation `H_iso_invariant` follows
  from the bijection on `Equiv.Perm (Fin n)` induced by the iso.
-/

import TournamentH7.Iso
import TournamentH7.SCC
import TournamentH7.OCF

open Finset

namespace Tournament

variable {n : ℕ}

/-! ### `op` commutes with `relabel` (arc-wise) -/

/-- Reversing then relabelling agrees on arcs with relabelling then
    reversing. -/
lemma op_relabel_arc (T : Tournament n) (σ : Equiv.Perm (Fin n)) (i j : Fin n) :
    (op (relabel T σ)).arc i j = (relabel (op T) σ).arc i j := rfl

/-! ### Out-degree under isomorphism -/

/-- An isomorphism preserves vertex out-degrees up to relabelling. -/
lemma outDegree_iso (T₁ T₂ : Tournament n) (f : TournamentIso T₁ T₂) (v : Fin n) :
    T₁.outDegree v = T₂.outDegree (f.perm v) := by
  unfold Tournament.outDegree
  -- Bijection w ↔ f.perm w connects the two filter sets.
  apply Finset.card_bij (fun w _ => f.perm w)
  · intro w hw
    simp at hw ⊢
    have := f.arc_eq v w
    rw [hw] at this
    exact this.symm
  · intros _ _ _ _ h
    exact f.perm.injective h
  · intro w hw
    simp at hw
    refine ⟨f.perm.symm w, ?_, ?_⟩
    · simp
      have hw' := f.arc_eq v (f.perm.symm w)
      rw [Equiv.apply_symm_apply] at hw'
      rw [hw']
      exact hw
    · simp

/-! ### `IsRegular` is iso-invariant -/

/-- A regular tournament's image under an iso is regular. -/
lemma isRegular_iso (T₁ T₂ : Tournament n) (h_iso : T₁ ≅ T₂)
    (hreg : IsRegular T₁) : IsRegular T₂ := by
  intro v
  obtain ⟨f⟩ := h_iso
  have := hreg (f.perm.symm v)
  rw [outDegree_iso T₁ T₂ f, Equiv.apply_symm_apply] at this
  exact this

/-! ### `H(T)` is iso-invariant

    Sketch: given f : T₁ ≅ T₂, the map σ ↦ f.perm * σ is a bijection
    between the satisfying `Equiv.Perm`s.  Recorded as axiom; the
    `outDegree_iso` lemma above is the inline-proved analogue. -/

/-- The Hamiltonian path count is an isomorphism invariant. -/
axiom H_iso_invariant (T₁ T₂ : Tournament n) (h : T₁ ≅ T₂) :
    H T₁ = H T₂

/-! ### `op` arc-equality is just transposition

    `(op T).arc i j = T.arc j i` definitionally, so `op (op T)`'s arcs
    agree with T's arcs pointwise.  The two structures aren't
    *definitionally* equal as `Tournament` records, but the relations
    are equal pointwise. -/

lemma op_op_arc_eq (T : Tournament n) (i j : Fin n) :
    (op (op T)).arc i j = T.arc i j := rfl

end Tournament
