/-
  TournamentH7.Iso — Tournament isomorphisms

  A *tournament isomorphism* between `T₁ T₂ : Tournament n` is a bijection
  σ : Fin n ≃ Fin n that carries the arc relation of T₁ to that of T₂.

  This module provides the basic infrastructure:
   • `TournamentIso T₁ T₂` — the type of isomorphisms.
   • Refl, symm, trans (groupoid structure).
   • Connections to `op T` (the reverse tournament).
   • The predicate `T₁ ≅ T₂`.

  Used by `SelfComplementary.lean` (now: an SC tournament is one with
  `T ≅ op T`) and by future modules formalising the iso-class graph G_n.
-/

import TournamentH7.Basic
import TournamentH7.GridReflection

namespace Tournament

variable {n : ℕ}

/-! ### Isomorphisms between tournaments -/

/-- A tournament isomorphism: a permutation that intertwines two
    tournament structures. -/
structure TournamentIso (T₁ T₂ : Tournament n) where
  /-- The underlying permutation. -/
  perm    : Equiv.Perm (Fin n)
  /-- σ takes T₁'s arcs to T₂'s arcs. -/
  arc_eq  : ∀ i j : Fin n, T₁.arc i j = T₂.arc (perm i) (perm j)

namespace TournamentIso

variable {T₁ T₂ T₃ : Tournament n}

/-- The identity is an isomorphism from T to T. -/
@[refl] def refl (T : Tournament n) : TournamentIso T T where
  perm   := Equiv.refl _
  arc_eq := fun _ _ => rfl

/-- Symmetry: invert the permutation. -/
@[symm] def symm (f : TournamentIso T₁ T₂) : TournamentIso T₂ T₁ where
  perm   := f.perm.symm
  arc_eq := fun i j => by
    have h := f.arc_eq (f.perm.symm i) (f.perm.symm j)
    rw [Equiv.apply_symm_apply, Equiv.apply_symm_apply] at h
    exact h.symm

/-- Transitivity: compose two isomorphisms. -/
@[trans] def trans (f : TournamentIso T₁ T₂) (g : TournamentIso T₂ T₃) :
    TournamentIso T₁ T₃ where
  perm   := f.perm.trans g.perm
  arc_eq := fun i j => by
    rw [f.arc_eq, g.arc_eq]; rfl

end TournamentIso

/-- The relation "T₁ is isomorphic to T₂". -/
def IsomorphicTo (T₁ T₂ : Tournament n) : Prop :=
  Nonempty (TournamentIso T₁ T₂)

scoped infix:50 " ≅ " => Tournament.IsomorphicTo

namespace IsomorphicTo

variable {T₁ T₂ T₃ : Tournament n}

@[refl] lemma refl (T : Tournament n) : T ≅ T := ⟨TournamentIso.refl T⟩

@[symm] lemma symm (h : T₁ ≅ T₂) : T₂ ≅ T₁ :=
  h.elim (fun f => ⟨f.symm⟩)

@[trans] lemma trans (h₁ : T₁ ≅ T₂) (h₂ : T₂ ≅ T₃) : T₁ ≅ T₃ :=
  h₁.elim (fun f => h₂.elim (fun g => ⟨f.trans g⟩))

end IsomorphicTo

/-! ### Connection to self-complementarity (defined in GridReflection.lean)

    `IsSelfComplementary T` (from `GridReflection.lean`) is defined as
    `∃ σ : Equiv.Perm (Fin n), ∀ i j, T.arc i j = (op T).arc (σ i) (σ j)`.
    This module's `IsomorphicTo` gives the cleaner statement:
        IsSelfComplementary T ↔ T ≅ op T. -/

lemma isSelfComplementary_iff_iso_op (T : Tournament n) :
    IsSelfComplementary T ↔ T ≅ op T := by
  constructor
  · rintro ⟨σ, h⟩
    exact ⟨{ perm := σ, arc_eq := h }⟩
  · rintro ⟨f⟩
    exact ⟨f.perm, f.arc_eq⟩

end Tournament
