/-
  TournamentH7.Examples — Worked examples of the formalisation.

  Demonstrates the formalisation by computing properties of concrete
  small tournaments.  All examples here are FULLY PROVED in Lean (no
  project axioms beyond Lean foundations).
-/

import TournamentH7.SmallTournaments
import TournamentH7.Iso
import TournamentH7.IsoProperties
import TournamentH7.StaircaseModel

namespace Tournament

/-! ### Example 1: 3-cycle properties -/

/-- The 3-cycle 0 → 1 → 2 → 0 is regular. -/
example : IsRegular threeCycle := threeCycle_isRegular

/-! ### Example 2: Transitive tournament properties -/

/-- The transitive tournament on n vertices has the base path. -/
example (n : ℕ) : HasBasePath (transitiveTournament n) := transitive_hasBasePath n

/-- The transitive tournament on 3 vertices has score sequence (0, 1, 2). -/
example : (transitiveTournament 3).outDegree ⟨0, by omega⟩ = 0 := by
  show (Finset.univ.filter
    (fun w : Fin 3 => (transitiveTournament 3).arc ⟨0, by omega⟩ w = true)).card = 0
  decide

example : (transitiveTournament 3).outDegree ⟨2, by omega⟩ = 2 := by
  show (Finset.univ.filter
    (fun w : Fin 3 => (transitiveTournament 3).arc ⟨2, by omega⟩ w = true)).card = 2
  decide

/-- The transitive tournament on 4 vertices is not regular. -/
example : ¬ IsRegular (transitiveTournament 4) :=
  transitive_not_regular 4 (by omega)

/-! ### Example 3: Iso identity is iso -/

/-- The identity is a tournament isomorphism. -/
example (n : ℕ) (T : Tournament n) : T ≅ T := IsomorphicTo.refl T

/-- Iso is symmetric. -/
example (n : ℕ) (T₁ T₂ : Tournament n) (h : T₁ ≅ T₂) : T₂ ≅ T₁ :=
  IsomorphicTo.symm h

/-- Iso preserves Hamiltonian path count. -/
example (n : ℕ) (T₁ T₂ : Tournament n) (h : T₁ ≅ T₂) : H T₁ = H T₂ :=
  H_iso_invariant T₁ T₂ h

/-! ### Example 4: Out-degree preservation under iso -/

/-- Concrete check: out-degree of vertex `f.perm v` in T₂ equals out-degree
    of `v` in T₁ when `f : T₁ ≅ T₂`. -/
example (n : ℕ) (T₁ T₂ : Tournament n) (f : TournamentIso T₁ T₂) (v : Fin n) :
    T₁.outDegree v = T₂.outDegree (f.perm v) :=
  outDegree_iso T₁ T₂ f v

/-! ### Example 5: Reachability composition -/

example (n : ℕ) (T : Tournament n) (a b c : Fin n)
    (r₁ : Reaches T a b) (r₂ : Reaches T b c) : Reaches T a c :=
  r₁.trans T r₂

/-! ### Example 6: Base-path descent -/

/-- Concrete: in any tournament with the base path, vertex 5 reaches vertex 0
    (when n ≥ 6). -/
example (n : ℕ) (T : Tournament n) (hbp : HasBasePath T) (h : 6 ≤ n) :
    Reaches T ⟨5, by omega⟩ ⟨0, by omega⟩ := by
  exact reaches_zero T hbp (by omega) ⟨5, by omega⟩

/-! ### Example 7: THM-330 — apex tile -/

/-- Concrete: if the apex arc 0 → (n-1) is present (along with the base path),
    then T is SC. -/
example (n : ℕ) (T : Tournament n) (hbp : HasBasePath T) (hn : 3 ≤ n)
    (h_apex : T.arc ⟨0, by omega⟩ ⟨n - 1, by omega⟩ = true) :
    IsStronglyConnected T :=
  apex_implies_SC T hbp hn (by omega) (by omega) h_apex

end Tournament
