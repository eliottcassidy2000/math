/-
  TournamentH7.Basic — Tournament definition and basic notions

  Defines:
   • `Tournament n`            : an n-vertex tournament as an asymmetric arc
                                  relation on `Fin n`.
   • `Tournament.beats`        : the strict beats relation.
   • `Tournament.outDegree`    : `s(v)` = number of vertices v beats.

  Reference: project canon definitions
    `/home/ubuntu/math/01-canon/definitions.md`.
-/

import Mathlib.Combinatorics.SimpleGraph.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Data.Finset.Card
import Mathlib.Data.Nat.Basic

open Finset

/-- An n-vertex tournament: an asymmetric, irreflexive arc relation on
    `Fin n`. For any two distinct vertices i, j exactly one of `arc i j` or
    `arc j i` holds. -/
structure Tournament (n : ℕ) where
  arc        : Fin n → Fin n → Bool
  irrefl     : ∀ i, arc i i = false
  total      : ∀ i j, i ≠ j → (arc i j = true ∨ arc j i = true)
  asym       : ∀ i j, ¬ (arc i j = true ∧ arc j i = true)

namespace Tournament

variable {n : ℕ}

/-- `beats T i j` says vertex i beats vertex j in the tournament T. -/
def beats (T : Tournament n) (i j : Fin n) : Prop := T.arc i j = true

/-- The out-degree (score) of a vertex. -/
def outDegree (T : Tournament n) (v : Fin n) : ℕ :=
  (Finset.univ.filter (fun w => T.arc v w = true)).card

/-- Notation `i →[T] j` for `beats T i j`. -/
scoped notation:50 i " →[" T "] " j => Tournament.beats T i j

end Tournament
