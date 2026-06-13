/-
  TournamentH7.Cycles — Directed cycles in tournaments

  Defines:
   • `Tournament.DirectedCycle T k`  : a directed k-cycle in T as an
                                       injection `Fin k → Fin n` with cyclic
                                       beats relation (using `% k` to step).
   • Helper: a cycle is *odd* iff its length is odd.

  Cycles are represented as `Fin k → Fin n` injections satisfying the
  cyclic arc condition; equality is on the underlying tuple, so each cyclic
  equivalence class is represented by `k` distinct tuples. The concrete
  representation does not matter downstream — `OCF` axiomatises everything
  we need at the count level.
-/

import TournamentH7.Basic
import Mathlib.Algebra.Ring.Parity

namespace Tournament

variable {n : ℕ}

/-- The cyclic successor in `Fin k` (uses `% k`; requires k ≥ 1). -/
private def cyclicSucc {k : ℕ} (h : 0 < k) (i : Fin k) : Fin k :=
  ⟨(i.val + 1) % k, Nat.mod_lt _ h⟩

/-- A directed k-cycle in T (k ≥ 3): an injection `Fin k → Fin n` such
    that consecutive vertices (cyclically, modulo k) are connected by arcs
    of T. -/
structure DirectedCycle (T : Tournament n) (k : ℕ) where
  toFun       : Fin k → Fin n
  inj         : Function.Injective toFun
  three_le    : 3 ≤ k
  arc_step    : ∀ i : Fin k,
                  T.arc (toFun i) (toFun (cyclicSucc (by omega) i)) = true

/-- A directed cycle has *odd length* iff `k` is odd. -/
def DirectedCycle.isOdd {T : Tournament n} {k : ℕ}
    (_ : DirectedCycle T k) : Prop := Odd k

end Tournament
