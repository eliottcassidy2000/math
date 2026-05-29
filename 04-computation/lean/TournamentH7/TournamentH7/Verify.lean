/-
  TournamentH7.Verify — Audit: print the axiom dependencies of the main
  theorem `Tournament.H_ne_seven`.

  This file should be inspected after `lake build` — the `#print axioms`
  output should list exactly the 7 axioms documented in `OCF.lean`
  (plus Lean's foundational `propext`, `Classical.choice`, `Quot.sound`).
-/

import TournamentH7.H7

open Tournament

theorem H_ne_seven_audit {n : ℕ} (T : Tournament n) : H T ≠ 7 := H_ne_seven T

#print axioms H_ne_seven_audit
