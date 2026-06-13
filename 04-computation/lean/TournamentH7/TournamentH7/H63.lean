/-
  TournamentH7.H63 — n ≤ 7 evidence: H(T) ≠ 63

  By the extended OCF, H = 63 forces
        1 + 2α₁ + 4α₂ + 8α₃ + 16α₄ + 32α₅ + 64α₆ = 63,
  i.e.
        α₁ + 2α₂ + 4α₃ + 8α₄ + 16α₅ + 32α₆ = 31.

  This admits MANY non-negative integer solutions (compared to H = 21 with
  only 4). Exhaustive computation proves absence through n = 7:

      H = 63 has zero occurrences among all 2,097,152 n = 7 tournaments
      (file `05-knowledge/results/audit_n7_exhaustive_s6.out`).

  IMPORTANT CORRECTION (opus-2026-05-29-S8):
  H = 63 is ACHIEVABLE at n = 8.  See
  `04-computation/h63_counterexample_audit_s8.py` and
  `05-knowledge/results/h63_counterexample_audit_s8.out`.
  Therefore this module must not export a universal theorem H(T) ≠ 63.

  ## Pattern observation

  Universally forbidden small H values currently supported by canon are
  7 = 7·3⁰ and 21 = 7·3¹.  The apparent continuation to 63 was a finite-n
  mirage: the n = 8 counterexample has Ω(T) = K31, so
  H(T) = I(K31, 2) = 1 + 2·31 = 63.
-/

import TournamentH7.OCF
import TournamentH7.Forbidden

namespace Tournament

variable {n : ℕ}

/-- **Axiom (finite verification).** For every tournament T on at most seven
    vertices, H(T) ≠ 63.

    *Justification.* Exhaustively verified at n ≤ 7 (2,097,152 tournaments
    enumerated by `04-computation/audit_n7_exhaustive_s6.py`, zero
    occurrences).  The corresponding universal statement is false at n = 8. -/
axiom H_ne_sixtythree_le_seven_axiom (hn : n ≤ 7) (T : Tournament n) : H T ≠ 63

/-- Exhaustive finite result: if n ≤ 7, then H(T) ≠ 63. -/
theorem H_ne_sixtythree_le_seven (hn : n ≤ 7) (T : Tournament n) : H T ≠ 63 :=
  H_ne_sixtythree_le_seven_axiom hn T

end Tournament
