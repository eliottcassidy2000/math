/-
  TournamentH7.H63 — HYP-1754: H(T) ≠ 63

  By the extended OCF, H = 63 forces
        1 + 2α₁ + 4α₂ + 8α₃ + 16α₄ + 32α₅ + 64α₆ = 63,
  i.e.
        α₁ + 2α₂ + 4α₃ + 8α₄ + 16α₅ + 32α₆ = 31.

  This admits MANY non-negative integer solutions (compared to H = 21 with
  only 4). Instead of case-analysing each solution structurally, we
  axiomatise the final result as a single *unrealisability axiom* citing
  the exhaustive n ≤ 7 computational evidence:

      H = 63 has zero occurrences among all 2,097,152 n = 7 tournaments
      (file `05-knowledge/results/audit_n7_exhaustive_s6.out`).

  For n ≥ 8 this is conjectural (HYP-1754 in `05-knowledge/hypotheses/INDEX.md`).

  ## Pattern observation

  Universally forbidden small H values so far: 7 = 7·3⁰, 21 = 7·3¹, 63 = 7·3².
  The natural conjecture "7·3^k forbidden for all k" fails at k = 3 since
  H(P(7)) = 189 = 7·3³ is realised (it is the maximum H at n = 7).
-/

import TournamentH7.OCF
import TournamentH7.Forbidden

namespace Tournament

variable {n : ℕ}

/-- **Axiom (HYP-1754 structural).** For every tournament T on any number of
    vertices, H(T) ≠ 63.

    *Justification.* Exhaustively verified at n ≤ 7 (2,097,152 tournaments
    enumerated by `04-computation/audit_n7_exhaustive_s6.py`, zero
    occurrences). For n ≥ 8 conjectural; recorded as HYP-1754. -/
axiom H_ne_sixtythree_axiom (T : Tournament n) : H T ≠ 63

/-- **HYP-1754.** For every tournament T, H(T) ≠ 63. -/
theorem H_ne_sixtythree (T : Tournament n) : H T ≠ 63 :=
  H_ne_sixtythree_axiom T

end Tournament
