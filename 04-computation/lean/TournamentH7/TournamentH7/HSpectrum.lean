/-
  TournamentH7.HSpectrum — Forbidden H-values

  Collects the universal forbidden-H theorems proven in this project, the
  finite n ≤ 7 absence of H = 63, and the parity result of Rédei.

  ## Forbidden H-values (universally, all n)

    · H(T) ≠ 7      — THM-343  (`H_ne_seven`,  H7.lean)
    · H(T) ≠ 21     — HYP-1753 (`H_ne_twentyone`, H21.lean)
    · H(T) even ≠ → False  (`H_ne_even`, Redei.lean)

  ## Finite verification only

    · If n ≤ 7, H(T) ≠ 63 (`H_ne_sixtythree_le_seven`, H63.lean).

  ## Pattern observation (kept here for documentation)

  Universally forbidden small odd H values currently supported by canon are
  {7, 21} = {7·3⁰, 7·3¹}.  The apparent continuation to 63 is false:
  H = 63 is realised at n = 8 by a tournament with Ω(T) = K31.
  References: `05-knowledge/results/audit_n7_exhaustive_s6.out`,
  `05-knowledge/results/h63_counterexample_audit_s8.out`.
-/

import TournamentH7.H7
import TournamentH7.H21
import TournamentH7.H63
import TournamentH7.Redei

namespace Tournament

variable {n : ℕ}

/-- Disjunction of the small universally-forbidden H-values:
    H(T) ∉ {7, 21} for every tournament T. -/
theorem H_not_in_forbidden_pair (T : Tournament n) :
    H T ≠ 7 ∧ H T ≠ 21 :=
  ⟨H_ne_seven T, H_ne_twentyone T⟩

/-- Finite verification bundle: if n ≤ 7, then H(T) ∉ {7, 21, 63}. -/
theorem H_not_in_forbidden_trio_le_seven (hn : n ≤ 7) (T : Tournament n) :
    H T ≠ 7 ∧ H T ≠ 21 ∧ H T ≠ 63 :=
  ⟨H_ne_seven T, H_ne_twentyone T, H_ne_sixtythree_le_seven hn T⟩

end Tournament
