/-
  TournamentH7.HSpectrum — Universally forbidden H-values

  Collects the three forbidden-H theorems proven in this project and the
  parity result of Rédei.

  ## Forbidden H-values (universally, all n)

    · H(T) ≠ 7      — THM-343  (`H_ne_seven`,  H7.lean)
    · H(T) ≠ 21     — HYP-1753 (`H_ne_twentyone`, H21.lean)
    · H(T) ≠ 63     — HYP-1754 (`H_ne_sixtythree`, H63.lean)
    · H(T) even ≠ → False  (`H_ne_even`, Redei.lean)

  ## Pattern observation (kept here for documentation)

  Universally forbidden small odd H values: {7, 21, 63} = {7·3⁰, 7·3¹, 7·3²}.
  The pattern fails at k = 3: H = 7·3³ = 189 is realised (it is the
  maximum H at n = 7, achieved by the Paley tournament on 7 vertices).
  Also 7·5 = 35 is realised (at n = 7).
  Reference: `05-knowledge/results/audit_n7_exhaustive_s6.out`.
-/

import TournamentH7.H7
import TournamentH7.H21
import TournamentH7.H63
import TournamentH7.Redei

namespace Tournament

variable {n : ℕ}

/-- Disjunction of the small universally-forbidden H-values:
    H(T) ∉ {7, 21, 63} for every tournament T. -/
theorem H_not_in_forbidden_trio (T : Tournament n) :
    H T ≠ 7 ∧ H T ≠ 21 ∧ H T ≠ 63 :=
  ⟨H_ne_seven T, H_ne_twentyone T, H_ne_sixtythree T⟩

end Tournament
