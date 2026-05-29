/-
  TournamentH7 — Lean 4 formalization of forbidden H-values in tournaments.

  ## Theorems formalised

  - `Tournament.H_ne_seven`         (THM-343)
  - `Tournament.H_ne_twentyone`     (HYP-1753)
  - `Tournament.H_ne_sixtythree`    (HYP-1754)
  - `Tournament.redei_existence`    (Rédei 1934, axiom)
  - `Tournament.redei_parity`       (Rédei 1934, axiom)
  - `Tournament.H_pos`              (corollary: H ≥ 1)
  - `Tournament.H_ne_two`           (corollary of parity)
  - `Tournament.H_ne_even`          (corollary of parity)
  - `Tournament.H_not_in_forbidden_trio`

  ## Module layout

    · TournamentH7.Basic        — Tournament structure.
    · TournamentH7.Cycles       — DirectedCycle, isOdd.
    · TournamentH7.SCC          — Reachability, IsSCC, HamPath, H(T).
    · TournamentH7.OCF          — External axioms (OCF, Moon, etc.).
    · TournamentH7.H7           — Theorem H_ne_seven (THM-343).
    · TournamentH7.Forbidden    — Chain step, binomial bound, extended OCF.
    · TournamentH7.H21          — Theorem H_ne_twentyone (HYP-1753).
    · TournamentH7.H63          — Theorem H_ne_sixtythree (HYP-1754).
    · TournamentH7.Redei        — Rédei 1934 axioms + corollaries.
    · TournamentH7.HSpectrum    — Forbidden-trio bundle.
    · TournamentH7.Tilings      — Base path, tile-complement T̃, score formula
                                  (project-novel; oracle-2026-05-11-S1).
    · TournamentH7.GridReflection — THM-280 (grid reflection ↔ complement;
                                  project-novel; opus-2026-04-03-S27).
    · TournamentH7.Verify       — Axiom audit (#print axioms).
-/

import TournamentH7.Basic
import TournamentH7.Cycles
import TournamentH7.SCC
import TournamentH7.OCF
import TournamentH7.H7
import TournamentH7.Forbidden
import TournamentH7.H21
import TournamentH7.H63
import TournamentH7.Redei
import TournamentH7.HSpectrum
import TournamentH7.Tilings
import TournamentH7.GridReflection
import TournamentH7.StaircaseModel
import TournamentH7.SelfComplementary
import TournamentH7.AntiAutomorphism
import TournamentH7.HPIPIdentity
import TournamentH7.Iso
import TournamentH7.IsoProperties
import TournamentH7.SCCounts
import TournamentH7.Verify
