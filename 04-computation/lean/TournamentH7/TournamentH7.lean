/-
  TournamentH7 — Lean 4 formalization of THM-343 (H(T) ≠ 7).

  Module layout:
    · TournamentH7.Basic   — Tournament structure.
    · TournamentH7.Cycles  — Directed cycles, odd cycles.
    · TournamentH7.SCC     — Reachability, SCC, Hamiltonian paths, H(T).
    · TournamentH7.OCF     — External axioms (OCF, Moon–Moser, Moon–Camion).
    · TournamentH7.H7      — Main theorem `Tournament.H_ne_seven`.
-/

import TournamentH7.Basic
import TournamentH7.Cycles
import TournamentH7.SCC
import TournamentH7.OCF
import TournamentH7.H7
import TournamentH7.Verify
