/-
klein-2026-07-02-S113 — AXIOM-FOOTPRINT AUDIT of the LRC(14) endgame surface.

Goal: state, unambiguously and machine-checked, EXACTLY what the LRC(14) proof depends on.
`#print axioms` on the top surface theorem reveals its complete axiom footprint. The target
"unconditional modulo the owner-sanctioned LRC(≤13) citation" reads, in the printout, as:

    lrc14_of_covering_far_22 depends on:
      [propext, Classical.choice, Quot.sound]   -- Lean's foundational axioms (always allowed)
      [Lean.ofReduceBool]                        -- the native_decide kernel-reduction axiom (the
                                                    machine-checked winData22 band gate/completeness)
    and takes as HYPOTHESES (not axioms):
      cite : LRCUpTo13          -- the owner-sanctioned LRC(<=13) citation node
      hcovfar : CoveringFarLonely 22   -- the SINGLE remaining analytic hypothesis

There must be NO `sorryAx`.  This file is audit-only: no new mathematical content.
-/
import TournamentH7.LRC14CoveringFarSurface
import TournamentH7.LRCWindowData22

open LonelyRunner LonelyRunner.LRC14 WindowData

-- The top endgame surface (S112): LRC(14) from {citation} + {CoveringFarLonely 22}.
#print axioms lrc14_of_covering_far_22

-- The window-census leg (should be: foundations + native_decide, NO sorry).
#print axioms hwindow22_closed

-- A certified infinite slice of the remaining hypothesis (block family, V ≥ 21).
#print axioms coveringFar_blockFamily

-- A certified deep-well instance of the remaining hypothesis (![1..12,182]).
#print axioms coveringFar_deepWell

-- The generic-window assembly (cut-parametric).
#print axioms lrc14_of_covering_far_of_window
