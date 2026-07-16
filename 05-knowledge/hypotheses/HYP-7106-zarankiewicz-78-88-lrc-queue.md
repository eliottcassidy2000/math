# HYP-7106 — K_{7,8}/K_{8,8} vs open Zarankiewicz + LRC(14) queue work

**Status:** RESOLVED (death-star-2026-07-16-S30). (1) THE OPEN ZARANKIEWICZ CASES MET
EXACTLY: K_{7,8} class-min = 108 = Z(7,8), K_{8,8} = 144 = Z(8,8) (full class-coloring
enumeration; free per-edge annealing found nothing lower); controls K_{6,7} = 54, K_{6,8}
= 72 (proved-optimal range, Kleitman) — the parallel-class book constructs
Z(m,n)-achieving drawings at every tested (m,n), incl. balanced-necklace layouts for
unequal parts and mixed class sizes 1..7 (within-class = 0 always). (2) LRC QUEUE (Lean
absorption): FragmentationLemma.lean moved from loose draft into the TournamentH7 project
and BUILDS GREEN — badArcs_periodic PROVED sorry-free; the λ > ½ branch PROVED; the
draft's arc-counting plan CORRECTED (it can hit ⌊Lw⌋+2 arcs — the periodicity/windowing
architecture documented and half-built); three localized sorries remain with exact
Mathlib pointers (window_bound, main branch, killer_budget).

(1) The cyclic bipartite book at mixed/even parts: K_{8,8} on Z₁₆ (parity parts, odd-sum
classes); K_{7,8} on Z₁₅ (parity-position parts; classes reindex by +4 under the
part-preserving rotation-by-2 — still circulant-tractable); controls K_{6,7}, K_{6,8}
(Zarankiewicz PROVED there, Kleitman). Compare class-coloring minima vs Z(7,8) = 108,
Z(8,8) = 144, Z(6,7) = 54, Z(6,8) = 72.
(2) LRC(14) queue (boxeph-S39): take an open item executable solo (band audit / Lean
absorption prep), honest scope.

-> THM-922, LEM-030, THM-913, boxeph-S39 queue; death-star-S30.
