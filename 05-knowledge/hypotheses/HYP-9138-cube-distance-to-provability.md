---
id: HYP-9138
title: "In the strategy cube, Collatz's Haar distance to the provable class tends to zero: delta_k / 2^(k-1) -> 0, where delta_k is the least number of sign flips of Collatz at level k reaching class (i)"
status: >
  OPEN HYPOTHESIS; FINITE-EXACT for k <= 9 (MaxSAT, re-verified). The
  exact minima are delta_k = 1, 2, 2, 4, 5, 9, 14, 23 for k = 2..9, i.e.
  Haar fractions 0.5, 0.5, 0.25, 0.25, 0.156, 0.141, 0.109, 0.090,
  roughly 2^(-k/3). By THM-4474 Theorem A, class (i) means that every cycle
  of the parity graph has odd density < log_3 2, so delta_k is a minimal
  "expanding-cycle transversal" of the de Bruijn graph B(2,k), under the
  rewiring induced by sign flips.
source: collatz-procgen-20260922 session, strategy-cube lane (2026-09-25)
depends_on:
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md
related:
  - 01-canon/theorems/THM-4475-price-of-provable-descent-tends-to-zero.md (the pairing-family analogue, proved)
---

# HYP-9138 -- Collatz's distance to provability in the strategy cube

**Statement.** `delta_k / 2^(k-1) -> 0`.

**Contrast with THM-4475.** In the pairing family, *no* periodic member is
provable (THM-4470 §5, the pair-0 obstruction), and aperiodic provable
members exist at price `2^(-0.05L)` (THM-4475). In the strategy cube,
periodic provable strategies exist at every level `k >= 2`, since the sign
pair `(sigma(1), sigma(-1)) = (+, -)` avoids both expanding fixed points
(THM-4474, Proposition F). The question is whether they approach Collatz in
Haar measure.

**Evidence and route.**
* The optimal flip sets are mostly `u -> d` conversions at residues
  `3 (mod 4)`. They always include `-1`, forced by the loop at `-1`, and
  they break the `-5` cycle.
* The nearest provable strategies are near-critical: `rho_max = 5/8` at
  `k = 6..9`.
* A THM-4475-style greedy on residues, flipping exactly the residues whose
  classes carry expanding cycles, is the natural route. The obstacle is
  that a residue flip rewires the parity graph globally, whereas a pair flip
  changes one pair.
