---
id: THM-3330
title: "AMM 12592: golden attainment to n <= 2047; the plain-rule slack law is linear (Estimate E refuted for rule A); Y-box reduction, Lucas diagonal of G_R, and the R/3 cancellation law"
status: >
  PROVED + VERIFIED-EXACT / AWAITING CROSS-SESSION HOSTILE AUDIT.
  (A) ATTAINMENT. Epochs close at the gamma* floor + D0 with exact new
  records: D0*(R) = 0 (R <= 128), 1 (256), 5 (512, D0 <= 4 dies), 15 (1024),
  38 (2048), 89 (4096), 192 (8192), in [401, 416] (16384) — each pinned by
  adjacent die/close pairs on a fast junk-flow engine certified bit-identical
  to the slow polynomial path across three implementations, with 16/16
  witnesses passing a fresh from-statement referee (new O(M^2) integer
  identity engine via t = (1-x)/x). Consequence:
  T(n) <= n + 1 + floor(gamma* n) + 15 for ALL critical n <= 2047.
  (B) GROWTH LAW / ESTIMATE E STATUS. Threshold ratios D0*(2R)/D0*(R) =
  5, 3, 2.53, 2.34, 2.16, ~2.1 decline to 2 from above: PLAIN rule A needs
  LINEAR slack, D0*(R)/R -> eps_inf in [0.0245, 0.028] (CONJECTURED limit;
  refutations of both pre-registered sublinear laws are exact: the polylog
  fit C(n+3,4) dies at R = 2048 (35/36/37 die), the 25/1024 saturation dies
  at R = 16384 (D0 = 400 dies at row 4055, predicted ~700 rows in advance by
  the two-phase law)). So Estimate E with o(R) slack is dead FOR RULE A;
  this says nothing about epoch feasibility (hazard discipline) and the
  bulk-rule route is open — see (D).
  (C) PROVED STRUCTURE. T4b exact edge lemma (initial junk = cells
  [0, R-2-d_0] exactly; tau* = (1-gamma*)/gamma*, death-delay bound
  2d_0 - R + 2 = (2gamma*-1)R + 2D0 + O(1) closed-form); T9a edge transport;
  T9b golden march threshold (cap decay along march diagonals iff t/d >
  1/phi — the golden ratio reappears as the march stability threshold); T9'
  two-phase death law (dawdle then exact 1/row march; death row = L + gap(L)
  exactly at six scales). Y-BOX REDUCTION (parity eliminated): Delta_i =
  (p-q) + 2 gamma_i admissible iff Bernstein cells y_k of gamma_i lie in
  [-C(d-1,k), +C(d-1,k-1)]; the epoch problem is sum x^i gamma_i = G_R in
  these parity-free boxes. LUCAS DIAGONAL: G_R = 1 + sum_k m_k p^k q^{k-1}
  is the UNIQUE such expansion, m_k explicit alternating with |m_k| ~ phi^R
  peaking at k/R ~ 0.276; step law 2(G_{R+1} - G_R) = -p L_{R-1}; doubling
  law G_{2R} = 2q G_R^2 + (p-q)(2G_R - 1) - p^R q^{R-1}. First closed-form
  epoch witness: gamma = (q^2, -p^2, 0) at R = 4.
  THE R/3 LAW (exact capacity computation, NOT a search negative):
  sign-coherent cellwise non-cancelling splittings of the Lucas-diagonal
  atoms are IMPOSSIBLE for 2 <= k < k*(R), k*/R -> 1/3, with deficit up to
  ~0.26R bits — cancellation is structurally mandatory in every all-R
  construction; winning witnesses confine deep cancellation to the top
  ~R/8 coefficient band (max ~0.21R bits).
  (D) THE DIE IS A BULK PHENOMENON. Steering-invariance: modifying steered
  cells k <= K changes later residuals only at absolute positions
  >= i + d_i - K, while every death reads position ~ (2gamma*-1)R + 2D0
  << d_0 — so no ballot/low-cell schedule can prevent it. The initial junk
  profile alternates in sign and the Pascal kernel (1,2,1) annihilates
  alternating profiles; the march is driven by the non-alternating residue.
  The golden route is therefore a BULK alternation-shaping clamp
  (in progress, next-wave lanes D1-D3), with (E-lin) as the rigorous
  fallback: closure at D0 = eps R for each eps > eps_inf would give
  C* <= 1 + gamma* + eps_inf unconditionally.
source: boxeph-2026-08-03-multifront-allR
depends_on:
  - THM-3329
  - THM-3302
  - THM-3002
  - THM-3026
related:
  - HYP-9061
  - THM-2966
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/amm12592_transient_fast_junkflow_boxeph.py
referee: 04-computation/amm12592_allR_referee_boxeph.py
notes:
  - 05-knowledge/results/amm12592-allR-transient-theorem-boxeph.md
  - 05-knowledge/results/amm12592-allR-GR-representation-boxeph.md
output: 05-knowledge/results/amm12592_allR_referee_boxeph.out
---

# THM-3330 -- the growth law, the Y-box/Lucas structure, and attainment to 2047

Statement, proofs, engines, growth tables, refutation ledger, and the
death-mechanism anatomy are in the two companion notes listed in the
frontmatter (Part I/II of the transient-theorem note; the G_R representation
note), with every lemma labeled PROVED / VERIFIED-finite / CONJECTURED and
every computation exact. Witnesses: 16 files under
`04-computation/amm12592_witness_*` and `amm12592_floor_witness_*`
(R = 1024 rule-B witness stored as SHA-256 + deterministic reproduction
record; sparse formats documented in-file).

The remaining program for `C* = log_5(5 phi^2)`:

1. **Upper (golden):** bulk alternation-shaping rule + its transient bound
   (Estimate E', next wave), or any all-R o(R)-slack construction; the
   rigorous fallback (E-lin) yields `C* <= 1 + gamma* + eps_inf` already if
   proved.
2. **Lower:** the deadline-bounded routing window (MISTAKE-361 repair).
