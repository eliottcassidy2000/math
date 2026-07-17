# THM-947 — The arc wire: the dictionary, the 7-overlap constraints, and the finish-up batch (death-star-2026-07-17-S40)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCArcWire.lean` + extensions to
`LRCDeviationPairs.lean`/`LRCMomentCertificates.lean` + `S40AxiomAudit.lean`; verify
the build report in the session log). Source: HYP-7186. (THM-946 reserved for
kind-pasteur's T_s renumber per the collision resolution.)

## Statement

1. `bad_at_witness` (**the dictionary's engine**): a band-failing runner at
   multiplier p carries a canonical integer witness with `14·|v·p − n·q| < q` — the
   discrete band IS the killer arc's `‖v·(p/q)‖ < 1/14` in exact integers; hence
   `bandCount` counts bad-arc memberships and `CoverageCapped(6)` = "no rational
   p/q meets seven bad arcs" (`coverageCapped_iff_no_seven`). The discrete-Bonferroni
   and killer/fragmentation arcs now speak one language.
2. `seven_overlap_pair_constraint` (**the interface**): two simultaneous bads force
   the integer near-proportionality `14·q·|v_i·n_j − v_j·n_i| < (|v_i|+|v_j|)·q`.
   Seven simultaneous bads yield 21 such constraints; on trapped strata (THM-939)
   the small-witness instances are forbidden low-mass relations — the formal hook
   through which trapped families resist deep coverage. Supplying
   `CoverageCapped(6)` on explicit strata reduces to ruling out 7-solvable
   constraint systems — a rigidity/counting question, no longer an analytic one.
3. **Finish-up batch**: `dilateTripleCount_eq` — the c = 3 dilate count
   `N₃(Q) = 2⌊(Q−1)/3⌋` at 14 | q (the ratio-3 ladder's own exact step price;
   conjecture `N_c = 2⌊(Q−1)/c⌋` now proved at c = 2, 3);
   `B5_eq_live_capped5` — B5 = liveCount identically on cap-5 strata;
   `S40AxiomAudit` — the S37–S40 arc under one audit.

## Referee

`arcwire_referee_deathstar_S40.py`: witness bound, pair constraints, c = 3 formula
(planted), cap identities — see 05-knowledge/results/.
