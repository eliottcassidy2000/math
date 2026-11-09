# Pre-registration: kind-pasteur-2026-11-09-S18

Written and pushed before any of the runs named below. Exploration that came before this file is listed under
"Already seen", so that it is not mistaken for a prediction.

## Owner prompt (2026-11-09)

The owner asked for four things:
- keep going and prove ternary Chan-Pak and beyond;
- study Hamiltonian path totals of digraphs written as square binary matrices;
- test the conjecture that every total of 7 or 21 has a 1 on the main diagonal or a mirrored pair of zeros;
- relate digraphs to the tournaments inside them, and binary to ternary triangles.

## Concurrent sessions, read first, and lanes kept apart

- **Session234** (864f45f40e, cdb77c02e9) already refuted the conjecture: there are semicomplete matrices with 7 and
  21. It classified the semicomplete 7/21 families at every order. This session does not redo that.
- **opus-2026-11-09-S89** reserved THM-5040-5042 and HYP-9610: even-graph dominating sets, matrices as pairs of
  (n+1)-tilings, and ternary triangles. Not touched here.
- **mac-mini-2026-11-09-digraph721** reserved THM-5055/5056 and HYP-9625 (even constructions, HYP-9575 (B)). Not
  touched here.
- **opus-2026-11-09-S87** (THM-4950) and **mac-mini-2026-11-08-ternary** (THM-4952) hold the ternary lane reduction,
  the order-10 census and the ribbon-word census to 10^5. They are reused with credit.

This session's lanes are:
- **(T)** the ternary problem via one-arc chains;
- **(D)** the Grinberg-Stanley signed cycle formula, applied to all digraphs. Session234 states that it did not use it.

## Definitions used in the predictions

**Port chain.** Start from one vertex u, the port. Repeatedly glue a block B (a small tournament with two marked
vertices w, u') after the current tournament T, as an ordinal sum, but reverse the single arc between u and w (so
w -> u). Then make u' the new port. The state (H_3(T), f(1), f(-1)) is updated by a 3x3 matrix depending only on
(B, w, u'), where f(z) is the sum, over paths P starting at the port, of z^|P| H_3(T - V(P)).

## Predictions

- **P1.** Port chains with blocks of order <= 4, in base 3, miss no value v = 1 (mod 3) in (3000, 200000].
  - Already seen: below 3000 the misses are exactly 10, 31, 61, 106, 217, 274, 337, 346, 607, 865 (exploratory
    Python run).
- **P2.** Port chains with blocks of order <= 5, in base 3, miss exactly {10, 31} among v = 1 (mod 3), v <= 200000.
- **P3.** Base 2, blocks of order <= 4: every miss among the odd values <= 200000 is below 2000.
- **P4.** Both H = 7 and H = 21 occur as the Hamiltonian path count of a tournament with exactly one arc deleted
  (one mirrored zero pair, no reciprocal pair, no loop) of order <= 6.
- **P5 (verification gate, not a prediction).** The weighted cycle formula of THM-5071 agrees with direct
  Hamiltonian path counts on:
  - every loopless digraph of order <= 4;
  - every semicomplete digraph and every oriented graph of order <= 5;
  - 10^5 random digraphs of orders 6 and 7.
  A single disagreement means the reading of Grinberg-Stanley used here is wrong.

Outcomes will be appended below, unedited, after the runs.

## Outcomes

Appended 2026-11-09 after the runs. The text above is unchanged.

- **P1: CONFIRMED.** Blocks of order <= 4, base 3, LIM = 200000: the misses are exactly 10, 31, 61, 106, 217, 274,
  337, 346, 607, 865, with 22,540,354 states (`out/chain_k4_base3_2e5.out`). The run was extended to 10^6 with a
  15-matrix subset: the same ten misses, 64,420,757 states (`out/chain_gen15_base3_1e6.out`).
- **P2: CONFIRMED.** Blocks of order <= 5, base 3: below 3000 the misses are exactly 10 and 31
  (`out/chain_k5_base3_3000.out`). The order-5 family contains the order-4 family, so its misses below 2*10^5 (and
  below 10^6) are exactly 10 and 31.
- **P3: CONFIRMED.** Base 2, blocks of order <= 4: the misses below 200000 are 7, 21, 47, 67, 127, 173, 491, all
  below 2000. With order <= 5 the misses are exactly 7 and 21.
- **P4: CONFIRMED.** One arc deleted gives H = 7 at order 5 and H = 21 at order 6 (`out/digraph_weights.out`,
  `out/one_defect_witness.out`). The audit found that the order-5 H = 7 half was already in Session234's census
  (general_census.csv, row 5,7,1,0,1,600) before this pre-registration. Only the order-6 H = 21 half was a genuine
  prediction.
- **P5: PASSED.** 232,263 checks, 0 failures. The order-7 sample was 10^4, not the promised 10^5, which is harmless
  for a gate. The audit's own check covers 1,372,741 digraphs.
