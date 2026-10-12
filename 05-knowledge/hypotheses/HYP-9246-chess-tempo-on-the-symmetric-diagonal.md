---
id: HYP-9246
title: "The tempo on the symmetric diagonal of chess: on a board fixed by a value-negating involution (board isometry composed with the colour swap) the mover's value is a sign t(B) in {+1, 0, -1}; the owner's axiom V(s0) >= 0 says t(s0) >= 0; in chess-with-passing t >= 0 always (strategy stealing); conjecture C1: with pawns only, facing each other on the same files (rank-reflection diagonal), t(B) >= 0 for every legal board, i.e. there is no rank-symmetric mutual zugzwang with kings and pawns; exhaustive in KPvKP (0 of 2,184 boards), untested beyond"
status: >
  OPEN. PROVED: the passing-chess statement (strategy stealing) and the symmetry facts used (THM-4612 statement 6).
  FINITE-EXACT relative to CITED Syzygy values (served by tablebase.lichess.ovh, every response archived in
  lichess_cache.json): the KPvKP census on the rank-reflection diagonal, 2,184 legal boards = 1,092 representatives
  under the file reflection, mover wins 450, draws 642, mutual zugzwangs 0. NUMERICAL (seeded samples of the
  diagonal): KPvKP rotation diagonal, 300 of 1,358 representatives: mover wins 83, draws 217, mutual zugzwangs 0; KRvKR rank diagonal, 200 of 1,248: 167 / 33 / 0; KQvKQ rank diagonal, 200 of 1,048: 175 / 25 / 0 (seed 2026; lichess throttles the API, so the secondary classes were sampled). The conjecture C1 is the sign statement at the pawn-only diagonal; it neither
  implies nor follows from the owner's axiom A0, which is the same sign statement at the full-board root.
session: opus-2026-10-12-S23
source: 05-knowledge/results/chess_proof_system_20261012.md (sections 2.3 and 8)
scripts:
  - 04-computation/experiments/chess_20261012/lichess_tools.py (tempo census; + lichess_tempo.out, lichess_cache.json, lichess_cache_<class>_<sym>.txt)
related:
  - THM-4612 (value symmetries; the rank reflection with colour swap fixes the initial board)
  - THM-4613 (traps; A0 says the maximal White trap contains the root)
---

# HYP-9246 — the tempo on the symmetric diagonal

## Setting

`tau` is a value-negating involution of legal states: a board isometry composed with the colour swap, the side swap, the swap of castling rights and the reflection of the en-passant square (THM-4612 statement 6: the rank reflection always; the file reflection and the rotation without castling rights; all eight isometries without pawns). A board `B` with `tau(B) = B` has the **tempo** `t(B) = V(B, White to move) = -V(B, Black to move)` (White's view), which is the mover's value `v` whichever side is to move, in `{+1, 0, -1}`; `t = -1` is a mutual zugzwang.

## Statements

1. **(Passing chess; PROVED.)** If either side may pass, `t(B) >= 0` for every `tau`-fixed board: a mover who lost would pass and hand the opponent the same symmetric board to move in, which the opponent then loses by symmetry.
2. **(The axiom.)** The owner's A0 is `t(s0) >= 0` for the initial board and the rank reflection: there is no zugzwang at the root of real chess. It is not a consequence of 1, because deleting the pass changes both players' options.
3. **(Conjecture C1, OPEN.)** For material consisting of kings and pawns only and the rank reflection (every pawn faces its mirror image on the same file, the kings on one file at ranks `r` and `9-r`), every legal board has `t(B) >= 0`.
   *Evidence.* KPvKP, exhaustive: 2,184 legal rank-symmetric boards (1,092 up to the file reflection, which preserves values); mover wins 450, draws 642, mutual zugzwangs 0 (Syzygy WDL, archived). Not tested with two or more pawns per side.
   *Why it might be true.* On the rank diagonal the pawns on a file either block each other or have already passed; the mover can always spend a tempo on the king, and the only way a king move loses is to release a pawn or walk into the opponent's king, which the symmetric opponent must then also face one move later. *Why it might fail.* With several pawns per side, blocked chains create reciprocal zugzwangs in asymmetric endgames (the trebuchet, which is point-symmetric about the midpoint of two adjacent squares and not about the board's centre, so it is not on any diagonal); whether some reciprocal zugzwang can be made rank-symmetric is exactly the question.
4. **(Observations on the other diagonals; NUMERICAL, seeded samples.)** On the rotation diagonal of KPvKP (pawns on mirrored files; a sample of 300 of the 1,358 representatives) the mover wins 83 and draws 217; on the rank diagonals of KRvKR (200 of 1,248) the mover wins 167 and draws 33, of KQvKQ (200 of 1,048) wins 175 and draws 25; no mutual zugzwang in any sample. The mover's winning share grows with the range of the material (28%, 84%, 88% in the samples against 41% exhaustive for KPvKP on the rank diagonal): the tempo converts to a win when the first check or skewer decides, and mostly to a draw with pawns. A seeded sample cannot exclude rare mutual zugzwangs in these classes; only the KPvKP rank row is exhaustive.

## Cheapest decisive tests

- KPPvKPP on the rank diagonal: about 54,000 legal boards; a seeded sample of a few hundred through the same tool decides C1 numerically, a full run exhausts it (6-man Syzygy exists).
- A mutual zugzwang on any pawn-only rank-symmetric board refutes C1 outright.

## Mechanism relation to the thread

A0 asserts the sign of the tempo at the unique root of the full game; C1 asserts it on the pawn-only slice of the diagonal that tablebases reach. Both are statements that symmetry alone does not give (statement 1 shows what symmetry does give, once passing is allowed), and both are about the one asymmetry of a symmetric board: who moves.
