---
id: THM-4550
title: "The owner's 8x8 board as a lonely-runner instrument, and the weave law of closed knight tours. (B) A bishop on ring k of the N x N board reaches 2N-3-2k squares, i.e. (N-2) + 2N lambda where lambda = (N-1-2k)/(2N) is the two-runner loneliness of the square's centre; the rook's reach is constant. (K) On the 2m x 2m board the knight's moves split as 8(m-1) inside a ring, 16 C(m,2) to the next ring, 16 C(m-1,2) across a ring. (D) For coprime a, b the two-speed loneliness is floor((a+b)/2)/(a+b): it lives on the anti-diagonal, equals 1/2 exactly for colour-preserving leapers, and its unique minimum 1/3 is the knight; on the 8x8 subdivision the knight is the only rider that never reaches the central four squares. (T) The first lonely time of two speeds is at most 4/9, attained only by the camel (1,3). (S) The odd lengths {1,3,...,2m-1} are lonely at 1/2, the even lengths are the tight progression 2{1..m}, their weave {1..2m} is tight. (W) In every closed knight tour of the 8x8 board the moves inside the outer ring equal the moves inside rings 1-2 in number, and there is at least one"
status: >
  PROVED (elementary) + FINITE-EXACT checks + INDEPENDENTLY AUDITED (2026-10-06, the session's adversarial
  audit subagent, own code, no claim refuted; corrections listed in the note's section 10 and MISTAKE-564).
  (D) and (S) are classical lonely-runner facts, re-proved; (W)'s colour lock is Posa's classical 4 x n
  argument (Schwenk 1991) transplanted from rows to rings.
  (B) bishop reach = 2N - 3 - 2 ring = (N-2) + 2N lambda(centre); queen 4N-5-2 ring; rook 2(N-1).
  (K) knight edges 4(2m-1)(2m-2) = 8(m-1) [within a ring; 8 per ring k >= 1, at its corners;
      two 4-cycles in ring 1, a matching of 16 squares in each ring k >= 2]
      + 16 C(m,2) [adjacent rings] + 16 C(m-1,2) [two rings]; no knight move crosses three rings.
  (D) delta(a,b) = max_t min(||ta||,||tb||) = floor(s/2)/s, s = a+b, gcd(a,b) = 1.
      Corollary: on the 2m x 2m subdivision of the unit cell the rider t(a,b) misses the closed
      central 2x2 block iff a+b is odd and < m; 8x8 is the smallest even board with an exception
      (closed block), the knight. (Open block: 6x6 already excludes the knight.)
  (T) tau(a,b) = min{t>0 : ||ta||,||tb|| >= 1/3}: tau(1,b) = 1/3 (3 not | b), tau(1,3k) = 1/3 + 1/(9k),
      tau(a,b) = 1/(3a) for a < b < 2a, tau(a,b) <= 2/(3a) for a >= 2, b > 2a. sup tau = 4/9, only (1,3).
      Three speeds: sup tau = 7/16, only (1,3,12) (HYP-9210 records the proof and the open n >= 4).
  (S) delta(1,3,...,2m-1) = 1/2; delta(2,4,...,2m) = 1/(m+1); delta(1,...,2m) = 1/(2m+1); the progression's
      lonely times at level 1/(n+1) are exactly k/(n+1), gcd(k, n+1) = 1.
  (W) Weave law. For a closed knight tour of the 8x8 board, #moves inside ring 3 = #moves inside rings 1-2 = k,
      and k >= 1 (every k = 1..8 occurs). Hall form: after deleting the 24 within-ring knight moves, the 16 squares
      X = W cap (ring 0 cup ring 3) of one colour have exactly the 16 neighbours B cap (ring 1 cup ring 2).
      6x6: #moves inside ring 2 = 4 + #moves inside the central 4x4 (>= 4).
  FINITE-EXACT companions: CP-SAT optimal ring profile of closed 8x8 tours (within-ring moves 1..16,
  two-ring jumps 2..32, m_23 17..40, m_13 2..24); the 2-wide frame (rings 2,3) has no knight tour, open or closed;
  the 8x8 minus its central 2x2 has closed tours; the 6x6 per-edge tour counts (9862 tours).
  NOT CLAIMED: (D) is classical folklore for two runners (re-derived); no LRC(14) consequence.
source: mac-mini-2026-10-06-chessboard (owner's 8x8 rings/scaffolds/sliders/knight-tour prompt)
depends_on: []
related:
  - 05-knowledge/results/chessboard_weave_20261006.md
  - 01-canon/theorems/THM-386-lrc-lonely-central-box-grounding.md
  - 01-canon/theorems/THM-3043-lrc-tight-instances-are-not-only-APs-and-what-FC2-does-not-transfer.md
  - 05-knowledge/hypotheses/HYP-9210-latest-first-lonely-time.md
  - 05-knowledge/hypotheses/HYP-9168-hp-blocking-number-equals-hall-bound.md (the Hall-obstruction genus)
scripts:
  - 04-computation/experiments/chessboard_weave_20261006_census.py
  - 04-computation/experiments/chessboard_weave_20261006_general_boards.py
  - 04-computation/experiments/chessboard_weave_20261006_rider_rings.py
  - 04-computation/experiments/chessboard_weave_20261006_tau2_check.py
  - 04-computation/experiments/chessboard_weave_20261006_weave_law.py
  - 04-computation/experiments/chessboard_weave_20261006_weave_values.py
  - 04-computation/experiments/chessboard_weave_20261006_nowithin_check.py
  - 04-computation/experiments/chessboard_weave_20261006_nowithin_witness.py
  - 04-computation/experiments/chessboard_weave_20261006_knight_cpsat8.py
---

# THM-4550 -- the owner's 8x8 board as a lonely-runner instrument, and the weave law

Full statements, proofs, censuses and the audit trail are in
[`chessboard_weave_20261006.md`](../../05-knowledge/results/chessboard_weave_20261006.md) (sections 1, 2, 5, 6, 10).

**(B) Rings are what the bishop sees.** On the `N x N` board, a bishop on ring `k` (Chebyshev shell about
the centre) reaches `2N - 3 - 2k` squares; the rook reaches `2(N-1)` from everywhere. With
`lambda = (N - 1 - 2k)/(2N)` the two-runner loneliness `min(||x||, ||y||)` of the square's centre (the board
read as one cell of the torus), `bishop reach = (N - 2) + 2N lambda`. So the rings are exactly the level
sets of the lonely-runner function and of the diagonal sliding range.

**(K) The knight's ring law.** Every knight move changes the ring by 0, 1 or 2; on `2m x 2m` the counts are
`8(m-1)`, `16 C(m,2)`, `16 C(m-1,2)` (eight inside each ring `k >= 1`, at its corners).

**(D) Two speeds.** For coprime `a, b`, `delta(a,b) = floor((a+b)/2)/(a+b)`: constant on anti-diagonals,
`1/2` exactly for colour-preserving leapers (`a + b` even), `1/2 - 1/(2(a+b))` otherwise; the knight is the
unique tight rider and, on the 8x8 subdivision, the only rider whose line never meets the closed central four
squares.

**(T) First passage.** The latest first lonely time of two runners is `4/9`, attained only by the camel `(1,3)`.

**(S) The scaffolds are the extremes.** The odd lengths are lonely at `1/2`; the even lengths are the tight
progression `2{1..m}`; their weave `{1..2m}` is tight. The progression's lonely times at level `1/(n+1)` are
exactly `k/(n+1)`, `gcd(k, n+1) = 1`.

**(W) The weave law.** Rings 0 and 3 hold `4 + 28 = 32` squares, rings 1 and 2 hold `12 + 20 = 32`; the only
knight moves inside rings 0 ∪ 3 are the eight corner moves of ring 3. Degree counting gives
`#(moves inside ring 3) = #(moves inside rings 1-2) = k` in every closed tour, and `k >= 1` because with `k = 0`
every move would flip both the colour and the ring pair, locking colour to ring pair on the whole board. Every
`k = 1..8` occurs. Hall form: after deleting the 24 within-ring knight moves, the 16 squares of one colour in
rings 0 ∪ 3 have exactly the 16 neighbours of the other colour in rings 1 ∪ 2.

**Audit.** The session's adversarial auditor re-derived (B), (K), (D), (S), (T), (W), the Hall form, the 6x6
analogue and the CP-SAT ring profile with independent code; no claim was refuted; the corrections are listed
in the note's section 10 and in MISTAKE-564.
