---
id: THM-900  # renumbered from 899 (lattice-law-boxhit first-pushed)
title: THE GUY / SQUARE-TRIANGULAR WEAVE AND THE BOTH-CLEAN TIGHTNESS — (G) Guy's crossing-number formula at odd n is the SQUARE OF A TRIANGULAR NUMBER: Z(2k+1) = T_{k−1}² = 1³ + ⋯ + (k−1)³ (Nicomachus), verified n = 5..21 (1, 9, 36, 100, 225, 441, 784, 1296, 2025); Z(9) = 36 is the second nonzero SQUARE-TRIANGULAR number (T₈ = 6²) — the point where squares-of-triangulars and triangular-squares meet, and simultaneously Σ of the first three cubes; (Q) THE MIDPOINT-QUADRUPLE CLOSED FORM: the number of unordered additive quadruples ({a,b} ≠ {c,d}, a+b = c+d) in [n] is ⌊n(n−2)(2n−5)/24⌋ (95, 125, 161 at n = 12, 13, 14 — exact), and E(AP_n) = n(n+1)(2n+1)/3 − n² (1469 at n = 13, the S321 constant); the EXTREMAL ANALOGY: Guy's problem minimizes cyclically-interleaved quadruples (optimal = T²), the Sidon problem minimizes additive-midpoint quadruples (optimal = 0), and these are the quadruple-extremal problems on the program's two axes (cyclic/multiplicative vs additive); (S) THE BOTH-CLEAN TIGHTNESS: greedy search for 13 speeds jointly Sidon + ratio-clean (Y* ≥ 40, q,p ≤ 13) + no 13-multiple FAILS below 25000 — the constraints FIGHT (Sidon forces spread; spread makes the 78 pairwise ratios dense near small rationals; the clean ratio channels are the thin Farey-14 gaps): THE EMPTINESS CONJECTURE — the joint level-5 admissible regime is empty at accessible scales — with the FAREY-CHANNEL PACKING argument named as the proof route (13 log-positions, 78 differences, a small-measure clean set: pigeonhole); if proved, S328's "the prefix recursions are essential" becomes a theorem
status: G/Q PROVED-EXACT (identities verified n = 5..21 and n = 12..14; the closed forms exact); S = an honest search failure with the structural diagnosis and the named conjecture + proof route (greedy is not exhaustive: a backtracking/SAT search or the packing proof is the next step)
source: opus-2026-07-16-S329 (owner: the both-clean search; near-AP quadruples vs square-triangular numbers and Guy's conjecture formula)
depends_on:
  - THM-897 (the admissibility discovery this quantifies)
  - S321/THM-882 (the additive axis), THM-863 (the ratio axis), S315/THM-873 (the Farey-14 channels)
related: [Guy's conjecture (cr(K_n) = Z(n): proved for n ≤ 12), Nicomachus, Pell/square-triangulars (A001110), the S327 fertility/toothpick thread]
verification: 05-knowledge/results/bothclean_search_guy_weave_opus_S329.out
---

# THM-900 — the Guy / square-triangular weave and the both-clean tightness

## (G) Guy's numbers are squares of triangulars

Z(n) = ¼⌊n/2⌋⌊(n−1)/2⌋⌊(n−2)/2⌋⌊(n−3)/2⌋ (Guy's conjectured cr(K_n), exact
for n ≤ 12). At odd n = 2k+1 the four factors pair into k(k−1)·k(k−1)/4:

> **Z(2k+1) = (k(k−1)/2)² = T_{k−1}² = 1³ + 2³ + ⋯ + (k−1)³** (Nicomachus),

verified n = 5..21: 1, 9, 36, 100, 225, 441, 784, 1296, 2025. The crossing
count of the odd complete graph is a sum of cubes. **Z(9) = 36 = T₈ = 6²**:
the second nonzero square-triangular number — the unique small point where
"squares of triangulars" (Guy's sequence) and "triangular squares" (the Pell
sequence 1, 36, 1225, 41616, …) coincide, and K₉ is the smallest complete
graph whose crossing number is simultaneously square, triangular, and a
cube-sum of both parities. (n = 13, the LRC vertex count: Z(13) = 225 = 15²
= T₅² — the fifth triangular squared, sitting beside the S313 knife-edge
1001 = C(14,4) in the ledger.)

## (Q) The quadruple closed forms (the S328 blockers quantified)

- E(AP_n) = n(n+1)(2n+1)/3 − n² (= 1469 at n = 13: S321's constant).
- Unordered midpoint quadruples of [n]: **⌊n(n−2)(2n−5)/24⌋** (95, 125, 161
  at n = 12, 13, 14): the exact count of the S₄-blockers a near-AP packet
  carries.
- The extremal analogy: Guy minimizes CYCLIC interleaving of quadruples
  (optimum T_{k−1}²); Sidon minimizes ADDITIVE midpoint quadruples (optimum
  0); these are the quadruple-extremal problems on the program's two axes,
  and the level-5 admissibility (THM-897) demands the Sidon side while the
  cyclic world of tournaments runs on Guy's side.

## (S) The both-clean tightness and the emptiness conjecture

Greedy construction of 13 speeds that are simultaneously Sidon, ratio-clean
(Y* ≥ 40 for q, p ≤ 13), and 13-multiple-free FAILS below 25000 (two step
pools). Diagnosis: the constraints fight — Sidon forces spread, spread
pushes the 78 pairwise ratios toward the dense small-rational web, and the
clean channels are only the thin Farey-14 gaps (S315). **Conjecture (the
emptiness of the joint regime at accessible scales):** no 13-packet below
[explicit scale] satisfies both cleanliness conditions; proof route = the
Farey-channel packing pigeonhole (13 log-positions ⟹ 78 differences vs a
clean set of small measure). If proved, THM-897's admissibility discovery
upgrades to: **the full-problem level-5 certificate is unreachable at small
scales, and the prefix decompositions are PROVABLY essential** — the ladder
does not collapse the recursion; it explains it.

## Named next

- The packing proof (or a SAT/backtracking search settling nonemptiness).
- The Y*-threshold trade-off curve (emptiness scale as a function of the
  cleanliness demands — the quantitative version).
- The Guy-side probe: tournaments on 2k+1 vertices vs the cylindrical
  optimal drawings (does the metagraph see Z(n)? the cyclic-quadruple count
  as a class invariant — one exact computation at n = 7).
