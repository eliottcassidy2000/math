---
id: HYP-9211
title: "Knight-torus tour blocking: for every n >= 5, the fewest moves whose removal leaves the n x n knight torus without a closed knight tour is 7, and the minimum sets are exactly the 8n^2 stars (all but one move at a square)"
status: >
  OPEN. FINITE-EXACT (stars only) for n = 5, 6, 7, 8 (census in
  04-computation/experiments/sixseven_20261006_knight_torusblock_n.out). Upper bound 7 trivial (a star);
  the 2-factor (Hall) bound is 7 for even n by Ore's bipartite f-factor count (d >= 6(|S|-|T|)+1); its
  obvious tight witnesses (T empty, or T a colour class minus a square) both give the star, and with restricted
  edge connectivity 14 (n = 5..8, 10, 12) no other 7-set kills every 2-factor for even n in {6, 8, 10, 12}.
source: mac-mini-2026-10-06-sixseven, 05-knowledge/results/sixes_and_sevens_20261006.md, section 1
related:
  - 01-canon/theorems/THM-4552-exceptional-knight-tori-six-is-a-two-place-product-seven-is-a-paley-torus.md
  - 05-knowledge/hypotheses/HYP-9168-hp-blocking-number-equals-hall-bound.md
  - 05-knowledge/results/chessboard_weave_20261006.md
---

# HYP-9211 — the knight torus is killed only by stripping a square

Let `G_n = Cay(Z_n^2, {(±1, ±2), (±2, ±1)})`. Let `β(n)` be the fewest edges whose deletion leaves no
Hamiltonian cycle.

**Conjecture.** For every `n >= 5`, `β(n) = 7`, and every 7-edge blocking set is a star: 7 of the 8 moves at
one square.

## Evidence

* **Exhaustive, through a fixed move `e0` (edge-transitivity).**

  | `n` | 6-sets leaving a Hamiltonian cycle | blocking 7-sets | of which stars |
  |---|---|---|---|
  | 5 | all 71,523,144 | 14 | 14 |
  | 6 | all 464,306,843 | 14 | 14 (chessboard session) |
  | 7 | all 2,231,243,664 | 14 | 14 |
  | 8 | all 8,637,487,551 | 14 (of 359,895,314,625) | 14 |
* **Independent method.** A CP-SAT lazy-cut hitting-set model gives the same answer at `n = 5`.

## Why it matters

The owner's sentence ("on the 6x6 torus ... 7, and the only way is to strip one square") is the case `n = 6`. The evidence says
the sentence is about the knight's degree (`8 - 1`), not about the board size. It is the knight's member of the
Hall genus (HYP-9168 for tournaments: the cheapest way to kill every Hamiltonian object is to starve a set).

## Proof routes

* **Robustness.** Show that `G_n` stays Hamiltonian after any 6 deletions.
  * This is the "(k − 2)-edge-fault-tolerant Hamiltonian" pattern for a k-regular graph (`k = 8`). It is known for hypercubes `Q_k` (Latifi–Zheng–Bagherzadeh 1992, CITED).
  * Connected Cayley graphs on abelian groups are Hamilton-connected or Hamilton-laceable (Chen–Quimpo 1981, CITED).
  * We have not checked whether fault tolerance `k − 2` is known for all abelian Cayley graphs. If it is, the lower half of HYP-9211 follows.
* **Exclusivity of stars.** Show that a 7-set not concentrated at a square leaves both a 2-factor and enough connectivity to merge its cycles. This is the knight analogue of the merge count in the HYP-9168 partial theorem `β >= min(hall, 3)`.
* **Odd `n`.** `G_n` is not bipartite, so Ore's count must be replaced by Tutte's `f`-factor theorem.
