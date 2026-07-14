---
id: HYP-6840
title: The r=8 deck no-go hunt — three-stage honest negative, the universal 1/7 buffering law (referee-corrected), and the flipped conclusion: full deck blocking at r=8 is as rigid as the r=7 tilings
status: RESOLVED-NEGATIVE for the hunted object (no full-coverage certificate found; three independent obstruction mechanisms identified and quantified); the buffering law VERIFIED; the sharp existence question OPEN and named
source: opus-2026-07-14-S301 (owner directive: progress or rigorous no-gos; separate signal from noise)
depends_on:
  - THM-767   # zero-variance, chamber locking (the mechanisms that resist covers)
  - THM-771   # the exact seven-owner defect frame
related: [THM-761, THM-777, HYP-6835, MISTAKE-146]
verification: 04-computation/lrc14_r8_single_lens_nogo_opus_S301.py,
  04-computation/lrc14_r8_static_deck_cover_search_opus_S301.py (+ .outs)
---

# HYP-6840 — the r=8 deck no-go hunt

## What was hunted, and why

THM-767's residual said the single-event pierce "fails structurally" at r = 8
(counts stay ≥ c), and a one-moment surviving cover was realized. The natural
no-go conjecture: there exist covering families V = 7P ⊔ W, |W| = 8, whose deck
is FULLY blocked over the whole closed core-safe set — an exact certificate that
the lens-7 method can never fire, bounding the method rigorously. This session
hunted that certificate three ways. All three failed, each for a QUANTIFIED
structural reason — and the failures flip the conclusion.

## The three-stage honest negative

1. **Measure-greedy, overlapping core** ({1,2,3,4,6}: |Ḡ_P| = 4/7): six free
   slots stall at 0.088 uncovered. Diagnosis: target too large (core overlaps
   inflate the safe set).
2. **Measure-greedy, spread core** ({2,5,9,11,13}: |Ḡ_P| = 12073/30030 ≈ 0.402,
   22 components): seven free slots + 4000 pair repairs stall at 0.067 uncovered
   in 662 dust pieces. Diagnosis: the endgame is FRAGMENTED DUST — each new
   exception's arcs land ~6/7 on covered ground; finishing greedily would need
   ~40 slots, not 8. Greedy/random construction CANNOT reach full coverage.
3. **Exact combinatorial search (the right formulation):** at c = 7 with coprime
   exceptions, zero-variance (THM-767) pins each exception to EXACTLY ONE sheet
   per chamber piece, so full coverage over Ḡ_P is equivalent to a finite
   condition: on every chamber piece (1,164 pieces over 22 safe components, 1,142
   interior walls, candidates w ≤ 90) the eight assigned sheets hit all of Z_7,
   AND at every interior wall the remaining seven form a PERFECT RAINBOW.
   Annealed search (24 restarts × 1200 steps): best subset violates 896
   constraints (832 piece + 64 wall). Nowhere near. Diagnosis: with 8 exceptions
   on 7 sheets, expected pairwise collisions per piece = C(8,2)/7 = 4; demanding
   ≤ 1 collision on EVERY piece simultaneously requires globally anti-correlated
   collision patterns — an exact-covering-system-grade rigidity.

## The universal 1/7 buffering law (VERIFIED; referee corrected the derivation)

For exceptions a, b at lens c (7 ∤ w's), partner b "buffers" an exit of a (keeps
the eventing sheet covered across a's wall) iff b's open bad set contains the
exit's u-point. Exact count over one window, verified on random pairs w ≤ 400:

> **buffered(a ← b) = w_a/7 + O(gcd(w_a, w_b)) — the fraction is 1/7 for EVERY
> partner, regardless of gcd.**

(The session's first derivation predicted w_a/(7d) — the exact referee corrected
it: the step-14w_b AP takes w_a/d values in the window but each with multiplicity
d; the d cancels. Caught before canon; recorded here per the
arithmetic-vs-governance lesson of MISTAKE-146.)

Consequence: seven partners buffer at most 7·(w_a/7 + O(d)) = w_a + O(Σd) of a's
w_a exits — **full buffering forces the seven buffer sets to be NEAR-EXACTLY
DISJOINT, simultaneously for all eight exceptions**: eight simultaneous
near-exact 7-fold AP-partitions of exit streams. This is the deck's version of an
exact covering system (DMNR territory), i.e. the same rigidity class as the r = 7
maintained tilings that chamber-locking forbids.

## The flipped conclusion (the signal)

The hunted no-go — "the deck method is bounded at r = 8" — is NOT supported.
Instead: **full deck blocking at r = 8 demands a rigid arithmetic object (the
eight-fold buffer partition) that generic, greedy, and annealed constructions
provably-in-practice cannot produce, and whose existence is the sharp open
question.** Three consequences:

1. The deck route plausibly EXTENDS past r = 7: for families outside the rigid
   packet class, free sheets over core-safe times should exist at r = 8 — the
   quantitative target is a union-of-buffers argument giving a positive fraction
   of unbuffered exits, each an almost-pierce.
2. The rigid packet class, if nonempty, is the deck's AP/GW corner — classify it
   (HYP-6835's chamber census is the right machinery; DMNR/covering-system tools
   apply to the buffer partitions).
3. NOISE removed: "r ≥ 8 is where the deck method dies" (my own S300 framing of
   the residual) was premature — the realized one-moment survivor does NOT
   extend to a full blocking; no full blocking was constructible under three
   escalating attacks.

## Sharp open questions (named for the fleet)

- **Q1 (existence):** does ANY covering 13-family have its lens-c deck fully
  blocked over the closed core-safe set with r = 8? Equivalent finite form at
  c = 7, small exceptions: an 8-subset whose sheet table is rainbow at every
  wall and 7-covering on every piece. SAT/CP over larger candidate pools (w up
  to a few hundred) would settle small strata exhaustively.
- **Q2 (the positive theorem):** prove that ≥ (1 − 7·(1/7))-type buffer counting
  plus disjointness defect δ > 0 forces an unbuffered exit in every event window
  — an r = 8 pierce for all families outside an explicitly-described packet
  class.
- **Q3:** connect the buffer partitions to DMNR: an exact 7-fold partition of an
  exit stream by partner-APs has moduli 14·d_ab-structured; do the two-largest-
  moduli-equal rigidity arguments transfer?

## Signal/noise ledger for this session (S301)

- SIGNAL: THM-777's floor (7/858 at {1..13}∖{6}, the detuning extremal — the
  double-extremal phenomenon); the universal 1/7 buffering law; the rainbow
  obstruction quantified; regime-2's normalized boundedness on floored strata.
- NOISE (removed): the raw r_P complementarity (MISTAKE-145, twice-refuted); the
  q ≤ 25 uniform finish (codex refutation, S300); "r=8 kills the deck method"
  (this file); my w_a/(7d) buffering guess (referee-corrected same-session).
