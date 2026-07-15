---
id: HYP-6840
title: The r=8 deck frontier — exact token/event-word criterion, persistent seven-stalk refuter to raw wall bounds, and the corrected owner-switch problem
status: UPDATED after MISTAKE-147 — full blocking on arbitrarily long core-safe wall runs is CONSTRUCTED; full blocking over an entire core-safe set remains OPEN; the global 1/7 law is VERIFIED but insufficient after core restriction
source: opus-2026-07-14-S301 (owner directive: progress or rigorous no-gos; separate signal from noise)
depends_on:
  - THM-767   # zero-variance, chamber locking (the mechanisms that resist covers)
  - THM-771   # the exact seven-owner defect frame
related: [THM-761, THM-777, THM-778, THM-779, THM-783, THM-784, THM-786, HYP-6835, MISTAKE-146, MISTAKE-147, MISTAKE-148]
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

## Sharp open questions — corrected after THM-779/783/784/786

- **Q0 (local existence): RESOLVED POSITIVELY.** THM-784 gives a full blocking
  interval `J=(5/16,7/20)` for the eight owners
  `{1,2,3,4,5,8,10,560N+1}` and `21N` consecutive covered walls. The seven
  slow owners alone form a fixed rainbow there. This does not yet make `J` a
  safe component of a five-speed LRC core, so the original **core-incidence**
  question remains open.  The stronger divisor-complete family
  `P={1,2,11,12,13}`, `W={1,4,5,6,8,9,10,182m+1}` places `2m` covered walls in
  a genuine core-safe interval, but still does not cover the whole safe set.
- **Q1 (whole-core existence): OPEN.**  THM-779 makes full blocking over any
  fixed `J` integer-decidable (piece surjectivity + wall rainbow + no
  simultaneity), but the empirical ceiling `K0=5` was false as a universal
  statement.  There are now exact fully blocked core-safe intervals with
  arbitrarily many walls: `P={1,2,11,12,13}`,
  `W={1,4,5,6,8,9,10,182m+1}`, and
  `J=[25/182,27/182]` give `2m` covered walls.  Whether one packet blocks the
  **entire** closed core-safe set is still open.
- **Q2 (positive theorem): REFRAMED AND PARTLY PROVED.**  Factor out intervals on which seven
  owners already form a persistent exact stalk, and bound genuine owner
  switches or centered-word complexity.  THM-779's exact redundancy cocycle
  gives the supportability equation; THM-778 supplies the centered-Christoffel
  event word and proves a simultaneous-wall pierce when equal 2-adic valuation
  has sufficiently large gcd.  THM-783 supplies the balanced-cluster laws;
  THM-786 proves the metric-extent pierce on its no-co-landing and sparse
  classes.  Its dense balanced-swap regime is the remaining universal case.
  THM-784 independently supplies a simpler unbounded raw-wall
  family, so this failure is not tied to divisor completeness.
  Raw wall density is not a proof coordinate.
- **Q3 (finite group / DMNR):** the rooted covered-state space is an `A8`
  torsor of size `20,160`, with one strongly connected component and
  monochromatic seven-cycles.  Connect the owner-switch word, after stalk
  factoring, to exact-covering-system rigidity.  The global `1/7` buffering
  count must be restricted to the core-safe exit set before it can govern this
  question.

The information-boundary lesson is now exact: runner or wall-event tournaments
without metric and token-fibre labels do not preserve blocking. Fast refinement
changes only event resolution. The appropriate quotient must retain the
owner-labelled chronological path, slow chamber, mesh ratios, and core-safe
base interval.

## Signal/noise ledger for this session (S301)

- SIGNAL: THM-777's floor (7/858 at {1..13}∖{6}, the detuning extremal — the
  double-extremal phenomenon); the universal 1/7 buffering law; the rainbow
  obstruction quantified; regime-2's normalized boundedness on floored strata.
- NOISE (removed): the raw r_P complementarity (MISTAKE-145, twice-refuted); the
  q ≤ 25 uniform finish (codex refutation, S300); "r=8 kills the deck method"
  (this file); my w_a/(7d) buffering guess (referee-corrected same-session).
