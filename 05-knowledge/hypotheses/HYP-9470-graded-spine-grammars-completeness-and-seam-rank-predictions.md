---
id: HYP-9470
title: "Graded spine grammars G(k, w). (A) Pre-registered census predictions P1-P5 (commit 756414ff80): G(0,2) misses between 10 and 60 odd numbers <= 300 besides 7, 21; G(0,3) misses at most 5 below 300 and at least one below 2000; G(1,1) contains 13 and misses at most 10 percent below 2000; G(1,2) misses nothing below 2000 besides 7, 21; 13 has minimal grades (0,2) and (2,0). (B) Conjecture: G(2,1) (two apexes over an interval chain, which contains the one-apex one-arc joins) attains every odd number except 7 and 21; equivalently the explicit N-rational series of THM-4878 on G(2,1) has range the odd numbers minus {7, 21}; this would prove Chan-Pak 6.2. (C) Conjecture: no grade G(0, w) without apexes is onto. (D) The two-reversed-arc seam rank 13 and the three-arc seam rank 62 are universal over all pairs of tournaments joined across two (three) reversed arcs."
status: PENDING (census outcomes recorded in the source note, section 6; (B) and (C) OPEN; (D) FINITE-EXACT within cutwidth <= 3 ambients, universality OPEN)
session: opus-2026-11-04-S79
date: 2026-11-04
related: [THM-4878, THM-4875, THM-4748, THM-4744, THM-4865, THM-4876, THM-1370, HYP-9468]
script: 04-computation/experiments/spine_grammar_census_20261104.cpp; 04-computation/experiments/spine_grammar_ranks_20261104.py; 04-computation/experiments/spine_grammar_prereg_20261104.md
---

# HYP-9470 - graded spine grammars: pre-registered census predictions, completeness conjectures, seam ranks

## Setting

See THM-4878 and the source note `05-knowledge/results/spine_grammars_cutwidth_apexes_hankel_seams_20261104.md`.
G(k, w) is the class of tournaments with k <= 2 apexes (x_1 -> x_2, each beating an arbitrary subset of the chain)
over a transitive chain with reversed arcs of cutwidth <= w. H on G(k, w) is an N-rational series over the finite
cut-state automaton of THM-4878 (1). "Attained by G(k, w)" includes ordered joins of members (products of counts).

## (A) Pre-registered predictions (written before the runs; `spine_grammar_prereg_20261104.md`)

- P1. G(0,2): besides 7 and 21, between 10 and 60 of the 150 odd numbers <= 300 are missed, and at least 3 percent of
  the odd numbers in [2^11, 2^12) are missed.
- P2. G(0,3): at most 5 odd numbers <= 300 are missed besides 7 and 21, and at least one odd number <= 2000 is
  missed besides 7 and 21.
- P3. G(1,1): 13 is attained, and at most 10 percent of the odd numbers <= 2000 are missed besides 7 and 21.
- P4. G(1,2): no odd number <= 2000 is missed besides 7 and 21.
- P5. 13 lies in G(0,2) and G(2,0) but in neither G(1,0) nor G(0,1).

Outcomes: PENDING (note, section 6).

## (B) Completeness of G(2,1) (OPEN; implies Chan-Pak)

Every odd number other than 7 and 21 is H of a tournament in G(2,1). Evidence: the one-apex one-arc joins are G(2,1)
words (THM-4878 (2)) and attain every odd number <= 2^32 except 7 and 21 (THM-4875 (D)); G(2,0) alone reaches 230,519
(Session 191). The statement is the range question for one explicit N-rational series (2,544 cut states; Hankel rank
in the note). A proof needs an all-N construction inside the grammar; Session 194's periodic-counter theorem shows
that fixed repeated blocks of length <= 6 in G(2,0) cannot be programmed to an arbitrary target, so the construction
must vary its blocks with the target (as the one-apex joins do through the Lorentzian law).

## (C) No apex-free grade is onto (OPEN)

For every fixed w, G(0, w) misses infinitely many odd numbers. Evidence: the censuses of G(0,1), G(0,2), G(0,3)
(note, section 6). A counting argument does not decide this (reduced words with H <= X number X^(O(1)) for every w),
so the obstruction, if real, is arithmetic.

## (D) Seam ranks (FINITE-EXACT within cutwidth <= 3; universality OPEN)

The composition law across a cut crossed by c reversed arcs has minimal rank 3, 13, 62 for c = 1, 2, 3 (apex-free
sides of cutwidth <= 3). Prediction: these ranks do not depend on the ambient cutwidth (as 3 and 13 did not between
G(0,2) and G(0,3)), i.e. 13 is the rank of the two-reversed-arc composition law over all tournament pairs, and the
G(0,4) computation gives the c = 4 term of the sequence 3, 13, 62, ... (result in the note if it finished).
