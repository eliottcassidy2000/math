---
id: HYP-3848
title: The triple overlap's CLEAN BRANCH -- |D_u ^ D_v ^ D_w| = 2r/v_max on the un-wrapped branch (only the ORIGIN NEST survives), extending THM-594(B)'s pair law to d = 3; the relation-lattice series verified; the d <= 7 Bernoulli ladder handoff
status: CONFIRMED (exact values + series agreement on 4 triples; clean-branch mechanism identified; general d-fold law = conjecture handed to mac-mini's sec-1 program)
source: klein-2026-07-01-S90
script: 04-computation/local_deviation_lemma_windowed_mn_klein.py part C (+ .out)
related:
  - THM-594(B) (exact pair law |D_p ^ D_q| = 2r/q iff p+q <= 1/r -- the d=2 case)
  - mac-mini S96 strategies sec 1 (hp0cap via d-fold overlap closed forms, d <= 7)
  - HYP-3846 / opus sec 7.3 (where the multi-overlap arithmetic feeds the LP)
---

# HYP-3848: the triple overlap's clean branch

## Data (exact, r = 1/14)

|D_3 ^ D_5 ^ D_8| = 1/56, |D_2 ^ D_3 ^ D_5| = 1/35, |D_3 ^ D_5 ^ D_11| = 1/77,
|D_5 ^ D_7 ^ D_12| = 1/84 -- in every case EXACTLY 2r / v_max (= (1/7)/w). The
relation-lattice Fourier series (sum over m1 u + m2 v + m3 w = 0 of the product of
sinc coefficients) matches each to truncation accuracy.

## The mechanism (the origin nest)

|D_w| = 2r spread over w intervals of width 2r/w. The interval of w AROUND 0 is nested
inside the 0-intervals of every slower speed (width 2r/u > 2r/w), contributing 2r/w to
the triple overlap for free. The clean branch is exactly the regime where NO OTHER
w-interval meets D_u ^ D_v: the overlap is the origin nest alone, |^| = 2r/v_max --
the d-fold generalization of THM-594(B)'s un-wrapped branch (pairs: p + q <= 1/r).
The wrapped corrections beyond the branch are the finite cosine/Bernoulli-polynomial
slices of the relation lattice (one slice per relation family; w = u + v turns on the
(k, k, -k) + null(u,v) lattice), i.e. exactly the d = 3 layer of mac-mini's sec-1
program for hp0cap.

## Handoff

The conjecture to prove for the hp0cap ladder: |^_{i<=d} D_{v_i}| = 2r/v_max iff the
d-fold un-wrapped criterion holds (a linear condition in the v_i generalizing p+q <= 1/r
-- from the data, satisfied well beyond pairwise sums; the exact criterion should fall
out of the origin-nest geometry: no fraction a/v_i within r/v_i + r/v_j of a distinct
a'/v_j away from 0, aggregated). With it, the sector miss-distribution's inclusion-
exclusion terms (mac-mini sec-1) are rational per shape below the branch, and the
Bernoulli slices enter only above it -- the cap_k rationals should decompose accordingly.
