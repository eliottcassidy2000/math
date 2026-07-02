---
id: HYP-3838
title: The GCD-NEST LAW + the MEDIANT-IN-THIRD-DANGER criterion -- |^D_Q| = 2r/max(Q/gcd Q) except when a pair-mediant point lands inside another speed's danger (onset (1,14,15) via 2/29); THE CAP-UNIVERSE NEST LEMMA: zero violations on ALL 8100 d>=3 subsets of {1..13}
status: CONFIRMED (exact; exhaustive on the cap universe; mechanism identified at every boundary case)
source: klein-2026-07-01-S91
script: 04-computation/exact_cap_ladder_decomposition_klein.py parts C/C' (+ .out)
related:
  - THM-594(B) (the d=2 law this extends), HYP-3848 (the origin-nest mechanism, S90)
  - HYP-3837 (the ladder this powers), kps HYP-3954 (c-averaged subtorus ledger -- homogeneous complement)
  - MISTAKE-093 / THM-596 bands (the SAME mediant 2/29 governs the final window and the d=3 onset)
---

# HYP-3838: the gcd-nest law and the mediant criterion

## The law

For a speed set Q with g = gcd(Q): |^_{v in Q} D_v(r)| = 2r / max(Q/g) -- the origin
nest after gcd reduction (dilation invariance folds the gcd; the fastest reduced speed's
0-interval nests inside all others) -- UNLESS a wrapped coincidence survives. Census at
r = 1/14, all 1140 triples <= 20: the law holds in 1061; all 79 violations are primitive
with the identified mechanism:

## The mediant criterion (every boundary case verified)

A triple wraps beyond the nest iff a COINCIDENCE POINT of one pair (the (v,w)-mediant
family, |aw - bv| small) lies INSIDE the third speed's danger. Onset cases: (1,14,15)
wraps because the (14,15)-mediant 2/29 = 0.0690 sits inside D_1 = (-1/14, 1/14) --
**the same mediant 2/29 that governs the final window** (MISTAKE-093, THM-596 bands);
(2,13,15) wraps via the (13,15) coincidence at ~0.4948 inside D_2's arc at 1/2. For the
family (1, v, v+1): wraps iff 2/(2v+1) < 1/14 iff v >= 14 -- exactly the observed onset.

## THE CAP-UNIVERSE NEST LEMMA (the payoff)

Within {1..13} at r = 1/14, pairwise sums are <= 25, below every mediant-in-danger
threshold: **all 8100 subsets Q with |Q| >= 3 satisfy the gcd-nest law -- zero
violations (exhaustive).** Consequence (with THM-594(B)): Lambda_P(1/14) is closed-form
rational for every P in the cap universe; the caps, the m_P floor, and the pentagon
census minimum are pair + nest arithmetic. mac-mini S96 sec-1's d <= 7 Bernoulli ladder
is needed only ABOVE the universe (speeds > 13 / other radii) -- below it, the
inclusion-exclusion truncates exactly at d = 2. Proof route for the lemma (not yet
written formally): the mediant criterion + the threshold arithmetic (pair sums <= 25 <
28) -- finite and elementary, a THM candidate.
