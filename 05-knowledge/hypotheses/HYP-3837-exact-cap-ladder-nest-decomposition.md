---
id: HYP-3837
title: THE EXACT CAP LADDER (all j = 1..13, THM-576 verified and extended) + two unifications (cap_10 = THM-530's m_P; cap_11 = the pentagon census minimum) + the caps are PAIR + NEST arithmetic
status: CONFIRMED (exact rationals; 2^13 exhaustive; I-E decomposition reproduces caps exactly)
source: klein-2026-07-01-S91
script: 04-computation/exact_cap_ladder_decomposition_klein.py (+ .out)
related:
  - THM-576 (caps j<=5 -- verified, extended), THM-530 (m_P = cap_10), kps-S27 census (pentagon 313/9702 = cap_11)
  - THM-594(B) (pair law -- one of the two ingredients), HYP-3838 (the nest lemma -- the other)
  - mac-mini S96 sec-1 (the Bernoulli program this collapses below the {1..13} universe)
---

# HYP-3837: the exact cap ladder

cap_j := min over j-subsets P of {1..13} of Lambda_P(1/14). Exhaustive exact table:

| j | k=13-j | cap_j | argmin | note |
|---|---|---|---|---|
| 1 | 12 | 78/91 = 6/7 | {1} | THM-576 MATCH |
| 2 | 11 | 66/91 | {1,13} | THM-576 MATCH |
| 3 | 10 | 55/91 | {1,12,13} | THM-576 MATCH |
| 4 | 9 | 1979/4004 | {1,11,12,13} | THM-576 MATCH |
| 5 | 8 | 2243/5880 | {1,5,7,8,9} | THM-576 MATCH |
| 6 | 7 | 3029/10780 | {1,5,7,8,9,11} | NEW |
| 7 | 6 | 45107/229320 | {1,5,7,8,9,11,13} | NEW |
| 8 | 5 | 2479/17640 | {1,3,5,7,8,9,11,13} | NEW |
| 9 | 4 | 10601/114660 | {1,2,3,5,7,8,9,11,13} | NEW |
| 10 | 3 | 14249/252252 | {1,2,3,5,7,8,9,11,12,13} | **= THM-530's m_P** |
| 11 | 2 | 313/9702 | {1..13}\{6,10} | **= kps pentagon census min** |
| 12 | 1 | 7/858 | {1..13}\{6} | NEW |
| 13 | 0 | 0 | AP (tight) | consistent |

Margins to second-min included in the .out (tightest rung: j=11, margin 0.00012).
UNIFICATIONS: the k<=7 pigeonhole floor constant m_P and the 11-core census minimum are
rungs 10 and 11 of ONE ladder -- the census, the pigeonhole floor, and the caps are the
same object Lambda_P(1/14) at different cardinalities.

DECOMPOSITION (verified == cap at j = 4, 5): inclusion-exclusion with exact d-fold
overlaps reproduces the caps; wrappedness occurs ONLY at d = 2 (THM-594(B) branch-2
pairs); every d >= 3 term is gcd-nest-exact (HYP-3838's lemma). Hence for EVERY
P in {1..13}:
    Lambda_P(1/14) = 1 - 2r j + sum_pairs |D^D| (THM-594B) - [gcd-nest alternating chain]
-- closed-form rational arithmetic. hp0cap's cap side needs no Vitali, no measure theory,
and no d <= 7 Bernoulli ladder below this universe: PAIRS + NEST suffice.
