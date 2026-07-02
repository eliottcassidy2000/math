---
id: HYP-3836
title: Profile crossing = the r-side scale-matching law; the envelope carrier ladder {1, m, m+1}; mean loneliness E(S) does not select the AP (GW beats it)
status: CONFIRMED (exact computation; exhaustive small-n envelopes)
source: klein-2026-07-01-S87
script: 04-computation/lonely_profile_farey_slope_klein.py (+ .out)
related:
  - HYP-3834, HYP-3835 (distribution function, collapse-rate law)
  - HYP-3763 (scale-matching law -- this is its r-side/measure avatar)
  - HYP-3822 (facility game; (6/7)^13 equidistributed value -- random sets land on it)
  - kind-pasteur S25-S27 (sub-critical inf meas >= 1/36: a point on this envelope)
---

# HYP-3836: profile crossing and the anatomy of the lonely envelope -- CONFIRMED

## 1. Profiles cross (n=14, exact table in .out)

At r = 0.001: PRIMES>13 (0.9747) < BLOCK 100..112 (0.9750) < RAND (0.9752) < GW (0.98289)
< CONSTR (0.98316) < AP (0.98322) -- spread sets cover best; union bound 1 - 26r = 0.974 is
asymptotically tight for them. At r = 0.05 the order has FLIPPED: AP (0.1752) < GW (0.1838)
< CONSTR (0.2316) < BLOCK < PRIMES. Crossover ~ r = 0.04. At r = 1/14: AP = GW = 0 <
CONSTR (4637/194040 = 0.0239) < BLOCK (0.1265) < RAND (0.129-0.139) < PRIMES (0.1519);
random sets bracket the equidistributed facility-game value (6/7)^13 = 0.1348 (HYP-3822).

Fine structure inside the tight pair: GW < AP for r <= 1/28 (the 24 spreads better
sub-critically), GW > AP on [0.05, 0.067], and both meet at slope 1666/6435 at the cusp --
the two tight sets EXCHANGE advantage twice before dying together.

## 2. The envelope carrier ladder (exhaustive small n)

n=4 (3 speeds <= 15): argmin Lambda walks {1,14,15} -> {1,10,11} -> {1,6,7} -> {1,3,4} ->
{1,2,3} (AP only from r ~ 9/40 to the cusp 1/4). n=5 (4 speeds <= 12): {1,9,10,11} ->
{1,5,6,7} -> {1,3,4,5} -> AP only AT the cusp. The carriers are {1} + consecutive blocks
sliding DOWN as r grows -- a mediant-like descent. At r = 19/100 (very near cusp 1/5) the
carrier {1,3,4,5} is LOOSE (M > 1/5): near-tight loose sets undercut the tight sets'
linear collapse until their isolation constant Lambda(1/n) > 0 bites -- exactly why the
envelope at the critical radius is THM-523's open inf-L lemma and not a tight-census fact.

## 3. Mean loneliness does not select the AP

E(S) = int_0^1 m_S dt = int_0^infty Lambda dr (layer-cake). Exact values (n=14):
E(GW) = 0.03046635 < E(AP) = 0.03051304 (**GW strictly beats AP**; difference -4.7e-5,
exact rationals in .out) < CONSTR (0.03262) < RAND (0.0354 ~ 1/28 = the iid value) <
BLOCK (0.04056) < PRIMES (0.04058). Exhaustive n=4 (speeds <= 15): global E-min is
{1,3,4} = 0.1095; AP ranks 35/409. E(AP_N) has the exact Farey closed form
(1/2) sum_{adjacent} 1/(qq'(q+q')) -- verified.

So the invariants SPLIT the tight locus in opposite ways: the collapse rate c is BLIND
within {AP, GW} at n=14 (both 1666/6435, HYP-3835) while E separates them (GW < AP);
conversely E cannot see tightness (loose {1,3,4} minimizes E at n=4). This is the
scale-matching law (HYP-3763) geometrized on the r-axis: coarse functionals (E, small-r
Lambda) prefer spread; only the cusp window prefers the tight sets. No single-radius or
averaged instrument ranks all of [0, 1/n] -- the ENVELOPE (a curve, not a number) is the
right object: kps's 1/36 (sub-critical), mac-mini's slope (cusp), THM-523's inf L (at the
critical radius) are three readings of it.
