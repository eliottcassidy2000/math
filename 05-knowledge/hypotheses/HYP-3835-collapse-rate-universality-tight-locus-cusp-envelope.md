---
id: HYP-3835
title: Collapse-rate universality over the tight locus -- c(GW) = c(AP); tight => c(S) >= 2·#maximizers/(n·v_max); the cusp envelope is linear with the universal constant
status: STUB -- CLAIMED klein-2026-07-01-S87, computation in progress this session
source: klein-2026-07-01-S87
script: 04-computation/lonely_profile_farey_slope_klein.py (to be written this session)
related:
  - HYP-3834 (the distribution function + AP collapse constant)
  - THM-523 (tight census {AP, GW T5}; GW = {1..11,13,24})
  - HYP-3824 (inf meas at r=1/14 is 0; floor = the linear slope)
  - HYP-2893 (GW = Jacobsthal acceleration)
---

# HYP-3835 (STUB): collapse-rate universality over the tight locus

## Claim (reserved, being tested)

1. **Universality**: the Goddyn-Wong tight set GW = {1..11, 13, 24} has the SAME collapse
   rate as the AP: c(GW) = c(AP) = 1666/6435. Hand-derivation: both tight sets have exactly
   the six witnesses t = k/14 (gcd(k,14)=1), and at witness k/14 the binding speeds are the
   unit pair {k^-1 mod 14, 14 - k^-1} -- all in {1,3,5,9,11,13}, contained in BOTH sets
   (12, the only element they disagree on below 14, is even, hence never binding; 24 is
   never binding at any k/14). The cusp slope sees only the mediant unit pairs -- it is
   BLIND to which tight set realizes them.
2. **Tight floor**: every tight 13-set S has c(S) >= (#maximizers)·(2)/(14·v_max(S)) > 0
   (each maximizer contributes a lonely interval of length >= eps·(1/v+ + 1/v-)/14).
   With any bounded-speed reduction (Tao arXiv:1701.02048), c is uniformly bounded below
   on the tight locus.
3. **Cusp envelope**: conditional on the THM-523 census (tight locus = {AP, GW} up to
   dilation; dilation preserves Lambda), the lonely envelope satisfies
   inf_S Lambda_S(r) = 1666/6435 · (1 - 14r) + O((1-14r)^2) as r -> 1/14^-,
   and near-cusp it is carried by BOTH tight families simultaneously (doubly attained).

## Why it matters

"The atom has no measure of its own" (HYP-3824): the content at the critical radius is the
APPROACH RATE. This hypothesis pins the approach rate as a single universal arithmetic
constant shared by all known tight sets -- the cusp cannot tell AP from GW. Any future
tight set must either contain the six unit pairs (same constant) or have witnesses at
deeper moduli 14d (different constant, detectable).

## Status

STUB. Hand-derivation above; exact rational verification (witness enumeration for GW at
moduli 14d, exact slopes) in progress this session.
