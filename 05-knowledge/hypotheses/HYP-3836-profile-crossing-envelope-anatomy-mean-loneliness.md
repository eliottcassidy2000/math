---
id: HYP-3836
title: Profile crossing = the r-side scale-matching law; spread sets carry the envelope at small r, tight sets at the cusp; mean loneliness E(S) does not select the AP
status: STUB -- CLAIMED klein-2026-07-01-S87, computation in progress this session
source: klein-2026-07-01-S87
script: 04-computation/lonely_profile_farey_slope_klein.py (to be written this session)
related:
  - HYP-3834, HYP-3835 (the distribution function, the collapse constant)
  - HYP-3763 (scale-matching law: coarse cap -> additive energy, fine cap -> entropy)
  - HYP-3822 (facility game; (6/7)^13 equidistributed optimum)
  - kind-pasteur S25-S27 (sub-critical inf meas >= 1/36 census)
---

# HYP-3836 (STUB): profile crossing and the anatomy of the lonely envelope

## Claim (reserved, being tested)

1. **Crossing**: the profiles Lambda_S(r) of different S CROSS: near r=0 spread sets
   (pairwise-coprime, large speeds) have SMALLER Lambda (union bound 1-2kr asymptotically
   tight; shared fractions force overlap for structured sets), while near r=1/14 the tight
   AP/GW are the minimizers. The lonely envelope inf_S Lambda_S(r) is therefore carried by
   DIFFERENT families at different radii -- the r-side geometrization of the scale-matching
   law (HYP-3763): no single S, and no single-radius functional, is extremal at all scales.
2. **Mean loneliness**: E(S) = int_0^1 m_S(t) dt = int_0^infty Lambda_S(r) dr (layer-cake).
   The AP does NOT minimize E among 13-sets (E is a coarse/averaged functional; predicted
   by the scale-matching law). Identify the actual small-E families.
3. **Envelope anatomy at small n**: exhaustive small-n (n=4,5,6), bounded speeds:
   the envelope is a piecewise-linear function whose carrier family changes finitely many
   times between r=0 and the cusp 1/n; locate the crossover radii.

## Why it matters

kind-pasteur's sub-critical floor (1/36 at r ~ 0.906/14) and mac-mini's cusp slope
(0.26(1-14r)) are two POINTS on one object -- the envelope. Making the envelope explicit
turns "inf meas" claims into claims about a single curve and its carrier families, and
explains why coarse instruments (moments, energy, E) rank differently from fine ones.

## Status

STUB. Computation in progress this session.
