---
id: HYP-3841
title: The overtaking-defect sieve -- the tangent-ladder certificate proves M >= 1/14 for covering sets from one or two sub-critical anchors; the defect concentrates on the deep-well family; the final window is overtaking-free
status: CONFIRMED as an instrument (exact computation; 9/11 single-anchor, 11/11 two-anchor ladder incl. the deep well)
source: klein-2026-07-01-S88
script: 04-computation/overtaking_defect_sieve_final_window_klein.py (+ .out)
related:
  - THM-592 (mac-mini radius-derivative structure; the ladder bound made per-set here)
  - HYP-3834/3835 (distribution function, collapse law)
  - HYP-3900 (opus simultaneous-peel; complementary residual instrument)
  - definitions.md deep-well isolation (the K-defect concentrates exactly there)
---

# HYP-3841: the overtaking-defect sieve and the tangent-ladder certificate

## The instrument

Lambda_S(r) is piecewise linear (HYP-3834); kinks are gap-deaths at r = d/(v+w) (convex,
slope jumps UP) and overtakings at r = d/(w-v) (concave, slope jumps DOWN; THM-592).
Define the exposed concave defect K(S; [r0, r1]) = total downward slope jump. Then

    TANGENT-LADDER:  Lambda(r) >= Lambda(r0) + (Lambda'(r0+) - K)(r - r0)  on [r0, r1].

A per-set FINITE certificate for M(S) >= 1/14 from a sub-critical anchor: the triple
(Lambda(r0), Lambda'(r0+), K). Covering sets always contain a multiple of 14 (q-witness),
so all their shallow points k/14 are dead and a sub-critical anchor is exactly what is
available -- the sieve is the covering-side counterpart of the tight-side unit-residue
rigidity.

## Verified (exact rational arithmetic)

- Anchor r0 = 1/16, target 1/14: certificate succeeds on 9/11 test sets (8/8 random
  covering sets, K in [0.045, 2.77]; the CLUSTER set).
- The 2 failures are exactly the DEEP-WELL family: CONSTR {1..12,182} (cert -0.0067,
  actual 0.0239, K = 1.989 from ONE overtaking kink of the far element) and TWO-SCALE
  {...,28,182} (cert -0.0138, actual 0.0300, K = 2.200). The defect concentrates
  precisely on the known danger-zone family (deep-well isolation).
- CHAINED LADDER (anchors 1/16 -> 1/15 -> 1/14), pure certificate (step 2 starts from
  step 1's certified value, not the actual): **11/11**, CONSTR +0.00875, TWO-SCALE
  +0.01185.
- **The final window [1/15, 1/14] is overtaking-free (K2 = 0) in every tested set** --
  overtaking there needs w - v in (14d, 15d], a one-unit-per-d difference band.

## The reduction this opens

Uniform-ize three quantities over covering sets: (i) anchor floor Lambda(1/15) >= a > 0;
(ii) slope bound |Lambda'(1/15+)| <= b (co-area: sum over surviving gaps of 1/v_L + 1/v_R);
(iii) final-window defect K = 0 (band emptiness -- CANDIDATE LEMMA: no covering set has an
overtaking kink in (1/15, 1/14): needs pairs with w - v in (14d, 15d] AND the d-th gap
alive there). Then M >= 1/14 for ALL covering sets follows from a - b/210 >= 0. The
three-quantity uniformization is a concrete finite program; (iii) looks provable outright.

## Honest limits

Test battery is 11 sets (8 random + 3 adversarial), not a census; the uniform (i)/(ii)
bounds are open (they are close to kps's HYP-3950 arc-count floor and the THM-579 CV
frame). The instrument's value: it converts "covering set" into a 3-number certificate
and isolates WHERE hardness lives (the deep well's mid-band overtakings, all BELOW 1/15).
