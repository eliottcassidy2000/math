---
id: HYP-3835
title: The COLLAPSE-RATE LAW -- c(S) = class-max harmonic sum; the AP MAXIMIZES the collapse rate; universality at n=14 is an even-class accident of the GW family; the cusp envelope is carried by the MINIMAL-c tight set
status: PARTIALLY-TRUE as originally stated (universality REFUTED as a general law at n=5,8) -- REFINED into the class-max law, VERIFIED exact on the full small-n tight censuses (n=4..8) + n=14, 20
source: klein-2026-07-01-S87
script: 04-computation/tight_census_collapse_rates_klein.py (+ .out); 04-computation/lonely_profile_farey_slope_klein.py
related:
  - HYP-3834 (distribution function, AP collapse constant)
  - HYP-3750 (klein-S48 tight classification: diff-closed / GW-type / CROSS-type -- the collapse rate is the HARMONIC DUAL of its lift data)
  - THM-523 (tight census {AP, GW} at n=14; the open "inf L > 0" lemma = the envelope AT r=1/14)
  - HYP-3824 ("AP is the minimizer" -- refined: true at n=14 only because GW ties it)
  - HYP-2893 (GW = Jacobsthal accelerations)
---

# HYP-3835: the collapse-rate law -- REFINED AND CONFIRMED

## The law (mini-theorem chain; elementary proofs sketched, every instance verified exactly)

Let S be TIGHT (M(S) = 1/n) with **no element divisible by n**. Then:

1. **(Shallow-witness rigidity.)** At t = k/n, m_S is a multiple of 1/n; since M = 1/n caps
   it and no v = 0 mod n, m_S(k/n) = 1/n EXACTLY for every k with gcd(k,n)=1. So ALL phi(n)
   primitive points k/n are maximizers, and for each k, S must contain speeds in BOTH
   residue classes +k^{-1} and -k^{-1} mod n (else m(k/n) >= 2/n > M, contradiction).
   [Dual of the THM-523 q-witness face: counterexamples must KILL all shallow witnesses;
   tight sets must KEEP every one of them.]

2. **(Class-max formula.)** Near the witness k/n all speeds in the binding classes bind at
   distance exactly 1/n, and the FASTEST one on each side controls the lonely interval:
       **c(S) = (1/n) sum_{k in (Z/n)*} [ 1/max(S cap class(+k^{-1})) + 1/max(S cap class(-k^{-1})) ].**

3. **(AP maximality.)** Since max(S cap class(u)) >= u (least positive residue),
       **c(S) <= c(AP_n) = (2/n) H×(n)**, with equality iff every unit class is realized by
   its least positive residue. The AP has the FASTEST possible collapse; sporadic tight
   sets can only collapse slower or equal. Every lift r -> r + jn of a unit-class element
   (HYP-3750's duplication+drop data) depresses c by exactly (1/n)(1/r - 1/(r+jn)).

## Verified census (exhaustive, exact; speeds <= 40/75, 40, 30, 24, 26)

| n | tight sets found | c | universal? |
|---|---|---|---|
| 4 | {1,2,3} | 2/3 | yes (only AP) |
| 5 | {1,2,3,4}; **{1,3,4,7}** | 5/6; **29/42** | **NO -- law refuted as universality** |
| 6 | {1,2,3,4,5}; {1,3,4,5,9} | 2/5; 2/5 | yes (cross-type lifts only even classes) |
| 7 | {1,2,3,4,5,6} | 7/10 | yes (only AP) |
| 8 | {1..7}; {1,2,3,4,5,7,12}; **{1,4,5,6,7,11,13}** | 44/105; 44/105; **328/1001** | **NO** |
| 14 | {1..13}; {1..11,13,24} (THM-523 census) | 1666/6435; 1666/6435 | yes |
| 20 | {1..19}; {1..17,19,36} | (2/20)H×(20); same = 543160/2909907 | yes |

Class-max formula MATCHES every tight set exactly (column checked in-script). All tight
sets found have all phi(n) shallow witnesses; **no tight set containing a multiple of n
was found at any tested n** (the deep-witness branch of the dichotomy is empty so far --
observed, not proved).

## Why n=14 is universal (and n=5, 8 are not)

The GW family swap n-2 -> 2n-4 lives in the classes -2, -4 mod n: for n EVEN these are
non-unit (even) classes, invisible to the collapse rate. Verified for the whole family
m = 7, 13, 19 (n = 8, 14, 20): c(GW_n) = c(AP_n) exactly. At n=5, {1,3,4,7} lifts the
UNIT class 2 (7 = 2 mod 5): c drops to (1/5)·2·(1 + 1/4 + 1/3 + 1/7)·... = 29/42 < 5/6.
At n=8, {1,4,5,6,7,11,13} lifts unit classes 3 (11) and 5 (13): c = (1/4)(1 + 1/7 + 1/11
+ 1/13) = 328/1001 < 44/105. **Universality at n=14 is an even-class accident, not a law.**

## The cusp envelope (conditional synthesis)

Lambda is decreasing in r, so for any loose set Lambda_S(r) >= Lambda_S(1/n) >= inf-isolation.
Hence near the cusp:
    inf_S Lambda_S(r) = min( c_min(n)·(1-nr) + O((1-nr)^2),  inf over loose S of Lambda_S(1/n) ),
where c_min(n) = min over tight sets of c(S). **The envelope AT r = 1/14 is exactly
THM-523's open lemma inf L > 0** (uniform lower bound on the gap-1/14 lonely measure);
conditional on it, the final window of the envelope is the line c_min·(1-14r) with
c_min(14) = 1666/6435 (census {AP, GW}, both attaining). At n=5 and n=8 the envelope's
cusp line is carried by the SPORADIC (29/42, 328/1001), NOT the AP -- refining HYP-3824's
"the AP is the minimizer": true at n=14 only because GW ties the AP there.

## Consequence for proof strategy (honest)

Measure-based cusp arguments must be aimed at the MINIMAL-c tight sets -- at n=14 the tie
makes AP and GW equally binding, but at general n the GW/cross-type sporadics are the
extremal objects for the lonely-measure envelope, not the AP. Any "certify the AP"
strategy for a measure floor is aiming at the wrong extremal whenever a unit-class lift
exists. The collapse rate is finite data mod n (which classes, which maxima) -- a
first-layer congruence invariant in the same spirit as the Gamma_0(N) localization
(HYP-3833): the cusp derivative is readable from residue classes alone.
