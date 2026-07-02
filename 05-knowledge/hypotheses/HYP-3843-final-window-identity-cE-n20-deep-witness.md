---
id: HYP-3843
title: The FINAL-WINDOW IDENTITY -- Lambda_AP == Lambda_GW as identical functions on (1/15, 1/14], r* = 1/15 exactly; the AP/GW tie breaks AT the 1/15 breakpoint, at no interior order; E separates the pair at n=20 too; deep-witness emptiness extends to n=9,10
status: CONFIRMED (exact)
source: klein-2026-07-01-S88
script: 04-computation/overtaking_defect_sieve_final_window_klein.py (+ .out)
related:
  - HYP-3834/3835/3836 (distribution function, collapse law, crossing)
  - mac-mini HYP-3850(d) (2nd-Farey-layer curvature -- answered sharply here)
  - THM-592/593 (radius-derivative frame)
---

# HYP-3843: the final-window identity and the n=20 / deep-witness extensions

## 1. Final-window identity (answers the curvature question)

Lambda_AP and Lambda_GW at n=14 are EQUAL AS FUNCTIONS on (1/15, 1/14] -- verified exactly
at every candidate breakpoint of both profiles and at midpoints; the LAST radius where
they differ is r = 1/15 itself. So the AP/GW tie does NOT break at second order, third
order, or any order inside the final window: on the last Farey window the two tight sets
have literally the same lonely measure, because only the six mediant gaps survive there
and both sets bound them with the same unit pairs {k^-1, 14-k^-1} (HYP-3835). The tie
breaks exactly AT the q+q' = 15 layer (r* = 1/15), where the profiles' 11 differing
radii end. Consequence for the envelope: on (1/15, 1/14] the tight locus contributes ONE
function, not two -- the census's "last mile" is single-valued; curvature separates AP
from GW only at and below 1/15.

## 2. (c, E) at n=20

E(AP_19) = 0.02121322, E(GW_19) = 0.02120207 (exact rationals in .out): **GW < AP again**
(diff -1.1e-5). Same pattern as n=14: collapse rate ties (HYP-3835 family universality),
mean loneliness separates, GW wins the mean. Second data point for "the 2(n-1) element
always wins the mean" and for (c, E) as a separating pair on the tight locus.

## 3. Deep-witness emptiness at n=9, 10

Exhaustive tight censuses (n=9 speeds <= 20; n=10 speeds <= 18): only the APs; no tight
set contains a multiple of n. The deep-witness branch of the HYP-3835 dichotomy stays
empty through n=10 (previously n<=8). The "tight => no multiple of n" proof target
(mac-mini HYP-3850(e)) gains two more supporting censuses.
