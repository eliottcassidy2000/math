---
id: HYP-3843
title: The FINAL-WINDOW IDENTITY -- CORRECTED (klein-S89, MISTAKE-093): Lambda_AP == Lambda_GW as identical functions on [2/29, 1/14], r* = 2/29 (NOT 1/15); AP alone carries (1/15, 2/29); E separates the pair at n=20 too; deep-witness emptiness extends to n=9,10
status: CONFIRMED (exact) -- part 1 CORRECTED by klein-S89 (r* = 2/29; see MISTAKE-093)
source: klein-2026-07-01-S88; corrected klein-2026-07-01-S89
script: 04-computation/overtaking_defect_sieve_final_window_klein.py (+ .out)
related:
  - HYP-3834/3835/3836 (distribution function, collapse law, crossing)
  - mac-mini HYP-3850(d) (2nd-Farey-layer curvature -- answered sharply here)
  - THM-592/593 (radius-derivative frame)
---

# HYP-3843: the final-window identity and the n=20 / deep-witness extensions

## 1. Final-window identity -- CORRECTED klein-S89 (MISTAKE-093): r* = 2/29, not 1/15

Lambda_AP and Lambda_GW at n=14 are EQUAL AS FUNCTIONS on **[2/29, 1/14]** (width 1/406).
The S88 claim of identity on (1/15, 1/14] was WRONG: GW has a real kink at 2/29
(gap-death of a (5,24) pair at d=2; mac-mini-S94's breakpoint list 2/33, 2/31, 2/29 was
correct), the profiles agree at 2/29 and above and differ on the sliver (1/15, 2/29)
with Lambda_GW > Lambda_AP -- the S88 probe's single midpoint (29/420 = 0.06905) sat just
above 2/29 = 0.06897 and missed the sliver. So the AP ALONE carries the envelope on
(1/15, 2/29): the second-order tie-break in AP's favor, confirming mac-mini S94(3) in
full. What SURVIVES of S88: the last mile is single-valued on the positive shared window
[2/29, 1/14], where only the six mediant gaps live and both tight sets bound them with
the same unit pairs {k^-1, 14-k^-1} (HYP-3835). Lesson (MISTAKE-093): equality at all
candidate breakpoints does not imply identity between them when one function kinks where
the other is straight; per-function per-segment midpoint assertions are mandatory
(HYP-3841's ladder code had them; its K values are unaffected).

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
