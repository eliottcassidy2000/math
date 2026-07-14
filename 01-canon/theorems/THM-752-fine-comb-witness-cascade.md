---
id: THM-752
title: The FINE-COMB WITNESS LEMMA (tooth threshold) + the ratio-13 cascade -- if the peeled core's 1/14-safe set contains an interval longer than ONE danger tooth (|I| > 1/(7w)), then M(V) >= 1/14 with an explicit rational witness inside I; LRC(<=13) gives a guaranteed core interval 1/(91 p_max), so w > 13 p_max closes unconditionally (all-consecutive-ratios > 13 recursively); closes klein-S304's loose-branch STALL with a verified exact witness; census 60/60 on spread covering bodies
status: PROVED (3 lines) + VERIFIED-EXACT (witnesses checked against all 13 speeds as Fractions: klein's stall {1,10,...,390} closes at t = 233/2912 with clearance 29/364; two more exemplars likewise; 60/60 census)
source: opus-2026-07-14-S285 (owner prompt: synthesize the remaining frontier and push it)
depends_on: []   # LRC(<=13) named citation for the cascade corollary only; the per-body lemma is citation-free
related:
  - klein-S304 (HYP-6670: the loose branch = iterated far-peel; THIS closes its crude-stall residual)
  - mac-mini S98/S102 (the dichotomy; THM-751 aligned monotonicity -- the binding side's twin)
  - spread13_lonely (kps/death-star: TOTAL ratio <= 13 closes; the cascade closes ALL-consecutive > 13: the two pincer the ratio spectrum)
  - THM-739/740 (the shaped-family clocks; this is the shapeless witness form)
---

# THM-752 -- the fine-comb witness lemma and the ratio-13 cascade

## The lemma (3 lines)

Let V be a 13-speed body, w in V, C = V minus {w}.  Suppose the 1/14-safe set of C contains an
interval I with **|I| > 1/(7w)** (one w-danger tooth).  The w-danger set is a union of disjoint
OPEN teeth of length exactly 1/(7w); a connected I longer than a tooth cannot sit inside one, so
I contains a point t* of the CLOSED w-safe set.  At t*: core clearances >= 1/14 (t* in I) and
||w t*|| >= 1/14: **M(V) >= 1/14**.  QED.  (The naive full-period threshold 1/w is 7x too
strong -- and empirically misses the loose stratum by exactly ~1.15x; the tooth threshold fires
with ~6x margin.  The margin IS the duty cycle 6/7 of the perspective frame's origin band.)

## The cascade (unconditional corollary, LRC(<=13) cited)

Around the core's LRC(<=13) witness, clearances degrade at rate p_max = max(C): the 1/14-safe
set contains an interval of length >= 2(1/13 - 1/14)/p_max = 1/(91 p_max).  This beats the
tooth 1/(7w) iff **w > 13 p_max**.  Recursing (peel largest, cite LRC(k) at each level):

> every body whose consecutive ratios all exceed 13 satisfies LRC(14), unconditionally.

Complementary to spread13_lonely (TOTAL ratio <= 13): the two theorems pincer the ratio
spectrum from both ends; the exact per-body interval test covers the middle.

## Verified exact witnesses (the loose-branch stall closes)

| body | fires at peel | witness t | min clearance over ALL 13 speeds |
|---|---|---|---|
| klein-S304 stall {1,10,21,24,56,65,77,135,219,265,335,367,390} | 390 | **233/2912** | **29/364 = 0.0797 >= 1/14** |
| {1,3,7,15,31,63,90,127,200,255,300,380,420} | 420 | 709/9800 | 709/9800 = 0.0723 |
| {11,13,29,41,73,97,143,190,240,280,330,360,390} | 390 | 901/120120 | 61/840 = 0.0726 |

CENSUS: 60/60 adversarial spread covering bodies (speeds to 450) fire.  With klein-S304's
iterated far-peel (rigorous on most of the stratum via the crude arc-count bound), the loose
branch is now DOUBLE-COVERED by two cheap exact certificates; the remaining uniform statement
("saturated-spread => some core interval > 1/(7w)") is the loose branch's single named lemma.

## Position

The witness is the S270 frozen-fan/sweeper mechanism in its terminal form: the core freezes an
interval; the largest speed is a fine comb whose duty cycle cannot blanket it.  Everything is
rational arithmetic per body (Lean shape: interval engine + one comparison), and the cascade is
a 6-line induction over the LRC(<=13) citations.

## Files

04-computation snippet in 05-knowledge/results/lrc14_finecomb_cascade_thm752_opus_S285.out
(lemma tests, witnesses, census); the frontier synthesis: 00-navigation/LRC14-FRONTIER-2026-07-14.md.
