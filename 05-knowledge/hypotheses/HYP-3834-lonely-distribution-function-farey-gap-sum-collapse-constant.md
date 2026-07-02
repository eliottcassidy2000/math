---
id: HYP-3834
title: The lonely distribution function is an exact finite gap-sum; the AP profile is the Farey gap-sum; the LRC14 collapse constant is exactly 1666/6435 = (2/n)·H×(n)
status: CONFIRMED (exact rational computation, n=4..16; all identities exact)
source: klein-2026-07-01-S87
script: 04-computation/lonely_profile_farey_slope_klein.py (+ .out)
related:
  - HYP-3824 (mac-mini numeric slope 0.26(1-14r) -- made exact here: 1666/6435)
  - HYP-3015/HYP-3018 (codex lonely-profile barcode, t-side; Lambda = its layer-cake)
  - THM-523 (tight census {AP, GW}; single-perturbation infimum 1/1260 -- reproduced exactly)
  - HYP-3763 (scale-matching), HYP-3762 (support rigidity)
  - 07-reflections/the-collapse-rate-law-and-the-lonely-envelope-klein.md
---

# HYP-3834: the lonely distribution function and the collapse constant -- CONFIRMED

## Objects

m_S(t) = min_{v in S} ||vt|| (the codex-S179 lonely profile, t-side).
**Lambda_S(r) = meas{t in [0,1): m_S(t) >= r}** -- the lonely DISTRIBUTION FUNCTION (r-side).
Layer-cake: int_0^infty Lambda_S(r) dr = int_0^1 m_S(t) dt = E(S) (mean loneliness).
Lambda is dilation-invariant (t -> ct measure-preserving mod 1): Lambda_{cS} = Lambda_S.

## Confirmed results (all EXACT rational arithmetic)

1. **Gap-sum formula** (elementary lemma): the danger set at radius r is the union of
   intervals (a/v - r/v, a/v + r/v); Lambda_S(r) = sum of positive gaps between consecutive
   danger intervals -- piecewise LINEAR in r, breakpoints where gaps die.

2. **AP = Farey** (exact equality verified at r = 1/100, 1/30, 1/20, 1/15, 9999/140000):
       Lambda_{1..13}(r) = sum over F_13-adjacent (a/q, a'/q') of (1 - r(q+q'))^+ / (qq').
   Mechanism: consecutive danger centers are the Farey fractions; gap between the intervals
   of an adjacent pair = 1/(qq') - r(1/q + 1/q') = (1 - r(q+q'))/(qq'). Tightness at r=1/14
   IS the mediant property: adjacent q+q' >= 14, equality exactly at the 6 primitive
   mediants k/14, gcd(k,14)=1. **The AP's tightness is the Farey mediant identity.**

3. **Last-mile exact linearity**: the only Farey-adjacent pairs with q+q' = 14 are the
   mediant parents of k/14; all others have q+q' >= 15. Hence Lambda_AP is EXACTLY linear
   on [1/15, 1/14]:
       Lambda_AP(r) = c(AP)·(1 - 14r),  c(AP) = sum_{6 mediants} 1/(qq')
                    = 2(1/13 + 1/33 + 1/45) = **1666/6435 = 0.2588966...**
   This is mac-mini's numeric 0.26 (HYP-3824), now exact. Verified three independent ways:
   direct Lambda computation, mediant-parent sum, and the witness formula below.

4. **Witness/binding structure at n=14**: AP has exactly phi(14)=6 maximizers t = k/14,
   binding pairs {k^{-1} mod 14, 14 - k^{-1}} = {1,13},{5,9},{3,11}. Contribution of each
   witness to the slope: (1/14)(1/v+ + 1/v-). Hence
       **c(AP_n) = (2/n)·H×(n)**,  H×(n) = sum_{u<n, gcd(u,n)=1} 1/u.
   Verified n = 4..16 (all exact matches; e.g. c = 2/3, 5/6, 2/5, 7/10, 44/105, 69/140,
   20/63, 671/1260, 92/385, 6617/13860, 1666/6435, 1205/4004, 11384/45045).
   Asymptotics: c(AP_n) ~ (2 phi(n)/n^2)(ln n + gamma + sum_{p|n} ln p/(p-1)) -> 0 like
   log n / n: **the atom hardens logarithmically**; even n (phi(n)/n small) collapse flatter.

5. **Mean loneliness closed form** (exact): E(AP_13) = (1/2) sum_{F_13-adjacent}
   1/(q q'(q+q')) = 62836087/2059318800 = 0.03051304.

6. **Cross-validation with canon**: Lambda_{1..12,182}(1/14) = 4637/194040 = 0.023897
   (mac-mini HYP-3824's 0.0239, exact). Lambda_{1..11,13,36}(1/14) = **1/1260 exactly** =
   THM-523's single-perturbation infimum, reproduced independently by this machinery.
   M({1..11,13,36}) = 3/41 with witnesses 17/41, 24/41, binding {5,36}.

## Why it matters

Fixed-radius instruments (global moments, Fourier) see the measure-0 atom and stall
(HYP-3791/3822). The distribution function works in the (r, measure) plane where the cusp
has a LINEAR approach with an exact arithmetic slope -- and the slope is a completely
explicit unit-harmonic sum. See HYP-3835 for the general law over the tight locus and
HYP-3836 for the envelope/crossing structure.
