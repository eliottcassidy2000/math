---
id: HYP-3834
title: The lonely distribution function is an exact finite gap-sum; the AP profile is the Farey gap-sum; the LRC14 collapse constant is exactly 1666/6435 = (2/n)·H×(n)
status: STUB -- CLAIMED klein-2026-07-01-S87, computation in progress this session
source: klein-2026-07-01-S87
script: 04-computation/lonely_profile_farey_slope_klein.py (to be written this session)
related:
  - HYP-3824 (mac-mini numeric slope 0.26(1-14r) -- made exact here)
  - HYP-3015/HYP-3018 (codex lonely-profile barcode, t-side object; Lambda is its layer-cake)
  - THM-523 (tight census {AP, GW})
  - HYP-3762, HYP-3763 (scale-matching)
  - 07-reflections/zeta2-governs-the-lonely-runner-floor.md
---

# HYP-3834 (STUB): the lonely distribution function and the collapse constant

## Claim (reserved, being tested)

For a speed set S, let m_S(t) = min_{v in S} ||vt|| (the codex-S179 "lonely profile") and

    Lambda_S(r) = meas{t in [0,1) : m_S(t) >= r}   (the lonely distribution function).

1. **Gap-sum formula**: Lambda_S(r) is an exact finite sum over consecutive pairs of
   danger-interval centers (fractions a/u, b/v with u,v in S): each pair contributes
   (spacing - r(1/u + 1/v))^+. Piecewise LINEAR in r with rational breakpoints.
2. **AP = Farey**: for S = {1..13}, the centers are exactly the Farey fractions F_13, and
   Lambda_AP(r) = sum over F_13-adjacent (a/q, a'/q') of (1 - r(q+q'))^+ / (qq').
   Tightness at r=1/14 IS the mediant identity (adjacent q+q' >= 14, equality exactly at
   the primitive mediants k/14, gcd(k,14)=1).
3. **Collapse constant**: c(AP) := lim_{r->1/14^-} Lambda_AP(r)/(1-14r)
   = sum over the 6 mediants k/14 of 1/(q q') = (2/14)·(1 + 1/3 + 1/5 + 1/9 + 1/11 + 1/13)
   = **1666/6435 = 0.258897...** -- the exact value of mac-mini's numeric 0.26 (HYP-3824).
4. **General n**: c(AP_n) = (2/n)·H×(n), H×(n) = sum_{u<n, gcd(u,n)=1} 1/u
   ~ (2·phi(n)/n^2)(ln n + gamma + sum_{p|n} ln p/(p-1)) -> the atom hardens like log n / n.

## Why it matters

Fixed-radius instruments (moments at r=1/14) hit the measure-0 atom; the distribution
function works in the (r, measure) PLANE, where the cusp has a linear slope with an exact
arithmetic value. Fubini: int_0^inf Lambda_S(r) dr = int_0^1 m_S(t) dt = mean loneliness.

## Status

STUB. Derivation by hand (this session): parts 2-4 verified by hand at n=14; script + exact
rational verification to follow today. Honest: part 1 is elementary; part 2 is the classical
Farey/mediant structure; the new content is the exact constant + the distribution-function
framing + the census consequences (HYP-3835).
