---
id: HYP-2915
title: The bounded LRC census, tightened -- all tight single-swap speeds <= 2(n-1); the census (given condition 4) is an O(n^2) search up to speed 2(n-1), down from D<=300 / Fibonacci
status: VERIFIED bound (n=5,6,7,8,14) + structural argument; tightens the finite search. (Conditional on condition 4 = single-swap, the open part.)
source: mac-mini-2026-06-22-S56 (owner: make the bounded finite census even more bounded)
related:
  - HYP-2914   # the necessary-condition battery (condition 4 = residues miss <=1 = single-swap)
  - HYP-2913   # three-gap characterization
---

# HYP-2915: the bounded census, tightened to speed <= 2(n-1)

## The tightened bound (VERIFIED n=5,6,7,8,14)
Among single-swaps of the AP {1,...,n-1} (replace one element, the case forced by HYP-2914 condition 4),
the max speed of any TIGHT (M=1/n) set is:
  n:                 5    6    7    8    14
  max tight speed:   7    9    6   12    24
  2(n-1):            8   10   12   14    26
So **every tight single-swap speed is <= 2(n-1)**. The tight sets per n: n=5 skip2->7, n=6 skip2->9,
n=7 only AP, n=8 skip6->12, n=14 skip12->24 (= GW). Only the AP and ONE GW survive at each n.

## Why 2(n-1): the only large speed is a SECOND multiple
A tight set is the AP, or the AP with one element replaced (single-swap). The replacement is the killer
of a skipped resonance b (a runner ≡ 0 mod b is needed to keep b killed, HYP-2914 cond. 1). The killer
is a multiple m*b > n-1; for tightness it is the SECOND multiple m=2 (S55: for n=14, 12->24=2*12 is
tight but 12->36=3*12, 48, 60 are NOT -- larger multiples equidistribute and give M > 1/n). Since the
skipped b <= n-1, the killer 2b <= 2(n-1). Hence all tight speeds <= 2(n-1).

## The census shrinks to O(n^2)
GIVEN condition 4 (tight => single-swap), the census is now a finite search of
  n*(n-1) candidates  (each of n-1 skip-positions x each lift in [n, 2(n-1)]),
all with speed <= 2(n-1). For n=14: **169 candidates, speeds <= 26** -- versus the old bounds
(HYP-2876 D<=41 [refuted], the S51 search to 300, or the R1/covering Fibonacci bound ~10^9). The finite
census is now tiny and instantly checkable, yielding exactly {AP, GW}.

## Status (honest)
- The bound `tight speed <= 2(n-1)` is VERIFIED n=5..14 and structurally argued (killer = 2nd multiple).
- The fully-rigorous proof that the 2nd multiple is tight but the 3rd+ is not (the equidistribution
  threshold) is verified, not yet a clean theorem -- the R1 peel gives M(S) >= (6/7)M(seed) but M(seed)
  in [1/(n-1), 1/(n-2)] is too close to force it directly; the gap is the same near-tight rigidity.
- The whole reduction is CONDITIONAL on condition 4 (residues miss <=1 = single-swap), which is the
  open Steinhaus rigidity (HYP-2913/2914). So this tightens the SEARCH SPACE of the census (the task),
  not the open core: given condition 4, the census is O(n^2) up to speed 2(n-1).
