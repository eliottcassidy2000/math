---
id: HYP-9436
title: "The crossing index of the odd harmonic sum is floor(exp(2n - gamma)/4) for every n >= 2 (the odd, half-shifted analogue of Comtet's harmonic-number problem); pre-registered exact test to n = 3000 and to the thresholds k*pi, k <= 500"
status: OPEN (PRE-REGISTERED; confirmatory run pending at the time of this commit)
session: opus-2026-10-28-S67
date: 2026-10-28
related: [THM-4837, THM-4838]
note: 05-knowledge/results/borwein_ladder_owner_numbers_vsc_twelve_20261028.md
script: 04-computation/experiments/borwein_ladder_vsc_twelve_20261028.py
---

# HYP-9436 - the odd harmonic crossing formula

## Statement

For a real threshold xi > 1 let m(xi) be the least m with

    O_m = 1 + 1/3 + 1/5 + ... + 1/(2m+1) > xi,

and put X(xi) = exp(2 xi - gamma)/4. **Conjecture:** m(n) = floor(X(n)) for every integer n >= 2.

The owner's 7, 56, 418 are m(2), m(3), m(4); 15 and 113 are 2m+1 at n = 2, 3; 80499 is 2m+1 at xi = 2 pi
(THM-4837). The formula is the one displayed in OEIS A092315. Its status there must be checked by the audit;
this file treats it as unproved in general.

## What is proved (THM-4837 (5))

Write M = m + 1. Then O_m = (1/2) ln(4M) + gamma/2 + r(M), with 0 < r(M) < 1/(48 M^2) for every M >= 1. The bound
follows from the Graham-Knuth-Patashnik envelope for H_N, and is checked directly for small M.

Hence **m(xi) = floor(X) whenever frac(X) >= X/(24 floor(X)^2)**, which is about 1/(24X). Otherwise
m(xi) is floor(X) - 1 or floor(X). The conjecture can therefore fail only when exp(2n - gamma)/4 lies just above an
integer, within about 1/(24 X). Deciding this for all n needs Diophantine information about exp(2n - gamma) that
nobody has. The classical analogue for H_k is Comtet's problem ("the least k with H_k > n is round(exp(n - gamma))",
Comtet 1967, Boas-Wrench 1971), which is open in general for the same reason.

## Exploratory evidence (before this pre-registration)

n = 2..60 passes the window test, so the formula is PROVED there. The tightest case is n = 2, with
frac(X)/window = 101.8. Under uniformly distributed fractional parts, the heuristic expected number of exceptions
over all n >= 2 is 0.0063.

## Pre-registered predictions

Confirmatory run: `python3 borwein_ladder_vsc_twelve_20261028.py 3000 500`.

- **(P1)** The window test has no exception for 61 <= n <= 3000, so the formula is PROVED for all 2 <= n <= 3000.
  Prediction: TRUE. Heuristically the expected number of hits for n >= 61 is about 10^-52.
- **(P2)** The fractional parts frac(X_n), 2 <= n <= 3000, pass a Kolmogorov-Smirnov uniformity test at the 1%
  level (D < 1.63/sqrt(N)). Prediction: TRUE.
- **(P3)** The owner's 2 pi generalised: for the thresholds xi = k pi with 1 <= k <= 500, the window test has no
  exception, so m(k pi) = floor(exp(2 k pi - gamma)/4). Prediction: TRUE. Values: m(pi) = 75 (denominator 151) and
  m(2 pi) = 40249 (denominator 80499).

## What would refute

The conjecture is refuted by any n with frac(X_n) inside the window together with an exact check (via the r(M)
bounds and, if needed, the next Euler-Maclaurin term) showing that O at M = floor(X) already exceeds n. A precision
failure (the guard is not met) is reported as UNDECIDED, not as a refutation.
