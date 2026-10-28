---
id: HYP-9436
title: "The crossing index of the odd harmonic sum is floor(exp(2n - gamma)/4) for every n >= 2 (the odd, half-shifted analogue of Comtet's harmonic-number problem); pre-registered exact test to n = 3000 and to the thresholds k*pi, k <= 500"
status: OPEN in general; PROVED for 2 <= n <= 3000; pre-registered P1, P2, P3 all CONFIRMED (2026-10-28)
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
(THM-4837). The formula is G. Kirchner's FORMULA line in OEIS A092315 (2020, untagged). Kirchner's comment gives an
explicit-remainder test (remainder -127t/(1920 m^4), 0 < t < 1, supported numerically) certifying n <= 1000. No
proof for all n is known, so it is treated as unproved in general.

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
  Prediction: TRUE. Heuristically the expected number of hits for n >= 61 is about 4 x 10^-54 (the audit computes 3.56e-54).
- **(P2)** The fractional parts frac(X_n), 2 <= n <= 3000, pass a Kolmogorov-Smirnov uniformity test at the 1%
  level (D < 1.63/sqrt(N)). Prediction: TRUE.
- **(P3)** The owner's 2 pi generalised: for the thresholds xi = k pi with 1 <= k <= 500, the window test has no
  exception, so m(k pi) = floor(exp(2 k pi - gamma)/4). Prediction: TRUE. Values: m(pi) = 75 (denominator 151) and
  m(2 pi) = 40249 (denominator 80499).

## What would refute

The conjecture is refuted by any n with frac(X_n) inside the window together with an exact check (via the r(M)
bounds and, if needed, the next Euler-Maclaurin term) showing that O at M = floor(X) already exceeds n. A precision
failure (the guard is not met) is reported as UNDECIDED, not as a refutation.

## Results (confirmatory run, 2026-10-28, after the pre-registration commit a14c71ec13)

Command: `python3 borwein_ladder_vsc_twelve_20261028.py 3000 500`.

- **(P1) CONFIRMED.** No window exception for 2 <= n <= 3000, so m(n) = floor(exp(2n - gamma)/4) is PROVED there. The
  tightest case is n = 2, with frac(X)/window = 101.8. n = 1 holds by direct check (m(1) = 1).
- **(P2) CONFIRMED.** Kolmogorov-Smirnov D = 0.0192 for N = 2999, against the 1% critical value 0.0298.
- **(P3) CONFIRMED.** No exception for the thresholds k pi, 1 <= k <= 500: m(pi) = 75 (denominator 151),
  m(2 pi) = 40249 (80499), m(3 pi) = 21553437 (43106875).

**OEIS A092315 status (read 2026-10-28).**
- G. Kirchner (2020) states the floor formula with an explicit-remainder test certifying n <= 1000. Its remainder
  bound is supported numerically; the GKP/DLMF envelope used here proves it.
- M. F. Hasler (2017) has the formula line floor(exp(2n - gamma)/4 + 1/8) (conjectured). It contradicts the entry's
  own DATA term a(17), so it is **false**: it differs from the proved
  value at n = 17, 31, 45, 57, 59, 63, 86, 118, ..., 375 of the n <= 3000 (THM-4837 (6)). The case n = 17 is also
  checked independently via digamma.

The rerun added the Hasler and digamma checks and the genus and Monster assertions of THM-4838. The pre-registered
quantities are unchanged.

**Audit (2026-10-28,
`04-computation/experiments/borwein_ladder_vsc_twelve_20261028_audit_S67/AUDIT.md`, SOUND WITH CORRECTIONS).**
- The pre-registration order is confirmed: a14c71ec13 precedes the confirmatory checkpoint 7a5cd669d0, and the P1-P3
  text is unchanged.
- The session's window test was high-precision floating point with a guard. The audit re-certified it in mpmath
  interval arithmetic for 2 <= n <= 3000 and k <= 500, finding no exception (tightest n = 2, frac/window = 101.84).
  "PROVED for 2 <= n <= 3000" stands.
- P1-P3 were near-certain predictions, so their confirmation adds little beyond the certificate.
