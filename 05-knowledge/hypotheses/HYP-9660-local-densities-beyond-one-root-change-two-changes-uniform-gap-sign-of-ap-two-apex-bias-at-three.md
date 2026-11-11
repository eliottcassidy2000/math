---
id: HYP-9660
title: "Local densities beyond THM-5128: (A) by-value limit of the two-root-change family, L(0) = 1 + (p-1)(sigma p - 1)/p^3 with sigma = 2 w_1 = 1.60718 instead of the letter-count 3/2; (B) a uniform spectral gap for the twisted digit operators, no rate above 0.93204 (the rate at 7) at any odd modulus; (C) the sign of the first-order constants A_p; (D) the two-apex grammar G(2,0) has no local bias except a slowly decaying one at 3; (E) the coprime states have density 8/pi^2"
status: >
  OPEN (opus-2026-11-10-S91, 2026-11-11). (A) is a derivation by a helper of this session whose exact ingredients
  (cascade identities X1-X4 of census_k2/cascade_exact.py) are FINITE-EXACT, but whose limit is NOT confirmed by the
  census: on [2^23, 2^24) the measured ratios exceed both candidate limits at p = 3, 5, 11, 17, as THM-5128 (D)
  predicts for primes with early divisible states. Unaudited; typed CONJECTURE. (B) VERIFIED numerically for the
  50 odd moduli m <= 101 and for 12 further moduli up to 1039 (floating point eigenvalues): no rate above 0.93204;
  the rule "a composite modulus has the rate of its worst prime power" is REFUTED (15, 45, 51, 65, 75, 85, 93). (C), (E) are
  questions with NUMERICAL data. (D) is FINITE-EXACT data for all 4^n words, n <= 16, with no theory. The
  root-change census law for general K is HYP-9506 (A) (kind-pasteur), OPEN; THM-5128 settles K = 1 only.
session: opus-2026-11-10-S91
date: 2026-11-11
source: 05-knowledge/results/s91_f1_local_densities_20261110/README.md
scripts:
  - 04-computation/experiments/s91_f1_local_densities_20261110/census_k2/ (census_k1.py, census_k2.py, cascade_exact.py, dict_k1.py and their .out)
  - 04-computation/experiments/s91_f1_local_densities_20261110/rates_exact/rates_scan.py, twisted_spectrum.py, bridges/b3_f1_data.py (m = 1039)
  - 04-computation/experiments/s91_f1_local_densities_20261110/ap_constants.py, divisible_states.py
  - 04-computation/experiments/s91_f1_local_densities_20261110/g20_local/full_hist.py, full_hist_report.py (+ .out)
related:
  - 01-canon/theorems/THM-5128-* (the K = 1 theory), THM-4916-*, 05-knowledge/hypotheses/HYP-9506-* (letter-count model, census OPEN)
  - 01-canon/theorems/THM-4878-* (grammars G(k, w)), THM-4885-*, THM-4886-*
---

# HYP-9660 - local densities beyond one root change

Notation of THM-5128: one-apex states `(x(h), y(h))`, digit chain on `(h mod m, h mod 2)`, twisted operators
`T_k`, rates `max_k rho(T_k)/2`, word-pair sheet shares `w_0 = 0.19640777`, `w_1 = 0.80359223`.

## (A) Two root changes, by value (CONJECTURE; helper derivation, not audited, not resolved numerically)

For the family with two root changes (three one-apex blocks; the `K = 2` words of the census of HYP-9506), ordered
by value, the limiting class factors at an odd prime `p` are conjectured to be

    L(0) = 1 + (p - 1)(sigma p - 1)/p^3,     L(c) = 1 - (sigma p - 1)/p^3  (c != 0),     sigma = 2 w_1 = 1.60718...

The letter-count model has the same shape with `sigma = 3/2` (THM-4916, HYP-9506 (A)). The mechanism is a cascade
(identities X1-X4 of cascade_exact.py, exact on 61,504 words): a `p`-divisible first block either stays confined
at the second join or resets to `(D_2, U_2) = (-1, 0)`, according to one binary feature of the first join whose
letter-count share is `3/4`; the conjecture is that by value this share becomes `w_1`, as in THM-5128 (D3) for
one change.

Evidence and its limit. The census of `[2^23, 2^24)` reproduces HYP-9506 (B)'s numbers and gives ratios
`L(0)/L(c != 0)` of `1.857, 1.614, 1.240, 1.374, 1.171, 1.534, 1.124, 1.102` at `p = 3, 5, 7, 11, 13, 17, 23, 31`;
the conjectured limits are `1.495, 1.298, 1.216, 1.140, 1.119, 1.092, 1.068, 1.051` (letter-count: `1.447, 1.274,
1.199, 1.130, 1.110, 1.085, 1.064, 1.047`). The excess is largest at `3, 5, 11, 17`, the primes with a balanced
divisible state `x = y = 2^m + 1` (THM-5128 (F)), and smallest at `7, 13, 23, 31`. This is the `1/ln X` term of
THM-5128 (D2) again, and at this size it hides the difference between the two candidates.

Decisive test: derive the first-order term for two changes (the same hyperbola argument, with the constant of the
first block), then pre-register the class factors of a window above `2^27`.

## (B) A uniform gap (CONJECTURE; VERIFIED numerically in the stated range)

For every odd `m > 1` and every `k` not divisible by `m`, `rho(T_k)/2 <= 0.93204`, the rate of `m = 7`.
Scanned completely: the 50 odd moduli `3 <= m <= 101` (rates_exact/scan_003_101.npz): the maximum `0.93204` is
attained exactly at the multiples of 7, and the largest rate at a modulus prime to 7 is `0.91374` (`m = 93`).
Isolated further values (interrupted runs): `0.8587, 0.8660, 0.9320, 0.8660, 0.8130, 0.9320, 0.8660, 0.9002,
0.8053, 0.9320` at `m = 121, 141, 161, 261, 281, 301, 321, 341, 361, 441`; `0.8754` at 127; `0.7636` at 1039.
The rate of a composite modulus is NOT the maximum over its prime-power divisors: `15, 45, 75` have `0.8922`
(divisors `0.8660, 0.7418, 0.8153`), `51` has `0.8690`, `65` `0.8482`, `85` `0.8503`, `93` `0.9137` (`31` has
`0.9002`), while the multiples of 7 stay at `0.93204`. So the gap is a statement about all characters, not about
primes. The slow prime moduli are the divisors of `2^j - 1` with small `j` (`7, 31, 127, 73`: `0.932, 0.900,
0.875, 0.852`); `3` has `sqrt(3)/2`. A uniform gap would make the remainder of THM-5128 (C2) uniform in the
modulus and is the input a level-of-distribution statement needs. Unsearched and able to refute: products of
3, 5 with 31 and 127 beyond 101 (`155, 381, 465`).

## (C) The first-order constants (QUESTION; NUMERICAL)

`A_p` (THM-5128 (D2)) is positive at `p = 3, 5, 7, 11, 17, 43` and negative at `p = 13, 31` (partial sums to
`2^26`; `p = 7` not converged). Is `A_p != 0` for every `p`? Is `A_p > 0` exactly when the first `p`-divisible
state is balanced or early enough to beat the deficit `-1/p^2` of the states before it? Does `sum_p A_p` converge?

## (D) The two-apex grammar (FINITE-EXACT data, no theory)

Over all `4^n` words of `G(2,0)` (THM-4878), `n <= 16`, the class factors of `H` are within `0.0004` of 1 at
`p = 5, 7, 11, 13` at `n = 16`, the odd classes modulo 8 and 16 within `0.0032`; at `p = 3` the class `0` keeps a
factor `1.0102, 1.0068, 1.0062, 1.0052, 1.0048` for `n = 12, ..., 16` (modulo 9: `1.0051` at `n = 16`). Question:
does the factor at 3 tend to 1, and at what rate? By THM-5128 (A) a rank-`r` seam with equidistributed states has
bias `p^(-r)`; the slow decay at 3 would then be a mixing rate, not a limit.

## (E) Coprime states (QUESTION; NUMERICAL)

The proportion of `h < 2^n` with `gcd(x(h), y(h)) = 1` is `0.8262, 0.7707, 0.8024, 0.8166, 0.8118` at
`n = 10, 14, 18, 22, 24`; the Euler product of THM-5128 (D3) predicts `8/pi^2 = 0.8106`. A proof needs the
equidistribution of THM-5128 (B) uniformly in `p` up to a power of `2^n`, i.e. (B) plus a large-sieve step.
