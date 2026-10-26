---
id: HYP-9410
title: "Collatz-line calibration of the no-sixth-Fermat-prime heuristic. The integers that Collatz attaches to 2^m + 1 and 2^p - 1 (the exact Syracuse runs 3^j 2^(m-2j) + 1 and 3^k 2^(p-k) - 1) are prime at the rate of the sieve-corrected random model (exact sieve to 2^16, then Mertens/ln x) once their algebraic and Aurifeuillean factorisations (w^r +- 1 when r divides both exponents, 4u^4 + 1) are removed; without that removal the model over-predicts by 3-6% (pre-registered L-P1 failed at z = -3.2). Genuine Fermat lines (m = 2^e) and F_m itself have no such factorisations, so the uncorrected model applies to them as is -- the same model that gives an expected 6.7e-10 (Boklan-Conway: < 2^-30) further Fermat primes."
status: >
  OPEN (opus-2026-10-26-S58).
  Pre-registered tests (PREREGISTRATION.md, PREREGISTRATION-2.md; each committed before its test):
  - L-P1 Fermat-type lines m = 8..1200: FAIL, 2548 observed against 2707.11 +- 49.51;
  - L-P2 Mersenne-type lines p = 8..400: pass, 1034 against 1087.58 +- 30.11;
  - L-P3 genuine Fermat lines e = 3..12: pass, 17 against 22.91 +- 4.43.
  Post-hoc correction (algebraic factorisations): z = -1.47, +0.51 and -1.33.
  Out-of-sample corrected tests:
  - L-P4 m = 1201..1500: pass, 677 against 671.98 +- 25.61 (the uncorrected model gave 697.27);
  - L-P5 p = 401..520: pass, 323 against 324.06 +- 17.55 (the uncorrected model gave 348.07).
  The immunity statements are PROVED (gcd argument, below). The Fermat-prime expectation is HEURISTIC (CITED:
  Boklan-Conway arXiv:1605.01371; the cruder 2/ln F_n model is computed here).
session: opus-2026-10-26-S58
source: 05-knowledge/results/s58_fermat_pepin_collatz_20261026/README.md
related:
  - 01-canon/theorems/THM-4793-the-fermat-line-under-collatz-valuation-two-descent-exact-merge-of-consecutive-lines-and-pepin-as-an-affine-cycle-type.md
  - 01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md
  - 01-canon/theorems/THM-4800-pepin-as-the-valuation-two-collatz-branch-trivial-cycle-half-turn-fermat-line-arc-mersenne-mirror.md
scripts:
  - 04-computation/fermat_pepin_s58/collatz_line_primes_20261026.py (+ prediction/result/posthoc/prediction2/result2 JSON)
---

# HYP-9410: Collatz-line calibration of the Fermat-prime heuristic

**Claim (heuristic, tested).** Consider the values on a Collatz line:
- the Fermat-type line `x_j = 3^j 2^(m-2j) + 1`, `1 <= j <= floor((m-1)/2)`;
- the Mersenne-type line `y_k = 3^k 2^(p-k) - 1`, `1 <= k <= p-1`.

After exact sieving by the primes up to `Q`, the probability that a surviving value `x` is prime is
`prod_{q <= Q} (1 - 1/q)^-1 / ln x`, provided `x` has no algebraic factorisation. That proviso excludes:
- `x_j = w^r + 1` when an odd `r > 1` divides `gcd(j, m)`;
- `x_j = 4u^4 + 1` when `j = 0 (mod 4)` and `m = 2 (mod 4)`;
- `y_k = w^r - 1` when `r = gcd(k, p) > 1`.

**Immunity (PROVED).**
- **Genuine Fermat lines.** For `m = 2^e`, `gcd(j, m)` is a power of 2 and `m = 0 (mod 4)` when `e >= 2`, so no value
  on the genuine Fermat line of `F_e` has these factorisations.
- **Fermat numbers.** `F_e` itself has none: `x^(2^e) + 1 = Phi_(2^(e+1))(x)` is irreducible, and the Aurifeuillean
  identity needs exponent `2 (mod 4)`.
- **Mersenne lines.** For prime `p`, `gcd(k, p) = 1` (`0 < k < p`), so the Mersenne lines of prime index are immune too.

**Consequence.** The random model that predicts no sixth Fermat prime contains no algebraic correction, because
Fermat numbers have no algebraic factorisations. The test shows that model is accurate on Collatz's own neighbours of
`F_e` exactly when they too are free of algebraic factorisations.
- With `P(F_n prime) ~ 2/ln F_n`, the expected number of Fermat primes beyond `F_32` is
  `sum_{n >= 33} 2/(2^n ln 2) = 6.7e-10`.
- Boklan–Conway's finer model gives `< 2^-30`.

**Predictions (falsifiable).**
- (a) For Fermat-type lines `m = 1501..2500` and Mersenne-type lines `p = 521..800`, the corrected model's prime
  counts hold within 2.5 sd, while the uncorrected model over-predicts by 3–6%.
- (b) Under the sieve-corrected model, a genuine Fermat line of `F_e` carries about `3.2–3.3` probable primes for
  `e = 9..12` (exact per-line expectations `3.20, 3.18, 3.32, 3.31`). Conjecturally the per-line expectation tends to
  a constant near `3.3`, the plain model's `2.43` times the lines' singular series. The observed counts on single
  lines should then be roughly Poisson around it.
  - The naive model with only the 2- and 3-corrections gives `3 ln(2 ln 2/ln 3)/(2 ln 2 - ln 3) = 2.43`. It misses a
    singular series of about `1.35`.
  - Example (PROVED): no point of a genuine Fermat line, and not `F_e` itself, is divisible by 7. Since
    `3/4 = -1 (mod 7)`, the line alternates mod 7 between `2^(2^e) + 1` and `1 - 2^(2^e)`, and `2^(2^e)` is never
    `+-1 (mod 7)`.
  - Against that, 5 divides a quarter of the points (random: a fifth) and 13 a sixth (random: a thirteenth).
  - Observed for `e = 3..12`: `0, 0, 0, 2, 3, 2, 0, 2, 2, 6`, total 17 against 22.9 predicted.

**Non-consequences.** Nothing here proves anything about Fermat primes. The test calibrates a heuristic on nearby
integers and does not touch the integers `F_n` themselves.
