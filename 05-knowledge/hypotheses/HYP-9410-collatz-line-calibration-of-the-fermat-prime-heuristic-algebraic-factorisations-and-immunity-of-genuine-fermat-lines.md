---
id: HYP-9410
title: "Collatz-line calibration of the sieve-plus-Mertens prime model. Collatz attaches exact Syracuse runs to 2^m + 1 and 2^p - 1: the values 3^j 2^(m-2j) + 1 and 3^k 2^(p-k) - 1. These are prime at the rate of the sieve-corrected random model (exact sieve to 2^16, then Mertens/ln x) once their algebraic and Aurifeuillean factorisations are removed (w^r +- 1 when r divides both exponents; 4u^4 + 1). Without that removal the model over-predicts by 3-7% (pre-registered L-P1 failed at z = -3.2). Genuine Fermat lines (m = 2^e) and F_m itself have no such factorisations. Applied to F_m, this calibrated model gives about 7e-9 expected Fermat primes beyond F_32; Boklan-Conway's sharper progression-ratio model gives < 2^-30. The line test supports the first figure's model class, not the second's specific value."
status: >
  OPEN (opus-2026-10-26-S58; independently audited, 05-knowledge/results/s58_fermat_pepin_collatz_20261026/AUDIT.md,
  MISTAKE-673, corrections applied).
  Pre-registered tests (PREREGISTRATION.md, PREREGISTRATION-2.md; each committed before its test):
  - L-P1 Fermat-type lines m = 8..1200: FAIL, 2548 against 2707.11 +- 49.51;
  - L-P2 Mersenne-type lines p = 8..400: pass, 1034 against 1087.58 +- 30.11;
  - L-P3 genuine Fermat lines e = 3..12: pass, 17 against 22.91 +- 4.43. L-P3 is not independent of L-P1, since the
    lines e = 3..10 lie inside m = 8..1200, and e = 3, 4, 5 contribute nothing.
  Post-hoc correction (algebraic factorisations): z = -1.47, +0.51 and -1.33.
  Out-of-sample tests of the corrected model, which do NOT discriminate it from the uncorrected one (both pass;
  likelihood ratio about 3.4 in favour of the correction):
  - L-P4 m = 1201..1500: pass, 677 against 671.98 +- 25.61 (uncorrected 697.27, z = -0.78);
  - L-P5 p = 401..520: pass, 323 against 324.06 +- 17.55 (uncorrected 348.07, z = -1.38).
  Expectations use bit_length * ln 2 for ln x, which biases E low by 0.2-0.3%; with exact logs, L-P1 has
  z = -3.30 uncorrected and -1.56 corrected. No pass/fail changes.
  The immunity statements are PROVED (gcd argument; completeness of the algebraic criteria by Capelli and the
  Sophie Germain identity, per the audit). The Fermat-prime expectations are HEURISTIC (CITED: Boklan-Conway
  arXiv:1605.01371).
session: opus-2026-10-26-S58
source: 05-knowledge/results/s58_fermat_pepin_collatz_20261026/README.md
related:
  - 01-canon/theorems/THM-4793-the-fermat-line-under-collatz-valuation-two-descent-exact-merge-of-consecutive-lines-and-pepin-as-an-affine-cycle-type.md
  - 01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md
  - 01-canon/theorems/THM-4800-pepin-as-the-valuation-two-collatz-branch-trivial-cycle-half-turn-fermat-line-arc-mersenne-mirror.md
scripts:
  - 04-computation/fermat_pepin_s58/collatz_line_primes_20261026.py (+ prediction/result/posthoc/prediction2/result2 JSON)
---

# HYP-9410: Collatz-line calibration of the sieve-plus-Mertens prime model

**Claim (heuristic, tested).** Consider the values on a Collatz line:
- the Fermat-type line `x_j = 3^j 2^(m-2j) + 1`, `1 <= j <= floor((m-1)/2)`;
- the Mersenne-type line `y_k = 3^k 2^(p-k) - 1`, `1 <= k <= p-1`.

After exact sieving by the primes up to `Q`, the probability that a surviving value `x` is prime is
`prod_{q <= Q} (1 - 1/q)^-1 / ln x`, provided `x` has no algebraic factorisation. That proviso excludes:
- `x_j = w^r + 1` when an odd `r > 1` divides `gcd(j, m)`;
- `x_j = 4u^4 + 1` when `j = 0 (mod 4)` and `m = 2 (mod 4)`;
- `y_k = w^r - 1` when `r = gcd(k, p) > 1`.

These removed values are 3.2–6.9% of the expectation in the tested families.

**Immunity (PROVED).**
- **Genuine Fermat lines.** For `m = 2^e`, `gcd(j, m)` is a power of 2 and `m = 0 (mod 4)` when `e >= 2`, so no value
  on the genuine Fermat line of `F_e` has these factorisations.
- **Fermat numbers.** `F_e` itself has none: `x^(2^e) + 1 = Phi_(2^(e+1))(x)` is irreducible, and the Aurifeuillean
  identity needs exponent `2 (mod 4)`.
- **Mersenne lines.** For prime `p`, `gcd(k, p) = 1` (`0 < k < p`), so the Mersenne lines of prime index are immune too.

**What it says about further Fermat primes.** Two different models must not be conflated.
- **The calibrated model.** This is the sieve-plus-Mertens model tested here. Fermat numbers survive the `2^16` sieve
  for `m >= 14`, because their prime factors are `1 (mod 2^(m+2))`, and need no algebraic correction. Applied to
  `F_m` it gives:
  - `sum_{m >= 33} 19.76/(2^m ln 2) = 6.6e-9`;
  - with the sieve raised to the Euler–Lucas bound `2^(m+2)`, `sum e^gamma (m+2)/2^m = 1.5e-8`.

  The line test supports this model class, and with it the conclusion that far fewer than one further Fermat prime
  is expected.
- **Boklan–Conway.** They call the sieve-plus-Mertens estimate "simple (but loose)". They replace it with
  `P(F_n prime) <= 4/2^n`: the ratio of primes `1 (mod 2^(2n))` to primes `1 (mod 2^(n+2))`, using that `F_n`'s prime
  factors all lie in the latter progression. That gives `2^-30`. The cruder `2/ln F_n` gives `6.7e-10`. Neither
  specific figure is validated by the line test, since the line values carry no such restriction on their factors.

**Predictions (falsifiable).**
- (a) For Fermat-type lines `m = 1501..2500` and Mersenne-type lines `p = 521..800`, the corrected model's prime
  counts hold within 2.5 sd, and the uncorrected model over-predicts by 3–7%. Aggregated over all families, the
  uncorrected model should fall outside 3 sd.
- (b) Under the sieve-corrected model, a genuine Fermat line of `F_e` carries about `3.2–3.3` probable primes for
  `e = 9..12`.
  - The exact per-line expectations are `3.20, 3.18, 3.32, 3.31`, and the audit's long-line limit at `e = 20–22` is
    `3.282`. Conjecturally the per-line expectation tends to a constant near `3.28`: the plain 2- and 3-corrected
    model's `3 ln(2 ln 2/ln 3)/(2 ln 2 - ln 3) = 2.43`, times a singular series near `1.35`.
  - Example (PROVED): no point of a genuine Fermat line, and not `F_e` itself, is divisible by 7. Since
    `3/4 = -1 (mod 7)`, the line alternates mod 7 between `2^(2^e) + 1` and `1 - 2^(2^e)`, and `2^(2^e)` is never
    `+-1 (mod 7)`.
  - Against that, 5 divides a quarter of the points and 13 a sixth.
  - Observed for `e = 6..12`: `2, 3, 2, 0, 2, 2, 6`.

**Non-consequences.** Nothing here proves anything about Fermat primes. The test calibrates a heuristic model class
on nearby integers and does not touch the integers `F_n` themselves.
