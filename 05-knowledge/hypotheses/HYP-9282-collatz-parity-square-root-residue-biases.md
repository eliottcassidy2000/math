---
id: HYP-9282
title: "Square-root residue-class biases of the Collatz stopping-time parity: the wavelet sums H_(q,r)(N) of t(n) = (-1)^tau(n) over residue classes modulo 3^q carry nonzero real main terms b_(q,r) sqrt(N) on many classes (75 of the 363 channels of levels 1..5 have octave-mean |t| > 5 over the 19 octaves 2^14..2^33, 20 have |t| > 8); the coefficients would form a fixed point of the critical 3-adic transfer operator of THM-4658 at s = 1/2"
status: >
  OPEN. FINITE-EXACT: all class sums modulo 243 at 1,473 grid points through 2^33 and all exact octave means of the
  level-1/2 channels through 2^34. NUMERICAL: octave-mean statistics, t-values and drifts (normalized by the
  independent-sign scale sqrt(N 6/3^q)). The asymptotic claim (a constant main term with exponent exactly 1/2), the
  fixed-point reading and any mechanism are CONJECTURAL. A brief literature search found no prior report.
session: opus-2026-10-15-S37
source: 05-knowledge/results/collatz_octave_memory_bias_energy_budgets_20261015.md
scripts:
  - 04-computation/experiments/collatz_residue_bias_deep_levels_20261015.py (output 05-knowledge/results/collatz_residue_bias_deep_levels_20261015.json)
  - 04-computation/experiments/collatz_octave_memory_census_20261015.py (output collatz_octave_memory_census_20261015.json/.npz)
related: [THM-4658, HYP-9283, THM-4629, THM-4581, HYP-9258]
---

# HYP-9282: square-root residue-class biases of the Collatz parity census

**Objects.** `t(n) = (−1)^τ(n)`, where `τ` is the number of shortcut steps `T(n) = n/2, (3n+1)/2` from `n` to 1. The wavelet channels are `H_(q,r)(N) = 3 Σ_(n<=N, n≡r (3^q)) t(n) − Σ_(n<=N, n≡r (3^(q−1))) t(n)`, and `H = H_(1,2)` is the mod-3 discrepancy of the amplitude programme. For independent random signs, `H_(q,r)(N)` has standard deviation `sqrt(6N/3^q)`; normalized values `g_(q,r)(N) = H_(q,r)(N)/sqrt(6N/3^q)` are reported.

## Conjecture

For many residue classes the channel has a nonzero real main term:

`H_(q,r)(N) = b_(q,r) sqrt(N) + F_(q,r)(N)`, with `b_(q,r) ≠ 0` and `F_(q,r)` a fluctuation of mean zero on logarithmic scales.

By THM-4658 statement 4(b), the coefficient distribution `β` on `Z_3` (`β(r + 3^q Z_3) ↔` the class excesses) would be a fixed point of the rank-one modified transfer operator `P̃_(1/2)`, at the line `Re s = 1/2` where the two branches of the 3-adic transfer operator are exactly balanced in `L²(Z_3)`.

## Evidence (FINITE-EXACT data, NUMERICAL statistics)

The octave means of `g_(q,r)` over the octaves `2^14..2^33` (19 octaves, 64 log-uniform samples each):

| channel | mean of `g` | sd over octaves | `t` | octaves with the sign of the mean | drift per octave |
|---|---:|---:|---:|---:|---:|
| `H_(1,2)` (= `H`) | +0.305 | 0.212 | +6.3 | 17/19 | 0.000 |
| `H_(1,1)` | −0.483 | 0.448 | −4.7 | 17/19 | +0.037 |
| `H_(2,3)` | −1.188 | 0.788 | −6.6 | 19/19 | +0.108 |
| `H_(2,6)` | +1.094 | 0.845 | +5.6 | 17/19 | −0.089 |
| `H_(3,9)` | +2.028 | 0.610 | +14.5 | 19/19 | +0.050 |
| `H_(3,18)` | −2.440 | 0.568 | −18.7 | 19/19 | −0.039 |
| `H_(3,16)` | −1.218 | 0.529 | −10.0 | 19/19 | +0.056 |
| `H_(4,51)` | −1.508 | 0.625 | −10.5 | 19/19 | +0.073 |
| `H_(4,55)` | −1.409 | 0.568 | −10.8 | 18/19 | −0.027 |
| `H_(4,22)` | +1.052 | 0.444 | +10.3 | 19/19 | +0.007 |
| `H_(5,78)` | +1.547 | 0.461 | +14.6 | 19/19 | −0.038 |
| `H_(5,155)` | −1.421 | 0.425 | −14.6 | 19/19 | −0.005 |
| `H_(5,59)` | +0.919 | 0.285 | +14.1 | 19/19 | +0.002 |
| `H_(5,79)` | +1.010 | 0.346 | +12.7 | 19/19 | −0.023 |

- **Overall.** 75 of the 363 channels of levels 1–5 have `|t| > 5`, and 20 have `|t| > 8`. The normalized profile of `H` decorrelates within about one octave (lag-one autocorrelation 0.06; HYP-9283), so the 19 octave means are close to independent. A zero-mean null would produce essentially no such channels.
- **Exponent.** The strongest channels drift by at most about ±2.5% per octave relative to their means, consistent with exponent `1/2 ± 0.04` on this range.
- **Splits of the multiples of `3^j`.** `D_j(N) = Σ_(3^j ∥ n) χ(n/3^j) t(n)`, with `χ` the nontrivial character mod 3, has mean `D_j/sqrt(N) ≈ −0.37, −0.62, +0.70, −0.02` for `j = 0, 1, 2, 3` (`t = −5.9, −6.3, +17.4, −0.8`). There is no simple sign alternation.
- **Consistency with the exact recursion.** The level-1 octave-mean averages satisfy the averaged refinement `m_(1,r) ≈ −(m_(1,e(r)) + m_(2,o(r)))/sqrt(2)` (THM-4658 statement 2) to within the octave scatter:
  - `+0.296` against the predicted `+0.281` for `r = 2`;
  - `−0.461` against `−0.489` for `r = 1`;
  - `+0.162` against `+0.211` for `r = 0` (census values, 2^14..2^34).

  This is necessary for a stationary main term, not evidence of one by itself.

## Pre-registered tests

- For the octaves `2^34..2^40` (outside both censuses), the octave means of `g_(3,18)`, `g_(3,9)`, `g_(5,155)`, `g_(5,78)`, `g_(5,59)` keep the signs `−, +, −, +, +`.
- Their averages stay within `±0.6` of `−2.44, +2.03, −1.42, +1.55, +0.92`.
- A sign change in a majority of these octaves would refute the constant-main-term reading for that channel.

## What would explain it, and what would not

- The white-sign model of `t` gives `b = 0` for every channel. So nonzero `b_(q,r)` measure a deterministic structure of the Collatz parity at exactly the `sqrt(N)` scale.
- The Mercator resonances (zeros of `D(s)` near `Re s = 1`) cancel exactly from every wavelet channel, as they must: `H(y) = −[S(y) + 3S(2y/3) + 3S(y/3)]` up to boundary terms. So the biases are not leakage of the density oscillation.
- A finite set of small sources contributes `O(1)`, not `sqrt(N)`.
- The analogy with Chebyshev's bias (a real `sqrt(x)` term from prime squares, S32) is suggestive; no mechanism of that type is identified here (ANALOGY).
