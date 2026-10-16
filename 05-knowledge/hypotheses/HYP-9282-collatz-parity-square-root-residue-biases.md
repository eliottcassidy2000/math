---
id: HYP-9282
title: "Square-root residue-class biases of the Collatz stopping-time parity: the wavelet sums H_(q,r)(N) of t(n) = (-1)^tau(n) over residue classes modulo 3^q carry nonzero real main terms b_(q,r) sqrt(N) on many classes; the coefficients would form a fixed point of the rank-one modified critical 3-adic transfer operator of THM-4658 at s = 1/2"
status: >
  REFUTED (2026-10-15, same session). The pre-registered out-of-sample test failed: in the octaves 2^34..2^36 of the
  exact census extended to 2^37, three of the five pre-registered channels changed sign in all or most octaves and
  missed their +-0.6 bands (H_(3,18): +0.04, +0.12, +0.15 against -2.44; H_(3,9): +0.01, -0.35, -0.46 against
  +2.03; H_(5,155): +0.45, +0.45, +0.23 against -1.42). The in-sample significance was an artefact of octave memory:
  adjacent octave means are correlated (pooled lag-one correlation 0.63-0.66 at every level), and a matched null
  with random odd signs and the exact halving relation t(2m) = -t(m) reproduces the out-of-sample behaviour of the
  Collatz channels. FINITE-EXACT data (class sums modulo 243 at 1,727 grid points through 2^37) are retained. The
  independent S37 referee found the same independence error (E2) and a false drift claim (E3).
session: opus-2026-10-15-S37
source: 05-knowledge/results/collatz_octave_memory_bias_energy_budgets_20261015.md
scripts:
  - 04-computation/experiments/collatz_residue_bias_deep_levels_20261015.py (in-sample, output collatz_residue_bias_deep_levels_20261015.json)
  - 04-computation/experiments/collatz_census_extension_2p37_20261015.py, collatz_census_extension_analysis_20261015.py (out-of-sample test)
  - 04-computation/experiments/collatz_bias_null_calibration_2p37_20261015.py (matched halving null at 2^37)
  - 04-computation/experiments/s37_audit_20261015/a4_bias_recompute.py, a4b_bias_sign_persistence.py (referee)
related: [THM-4658, HYP-9283, HYP-9284, THM-4629, THM-4581, HYP-9258]
---

# HYP-9282: square-root residue-class biases of the Collatz parity census (REFUTED)

**Objects.** `t(n) = (−1)^τ(n)`, where `τ` is the number of shortcut steps `T(n) = n/2, (3n+1)/2` from `n` to 1. The wavelet channels are `H_(q,r)(N) = 3 Σ_(n<=N, n≡r (3^q)) t(n) − Σ_(n<=N, n≡r (3^(q−1))) t(n)`, and `H = H_(1,2)` is the mod-3 discrepancy of the amplitude programme. Normalized values are `g_(q,r)(N) = H_(q,r)(N)/sqrt(6N/3^q)`, the independent-sign scale.

## The conjecture as proposed

For many residue classes `H_(q,r)(N) = b_(q,r) sqrt(N) + F_(q,r)(N)`, with `b_(q,r) ≠ 0` and `F_(q,r)` a fluctuation of mean zero on logarithmic scales. The coefficients would then be a fixed point of `P̃_(1/2)` (THM-4658 statement 4(b)). Note that a fluctuation of the same order as the main term violates the `o(N^σ)` hypothesis of statement 4(b); the fixed-point reading needed an unstated logarithmic-mean version (referee, I7).

## Why it was proposed (in-sample, octaves 2^14..2^32)

The octave means of `g_(q,r)` over 19 octaves (64 log-uniform samples each) gave naive t-values (mean over octaves divided by sd/sqrt(19)) with 75 of the 363 channels at `|t| > 5` and 20 at `|t| > 8`. 27 channels kept one sign in all 19 octaves. Examples: `H_(3,18)` −2.44 (t −18.7), `H_(3,9)` +2.03 (+14.5), `H_(5,78)` +1.55 (+14.6), `H_(5,155)` −1.42 (−14.6), `H_(5,59)` +0.92 (+14.1). All of these numbers are exact reproductions (referee).

## What went wrong

1. **Octave memory (E2).** The note justified independence of the 19 octave means by the one-octave point autocorrelation 0.06 of `H`'s profile. That is not the correlation of adjacent octave means (0.28-0.36 for `H`), and it says nothing about the other channels. Measured over 2^14..2^37, the pooled lag-one correlation of octave means is 0.63, 0.65, 0.66, 0.65, 0.65 at levels 1-5 (lag two about 0.30). THM-4658 statement 6 explains a large part: every halving eigen-combination inherits its value from the previous octave exactly. With this memory the naive t-values are inflated about 1.7-2.5 times; an AR(1) correction leaves 7-17 channels at `|t| > 5`. Only 242 of the 363 channels are free (sibling triples sum to zero).
2. **The null was wrong.** "A zero-mean null would produce essentially no such channels" holds for independent octave means only. The matched *halving null* keeps the exact relation `t(2m) = −t(m)` and draws the odd signs at random (splitmix64 hash). At 2^37, with the same statistics on the same octaves, its 12 replicates give 23-49 (median 33) channels at naive `|t| > 5`, 3-11 (median 6) at `|t| > 8`, and 5-15 (median 8) channels of constant sign over 19 octaves; the largest naive `|t|` in a null replicate is 32.1. Collatz has 75, 20 and 27: more in-sample persistence than this null, i.e. longer memory, not a demonstrated main term.
3. **Drift (E3).** "Drift at most about ±2.5% per octave" was false for 5 of the 14 table rows, all drifting toward zero (`H_(2,3)` +9.1%/octave, `H_(2,6)` −8.2%, `H_(1,1)` +7.6%, `H_(4,51)` +4.8%, `H_(3,16)` +4.6%).
4. **Automatic consistency.** The averaged refinement `m_(1,r) ≈ −(m_(1,e(r)) + m_(2,o(r)))/sqrt(2)` and the level-by-level fixed-point check of statement 4(b) (residuals 6-12% of the bias size) hold for every sequence obeying the exact recursion, up to end effects (THM-4658 remark on statement 6). They checked code, not the hypothesis.
5. **Mercator cancellation.** Exact for `H` only (`H(y) = −[S(y) + 3S(⌊2y/3⌋) + 3S(⌊y/3⌋)]` for every `y`); for the other channels the 26.5-cycle amplitude is numerically at most 8e-5 of its amplitude in `S` (referee, I6).

## The out-of-sample test (pre-registered, FAILED)

Census extended exactly to 2^37 (`collatz_census_extension_2p37_20261015.py`: 245 uint64-overflow words recomputed with Python integers; spot checks 0 mismatches; reproduces `S(2^34) = 24953018`, `H(2^34) = 26869` and every grid sample of the 2^34 census). Octave means of the pre-registered channels:

| channel | in-sample mean | predicted | octave 33 | 34 | 35 | 36 | verdict |
|---|---:|---|---:|---:|---:|---:|---|
| `H_(3,18)` | −2.440 | sign −, band ±0.6 | −1.042 | +0.044 | +0.124 | +0.154 | refuted |
| `H_(3,9)` | +2.028 | sign +, band ±0.6 | +0.953 | +0.006 | −0.351 | −0.460 | refuted |
| `H_(5,155)` | −1.421 | sign −, band ±0.6 | +0.186 | +0.454 | +0.445 | +0.226 | refuted |
| `H_(5,78)` | +1.547 | sign +, band ±0.6 | +1.265 | +1.318 | +1.478 | +1.475 | kept |
| `H_(5,59)` | +0.919 | sign +, band ±0.6 | +0.163 | +0.380 | +0.602 | +1.027 | kept |

- Of the 20 channels with naive `|t| > 8`, 7 keep the sign of their out-of-sample mean (octaves 34-36); the slope of out-of-sample on in-sample means is 0.05.
- For the 75 channels with `|t| > 5`: 48 keep sign, slope 0.19. For all 363 channels: slope 0.20, correlation 0.21.
- The halving null gives the same picture (12 replicates: top five 7-15 of 15 octave signs kept, median 10, and no replicate passes the pre-registered criteria; channels with `|t| > 5`: 50-76% keep sign; all-channel slope 0.03-0.46, median 0.23). So the out-of-sample persistence of the Collatz channels is at the level that octave memory alone produces.

## What survives

- The data: class sums modulo 243 at every grid point through 2^37 (FINITE-EXACT), and the in-sample statistics as exact numbers.
- A genuinely Collatz-specific signal is the *size*, not the sign, of the channels: see HYP-9284 (the mod-3 discrepancy is suppressed by an order of magnitude relative to the halving null, and the normalized channel energies at levels 3-5 declined between 2^28 and 2^37).
- The level-1 discrepancy `H` kept its positive sign in the octaves 33-36 (+0.11, +0.11, +0.09, +0.05) at about a third of its in-sample mean +0.305; no claim is attached.
