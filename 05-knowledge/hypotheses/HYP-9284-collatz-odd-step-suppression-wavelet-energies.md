---
id: HYP-9284
title: "Odd-step suppression of the Collatz wavelet discrepancies: relative to the halving null (random odd signs with t(2m) = -t(m)), the mod-3 discrepancy H of t(n) = (-1)^tau(n) carries 3-30 times less energy in sample and 25-318 times less on 2^33..2^36, and decorrelates within about one octave, and the deep wavelet energies, which match the null up to about 2^28, declined from there and lie below every null replicate on 2^30..2^37; conjecture: sup_X E_X/X^2 < infinity for H (a constant energy budget holds), with the suppression produced by the odd-step relation t(n) = -t((3n+1)/2)"
status: >
  OPEN (the conjecture). NUMERICAL against 12 halving-null replicates on the FINITE-EXACT census through 2^37
  (min-median-max quoted); pre-registered tests on 2^37..2^40: P2, P3, P5, P6 PASSED, P4 FAILED (the quiet stretch
  ended at 2^38). PROVED
  (THM-4658 statement 7): on densities of any finite 3-adic depth the excess operator has eigenvalue 1 only on
  Re s = 0. HEURISTIC: the link to the suppression (excess eigen-distributions of unbounded depth exist for every s;
  second referee E7).
session: opus-2026-10-15-S37
source: 05-knowledge/results/collatz_octave_memory_bias_energy_budgets_20261015.md
scripts:
  - 04-computation/experiments/collatz_bias_null_calibration_2p37_20261015.py (null replicates and identical Collatz statistics)
  - 04-computation/experiments/collatz_census_extension_2p37_20261015.py, collatz_census_extension_2p40_20261015.py (censuses)
  - 04-computation/experiments/collatz_census_2p40_evaluation_20261015.py (P2-P6)
  - 04-computation/experiments/collatz_excess_operator_resonance_scan_20261015.py (no finite-level excess modes)
related: [HYP-9282, HYP-9283, HYP-9285, THM-4658]
---

# HYP-9284: the odd-step relation suppresses the Collatz wavelet discrepancies

**Setting.** `t(n) = (−1)^τ(n)`, `H(N) = Σ_(n<=N) (3[n ≡ 2 (3)] − 1) t(n)`, `E_X = Σ_(X<=n<2X) H(n)²`, and the wavelet channels `g_(q,r) = H_(q,r)/sqrt(6N/3^q)` of HYP-9282 (mean square about 1 for independent signs). The **halving null** is `t_null(n) = (−1)^(v_2(n)) s(m)` with `m` the odd part of `n` and `s(m) = ±1` from a splitmix64 hash: it keeps `t(2m) = −t(m)` exactly, hence every halving identity of THM-4658 (statement 6), and has no odd-step structure.

## Conjecture

`sup_X E_X/X² < ∞` for the Collatz `H`, and the mechanism is the odd-step relation: any argument that uses only the halving structure and independence of the odd signs must fail, because it would apply to the halving null. A stronger question (not conjectured): `E_X/X² → 0`, i.e. square-root cancellation fails in favour of more cancellation in the wavelet channels.

## Evidence (NUMERICAL, same statistics on the same octaves)

| statistic | Collatz | halving null (12 replicates: min-median-max) |
|---|---:|---|
| profile of `H`: sd of `H/sqrt(2N)` | 0.27 | 0.63-1.27-1.95 |
| profile of `H`: one-octave autocorrelation | 0.13 | 0.65-0.84-0.91 |
| window energy of `H`, mean on `2^14..2^34` | 0.51 | 1.5-7.1-14.4 |
| window energy of `H`, maximum on `2^14..2^34` | 1.85 | 5.5-30.1-42.7 |
| window energy of `H`, mean on `(2^33, 2^36]` | 0.07 | 1.8-7.1-22.9 |
| level-1 channel energy (in sample) | 0.28 | 0.40-1.62-3.52 |
| level energies 3, 4, 5 (in sample) | 1.09, 0.95, 0.98 | 0.67-1.33, 0.88-1.11, 0.84-1.11 |
| level energies 3, 4, 5: octaves 34-36 over 14-32 | 0.20, 0.46, 0.65 | min 0.52, 0.77, 0.78 |
| level energies 3, 4, 5: octaves 30-36 over 14-29 | 0.34, 0.48, 0.58 | min 0.68, 0.77, 0.79 |

("Below every replicate" is a rank among 12, p ≈ 1/13 per level, on a split chosen after seeing the data; the size of the gap, e.g. 0.20 against a null minimum of 0.52 at level 3, is the stronger evidence. The pre-registered P5 is the clean test.)
| level 2: octaves 30-36 over 14-29 | 0.26 | min 0.22 (not outside the null) |

- Collatz's exact dyadic energies never exceeded 1.77 for `X >= 4` (maximum at 2^24) and were 0.04-0.10 on 2^33..2^36; the halving null's in-sample window means exceed 2 in 11 of 12 replicates.
- The leading-mode reading of HYP-9285 explains part of the level-1 picture: `H(N) = −[S(N) + 3S(⌊2N/3⌋) + 3S(⌊N/3⌋)]` cancels the leading D-zero power but retains an O_ρ(N^(Re ρ-1)) floor residue per mode (MISTAKE-620). The finite fit of Collatz S has about 2% holdout error; that does not prove a complete modal expansion or exact absence of modes in a null sequence.
- The deep-level decline started at about the same octave for levels 3-5 (all peak or plateau near 2^27 and fall from 2^28 on) and is weaker at deeper levels (octave 31: 0.62, 0.44, 0.60; octave 36: 0.24, 0.43, 0.62); whether this is the onset of sub-square-root growth or a slow modulation is open.

## Pre-registered tests (2^37..2^40; `collatz_census_2p40_preregistration_20261015.md`)

- P2: `e_X <= 2` for `X = 2^37, 2^38, 2^39`. P3: `e_X <= 1` for the same windows. P4 (about even odds): their mean below 0.25. P5: level energies 3, 4, 5 averaged over octaves 37-39 below 0.8. P6: the octave mean of `H/sqrt(2N)` positive in at least two of the octaves 37-39.
- Results: P2 PASSED and P3 PASSED (`e_X` = 0.069, 0.931, 0.078). P4 FAILED (mean 0.359 > 0.25). P5 PASSED (0.50, 0.71, 0.71) with the levels recovering (octave 39: 0.69, 0.81, 0.76). P6 PASSED (−0.01, +0.50, +0.14).
- Reading: the deep-level dip of 2^28..2^37 looks like a slow modulation, so the stronger question `E_X/X² → 0` gets no support; the level-1 suppression against the halving null (an order of magnitude) persists through 2^39.
