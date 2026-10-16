# Pre-registration: Collatz parity census 2^37..2^40 (S37)

2026-10-15, opus-2026-10-15-S37. Written and pushed while the streaming census
`04-computation/experiments/collatz_census_extension_2p40_20261015.py` was running and before any of its output
beyond the parity-table spot check had been read. The census extends the exact census of
`collatz_census_extension_2p37_20261015.npz` (through 2^37) to 2^40. Notation as in
[the S37 note](collatz_octave_memory_bias_energy_budgets_20261015.md): `t(n) = (−1)^τ(n)`, `S(N) = Σ_(n<=N) t(n)`,
`H(N) = Σ_(n<=N) (3[n ≡ 2 (3)] − 1) t(n)`, `e_X = E_X/X² = Σ_(X<=n<2X) H(n)²/X²`, `g = H/sqrt(2N)`.

## P1. Explicit-formula waveform of S (main test)

Model: `S(N) ≈ Re Σ_ρ c_ρ N^ρ` over the 28 zeros `ρ` of `D(s) = 1 + 2^(−s) + (3/2)^s/3` with `Re ρ > 0.65` and
`0 < Im ρ < 3000`, coefficients fitted by least squares on the grid points of `[2^22, 2^37)`.
The model was chosen on holdouts before this file was written: fitted on `2^22..2^34` it predicts `2^34..2^37` to
2.05% relative rms, against 19% for the single Mercator zero `s_0 ≈ 0.99821 + 240.248i` and 7.5% for the 12 zeros
with `Re ρ > 0.85`.
Predictions are stored in
[collatz_mercator_waveform_prediction_20261015.json](collatz_mercator_waveform_prediction_20261015.json)
(script `04-computation/experiments/collatz_mercator_waveform_prediction_20261015.py`). Among them:
`S(2^38) ≈ +2.738e8`, `S(2^39) ≈ −5.108e8`, `S(2^40) ≈ +9.618e8`.

**Pass criterion:** relative rms error of the predicted `S` over the grid points of `[2^37, 2^40)` at most 0.03.

## P2–P6. Energies and channels

- **P2 (budget):** `e_X <= 2` for the dyadic windows `X = 2^37, 2^38, 2^39`.
- **P3 (suppression against the halving null):** `e_X <= 1` for the same three windows. The halving-null replicates
  at 2^37 have window means 4.6–22.9 on `(2^33, 2^36]`.
- **P4 (quiet phase continues; genuinely uncertain, about even odds):** the mean of the three `e_X` is below 0.25.
- **P5 (multi-level decline persists):** the mean square of the normalized wavelet channels
  `g_(q,r) = H_(q,r)/sqrt(6N/3^q)`, averaged over the octaves 37–39 and over the channels of a level, is below 0.8 for
  each of the levels 3, 4, 5. The in-sample (2^14..2^33) level averages are about 1.1.
- **P6 (sign of the mod-3 discrepancy):** the octave mean of `g` is positive in at least two of the octaves 37, 38, 39.

Any failure will be reported as such in the S37 note and the hypothesis files.

## Results (added after the census finished; evaluation `collatz_census_2p40_evaluation_20261015.py`)

- **P1 PASSED:** relative rms error 0.0195 (192 grid points); `S(2^38)`, `S(2^39)`, `S(2^40)` = 281,675,896, −525,856,530, 989,143,994 against the predicted 273,757,038, −510,780,478, 961,770,189.
- **P2 PASSED, P3 PASSED:** `e_X` = 0.069, 0.931, 0.078.
- **P4 FAILED:** the mean is 0.359; the quiet stretch ended with a coarse-mean excursion at 2^38.
- **P5 PASSED:** 0.50, 0.71, 0.71 for levels 3, 4, 5 (recovering: octave 39 gives 0.69, 0.81, 0.76).
- **P6 PASSED:** octave means of `g` = −0.01, +0.50, +0.14.
