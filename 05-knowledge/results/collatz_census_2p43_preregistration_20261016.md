# Pre-registration: Collatz parity census 2^42..2^43, replication of 2^40..2^42, and the certificate at the 1054/665 zero (S38, second stage)

2026-10-16, opus-2026-10-16-S38. Written and pushed after the resumed census of `[2^40, 2^42)` (Run B) had finished and
been evaluated ([first-stage pre-registration](collatz_census_2p42_preregistration_20261016.md),
[results](collatz_census_2p42_evaluation_20261016.json): Q1, Q2, Q3, Q5, Q6 passed, Q4 failed), and BEFORE the checked
census of `[2^40, 2^43)` (Run C) was started. Notation as before. The evaluation grid is `N_j = round(2^(10 + j/1024))`
restricted to `[2^42, 2^43)` (1024 points).

## The run

`collatz_census_dense_resumable_checked_20261016.py 43 37 1024 <out> --resume dense40_state.json --check-SH 989143994 -55691`:
the checked version of the census script (checksum partners for the weighted moments and the energy sums, grid
Lipschitz checks, strict resume check; MISTAKE-623). It recomputes `[2^40, 2^42)` independently of Run B and extends to
`2^43`.

## Predictions

[collatz_explicit_formula_prediction_2p43_20261016.json](collatz_explicit_formula_prediction_2p43_20261016.json)
(script `collatz_explicit_formula_next_octave_prediction_20261016.py`): the first stage's model M2 (the 69 zeros of
`D(s) = 1 + 2^(−s) + (3/2)^s/3` with `1/2 < Re ρ < 1`, `Im ρ < 4640.1`) refitted on the dense grid of `[2^22, 2^42)`.
Its holdout (fit `2^22..2^40`, test `2^40..2^42`) is the first stage's Q2 value: relative rms 0.00115, residual
`1.73 sqrt(N)`.

## Tests

- **R1:** relative rms error of the M2 prediction over the 1024 grid points of `[2^42, 2^43)` at most **0.0035** (three
  times the holdout).
- **R2:** `e_X = E_X/X² <= 2` for `X = 2^42`.
- **R3:** the boundary residue estimate at `ρ* = 0.99999921 + 3014.0202i` before the future tail (S37 certificate code,
  moments through 2^43) has modulus in `[0.0098901, 0.0100899]` (within 1% of 0.00999; the values at
  `N = 2^32, ..., 2^42` lie in `[0.009969, 0.009995]`).
- **R4:** the conditional Arb certificate (`collatz_pole_certificate_other_zeros_20261015.py`) at `ρ*` with exact data
  through `N = 2^43` and the budget `E_X <= 2X²` for dyadic `X >= 2^43` (row `(a, b) = (2, 0)`) gives a positive lower
  bound for `|r|`. The expected radius is about 0.0092 against `|r| ≈ 0.00999`, a margin of about 8% (the radius at 2^42
  is 0.0130, so the certificate fails there: `[−0.00304, 0.0230]`).
- **R5:** replication: Run C reproduces Run B exactly on `[2^40, 2^42)` (all 2048 grid values `S`, `H`, the exact prefix
  sums of `H` and `H²`, and every block-moment row).

The evaluation script `collatz_census_2p43_evaluation_20261016.py` is committed with this file. Any failure will be
reported as such.

## Notes recorded before the new data

- Run B is a single run; its integrity checks (grid Lipschitz bounds, monotone prefix sums, moment bounds, moment total
  against `S`) all pass ([output](collatz_census_dense_2p42_integrity_20261016.out)). R5 tests it fully.
- R4 is conditional on the unproved budget; a positive lower bound would give `limsup log B(X)/log X >= 0.99999921`
  under that budget (THM-4659, THM-4660, HYP-9283), extending the S37 value 0.99821. It needs exact moments through
  `2^43`, which is why Run C recomputes `[2^40, 2^42)` with checksummed accumulators instead of resuming from Run B.
