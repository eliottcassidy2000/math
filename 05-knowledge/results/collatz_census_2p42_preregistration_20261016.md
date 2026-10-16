# Pre-registration: Collatz parity census 2^40..2^42 on a dense grid (S38)

2026-10-16, opus-2026-10-16-S38. Written and pushed after the dense census through 2^40 had finished and BEFORE the
resumed census of `[2^40, 2^42)` was started (`collatz_census_dense_resumable_20261016.py 42 37 1024 ... --resume`).
Notation as in [the S37 note](collatz_octave_memory_bias_energy_budgets_20261015.md): `t(n) = (−1)^τ(n)`,
`S(N) = Σ_(n<=N) t(n)`, `H(N) = Σ_(n<=N) (3[n ≡ 2 (3)] − 1) t(n)`, `e_X = E_X/X²`. The evaluation grid is
`N_j = round(2^(10 + j/1024))` restricted to `[2^40, 2^42)` (2048 points).

## Models (fitted before the new data exist)

Predictions: [collatz_explicit_formula_dense_preregistration_20261016.json](collatz_explicit_formula_dense_preregistration_20261016.json)
(script `04-computation/experiments/collatz_explicit_formula_dense_preregistration_20261016.py`).

- **M1:** S37's model (the 28 zeros of `D(s) = 1 + 2^(−s) + (3/2)^s/3` with `0.65 < Re ρ < 1`, `Im ρ < 3000`),
  refitted on the 64-point grid of `[2^22, 2^40)`. Holdout (fit `2^22..2^37`, test `2^37..2^40`): 0.0195 on the
  coarse grid, 0.0191 on the dense grid.
- **M2:** the 69 zeros with `1/2 < Re ρ < 1` and `Im ρ < 4640.1` (the Nyquist height of 1024 points per octave), chosen as the best of twelve candidate zero sets by the holdout below, refitted on the dense grid of `[2^22, 2^40)`. Holdout on the dense grid: 0.00167 (in-sample 0.00061) (residual
  `1.10 sqrt(N)` rms).

## Tests

- **Q1:** relative rms error of M1 over the 2048 grid points at most 0.03.
- **Q2:** relative rms error of M2 at most 0.005 (three times its holdout error).
- **Q3:** `e_X <= 2` for `X = 2^40` and `X = 2^41` (the budget).
- **Q4:** `e_X <= 1` for both.
- **Q5:** mean window energy of `H` over the window starts in `[2^40, 2^41)` (dense windows `[N_j, N_(j+1024))`)
  below 1 (every halving-null replicate of S37 had window means at least 1.5).
- **Q6:** rms of the M2 residual divided by `sqrt(N)` over the 2048 points at most 3.0.

The evaluation script `collatz_census_2p42_evaluation_20261016.py` was committed with this file. Any failure will be
reported as such.

## Notes recorded before the new data

- The holdout error drops from 0.0118 to 0.0017 when the zeros with `4200 < Im ρ < 4640` enter. On a log-uniform grid
  with an integer number of points per octave, every zero with `2^(−i Im ρ) ≈ −1` (the Mercator family near `Re s = 1`,
  frequencies `k + 1/2` cycles per octave) aliases onto another member of the family: on this grid `5307.3 → 3975.0`,
  `5787.8 → 3494.5`, `6268.3 → 3014.0`, `6748.8 → 2533.5`, `8561.6 → 720.7`, `9042.1 → 240.2`. The fitted
  coefficients of those zeros therefore mix several modes, and the near-Nyquist zeros partly absorb aliased power. The
  test is on the same grid formula, so the aliasing is the same in and out of sample; the individual coefficients of
  the family members above `Im ≈ 2500` are not residue estimates.
- Fitted against boundary residues (2^40 block moments) for the ten largest low modes: `s_0` 0.18427 / 0.18398,
  `0.8635 + 131.16i` 0.5081 / 0.5065, `0.9849 + 720.74i` 0.0648 / 0.0701, `0.9497 + 1572.56i` 0.2025 / 0.2036,
  `0.6939 + 22.25i` 2.579 / 2.570.
- The first run of the dense 2^40 census ended in an internal inconsistency (final class counts against the chained
  block increments, two bit-level discrepancies) and was discarded; the rerun, with per-block invariants, has 0
  inconsistent blocks, gives `S(2^40) = 989143994`, `H(2^40) = −55691` (as the S37 census) and agrees exactly with the
  S37 64-point census at all 1,920 shared grid points.

## Errata (added after the S38 audit, before the evaluation; criteria and predictions unchanged)

- The aliasing note above overstates the scope. A partner exists only when the density makes `2P log2 3` close to an
  even integer; 1024 is such a density (`1024 log2 3 = 1623.0016`), 1000 is not. The aliasing is approximate (common
  detuning 0.00386 cycles per octave). THM-4659 statement 5 is corrected accordingly.
- "the ten largest low modes" lists five.
- `Q2_bound` and `Q6_bound` were written into the prediction JSON by hand from the bounds stated in this file; the
  committed prediction script does not write them. The evaluation script reads them from the JSON.
- Run B was first launched at 11:43:48 and relaunched at 14:20:09 after the session restart, with the census script's
  only post-registration edit (the density argument parsed as a float). The evaluation grid is provably unchanged
  (integer density 1024, same rounding).
