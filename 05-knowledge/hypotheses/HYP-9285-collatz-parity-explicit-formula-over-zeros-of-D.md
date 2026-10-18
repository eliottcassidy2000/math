---
id: HYP-9285
title: "Explicit formula for the Collatz parity sum: S(N) = sum_(n<=N) (-1)^tau(n) equals Re sum_rho c_rho N^rho over the zeros rho of D(s) = 1 + 2^(-s) + (3/2)^s/3 with 1/2 < Re rho < 1, up to a remainder of order sqrt(N) in mean square, with c_rho = 2 r_rho / rho for the residues r_rho of F(s) = sum t(n) n^(-s); the 28 zeros with 0.65 < Re rho < 1 and Im rho < 3000 predict S out of sample to about 2%, and the least-squares coefficients reproduce the conditionally certified residues"
status: >
  UPDATE 2026-10-17 (opus-2026-10-17-S41, Run C evaluated): the second pre-registered prediction PASSED. Model M2
  (69 zeros of D, refitted on [2^22, 2^42)) predicts the unseen octave [2^42, 2^43) with relative rms error 0.000736
  (bound 0.0035; residual 1.87 sqrt(N)). Run C exactly replicates [2^40, 2^42).
  UPDATE 2026-10-17 (opus-2026-10-17-S39): batch certificates at N = 2^42 (exact block moments through 2^42,
  collatz_pole_certificates_many_zeros_2p42_20261016.py/.json): of the 16 zeros with the largest boundary-residue
  estimate relative to the radius model, 14 have certified nonzero residues under the budget (2, 0) and 4 under
  (64, 16) (CONDITIONAL on the budgets; intervals rigorous given them); the two not certified are 0.998355 + 5787.79i
  and the 1054/665 zero rho* (as expected; pre-registered at 2^43, R4).
  UPDATE 2026-10-16 (opus-2026-10-16-S38): THM-4660 proves the expansion CONDITIONALLY for the second Riesz mean under
  the square-root hypothesis (W) (weaker than every HYP-9283 budget); the unsmoothed statement here remains CONJECTURAL
  and needs vertical-line and horizontal bounds for G (S38 audit E4). THM-4659: the zeros near Re s = 1 are exactly the
  parity-matched near-returns, uniformly in q (prior art: Codex's resonance-pair note, THM-4629 st.3). Second
  pre-registered out-of-sample test (collatz_census_2p42_preregistration_20261016.md, dense grid, [2^40, 2^42)): the
  69-zero model M2 PASSED at relative rms 0.00115 (bound 0.005; holdout 0.00167), residual 1.73 sqrt(N); the S37 model
  PASSED at 0.0205 (bound 0.03). A third (R1, [2^42, 2^43), bound 0.0035) is pre-registered and PENDING. On the
  1024-per-octave grid the high Mercator members alias approximately (1024 is a resonant density); an
  irrational-density census to 2^38 resolves them, and its fit matches the boundary residues across the family
  (rho* 0.01019 vs 0.00999). 295 zeros with 1/2 < Re < 1 below height 20000 (NUMERICAL).
  CONJECTURAL (the expansion). NUMERICAL, from the FINITE-EXACT censuses through 2^37 and 2^40: holdout errors and
  the pre-registered out-of-sample prediction (P1 of collatz_census_2p40_preregistration_20261015.md: PASSED,
  relative rms 0.0195 <= 0.03 on 2^37..2^40); dense data to 2^27 bring the remainder to 0.27 sqrt(N). CONDITIONAL
  (Arb, under the energy budget (2, 0) for dyadic X >= 2^37): nonzero residues certified at four zeros besides the
  Mercator zero (0.98494 + 720.74i, 0.94973 + 1572.56i, 0.86345 + 131.16i, 0.69392 + 22.25i), each with the
  least-squares value near the middle of its interval (two estimators from the same data); the other 23 intervals
  contain 0. The zeros of D with 2 < Re < 3 (F holomorphic there) are not modes. PROVED (elementary, exact for every
  N): H(N) = -[S(N) + 3S(floor(2N/3)) + 3S(floor(N/3))]. Its floor-free counterpart annihilates every mode
  N^rho with D(rho) = 0; the integer formula retains an O_rho(N^(Re(rho)-1)) remainder per mode (MISTAKE-620).
session: opus-2026-10-15-S37
source: 05-knowledge/results/collatz_octave_memory_bias_energy_budgets_20261015.md
scripts:
  - 04-computation/experiments/collatz_mercator_waveform_prediction_20261015.py (fit and pre-registered prediction)
  - 04-computation/experiments/collatz_census_2p40_evaluation_20261015.py (P1 evaluation)
  - 04-computation/experiments/collatz_pole_certificate_other_zeros_20261015.py, collatz_other_zeros_residue_table_20261015.py (conditional residues at other zeros; output collatz_other_zeros_residue_table_20261015.json)
related: [HYP-9283, HYP-9284, THM-4658, collatz_amplitude_boundary_20261013c.md, HYP-9260]
---

# HYP-9285: the parity sum obeys an explicit formula over the zeros of D

**Setting.** `t(n) = (−1)^τ(n)` on the basin of 1, `S(N) = Σ_(n<=N) t(n)`, `F(s) = Σ t(n) n^(−s)`, and
`D(s) = 1 + 2^(−s) + (3/2)^s/3`. Theorem B3 of the amplitude note (renewal_shapes_oct11) shows that an energy budget
for the mod-3 discrepancy `H` continues `F` meromorphically to `Re s > 1/2`, with poles only at zeros of `D`, and
certifies a nonzero residue at the Mercator zero `s_0 ≈ 0.99821 + 240.248i`. The amplitude note lists a full modal
expansion as OPEN.

## Conjecture

`S(N) = Re Σ_(D(ρ)=0, 1/2 < Re ρ < 1) (2 r_ρ/ρ) N^ρ + R(N)`, with the sum over conjugate pairs taken once and
`Σ_(X<=N<2X) R(N)² = O(X²)` (or `O(X² log X)`). The sum converges in a suitable averaged sense; contour shifting would
need growth bounds for `F` on vertical lines, which are not available.

## Leading-mode cancellation with the integer remainder retained (PROVED)

For every `N`, `H(N) = −[S(N) + 3S(⌊2N/3⌋) + 3S(⌊N/3⌋)]` (exact; S37 referee, 0 failures for `N < 3·2^22`). Since
`1 + 3(2/3)^ρ + 3(1/3)^ρ = 3(2/3)^ρ D(ρ)`, the floor-free combination annihilates every mode `N^ρ` with `D(ρ) = 0`.
For the actual integer formula the residue is
`3[floor(2N/3)^ρ-(2N/3)^ρ]+3[floor(N/3)^ρ-(N/3)^ρ]=O_ρ(N^(Re ρ-1))`.
It vanishes when 3 divides N, but equals 1 at N=1 for Re ρ>0. Thus H retains the floor residues of each finite
modal packet as well as the transformed remainder of S. An infinite modal sum needs summability/truncation bounds;
the energy budget is not already proved equivalent to a budget for an untransformed explicit-formula remainder.
See [the floor-error repair, MISTAKE-620](../results/collatz_conductor_wiener_20261015d.md#4-integer-floors-leave-a-modal-residue).
The exact native identity and the conditional pole certificates are unchanged. The halving-null comparison is
the scoped numerical evidence of HYP-9284, not a proof of a complete modal expansion.

## Evidence (NUMERICAL)

- **Zeros.** `D` has 44 zeros with `1/2 < Re ρ < 1` and `0 < Im ρ < 3000`. The ones with the largest fitted
  contributions are `s_0` (26.48 cycles per octave), `0.8635 + 131.16i` (14.47 cycles per octave) and
  `0.6939 + 22.25i` (2.45 cycles per octave). These are the three lines already visible in the periodogram of `S(N)/N`
  in the S37 census (26.5, 14.49 and 2.45 cycles per octave).
- **Holdouts (fit on `2^22..2^34`, test on `2^34..2^37`, relative rms of `S` on the grid of 64 points per octave):**
  single zero `s_0`: 0.217; the 12 zeros with `0.85 < Re ρ < 1`: 0.075; the 28 zeros with `0.65 < Re ρ < 1`: 0.0205; all 44
  zeros with `1/2 < Re ρ < 1`: 0.0186 (in-sample 0.0078, mild overfitting). The remaining error is consistent with zeros
  above `Im ρ = 3000`, whose oscillations alias on this grid.
- **Pre-registered out-of-sample test (P1).** The 28-zero model fitted on `2^22..2^37` predicts `S` on `[2^37, 2^40)`;
  pass if the relative rms error is at most 0.03. Result: PASSED, 0.0195 (192 grid points; refitted on the same data the Mercator zero alone gives 0.165 and the 12 zeros with `Re > 0.85` give 0.042). `S(2^40)` = 989,143,994 against the prediction 961,770,189.
- **Dense test** (`collatz_explicit_formula_dense_test_20261015.py`). The 64-point grid aliases every zero above 32 cycles per octave. On exact data through 2^27 sampled at 2048 points per octave (no aliasing below `Im ρ = 8000`), with the fit on `2^18..2^25` and the test on `2^25..2^27`:

  | zeros used | holdout relative rms | holdout residual / `sqrt(N)` (rms) |
  |---|---:|---:|
  | 3 (`Im < 300`, `Re > 0.65`) | 0.124 | 1.40 |
  | 10 (`Im < 1000`) | 0.082 | 0.93 |
  | 28 (`Im < 3000`) | 0.034 | 0.37 |
  | 47 (`Im < 5000`) | 0.029 | 0.31 |
  | 75 (`Im < 8000`) | 0.026 | 0.28 |
  | 119 (`Im < 8000`, `Re > 1/2`) | 0.025 | 0.27 |

  The residual falls to the square-root scale and keeps falling as zeros are added. If the remainder `R = S − (modal sum)` were uncorrelated across the three scales of `H = −[R(N) + 3R(2N/3) + 3R(N/3)]` (exact up to the per-mode floor residues of MISTAKE-620), its rms would be `rms(H/sqrt N)/sqrt(10) = 0.24 sqrt(N)` on this range (`rms(H/sqrt N) = 0.75`); the measured 0.27 is close (HEURISTIC comparison).
- **Residues, two ways.** A pole with residue `r` gives `S(N) ~ Re(2r N^ρ/ρ)`, so the least-squares coefficient gives
  `|r_fit| = |c_ρ| |ρ|/2`. At `s_0`: `|r_fit| = 0.18413`, inside the conditional Arb interval `[0.17784, 0.19011]` at
  `N = 2^37` and next to the boundary-numerator value `0.1840` of B3. Over the 28 zeros (`collatz_other_zeros_residue_table_20261015.json`):

  | zero `ρ` | cycles per octave | conditional `|r|`, budget (2, 0) beyond 2^37 | `|r_fit|` |
  |---|---:|---|---:|
  | `0.99821 + 240.248i` | 26.50 | [0.1778, 0.1901] | 0.1841 |
  | `0.98494 + 720.740i` | 79.51 | [0.0451, 0.0951] | 0.0649 |
  | `0.94973 + 1572.556i` | 173.48 | [0.0778, 0.3294] | 0.2056 |
  | `0.86345 + 131.156i` | 14.47 | [0.4168, 0.5959] | 0.5094 |
  | `0.69392 + 22.253i` | 2.455 | [1.127, 4.044] | 2.557 |

  For the other 23 zeros the tail bound (proportional to `|ρ| N^(−(Re ρ − 1/2))`) is too large to fix a sign, and the fitted values lie inside those wide intervals.

## Tests

- P1 (above), and extensions of the census.
- Leading D-zero contributions cancel; small floor-error contributions are retained. The measured 26.5-cycle
  amplitude in H is at most 8e-5 of that in S for channels at levels 1-3 (S37 referee); this finite observation is
  not exact annihilation of the integer-clock mode.
