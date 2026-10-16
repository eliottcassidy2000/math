---
id: HYP-9285
title: "Explicit formula for the Collatz parity sum: S(N) = sum_(n<=N) (-1)^tau(n) equals Re sum_rho c_rho N^rho over the zeros rho of D(s) = 1 + 2^(-s) + (3/2)^s/3 with Re rho > 1/2, up to a remainder of order sqrt(N) in mean square, with c_rho = 2 r_rho / rho for the residues r_rho of F(s) = sum t(n) n^(-s); the 28 zeros with Re rho > 0.65 and Im rho < 3000 predict S out of sample to about 2%, and the least-squares coefficients reproduce the conditionally certified residues"
status: >
  CONJECTURAL (the expansion). NUMERICAL, from the FINITE-EXACT census through 2^37; the 2^40 extension is pending: holdout and pre-registered
  out-of-sample prediction errors (P1 of collatz_census_2p40_preregistration_20261015.md: PENDING-P1). CONDITIONAL
  (Arb, under the energy budget (2, 0) for dyadic X >= 2^37): nonzero residues certified at four zeros besides the
  Mercator zero (0.98494 + 720.74i, 0.94973 + 1572.56i, 0.86345 + 131.16i, 0.69392 + 22.25i), each with the
  least-squares value near the middle of its interval; all 28 least-squares values lie inside their (mostly wide)
  intervals. PROVED (elementary, exact for every
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

`S(N) = Re Σ_(D(ρ)=0, Re ρ > 1/2) (2 r_ρ/ρ) N^ρ + R(N)`, with the sum over conjugate pairs taken once and
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
  single zero `s_0`: 0.19; the 12 zeros with `Re ρ > 0.85`: 0.075; the 28 zeros with `Re ρ > 0.65`: 0.0205; all 44
  zeros with `Re ρ > 1/2`: 0.0186 (in-sample 0.006, mild overfitting). The remaining error is consistent with zeros
  above `Im ρ = 3000`, whose oscillations alias on this grid.
- **Pre-registered out-of-sample test (P1).** The 28-zero model fitted on `2^22..2^37` predicts `S` on `[2^37, 2^40)`;
  pass if the relative rms error is at most 0.03. Result: PENDING-P1.
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
